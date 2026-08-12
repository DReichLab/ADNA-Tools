package adnatools;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.help.HelpFormatter;
import org.apache.commons.lang3.mutable.MutableInt;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.cram.ref.CRAMReferenceSource;
import htsjdk.samtools.cram.ref.ReferenceSource;

/**
 * Mark or remove PCR/optical duplicates in coordinate-sorted SAM/BAM/CRAM files.
 * Duplicates share the same reference, alignment start, read length, strand, optional tag (typically XD),
 * and library (LB read group tag). Reads in different libraries are not marked as duplicates unless
 * those libraries are grouped together with --library-group.
 * The read with the highest average base quality is retained; others are marked duplicate.
 */
public class MarkDuplicates {
	private final String duplicateTag;
	private final boolean useTagInKey;
	private final boolean removeDuplicates;
	private final LibraryDeduplicationScope libraryScope;
	private long distinctReadGroups = 0;
	private long totalDuplicateReads = 0;
	private final Map<Integer, MutableInt> duplicateCountHistogram = new HashMap<Integer, MutableInt>();
	private int maximumDuplicateCount = 0;

	public MarkDuplicates(String duplicateTag, boolean useTagInKey, boolean removeDuplicates,
			LibraryDeduplicationScope libraryScope) {
		this.duplicateTag = duplicateTag;
		this.useTagInKey = useTagInKey;
		this.removeDuplicates = removeDuplicates;
		this.libraryScope = libraryScope;
	}

	/**
	 * Maps each read to a deduplication scope from its read group library (LB).
	 * By default each library is deduplicated independently; explicit groups combine libraries.
	 */
	static class LibraryDeduplicationScope {
		private final Map<String, String> libraryToScope;

		LibraryDeduplicationScope() {
			libraryToScope = new HashMap<String, String>();
		}

		LibraryDeduplicationScope(String[] groupSpecs) {
			libraryToScope = new HashMap<String, String>();
			if (groupSpecs == null) {
				return;
			}
			for (int groupIndex = 0; groupIndex < groupSpecs.length; groupIndex++) {
				String scopeId = Integer.toString(groupIndex);
				for (String library : parseLibrariesFromGroupSpec(groupSpecs[groupIndex])) {
					library = library.trim();
					if (library.isEmpty()) {
						continue;
					}
					if (libraryToScope.containsKey(library)) {
						throw new IllegalArgumentException(
								"Library '" + library + "' appears in more than one library group");
					}
					libraryToScope.put(library, scopeId);
				}
			}
		}

		String scopeFor(SAMRecord record) {
			String library = getLibrary(record);
			if (libraryToScope.containsKey(library)) {
				return libraryToScope.get(library);
			}
			return library;
		}

		// separate based on spaces or commas
		static String[] parseLibrariesFromGroupSpec(String groupSpec) {
			return groupSpec.split("[\\s,]+");
		}
	}

	static String getLibrary(SAMRecord record) {
		SAMReadGroupRecord readGroup = record.getReadGroup();
		if (readGroup == null) {
			return "";
		}
		String library = readGroup.getLibrary();
		return library != null ? library : "";
	}

	public static void main(String[] args) throws ParseException, IOException {
		CommandLineParser parser = new DefaultParser();
		Options options = new Options();
		options.addRequiredOption("i", "input", true, "Input SAM/BAM/CRAM filename");
		options.addRequiredOption("o", "output", true, "Output SAM/BAM/CRAM filename");
		options.addOption("r", "reference", true, "CRAM reference: required when any CRAM file is used");
		options.addOption("t", "tag", true, "Optional SAM tag to include in duplicate key (default: XD). Use with --no-tag to omit.");
		options.addOption(null, "no-tag", false, "Do not include a SAM tag in the duplicate key");
		options.addOption("s", "summary", true, "Optional summary file with duplicate statistics");
		options.addOption("R", "remove-duplicates", false, "Remove duplicate reads from output instead of marking them");
		options.addOption("b", "BAM", false, "Write BAM output (also inferred from .bam extension)");
		options.addOption("c", "compression", true, "HTSJDK compression level for BAM/CRAM (0-9, default 5)");
		options.addOption(Option.builder("g").longOpt("library-group").hasArg()
				.desc("Library names (LB) to deduplicate together, separated by commas or spaces; repeatable for multiple groups")
				.get());

		CommandLine commandLine;
		try {
			commandLine = parser.parse(options, args);
		} catch (ParseException e) {
			HelpFormatter formatter = HelpFormatter.builder().get();
			formatter.printHelp("MarkDuplicates", "Mark or remove duplicates in coordinate-sorted alignments.", options, Driver.submittingHelp, true);
			System.err.println(e.getLocalizedMessage());
			throw e;
		}

		String inputFilename = commandLine.getOptionValue("input");
		String outputFilename = commandLine.getOptionValue("output");
		boolean useTagInKey = !commandLine.hasOption("no-tag");
		String duplicateTag = DemultiplexSAM.duplicatesSAMTag;
		if (commandLine.hasOption("tag")) {
			duplicateTag = commandLine.getOptionValue("tag");
		}
		boolean removeDuplicates = commandLine.hasOption("remove-duplicates");
		String summaryFilename = commandLine.getOptionValue("summary");
		int compression = Integer.valueOf(commandLine.getOptionValue("compression", "5"));
		boolean useBAM = commandLine.hasOption("BAM") || Driver.isBAMFilename(outputFilename);

		boolean isCramPresent = AlignmentComparison.isCramFilename(inputFilename)
				|| AlignmentComparison.isCramFilename(outputFilename);
		String cramReference = commandLine.getOptionValue("reference");
		if (isCramPresent && cramReference == null) {
			throw new RuntimeException("Missing CRAM reference file");
		}

		LibraryDeduplicationScope libraryScope = new LibraryDeduplicationScope(
				commandLine.getOptionValues("library-group"));
		MarkDuplicates markDuplicates = new MarkDuplicates(duplicateTag, useTagInKey, removeDuplicates, libraryScope);
		markDuplicates.markDuplicates(inputFilename, outputFilename, cramReference, compression, useBAM);

		if (summaryFilename != null) {
			markDuplicates.writeSummary(summaryFilename);
		}
	}

	public void markDuplicates(String inputFilename, String outputFilename, String reference, int compression,
			boolean useBAM) throws IOException {
		SamReaderFactory samReaderFactory = SamReaderFactory.makeDefault();
		CRAMReferenceSource cramReference = null;
		if (reference != null) {
			Path referencePath = Paths.get(reference);
			cramReference = new ReferenceSource(referencePath);
			samReaderFactory.referenceSource(cramReference);
		}

		SamReader reader;
		if (AlignmentComparison.isCramFilename(inputFilename)) {
			reader = SamReaderFactory.makeDefault().referenceSource(cramReference)
					.open(SamInputResource.of(new BufferedInputStream(new FileInputStream(inputFilename))));
		} else {
			reader = samReaderFactory.open(SamInputResource.of(new BufferedInputStream(new FileInputStream(inputFilename))));
		}

		SAMFileHeader header = reader.getFileHeader();
		SAMFileWriterFactory outputFileFactory = new SAMFileWriterFactory();
		outputFileFactory.setCompressionLevel(compression);
		BufferedOutputStream outputFile = new BufferedOutputStream(new FileOutputStream(outputFilename));
		SAMFileWriter output;
		if (AlignmentComparison.isCramFilename(outputFilename)) {
			output = outputFileFactory.makeCRAMWriter(header, outputFile, new java.io.File(reference));
		} else if (useBAM || outputFilename.toLowerCase().endsWith(".bam")) {
			output = outputFileFactory.makeBAMWriter(header, true, outputFile);
		} else {
			output = outputFileFactory.makeSAMWriter(header, true, outputFile);
		}

		try (SAMRecordIterator iterator = reader.iterator()) {
			List<SAMRecord> recordsAtPosition = new ArrayList<SAMRecord>();
			String currentReference = null;
			int currentStart = -1;

			while (iterator.hasNext()) {
				SAMRecord record = iterator.next();
				if (isEligibleForDuplicateMarking(record)) {
					String referenceName = record.getReferenceName();
					int alignmentStart = record.getAlignmentStart();
					if (currentReference != null
							&& (!referenceName.equals(currentReference) || alignmentStart != currentStart)) {
						writeProcessedPosition(recordsAtPosition, output);
						recordsAtPosition.clear();
					}
					currentReference = referenceName;
					currentStart = alignmentStart;
					recordsAtPosition.add(record);
				} else {
					if (!recordsAtPosition.isEmpty()) {
						writeProcessedPosition(recordsAtPosition, output);
						recordsAtPosition.clear();
						currentReference = null;
						currentStart = -1;
					}
					output.addAlignment(record);
				}
			}
			if (!recordsAtPosition.isEmpty()) {
				writeProcessedPosition(recordsAtPosition, output);
			}
		} finally {
			output.close();
			reader.close();
		}
	}

	static boolean isEligibleForDuplicateMarking(SAMRecord record) {
		return !record.getReadUnmappedFlag()
				&& !record.isSecondaryOrSupplementary();
	}

	static String duplicateKey(SAMRecord record, String duplicateTag, boolean useTagInKey,
			LibraryDeduplicationScope libraryScope) {
		StringBuilder key = new StringBuilder();
		key.append(record.getReferenceName()).append('\t');
		key.append(record.getAlignmentStart()).append('\t');
		key.append(record.getReadLength()).append('\t');
		key.append(record.getReadNegativeStrandFlag());
		if (useTagInKey) {
			key.append('\t');
			String tagValue = record.getStringAttribute(duplicateTag);
			key.append(tagValue != null ? tagValue : "");
		}
		key.append('\t').append(libraryScope.scopeFor(record));
		return key.toString();
	}

	static double averageBaseQuality(SAMRecord record) {
		byte[] qualities = record.getBaseQualities();
		if (qualities == null || qualities.length == 0) {
			return 0.0;
		}
		long sum = 0;
		for (byte quality : qualities) {
			sum += (quality & 0xff);
		}
		return (double) sum / qualities.length;
	}

	private void writeProcessedPosition(List<SAMRecord> recordsAtPosition, SAMFileWriter output) {
		Map<String, List<SAMRecord>> groups = new HashMap<String, List<SAMRecord>>();
		for (SAMRecord record : recordsAtPosition) {
			String key = duplicateKey(record, duplicateTag, useTagInKey, libraryScope);
			if (!groups.containsKey(key)) {
				groups.put(key, new ArrayList<SAMRecord>());
			}
			groups.get(key).add(record);
		}

		for (List<SAMRecord> group : groups.values()) {
			markDuplicateGroup(group);
		}

		for (SAMRecord record : recordsAtPosition) {
			if (!removeDuplicates || !record.getDuplicateReadFlag()) {
				output.addAlignment(record);
			}
		}
	}

	/**
	 * Read with the highest average base quality is the group representative. 
	 * All other reads are flagged as duplicates. 
	 * @param group
	 */
	private void markDuplicateGroup(List<SAMRecord> group) {
		distinctReadGroups++;

		SAMRecord bestRecord = group.get(0);
		double bestAverageQuality = averageBaseQuality(bestRecord);
		for (int i = 1; i < group.size(); i++) {
			SAMRecord candidate = group.get(i);
			double candidateAverageQuality = averageBaseQuality(candidate);
			if (candidateAverageQuality > bestAverageQuality
					|| (candidateAverageQuality == bestAverageQuality
							&& candidate.getReadName().compareTo(bestRecord.getReadName()) < 0)) {
				bestRecord = candidate;
				bestAverageQuality = candidateAverageQuality;
			}
		}

		int duplicateCount = group.size() - 1;
		if (duplicateCount > 0) {
			totalDuplicateReads += duplicateCount;
		}
		maximumDuplicateCount = Math.max(maximumDuplicateCount, duplicateCount);
		if (!duplicateCountHistogram.containsKey(duplicateCount)) {
			duplicateCountHistogram.put(duplicateCount, new MutableInt(0));
		}
		duplicateCountHistogram.get(duplicateCount).increment();

		for (SAMRecord record : group) {
			record.setDuplicateReadFlag(record != bestRecord);
		}
	}

	public long getDistinctReadGroups() {
		return distinctReadGroups;
	}

	public long getTotalDuplicateReads() {
		return totalDuplicateReads;
	}

	public int[] getDuplicateCountHistogram() {
		int[] histogram = new int[maximumDuplicateCount + 1];
		for (int duplicateCount = 0; duplicateCount <= maximumDuplicateCount; duplicateCount++) {
			histogram[duplicateCount] = duplicateCountHistogram.getOrDefault(duplicateCount, new MutableInt(0)).intValue();
		}
		return histogram;
	}

	public void writeSummary(String summaryFilename) throws IOException {
		try (BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(summaryFilename)))) {
			writer.write("distinct_aligned_reads\t" + distinctReadGroups);
			writer.newLine();
			writer.write("total_duplicate_reads\t" + totalDuplicateReads);
			writer.newLine();
			writer.write("duplicate_count\tread_groups");
			writer.newLine();
			int[] histogram = getDuplicateCountHistogram();
			for (int duplicateCount = 0; duplicateCount < histogram.length; duplicateCount++) {
				writer.write(duplicateCount + "\t" + histogram[duplicateCount]);
				writer.newLine();
			}
		}
	}
}
