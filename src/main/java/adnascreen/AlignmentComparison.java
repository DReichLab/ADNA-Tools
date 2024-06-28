package adnascreen;

import java.io.BufferedInputStream;
import java.io.FileInputStream;
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import htsjdk.samtools.SAMFlag;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.cram.ref.CRAMReferenceSource;
import htsjdk.samtools.cram.ref.ReferenceSource;

/**
 * Compare the contents of a sorted alignment file to a series of sorted input source files. 
 * @author mmah
 *
 */
public class AlignmentComparison {
	public static final int MATCH = 0;
	public static final int NO_MATCH = 1;
	
	public static int main(String[] args) throws IOException, ParseException{
		CommandLineParser parser = new DefaultParser();
		Options options = new Options();
		options.addRequiredOption("c", "check", true, "SAM/BAM/CRAM file to check. This program checks that all reads are present in this file.");
		options.addRequiredOption("i", "input", true, "Input SAM/BAM/CRAM file(s). This program checks that all reads in the input files are present in the output file.");
		options.addOption("r", "reference", true, "CRAM reference: if any cram files are used, the reference must be specified and the same for all CRAM files.");
		CommandLine commandLine	= parser.parse(options, args);
		
		String checkFilename = commandLine.getOptionValue("check");
		String[] inputFilenames = commandLine.getOptionValues("input");
		
		// if there is a CRAM file, we need to make sure there is a reference
		boolean isCramPresent = false;
		for (String inputFilename : inputFilenames) {
			if (inputFilename.toLowerCase().endsWith(".cram")) {
				isCramPresent = true;
				break;
			}
		}
		if (!isCramPresent) {
			isCramPresent = checkFilename.toLowerCase().endsWith(".cram");
		}
		if (isCramPresent && !commandLine.hasOption("reference")) {
			throw new RuntimeException("Missing CRAM reference file");
		}
		CRAMReferenceSource cramReference = null;
		if (commandLine.hasOption("reference")) {
			String referenceFilename = commandLine.getOptionValue("reference");
			Path referencePath = Paths.get(referenceFilename);
			cramReference = new ReferenceSource(referencePath);
		}
		
		String[] tagsToCheck = commandLine.getArgs();
		boolean match = compareAlignmentFiles(checkFilename, inputFilenames, tagsToCheck, cramReference);
		return match ? MATCH : NO_MATCH;
	}
	
	/**
	 * 
	 * @param x
	 * @param y
	 * @return true if names match
	 */
	public static boolean compareNames(SAMRecord x, SAMRecord y) {
		if (x != null && y != null) {
			return x.getReadName().equals(y.getReadName());
		} else if (x == null && y == null) {
			throw new IllegalArgumentException("comparing null alignments");
		}
		return false; // exactly one SAM record is null
	}
	
	public static boolean checkTag(SAMRecord x, SAMRecord y, String tag) {
		Object xTag = x.getAttribute(tag);
		Object yTag = y.getAttribute(tag);
		if(xTag != null) {
			return xTag.equals(yTag);
		} else if (yTag != null) {
			return yTag.equals(xTag);
		}
		return true; // both null
	}

	/**
	 * Determine whether two alignments are the same. 
	 * Same is defined by sequence, quality values, flags (except duplicate), alignment, optional tagsToCheck
	 * @param x
	 * @param y
	 * @param tagsToCheck
	 * @return
	 */
	public static boolean compareFullAlignments(SAMRecord x, SAMRecord y, String[] tagsToCheck) {
		if (x != null && y != null) {
			byte[] xSequence = x.getReadBases();
			byte[] ySequence = y.getReadBases();
			byte[] xQuality = x.getBaseQualities();
			byte[] yQuality = y.getBaseQualities();
			
			if (xSequence.length != ySequence.length || xSequence.length != xQuality.length || xQuality.length != yQuality.length)
				return false;
			for (int i = 0; i < xSequence.length; i++) {
				if (xSequence[i] != ySequence[i] || xQuality[i] != yQuality[i])
					return false;
			}
			
			int flag_mask = 0xffff ^ SAMFlag.DUPLICATE_READ.intValue(); // exclude duplicate bit
			if (!x.getReadName().equals(y.getReadName()) ||
					((x.getFlags() & flag_mask) != (y.getFlags() & flag_mask)) ||
					!x.getReferenceName().equals(y.getReferenceName()) ||
					x.getAlignmentStart() != y.getAlignmentStart() ||
					x.getMappingQuality() != y.getMappingQuality() ||
					!x.getCigar().equals(y.getCigar())
					)
				return false;
			for (String tagToCheck : tagsToCheck) {
				if (!checkTag(x, y, tagToCheck))
					return false;
			}
			return true;
		} else if (x == null && y == null) {
			throw new IllegalArgumentException("comparing null alignments");
		}
		return false; // exactly one SAM record is null
	}
	
	public static boolean compareAlignmentFiles(String finalOutputFilename, String[] inputFilenames, String[] tagsToCheck, CRAMReferenceSource reference) throws IOException {
		SamReader outputReader = null;
		SamReader[] inputReaders = new SamReader[inputFilenames.length];
		SamReaderFactory samReaderFactory = SamReaderFactory.makeDefault();
		if (reference != null)
			samReaderFactory.referenceSource(reference);
		
		try {
			// open files
			outputReader = samReaderFactory.open(SamInputResource.of(new BufferedInputStream(new FileInputStream(finalOutputFilename))));
			SAMRecordIterator outputIterator = outputReader.iterator();
			
			SAMRecordIterator[] inputIterators = new SAMRecordIterator[inputFilenames.length];
			SAMRecord[] currentInputRecords = new SAMRecord[inputFilenames.length];
			for (int n = 0; n < inputFilenames.length; n++) {
				inputReaders[n] = samReaderFactory.open(SamInputResource.of(new BufferedInputStream(new FileInputStream(inputFilenames[n]))));
				inputIterators[n] = inputReaders[n].iterator();
				currentInputRecords[n] = inputIterators[n].next();
			}
			
			// each final output file record should appear in exactly one input file
			while(outputIterator.hasNext()) {
				SAMRecord toFind = outputIterator.next();
				boolean found = false;
				for (int n = 0; n < inputIterators.length; n++) {
					if(compareNames(toFind, currentInputRecords[n])) {
						if(compareFullAlignments(toFind, currentInputRecords[n], tagsToCheck)) {
							found = true;
							if (inputIterators[n].hasNext())
								currentInputRecords[n] = inputIterators[n].next();
							else
								currentInputRecords[n] = null;
							break;
						}
					}
				}
				if (!found) {
					return false;
				}
			}
			// there should be no input records remaining
			for (SAMRecord remaining : currentInputRecords) {
				if (remaining != null)
					return false;
			}
			return true;
		}
		finally { // cleanup open files
			if (outputReader != null) {
				outputReader.close();
			}
			for (SamReader reader : inputReaders) {
				if (reader != null)
					reader.close();
			}
		}
	}
}
