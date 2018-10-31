package adnascreen;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import htsjdk.samtools.fastq.FastqReader;

/**
 * Make a pass through one lane of Illumina sequencer output. 
 * For reads that have Reichlab indices and barcodes that look like the input sets, 
 * output to a series of size-matched files and count them. 
 * This performs merging of forward and reverse reads, assuming a minimum overlap
 * and a minimum resulting length. 
 * Adapters and barcodes are trimmed. 
 */
public class IndexAndBarcodeScreener {
	
	public static final String RAW = "raw";
	public static final String MERGED = "merged";
	public static final String OLIGO = "oligo";

	public static void main(String []args) throws IOException, ParseException{
		CommandLineParser parser = new DefaultParser();
		Options options = new Options();
		options.addRequiredOption("i", "i5-indices", true, 
				"File containing one valid i5 index sets per line");
		options.addRequiredOption("j", "i7-indices", true, 
				"File containing one valid i7 index sets per line");
		options.addRequiredOption("b", "barcodes", true, 
				"File containing one valid barcodes sets per line with ':'-delimited elements");
		options.addOption("m", "mismatch-penalty-max", true, "Max allowable penalty for mismatch while aligning for merge");
		options.addOption("j", "mismatch-penalty-high", true, "Penalty for mismatch at high-quality bases while aligning for merge");
		options.addOption("k", "mismatch-penalty-low", true, "Penalty for mismatch at non-high-quality bases while aligning for merge");
		options.addOption("q", "mismatch-quality-threshold", true, "Threshold for determining penalty for mismatch while aligning merge");
		options.addOption("o", "minimum-overlap", true, "Minimum bases of overlap for paired reads to merge");
		options.addOption("l", "minimum-length", true, "Minimum length for output merged reads");
		options.addOption("n", "number-output-files", true, "Number of output files to divide merged reads between");
		options.addOption("h", "hamming-distance", true, "Max hamming distance for index or barcode match");
		options.addOption("r", "read-group-file", true, "Output file for read group");
		options.addOption("c", "barcode-count", true, "File containing prior pass's counts of keys with(out) barcodes by index pair to determine whether to demultiplex with barcodes");
		options.addOption("t", "barcode-threshold", true, "Threshold for count to use barcode length over no barcodes");
		options.addOption("x", "index-barcode-keys", true, "Index-barcode keys for setting explicit barcode lengths");
		options.addOption("z", "positive-oligo", true, "Provide count for reads matching provided positive oligo sequence");
		CommandLine commandLine	= parser.parse(options, args);
		
		BarcodeMatcher i5Indices = null, i7Indices = null;
		BarcodeMatcher barcodes = null;
		// We keep statistics for each 4-tuple of indices and barcodes
		SampleSetsCounter sampleSetCounter = new SampleSetsCounter();
		final int maxPenalty = Integer.valueOf(commandLine.getOptionValue("mismatch-penalty-max", "3"));
		final int mismatchPenaltyHigh = Integer.valueOf(commandLine.getOptionValue("mismatch-penalty-high", "3"));
		final int mismatchPenaltyLow = Integer.valueOf(commandLine.getOptionValue("mismatch-penalty-low", "1"));
		final int mismatchBaseQualityThreshold = Integer.valueOf(commandLine.getOptionValue("mismatch-penalty-threshold", "20"));
		final int minOverlap = Integer.valueOf(commandLine.getOptionValue('o', "15"));
		final int minMergedLength = Integer.valueOf(commandLine.getOptionValue('l', "30"));
		final int numOutputFiles = Integer.valueOf(commandLine.getOptionValue('n', "25"));
		final int maxHammingDistance = Integer.valueOf(commandLine.getOptionValue('h', "1"));
		final float barcodeToNoBarcodeThreshold = Float.valueOf(commandLine.getOptionValue("barcode-threshold", "0.05"));
		String readGroupFilename = commandLine.getOptionValue("read-group-file", "read_group");
		
		try{
			i5Indices = new BarcodeMatcher(commandLine.getOptionValue("i5-indices"), maxHammingDistance);
			i7Indices = new BarcodeMatcher(commandLine.getOptionValue("i7-indices"), maxHammingDistance);
			barcodes = new BarcodeMatcher(commandLine.getOptionValue('b'), maxHammingDistance);
		} catch(IOException e){
			System.exit(1);
		}
		
		// optional specification of barcode lengths from index-barcode key file
		Map<IndexAndBarcodeKey, Integer> barcodeLengthsFromSampleSheet = null;
		String explicitIndexFile = commandLine.getOptionValue("index-barcode-keys", null);
		if (explicitIndexFile != null) {
			barcodeLengthsFromSampleSheet = barcodeLengthsByIndexPair(explicitIndexFile, barcodes);
		}
		
		// A previous pass through the data is needed to count the number of paired reads
		// that demultiplex with barcodes
		final String barcodeCountStatisticsFilename = commandLine.getOptionValue("barcode-count", null);
		SampleSetsCounter barcodeCountStatistics = null;
		if(barcodeCountStatisticsFilename != null){
			File barcodeCountStatisticsFile = new File(barcodeCountStatisticsFilename);
			barcodeCountStatistics = new SampleSetsCounter(barcodeCountStatisticsFile);
		}
		Map<IndexAndBarcodeKey, Integer> barcodeLengthByIndexPairCache = new HashMap<IndexAndBarcodeKey, Integer>();
		
		// positive oligo
		// count the number of appearances of this sequence
		// This is for wetlab diagnostic purposes
		String positiveOligoSequence = commandLine.getOptionValue("positive-oligo", null);
		Read positiveOligo = null;
		Read positiveOligoReverseComplement = null;
		if(positiveOligoSequence != null) {
			String qualityString = String.join("", Collections.nCopies(positiveOligoSequence.length(), "I")); // oligo sequence is max quality
			positiveOligo = new Read("", positiveOligoSequence, qualityString);
			positiveOligoReverseComplement = positiveOligo.reverseComplement();
		}

		String[] remainingArgs = commandLine.getArgs();
		PrintWriter [] fileOutputs = new PrintWriter[numOutputFiles];
		try(
				FileInputStream r1File = new FileInputStream(remainingArgs[0]);
				FileInputStream r2File = new FileInputStream(remainingArgs[1]);
				FileInputStream i1File = new FileInputStream(remainingArgs[2]);
				FileInputStream i2File = new FileInputStream(remainingArgs[3]);

				FastqReader r1Reader = new FastqReader(new BufferedReader(new InputStreamReader(new GZIPInputStream(r1File))));
				FastqReader r2Reader = new FastqReader(new BufferedReader(new InputStreamReader(new GZIPInputStream(r2File))));
				FastqReader i1Reader = new FastqReader(new BufferedReader(new InputStreamReader(new GZIPInputStream(i1File))));
				FastqReader i2Reader = new FastqReader(new BufferedReader(new InputStreamReader(new GZIPInputStream(i2File))));
				){
			// prepare output files for multiple parallel processing jobs downstream
			// for load balancing purposes, these are not demultiplexed
			String outputFilenameRoot = remainingArgs[4];
			for(int i = 0; i < numOutputFiles; i++){
				// start counting from 1 for filenames
				String outputFilename = String.format("%s_%03d.fastq.gz", outputFilenameRoot, i + 1);
				fileOutputs[i] = new PrintWriter(new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outputFilename)))));
			}
			
			String readGroup = null;
			int pairedReadOutputCount = 0;
			while(r1Reader.hasNext() && r2Reader.hasNext() && i1Reader.hasNext() && i2Reader.hasNext()){
				Read r1 = new Read(r1Reader.next());
				Read r2 = new Read(r2Reader.next());
				Read i1 = new Read(i1Reader.next());
				Read i2 = new Read(i2Reader.next());
				
				// Lookup by index pair whether barcodes are used
				IndexAndBarcodeKey keyIndexOnly = MergedRead.findExperimentKey(r1, r2, i1, i2, 
						i5Indices, i7Indices, null, 0);
				int barcodeLength = -1;
				if(barcodeCountStatistics != null && keyIndexOnly != null){
					// We assume that for a given index pair, barcodes are all the same length
					// We first use the index-barcode key file, if given, then
					// We use the 4-tuple with maximum count to determine the barcode length for this index pair
					if(barcodeLengthByIndexPairCache.containsKey(keyIndexOnly)){
						barcodeLength = barcodeLengthByIndexPairCache.get(keyIndexOnly);
					} else{
						if (barcodeLengthsFromSampleSheet != null && barcodeLengthsFromSampleSheet.containsKey(keyIndexOnly)) {
							barcodeLength = barcodeLengthsFromSampleSheet.get(keyIndexOnly);
						} else {
							barcodeLength = barcodeLengthFromPriorPassCounts(barcodeCountStatistics, keyIndexOnly, barcodes, barcodeToNoBarcodeThreshold);
						}
						barcodeLengthByIndexPairCache.put(keyIndexOnly, barcodeLength);
					}
				}
				
				// update key if barcodes are used, otherwise reuse the index pair
				IndexAndBarcodeKey key = (barcodeLength > 0) ? MergedRead.findExperimentKey(r1, r2, i1, i2, 
						i5Indices, i7Indices, barcodes, barcodeLength) : keyIndexOnly;
				
				IndexAndBarcodeKey keyFlattened = null;
				sampleSetCounter.increment(); // statistics recording
				MergedRead merged = null;
				if(key != null){
					keyFlattened = key.flatten();
					int r1BarcodeLength = barcodes.getBarcodeLength(keyFlattened.getP5Label());
					int r2BarcodeLength = barcodes.getBarcodeLength(keyFlattened.getP7Label());
					sampleSetCounter.increment(keyFlattened, RAW);
					merged = MergedRead.mergePairedSequences(r1, r2, key, 
							r1BarcodeLength, r2BarcodeLength, maxPenalty, minOverlap, minMergedLength,
							mismatchPenaltyHigh, mismatchPenaltyLow, mismatchBaseQualityThreshold);
					// read group consistency
					String readGroupForThisRead = r1.getFASTQHeader().getReadGroupElements();
					if(readGroup == null){
						readGroup = readGroupForThisRead;
					} else { // read groups are expected to match for all reads in lane
						if(!readGroup.equals(readGroupForThisRead)){
							throw new IllegalStateException();
						}
					}
					// count occurrences of positive oligo
					if(positiveOligo != null && merged != null) {
						if(Read.alignmentAssessment(positiveOligo, merged, 0, 0, positiveOligo.length(), maxPenalty, mismatchPenaltyHigh, mismatchPenaltyLow, mismatchBaseQualityThreshold)
							|| Read.alignmentAssessment(positiveOligoReverseComplement, merged, 0, 0, positiveOligo.length(), maxPenalty, mismatchPenaltyHigh, mismatchPenaltyLow, mismatchBaseQualityThreshold)) {
							sampleSetCounter.increment(keyFlattened, OLIGO);
						}
					}
				}
				// output to file and more statistics recording
				// only reads that pass Illumina's pass filter (PF) [aka chastity filter]
				if(merged != null && !merged.getFASTQHeader().isFiltered()){
					// separate into different files
					fileOutputs[pairedReadOutputCount % numOutputFiles].println(merged.toString());
					pairedReadOutputCount++;
					sampleSetCounter.increment(keyFlattened, MERGED);
				}
			}
			// output map statistics
			PrintStream statisticsOutput = System.out;
			statisticsOutput.println(sampleSetCounter.toStringSorted(RAW));
			// output read group
			if(readGroupFilename != null){
				try(PrintWriter readGroupFile = new PrintWriter(readGroupFilename)){
					readGroupFile.println(readGroup);
				}
			}
		} catch(IOException e){
			System.err.println(e);
			System.exit(1);
		} finally {
			for(int i = 0; i < numOutputFiles; i++){
				if(fileOutputs[i] != null){
					fileOutputs[i].close();
				}
			}
		}
	}
	
	/**
	 * Find the barcode length using counts from a previous pass through the index pair data. 
	 * The top barcode length is chosen, unless that is 0 (non-barcoded). If 0, the top non-zero barcode is checked, 
	 * and if that barcode is a substantial fraction of the non-barcoded count, the non-zero barcode length is returned. 
	 * @param barcodeCountStatistics
	 * @param keyIndexOnly
	 * @param barcodes
	 * @param threshold Use non-zero barcode count if (non-zero barcode count > (non-barcoded count * threshold))
	 * @return 
	 */
	public static int barcodeLengthFromPriorPassCounts(SampleSetsCounter barcodeCountStatistics, 
			IndexAndBarcodeKey keyIndexOnly, BarcodeMatcher barcodes, float threshold){
		if(barcodeCountStatistics == null)
			throw new IllegalArgumentException("No counts");
		else if (keyIndexOnly == null)
			throw new IllegalArgumentException("No index pair");

		SampleCounter countsForKeyIndexOnly = barcodeCountStatistics.get(keyIndexOnly);
		List<String> barcodePairStrings = countsForKeyIndexOnly.getLabelList();
		// sort barcode pairs by count, tail is highest count
		SortedMap<Integer, String> barcodePairsByCount = new TreeMap<Integer, String>();
		for(String barcodePairString : barcodePairStrings){
			barcodePairsByCount.put(countsForKeyIndexOnly.get(barcodePairString), barcodePairString);
		}
		
		int barcodeLength = -1;
		Integer max;
		max = barcodePairsByCount.lastKey();
		String maxPairLabel = barcodePairsByCount.remove(max);
		barcodeLength = barcodes.getBarcodePairLength(maxPairLabel);
		
		if(barcodeLength == 0 && barcodePairsByCount.size() > 0){
			// if a minority of barcodes pass the barcode check, we can falsely assume there are no barcodes
			// check whether a substantial fraction of barcodes are for a single barcode pair
			Integer second = barcodePairsByCount.lastKey();
			if(second > threshold * max){
				String secondPairLabel = barcodePairsByCount.remove(second);
				int secondLength = barcodes.getBarcodePairLength(secondPairLabel);
				if(secondLength > 0){
					barcodeLength = secondLength;
				} else {
					throw new IllegalStateException("Two index-barcode keys without barcodes. There should be at most one.");
				}
			}
		}
		
		return barcodeLength;
	}
	
	/**
	 * Use the content from a file containing the index and barcode keys for samples to determine the barcode lengths for index pairs
	 * Each index pair should have only one barcode length associated with it. 
	 * This function assumes that barcode and index base pair sequences are directly in the input file, 
	 * and that the IndexBarcodeKey uses these directly without other labels. 
	 * @param explicitIndexFile
	 * @param barcodes
	 * @return
	 * @throws IOException
	 * @throws ParseException
	 */
	public static Map<IndexAndBarcodeKey, Integer> barcodeLengthsByIndexPair(String explicitIndexFile, BarcodeMatcher barcodes) throws IOException, ParseException{
		HashMap<IndexAndBarcodeKey, Integer>  barcodeLengths = new HashMap<IndexAndBarcodeKey, Integer>();
		try(BufferedReader reader = new BufferedReader(new FileReader(explicitIndexFile))){
			String entryLine;
			while((entryLine = reader.readLine()) != null && entryLine.length() > 0){
				String [] fields = entryLine.split("\t");
				String keyString = fields[0];
				String [] basePairSequenceGroups = keyString.split(String.valueOf(IndexAndBarcodeKey.FIELD_SEPARATOR));
				
				String i5 = basePairSequenceGroups[0];
				String i7 = basePairSequenceGroups[1];
				IndexAndBarcodeKey indexOnlyKey = new IndexAndBarcodeKey(i5, i7, null, null);
				
				String stringP5 = "";
				String stringP7 = "";
				if (basePairSequenceGroups.length > 2)
					stringP5 = basePairSequenceGroups[2].split(String.valueOf(BarcodeMatcher.BARCODE_DELIMITER))[0];
				if (basePairSequenceGroups.length > 3)
					stringP7 = basePairSequenceGroups[3].split(String.valueOf(BarcodeMatcher.BARCODE_DELIMITER))[0];
				DNASequence singleP5 = new DNASequence(stringP5);
				DNASequence singleP7 = new DNASequence(stringP7);
				int length1 = singleP5.length();
				int length2 = singleP7.length();
				if(length1 == length2) {
					if(barcodeLengths.containsKey(indexOnlyKey)){
						// check that length matches any existing
						int existingLength = barcodeLengths.get(indexOnlyKey);
						if(existingLength != length1) {
							throw new IllegalStateException("barcode length mismatch for multiple reads on same index pair");
						}
					}
					// check that barcodes are valid
					if(length1 > 0 && (barcodes.find(singleP5) == null || barcodes.find(singleP7) == null)){
						throw new IllegalStateException("barcode not found: " + singleP5.toString() + " " + singleP7.toString());
					}
					// store barcode length for this index pair
					barcodeLengths.put(indexOnlyKey, length1);
				}
				else
					throw new IllegalStateException("barcode length mismatch on single read");
			}
		}
		return barcodeLengths;
	}
}
