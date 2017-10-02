package adnascreen;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.io.PrintStream;
import java.io.PrintWriter;
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
		
		// A previous pass through the data is needed to count the number of paired reads
		// that demultiplex with barcodes
		final String barcodeCountStatisticsFilename = commandLine.getOptionValue("barcode-count", null);
		SampleSetsCounter barcodeCountStatistics = null;
		if(barcodeCountStatisticsFilename != null){
			File barcodeCountStatisticsFile = new File(barcodeCountStatisticsFilename);
			barcodeCountStatistics = new SampleSetsCounter(barcodeCountStatisticsFile);
		}
		Map<IndexAndBarcodeKey, Integer> barcodeLengthByIndexPairCache = new HashMap<IndexAndBarcodeKey, Integer>();

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
					// We use the 4-tuple with maximum count to determine the barcode length for this index pair
					if(barcodeLengthByIndexPairCache.containsKey(keyIndexOnly)){
						barcodeLength = barcodeLengthByIndexPairCache.get(keyIndexOnly);
					} else{
						barcodeLength = barcodeLengthFromPriorPassCounts(barcodeCountStatistics, keyIndexOnly, barcodes, barcodeToNoBarcodeThreshold);
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
}
