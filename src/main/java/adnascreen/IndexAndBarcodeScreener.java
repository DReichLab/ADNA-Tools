package adnascreen;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;
import java.util.concurrent.RejectedExecutionHandler;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
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
	
	private int maxPenalty = 3;
	private int mismatchPenaltyHigh = 3;
	private int mismatchPenaltyLow= 1;
	private int mismatchBaseQualityThreshold = 20;
	private int minOverlap = 15;
	private int minMergedLength = 30;
	
	private int numOutputFiles = 25;
	private float barcodeToNoBarcodeThreshold = 0.05f;
	private String readGroupFilename; 
	
	private BarcodeMatcher i5Indices;
	private BarcodeMatcher i7Indices;
	private BarcodeMatcher barcodes;
	
	boolean reverseComplementI5 = false;
	private Map<IndexAndBarcodeKey, Integer> barcodeLengthsFromSampleSheet = null;
	private SampleSetsCounter barcodeCountStatistics = null;
	
	private Read positiveOligo = null;
	private Read positiveOligoReverseComplement = null;
	
	private PrintStream printStream = System.out;
	
	public IndexAndBarcodeScreener(BarcodeMatcher i5Indices, BarcodeMatcher i7Indices, BarcodeMatcher barcodes){
		this.i5Indices = i5Indices;
		this.i7Indices = i7Indices;
		this.barcodes = barcodes;
	}
	
	class MergeResult{
		boolean isOligo;
		MergedRead merged;
		IndexAndBarcodeKey keyFlattened;
		String readGroup;
		
		MergeResult(IndexAndBarcodeKey keyFlattened, MergedRead merged, boolean isOligo, String readGroup){
			this.keyFlattened = keyFlattened;
			this.merged = merged;
			this.isOligo = isOligo;
			this.readGroup = readGroup;
		}
	}
	
	class SynchronizedOutput{
		PrintWriter [] fileOutputs;
		// We keep statistics for each 4-tuple of indices and barcodes
		SampleSetsCounter sampleSetCounter;
		String readGroup = null;
		
		private int pairedReadOutputCount = 0; // used only for distributing reads across output files
		
		public SynchronizedOutput(int numOutputFiles, String outputFilenameRoot) throws IOException {
			sampleSetCounter = new SampleSetsCounter();
			fileOutputs = new PrintWriter[numOutputFiles];
			
			// prepare output files for multiple parallel processing jobs downstream
			// for load balancing purposes, these are not demultiplexed
			for(int i = 0; i < numOutputFiles; i++){
				// start counting from 1 for filenames
				String outputFilename = String.format("%s_%03d.fastq.gz", outputFilenameRoot, i + 1);
				fileOutputs[i] = new PrintWriter(new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outputFilename)))));
			}
		}
		
		/**
		 * Update sample counters for raw, oligo, merged categories
		 * Write merged read to load balanced files
		 * @param mergeResult
		 */
		public synchronized void updateCountersAndWriteMergeToFile(MergeResult mergeResult) {
			sampleSetCounter.increment(); // statistics recording
			if(mergeResult.keyFlattened != null) {
				sampleSetCounter.increment(mergeResult.keyFlattened.toString(), RAW);
				if(mergeResult.isOligo)
					sampleSetCounter.increment(mergeResult.keyFlattened.toString(), OLIGO);

				if(readGroup == null){
					readGroup = mergeResult.readGroup;
				} else { // read groups are expected to match for all reads in lane
					if(!readGroup.equals(mergeResult.readGroup)){
						throw new IllegalStateException("FASTQ read group mismatch");
					}
				}
				// output to file and more statistics recording
				// only reads that pass Illumina's pass filter (PF) [aka chastity filter]
				if(mergeResult.merged != null && !mergeResult.merged.getFASTQHeader().isFiltered()){
					// separate into different files
					fileOutputs[pairedReadOutputCount % fileOutputs.length].println(mergeResult.merged.toString());
					pairedReadOutputCount++;
					if (pairedReadOutputCount >= fileOutputs.length)
						pairedReadOutputCount -= fileOutputs.length;
					sampleSetCounter.increment(mergeResult.keyFlattened.toString(), MERGED);
				}
			}
		}
		
		public synchronized void cleanup() {
			for(int i = 0; i < fileOutputs.length; i++){
				if(fileOutputs[i] != null){
					fileOutputs[i].close();
					fileOutputs[i] = null; 
				}
			}
		}
	}

	public static void main(String []args) throws IOException, ParseException, InterruptedException, ExecutionException{
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
		options.addOption("y", "reverse-complement-i5", false, "Whether i5 index should be reverse complemented (NextSeq is not)");
		options.addOption(null, "threads", true, "Number of threads for thread pool");
		
		options.addOption(null, "fixed-i5", true, "Assume all fragments have this i5 sequence label");
		options.addOption(null, "fixed-i7", true, "Assume all fragments have this i7 sequence label");
		CommandLine commandLine	= parser.parse(options, args);
		
		final int maxHammingDistance = Integer.valueOf(commandLine.getOptionValue('h', "1"));
		String readGroupFilename = commandLine.getOptionValue("read-group-file", "read_group");
		
		BarcodeMatcher i5Indices = new BarcodeMatcher(commandLine.getOptionValue("i5-indices"), maxHammingDistance);
		BarcodeMatcher i7Indices = new BarcodeMatcher(commandLine.getOptionValue("i7-indices"), maxHammingDistance);
		BarcodeMatcher barcodes = new BarcodeMatcher(commandLine.getOptionValue('b'), maxHammingDistance);
		
		IndexAndBarcodeScreener screener = new IndexAndBarcodeScreener(i5Indices, i7Indices, barcodes);
		screener.setMaxPenalty(Integer.valueOf(commandLine.getOptionValue("mismatch-penalty-max", "3")));
		screener.setMismatchPenaltyHigh(Integer.valueOf(commandLine.getOptionValue("mismatch-penalty-high", "3")));
		screener.setMismatchPenaltyLow(Integer.valueOf(commandLine.getOptionValue("mismatch-penalty-low", "1")));
		screener.setMismatchBaseQualityThreshold(Integer.valueOf(commandLine.getOptionValue("mismatch-penalty-threshold", "20")));
		screener.setMinOverlap(Integer.valueOf(commandLine.getOptionValue('o', "15")));
		screener.setMinMergedLength(Integer.valueOf(commandLine.getOptionValue('l', "30")));
		screener.setNumOutputFiles(Integer.valueOf(commandLine.getOptionValue('n', "25")));
		screener.setBarcodeToNoBarcodeThreshold(Float.valueOf(commandLine.getOptionValue("barcode-threshold", "0.05")));
		screener.setReadGroupFilename(readGroupFilename);
		
		screener.setReverseComplementI5(commandLine.hasOption('y'));
		
		int numThreads = Integer.valueOf(commandLine.getOptionValue("threads", "1"));
		
		// optional specification of barcode lengths from index-barcode key file
		String explicitIndexFile = commandLine.getOptionValue("index-barcode-keys", null);
		if (explicitIndexFile != null) {
			screener.setBarcodeLengthsFromSampleSheet(barcodeLengthsByIndexPair(explicitIndexFile, barcodes));
		}
		
		// A previous pass through the data is needed to count the number of paired reads
		// that demultiplex with barcodes
		final String barcodeCountStatisticsFilename = commandLine.getOptionValue("barcode-count", null);
		if(barcodeCountStatisticsFilename != null){
			File barcodeCountStatisticsFile = new File(barcodeCountStatisticsFilename);
			screener.setBarcodeCountStatistics(new SampleSetsCounter(barcodeCountStatisticsFile));
		}
		
		// positive oligo
		// count the number of appearances of this sequence
		// This is for wetlab diagnostic purposes
		String positiveOligoSequence = commandLine.getOptionValue("positive-oligo", null);
		if(positiveOligoSequence != null) {
			screener.setPositiveOligo(positiveOligoSequence);
		}
		
		String[] remainingArgs = commandLine.getArgs();
		String outputFilenameRoot = remainingArgs[remainingArgs.length-1]; // last argument determines output filenames
		String r1Filename = remainingArgs[0];
		String r2Filename = remainingArgs[1];
		String i1Filename = null;
		String i2Filename = null;
		String i5Label = null;
		String i7Label = null;

		// If there are fixed indices, we use those and do not have index reads
		if(commandLine.hasOption("fixed-i5") && commandLine.hasOption("fixed-i7")) {
			i5Label = commandLine.getOptionValue("fixed-i5");
			i7Label = commandLine.getOptionValue("fixed-i7");
			if(remainingArgs.length != 3) {
				throw new IllegalStateException("Unexpected number of remaining arguments when indices fixed: " + remainingArgs.length);
			}
		}
		else { // processing with index reads
			if(remainingArgs.length != 5) {
				throw new IllegalStateException("Unexpected number of remaining arguments with index reads: " + remainingArgs.length);
			}
			i1Filename = remainingArgs[2];
			i2Filename = remainingArgs[3];
		}
		// start processing
		screener.performScreeningMergeTrim(numThreads, outputFilenameRoot, r1Filename, r2Filename, i1Filename, i2Filename, i5Label, i7Label);
	}

	protected void performScreeningMergeTrim(int numThreads, String outputFilenameRoot, String r1Filename, String r2Filename, String i1Filename, String i2Filename,
			String i5Label, String i7Label) throws IOException, ParseException, InterruptedException, ExecutionException {
		BlockingQueue<Runnable> inputsQueue = new ArrayBlockingQueue<Runnable>(numOutputFiles * numThreads);
		BlockingQueue<Future<MergeResult>> resultsQueue = new ArrayBlockingQueue<Future<MergeResult>>(numOutputFiles * numThreads);
		SynchronizedOutput output = new SynchronizedOutput(numOutputFiles, outputFilenameRoot);
		// if input thread cannot submit new job, run that job in the input thread
		// It would be a better design to block (so there is always a producer thread), 
		// but that is not available from standard java libraries
		RejectedExecutionHandler rejectedExecutionHandler = new ThreadPoolExecutor.CallerRunsPolicy();
		ExecutorService pool = new ThreadPoolExecutor(numThreads, numThreads, 5, TimeUnit.MINUTES, inputsQueue, rejectedExecutionHandler);
				
		// Exactly one thread handles file input to preserve order
		Future<Boolean> inputFuture = pool.submit(new Callable<Boolean>() {
			public Boolean call() throws IOException, InterruptedException{
				try {
					// queue each paired read for merging
					enqueuePairedReads(pool, resultsQueue, r1Filename, r2Filename, i1Filename, i2Filename, i5Label, i7Label);
					return new Boolean(true);
				} catch (IOException | InterruptedException e) {
					throw(e);
				} finally {
					pool.shutdown();
				}
			}
		});
		
		// the main thread serves as the output thread
		while(!pool.isTerminated() || !resultsQueue.isEmpty()) {
			// fetch next result, with timeout
			Future<MergeResult> mergeResultFuture = resultsQueue.poll(1, TimeUnit.SECONDS);
			if(mergeResultFuture != null) {
				try {
					MergeResult mergeResult = mergeResultFuture.get();
					output.updateCountersAndWriteMergeToFile(mergeResult);
				} catch (ExecutionException e) {
					System.err.println(e);
					System.exit(1);
				}
			}
		}
		inputFuture.get(); //  check for input thread exception
		
		// output map statistics
		printStream.println(output.sampleSetCounter.toStringSorted(RAW));
		// output read group
		if(readGroupFilename != null){
			try(PrintWriter readGroupFile = new PrintWriter(readGroupFilename)){
				readGroupFile.println(output.readGroup);
			}
		}
		output.cleanup();
	}
	
	/**
	 * Find the index/barcodes for each read pair, then submit to thread pool to merge
	 * There are two separate cases:
	 * 1. index reads are set at command line (i5Label and i7Label non-null)
	 * 2. index reads are read from fastq
	 * @param pool
	 * @param resultsQueue
	 * @param r1Filename
	 * @param r2Filename
	 * @param i1Filename
	 * @param i2Filename
	 * @param i5Label
	 * @param i7Label
	 * @throws FileNotFoundException
	 * @throws IOException
	 * @throws InterruptedException
	 */
	protected void enqueuePairedReads(ExecutorService pool, BlockingQueue<Future<MergeResult>> resultsQueue, 
			String r1Filename, String r2Filename, String i1Filename, String i2Filename, String i5Label, String i7Label) throws FileNotFoundException, IOException, InterruptedException {
		Map<IndexAndBarcodeKey, Integer> barcodeLengthByIndexPairCache = new HashMap<IndexAndBarcodeKey, Integer>();
		
		if(i5Label != null && i7Label != null) {
			if(i5Indices.getBarcodeLength(i5Label) == 0) {
				throw new RuntimeException("Bad index label: " + i5Label);
			}
			if(i7Indices.getBarcodeLength(i7Label) == 0) {
				throw new RuntimeException("Bad index label: " + i7Label);
			}
			try(
					FileInputStream r1File = new FileInputStream(r1Filename);
					FileInputStream r2File = new FileInputStream(r2Filename);

					FastqReader r1Reader = new FastqReader(new BufferedReader(new InputStreamReader(new GZIPInputStream(r1File))));
					FastqReader r2Reader = new FastqReader(new BufferedReader(new InputStreamReader(new GZIPInputStream(r2File))));
					){
				while(r1Reader.hasNext() && r2Reader.hasNext() ){
					Read r1 = new Read(r1Reader.next());
					Read r2 = new Read(r2Reader.next());

					// Lookup by index pair whether barcodes are used
					IndexAndBarcodeKey keyIndexOnly = MergedRead.findExperimentKey(r1, r2, i5Label, i7Label, null, -1);
					int barcodeLength = findBarcodeLength(barcodeCountStatistics, keyIndexOnly, barcodeLengthByIndexPairCache, 
							barcodeLengthsFromSampleSheet, barcodes, barcodeToNoBarcodeThreshold);

					// update key if barcodes are used, otherwise reuse the index pair
					IndexAndBarcodeKey key = (barcodeLength > 0) ? MergedRead.findExperimentKey(r1, r2, i5Label, i7Label, barcodes, -1) : keyIndexOnly;

					Future<MergeResult> x = pool.submit(new Callable<MergeResult>() {
						public MergeResult call() {
							MergeResult mergeResult = merge(key, r1, r2);
							return mergeResult;
						}
					});
					resultsQueue.put(x);
				}
			}
		} else {
			try( 
					FileInputStream r1File = new FileInputStream(r1Filename);
					FileInputStream r2File = new FileInputStream(r2Filename);
					FileInputStream i1File = new FileInputStream(i1Filename);
					FileInputStream i2File = new FileInputStream(i2Filename);

					FastqReader r1Reader = new FastqReader(new BufferedReader(new InputStreamReader(new GZIPInputStream(r1File))));
					FastqReader r2Reader = new FastqReader(new BufferedReader(new InputStreamReader(new GZIPInputStream(r2File))));
					FastqReader i1Reader = new FastqReader(new BufferedReader(new InputStreamReader(new GZIPInputStream(i1File))));
					FastqReader i2Reader = new FastqReader(new BufferedReader(new InputStreamReader(new GZIPInputStream(i2File))));
					){

				while(r1Reader.hasNext() && r2Reader.hasNext() && i1Reader.hasNext() && i2Reader.hasNext()){
					Read r1 = new Read(r1Reader.next());
					Read r2 = new Read(r2Reader.next());
					Read i1 = new Read(i1Reader.next());
					Read i2 = new Read(i2Reader.next());

					if(reverseComplementI5) {
						i2 = i2.reverseComplement();
					}

					// Lookup by index pair whether barcodes are used
					IndexAndBarcodeKey keyIndexOnly = MergedRead.findExperimentKey(r1, r2, i1, i2, 
							i5Indices, i7Indices, null, 0);
					int barcodeLength = findBarcodeLength(barcodeCountStatistics, keyIndexOnly, barcodeLengthByIndexPairCache, 
							barcodeLengthsFromSampleSheet, barcodes, barcodeToNoBarcodeThreshold);

					// update key if barcodes are used, otherwise reuse the index pair
					IndexAndBarcodeKey key = (barcodeLength > 0) ? MergedRead.findExperimentKey(r1, r2, i1, i2, 
							i5Indices, i7Indices, barcodes, barcodeLength) : keyIndexOnly;

					Future<MergeResult> x = pool.submit(new Callable<MergeResult>() {
						public MergeResult call() {
							MergeResult mergeResult = merge(key, r1, r2);
							return mergeResult;
						}
					});
					resultsQueue.put(x);
				}
			}
		}
	}
	
	/**
	 * 
	 * @param key
	 * @param r1
	 * @param r2
	 * @return
	 */
	public MergeResult merge(IndexAndBarcodeKey key, Read r1, Read r2) {
		IndexAndBarcodeKey keyFlattened = null;
		MergedRead merged = null;
		boolean isOligo = false;
		String readGroup = null;
		if(key != null){
			keyFlattened = key.flatten();
			int r1BarcodeLength = barcodes.getBarcodeLength(keyFlattened.getP5Label());
			int r2BarcodeLength = barcodes.getBarcodeLength(keyFlattened.getP7Label());
			merged = MergedRead.mergePairedSequences(r1, r2, key, 
					r1BarcodeLength, r2BarcodeLength, maxPenalty, minOverlap, minMergedLength,
					mismatchPenaltyHigh, mismatchPenaltyLow, mismatchBaseQualityThreshold);
			// read group consistency
			readGroup = r1.getFASTQHeader().getReadGroupElements();
			// count occurrences of positive oligo
			if(positiveOligo != null && merged != null) {
				if(Read.alignmentAssessment(positiveOligo, merged, 0, 0, positiveOligo.length(), maxPenalty, mismatchPenaltyHigh, mismatchPenaltyLow, mismatchBaseQualityThreshold)
					|| Read.alignmentAssessment(positiveOligoReverseComplement, merged, 0, 0, positiveOligo.length(), maxPenalty, mismatchPenaltyHigh, mismatchPenaltyLow, mismatchBaseQualityThreshold)) {
					isOligo = true;
				}
			}
		}
		MergeResult result = new MergeResult(keyFlattened, merged, isOligo, readGroup);
		return result;
	}
	
	/**
	 * We assume that for a given index pair, barcodes are all the same length
	 * We first use the index-barcode key file, if given, then
	 * We use the 4-tuple with maximum count to determine the barcode length for this index pair
	 * @param barcodeCountStatistics
	 * @param keyIndexOnly
	 * @param barcodeLengthByIndexPairCache
	 * @param barcodeLengthsFromSampleSheet
	 * @param barcodes
	 * @param barcodeToNoBarcodeThreshold
	 * @return
	 */
	public static int findBarcodeLength(SampleSetsCounter barcodeCountStatistics, IndexAndBarcodeKey keyIndexOnly, 
			Map<IndexAndBarcodeKey, Integer> barcodeLengthByIndexPairCache, Map<IndexAndBarcodeKey, Integer> barcodeLengthsFromSampleSheet, 
			BarcodeMatcher barcodes, float barcodeToNoBarcodeThreshold) {
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
		return barcodeLength;
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

		SampleCounter countsForKeyIndexOnly = barcodeCountStatistics.get(keyIndexOnly.toString());
		List<String> barcodePairStrings = countsForKeyIndexOnly.getLabelList();
		// sort barcode pairs by count, tail is highest count
		SortedMap<Long, String> barcodePairsByCount = new TreeMap<Long, String>();
		for(String barcodePairString : barcodePairStrings){
			barcodePairsByCount.put(countsForKeyIndexOnly.get(barcodePairString), barcodePairString);
		}
		
		int barcodeLength = -1;
		Long max;
		max = barcodePairsByCount.lastKey();
		String maxPairLabel = barcodePairsByCount.remove(max);
		barcodeLength = barcodes.getBarcodePairLength(maxPairLabel);
		
		if(barcodeLength == 0 && barcodePairsByCount.size() > 0){
			// if a minority of barcodes pass the barcode check, we can falsely assume there are no barcodes
			// check whether a substantial fraction of barcodes are for a single barcode pair
			Long second = barcodePairsByCount.lastKey();
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

	public int getMaxPenalty() {
		return maxPenalty;
	}

	public void setMaxPenalty(int maxPenalty) {
		this.maxPenalty = maxPenalty;
	}

	public int getMismatchPenaltyHigh() {
		return mismatchPenaltyHigh;
	}

	public void setMismatchPenaltyHigh(int mismatchPenaltyHigh) {
		this.mismatchPenaltyHigh = mismatchPenaltyHigh;
	}

	public int getMismatchPenaltyLow() {
		return mismatchPenaltyLow;
	}

	public void setMismatchPenaltyLow(int mismatchPenaltyLow) {
		this.mismatchPenaltyLow = mismatchPenaltyLow;
	}

	public int getMismatchBaseQualityThreshold() {
		return mismatchBaseQualityThreshold;
	}

	public void setMismatchBaseQualityThreshold(int mismatchBaseQualityThreshold) {
		this.mismatchBaseQualityThreshold = mismatchBaseQualityThreshold;
	}

	public int getMinOverlap() {
		return minOverlap;
	}

	public void setMinOverlap(int minOverlap) {
		this.minOverlap = minOverlap;
	}

	public int getMinMergedLength() {
		return minMergedLength;
	}

	public void setMinMergedLength(int minMergedLength) {
		this.minMergedLength = minMergedLength;
	}

	public int getNumOutputFiles() {
		return numOutputFiles;
	}

	public void setNumOutputFiles(int numOutputFiles) {
		this.numOutputFiles = numOutputFiles;
	}

	public Read getPositiveOligo() {
		return positiveOligo;
	}

	public void setPositiveOligo(String positiveOligoSequence) {
		String qualityString = String.join("", Collections.nCopies(positiveOligoSequence.length(), "I")); // oligo sequence is max quality
		positiveOligo = new Read("", positiveOligoSequence, qualityString);
		positiveOligoReverseComplement = positiveOligo.reverseComplement();
	}
	
	public float getBarcodeToNoBarcodeThreshold() {
		return barcodeToNoBarcodeThreshold;
	}

	public void setBarcodeToNoBarcodeThreshold(float barcodeToNoBarcodeThreshold) {
		this.barcodeToNoBarcodeThreshold = barcodeToNoBarcodeThreshold;
	}

	public String getReadGroupFilename() {
		return readGroupFilename;
	}

	public void setReadGroupFilename(String readGroupFilename) {
		this.readGroupFilename = readGroupFilename;
	}

	public boolean isReverseComplementI5() {
		return reverseComplementI5;
	}

	public void setReverseComplementI5(boolean reverseComplementI5) {
		this.reverseComplementI5 = reverseComplementI5;
	}

	public Map<IndexAndBarcodeKey, Integer> getBarcodeLengthsFromSampleSheet() {
		return barcodeLengthsFromSampleSheet;
	}

	public void setBarcodeLengthsFromSampleSheet(Map<IndexAndBarcodeKey, Integer> barcodeLengthsFromSampleSheet) {
		this.barcodeLengthsFromSampleSheet = barcodeLengthsFromSampleSheet;
	}

	public SampleSetsCounter getBarcodeCountStatistics() {
		return barcodeCountStatistics;
	}

	public void setBarcodeCountStatistics(SampleSetsCounter barcodeCountStatistics) {
		this.barcodeCountStatistics = barcodeCountStatistics;
	}

	public PrintStream getPrintStream() {
		return printStream;
	}

	public void setPrintStream(OutputStream stream) {
		this.printStream = new PrintStream(stream);
	}
}
