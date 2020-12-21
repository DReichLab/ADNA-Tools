package adnascreen;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Queue;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMFormatException;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class DemultiplexSAM {
	public static final String ALIGNED = "aligned";
	public static final String DEMULTIPLEXED = "demultiplexed";
	public static final String duplicatesSAMTag = "XD";
	
	/**
	 * Add the numTopSamples from a statistics file to the queue of files to be opened for output. 
	 * Each of these must have at least minimumReads. 
	 * @param statisticsFilename
	 * @param numTopSamples
	 * @param minimumReads
	 * @param thresholdReads
	 * @throws FileNotFoundException
	 * @throws IOException
	 */
	public static Queue<IndexAndBarcodeKey> selectTopSamples(String statisticsFilename, int numTopSamples, int minimumReads, int thresholdReads) throws FileNotFoundException, IOException {
		Queue<IndexAndBarcodeKey> topOutputFiles = new LinkedList<IndexAndBarcodeKey>();
		File statisticsFile = new File(statisticsFilename);
		try(BufferedReader reader = new BufferedReader(new FileReader(statisticsFile))){
			Long.valueOf(reader.readLine());
			int n = 0;
			String entryLine = null;
			while((entryLine = reader.readLine()) != null) {
				String [] fields = entryLine.split("\t", 2);					
				String keyString = fields[0]; 
				String keyValuePairs = fields[1];
				
				SampleCounter counts = new SampleCounter(keyValuePairs);
				long rawCount = counts.get(IndexAndBarcodeScreener.RAW);
				if ( (n++ < numTopSamples || ((thresholdReads >= 0) && (rawCount >= thresholdReads) ) )
						&& (rawCount >= minimumReads)){
					IndexAndBarcodeKey key = new IndexAndBarcodeKey(keyString);
					topOutputFiles.add(key);
				}
			}
		}
		return topOutputFiles;
	}
	
	public static void main(String [] args) throws IOException, ParseException {
		CommandLineParser parser = new DefaultParser();
		Options options = new Options();
		options.addRequiredOption("s", "statisticsFilename", true, "Statistics file sorted in order of output");
		options.addOption("f", "barcodeFile", true, "Barcode file to mark for duplicates");
		options.addOption("n", "numSamples", true, "Number of top samples to output");
		options.addOption("m", "maximumConcurrentOpenFiles", true, "Maximum number of samples to demultiplex concurrently [restricted by OS]");
		options.addOption("r", "minimumReads", true, "Minimum number of reads to process");
		options.addOption("b", "BAM", false, "Use bam files for output");
		options.addOption("e", "explicit", true, "Explicit indices to demultiplex");
		options.addOption(null, "bufferSize", true, "Output file buffer size for performance");
		options.addOption("c", "compression", true, "BAM compression level for htsjdk 0-9 default 5");
		options.addOption("null", "async", false, "Use asynchronous threading for output");
		options.addOption("t", "thresholdReads", true, "Threshold number of reads to trigger processing.");
		options.addOption("o", "outputDirectory", true, "Directory to use for output files");
		
		options.addOption(null, "stdoutFile", true, "Use this file for stdout");
		options.addOption(null, "stderrFile", true, "Use this file for stderr");
		
		CommandLine commandLine	= parser.parse(options, args);
		
		int numTopSamples = Integer.valueOf(commandLine.getOptionValue('n', "1000"));
		int maximumConcurrentOpenFiles = Integer.valueOf(commandLine.getOptionValue('m', "1000"));
		int minimumReads = Integer.valueOf(commandLine.getOptionValue('r', "1"));
		int bufferSize = Integer.valueOf(commandLine.getOptionValue("bufferSize", "1048576")); // 1 MB
		int compressionLevel = Integer.valueOf(commandLine.getOptionValue("compression", "5"));
		boolean useBAM = commandLine.hasOption('b');
		boolean useAsyncThreads = commandLine.hasOption("async");
		String fileExtension = useBAM ? ".bam" : ".sam";
		String explicitIndexFile = commandLine.getOptionValue("explicit", null);
		String barcodeFilename = commandLine.getOptionValue("barcodeFile", null);
		int thresholdReads = Integer.valueOf(commandLine.getOptionValue("thresholdReads", "-1"));
		String outputDirectory = commandLine.getOptionValue("outputDirectory", ".");
		
		String stdoutFilename = commandLine.getOptionValue("stdoutFile", null);
		PrintStream stdout = (stdoutFilename == null) ? System.out : new PrintStream(stdoutFilename);
		String stderrFilename = commandLine.getOptionValue("stderrFile", null);
		PrintStream stderr = (stderrFilename == null) ? System.err : new PrintStream(stderrFilename);
		
		Queue<IndexAndBarcodeKey> outputFilesAll = new LinkedList<IndexAndBarcodeKey>();
		
		SAMSequenceDictionary alignmentReference = null;
		SAMFileWriterFactory outputFileFactory = new SAMFileWriterFactory();
		outputFileFactory.setBufferSize(bufferSize);
		outputFileFactory.setCompressionLevel(compressionLevel);
		outputFileFactory.setUseAsyncIo(useAsyncThreads);
		
		// allow explicit additions to list of samples to demultiplex
		// these will always be demultiplexed, independent of the top number of samples or number of raw reads
		if(explicitIndexFile != null){
			try(BufferedReader reader = new BufferedReader(new FileReader(explicitIndexFile))){
				String entryLine;
				while((entryLine = reader.readLine()) != null){
					try {
						String [] fields = entryLine.split("\t");
						String keyString = fields[0];
						IndexAndBarcodeKey key = new IndexAndBarcodeKey(keyString);
						outputFilesAll.add(key);
					}
					catch(Exception e) {
						stderr.println(e.toString());
					}
				}
			}
		}
		
		// read statistics file with top keys
		// open a SAM/BAM file for each key for demultiplexing its data
		String statisticsFilename = commandLine.getOptionValue('s');
		Queue<IndexAndBarcodeKey> topOutputFiles = selectTopSamples(statisticsFilename, numTopSamples, minimumReads, thresholdReads);
		outputFilesAll.addAll(topOutputFiles);
		
		stderr.println("Outputting " + outputFilesAll.size() + " files");
		
		// we write barcode sequences into reads to mark duplicates
		BarcodeMatcher barcodes = new BarcodeMatcher();
		if(barcodeFilename != null)
			barcodes.loadFile(barcodeFilename);
		
		File statisticsFile = new File(statisticsFilename);
		SampleSetsCounter statistics = new SampleSetsCounter(statisticsFile);
		
		while(outputFilesAll.size() > 0) {
			// We may need multiple passes through the input files due to concurrent open file limit
			Map<IndexAndBarcodeKey, SAMFileWriter> outputFilesConcurrent = new HashMap<IndexAndBarcodeKey, SAMFileWriter>(numTopSamples);
			// prepare as many output files as concurrently possible
			while(outputFilesConcurrent.size() < maximumConcurrentOpenFiles && outputFilesAll.size() > 0) {
				IndexAndBarcodeKey key = outputFilesAll.remove();
				outputFilesConcurrent.put(key, null); // mark this key for output in this pass
				// we delay opening SAM/BAM file writer until the SAM/BAM header is available
				// this is after we have opened the first SAM/BAM input file
			}		

			// iterate through input files
			List<String> samFilenamesToProcess = commandLine.getArgList();
			for(String filename : samFilenamesToProcess){
				SamInputResource bufferedSAMFile = SamInputResource.of(new BufferedInputStream(new FileInputStream(filename), bufferSize));
				try(
						SamReader reader = SamReaderFactory.makeDefault().open(bufferedSAMFile);
						){
					SAMFileHeader header = reader.getFileHeader();
					SAMSequenceDictionary currentAlignmentReference = header.getSequenceDictionary();
					if(alignmentReference == null){
						alignmentReference = currentAlignmentReference;
					} else if(!alignmentReference.equals(currentAlignmentReference)){
						throw new IllegalArgumentException("SAM references do not match");
					}

					SAMRecordIterator i = reader.iterator();
					while(i.hasNext()){
						// iterate through alignments
						try{
							SAMRecord record = i.next();
							// parse out key, which is a 4-tuple of indices and barcodes
							String readName = record.getReadName();
							String [] readNameParts = readName.split(String.valueOf(MergedRead.KEY_SEPARATOR));
							IndexAndBarcodeKey key = new IndexAndBarcodeKey(readNameParts[1]);

							// for deduplication, write barcodes and read length to tag
							int length = record.getReadLength();
							DNASequence p5Barcode = barcodes.getBarcode(key.getP5Label());
							DNASequence p7Barcode = barcodes.getBarcode(key.getP7Label());
							String duplicate_marker = 
									(p5Barcode != null ? p5Barcode.toString() : "") + IndexAndBarcodeKey.FIELD_SEPARATOR +
									(p7Barcode != null ? p7Barcode.toString() : "") + IndexAndBarcodeKey.FIELD_SEPARATOR +
									length;
							record.setAttribute(duplicatesSAMTag, duplicate_marker);
							// remove the key from the read name
							String readNameNoKey = readNameParts[0];
							record.setReadName(readNameNoKey);

							IndexAndBarcodeKey keyFlattened = key.flatten();
							// record statistics
							// count of demultiplexed reads is for checking consistency
							statistics.increment(keyFlattened.toString(), DEMULTIPLEXED);
							if(!record.getReadUnmappedFlag()){ // read is mapped
								statistics.increment(keyFlattened.toString(), ALIGNED);
							}

							// write only to open files for top keys
							if(outputFilesConcurrent.containsKey(keyFlattened)){
								// find file corresponding to this key
								SAMFileWriter output = outputFilesConcurrent.get(keyFlattened);
								if(output == null){ // open new file, if none exists for this key
									String outputFilename = (keyFlattened.toString() + fileExtension).replace(':', '-'); // Cromwell chokes on files with ':'
									String outputPath = outputDirectory + "/" + outputFilename;
									BufferedOutputStream outputFile = new BufferedOutputStream(new FileOutputStream(outputPath), bufferSize);
									if(useBAM){
										output = outputFileFactory.makeBAMWriter(header, false, outputFile);
									} else {
										output = outputFileFactory.makeSAMWriter(header, false, outputFile);
									}
									outputFilesConcurrent.put(keyFlattened, output); // 
								}
								// write alignment to file
								output.addAlignment(record);
							}
						} catch (SAMFormatException e){
							stderr.print(filename + "\t");
							stderr.println(e);
							// ignore this record and continue to the next
						} catch (Exception e){
							stderr.print(filename + "\t");
							stderr.println(e);
							e.printStackTrace(stderr);
						}
					}
				}
			}
			// cleanup, close all output files
			for(SAMFileWriter writer : outputFilesConcurrent.values()){
				if(writer != null){
					writer.close();
				}
			}
			outputFilesConcurrent.clear();
		}
		stdout.println(statistics.toStringSorted(IndexAndBarcodeScreener.RAW));
		stdout.close();
		stderr.close();
	}
}
