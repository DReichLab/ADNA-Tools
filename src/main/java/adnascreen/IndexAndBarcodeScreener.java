package adnascreen;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.io.PrintStream;
import java.io.PrintWriter;
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
		options.addOption("p", "mismatch-penalty", true, "Penalty for mismatch while aligning merge");
		options.addOption("o", "minimum-overlap", true, "Minimum bases of overlap for paired reads to merge");
		options.addOption("l", "minimum-length", true, "Minimum length for output merged reads");
		options.addOption("n", "number-output-files", true, "Number of output files to divide merged reads between");
		options.addOption("h", "hamming-distance", true, "Max hamming distance for index or barcode match");
		options.addOption("r", "read-group-file", true, "Output file for read group");
		CommandLine commandLine	= parser.parse(options, args);
		
		BarcodeMatcher i5Indices = null, i7Indices = null;
		BarcodeMatcher barcodes = null;
		// We keep statistics for each 4-tuple of indices and barcodes
		SampleSetsCounter sampleSetCounter = new SampleSetsCounter();
		final int maxPenalty = Integer.valueOf(commandLine.getOptionValue('p', "3"));
		final int minOverlap = Integer.valueOf(commandLine.getOptionValue('o', "15"));
		final int minMergedLength = Integer.valueOf(commandLine.getOptionValue('l', "30"));
		final int numOutputFiles = Integer.valueOf(commandLine.getOptionValue('n', "25"));
		final int maxHammingDistance = Integer.valueOf(commandLine.getOptionValue('h', "1"));
		String readGroupFilename = commandLine.getOptionValue("read-group-file", "read_group");
		
		try{
			i5Indices = new BarcodeMatcher(commandLine.getOptionValue("i5-indices"), maxHammingDistance);
			i7Indices = new BarcodeMatcher(commandLine.getOptionValue("i7-indices"), maxHammingDistance);
			barcodes = new BarcodeMatcher(commandLine.getOptionValue('b'), maxHammingDistance);
		} catch(IOException e){
			System.exit(1);
		}

		String[] remainingArgs = commandLine.getArgs();
		PrintWriter [] fileOutputs = new PrintWriter[numOutputFiles];;
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
				
				IndexAndBarcodeKey key = MergedRead.findExperimentKey(r1, r2, i1, i2, 
						i5Indices, i7Indices, barcodes);
				IndexAndBarcodeKey keyFlattened = null;
				sampleSetCounter.increment(); // statistics recording
				MergedRead merged = null;
				if(key != null){
					keyFlattened = key.flatten();
					sampleSetCounter.increment(keyFlattened, RAW);
					merged = MergedRead.mergePairedSequences(r1, r2, key, 
						barcodes.getBarcodeLength(), maxPenalty, minOverlap, minMergedLength);
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
				if(merged != null){
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
}
