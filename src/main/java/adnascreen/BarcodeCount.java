package adnascreen;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.util.zip.GZIPInputStream;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import htsjdk.samtools.fastq.FastqReader;

public class BarcodeCount {
	public static final String WITHOUT_BARCODES = "without_barcodes";
	public static final String WITH_BARCODES = "with_barcodes";

	public static void main(String []args) throws IOException, ParseException{
		CommandLineParser parser = new DefaultParser();
		Options options = new Options();
		options.addRequiredOption("i", "i5-indices", true, 
				"File containing one valid i5 index sets per line");
		options.addRequiredOption("j", "i7-indices", true, 
				"File containing one valid i7 index sets per line");
		options.addRequiredOption("b", "barcodes", true, 
				"File containing one valid barcodes sets per line with ':'-delimited elements");
		options.addOption("h", "hamming-distance", true, "Max hamming distance for index or barcode match");
		options.addOption("y", "reverse-complement-i5", false, "Whether i5 index should be reverse complemented (NextSeq is not)");
		
		options.addOption(null, "fixed-i5", true, "Assume all fragments have this i5 sequence label");
		options.addOption(null, "fixed-i7", true, "Assume all fragments have this i7 sequence label");
		CommandLine commandLine	= parser.parse(options, args);
		
		BarcodeMatcher i5Indices = null, i7Indices = null;
		BarcodeMatcher barcodes = null;
		// We keep statistics for each 4-tuple of indices and barcodes
		SampleSetsCounter sampleSetCounter = new SampleSetsCounter();
		final int maxHammingDistance = Integer.valueOf(commandLine.getOptionValue('h', "1"));
		final boolean reverseComplementI5 = commandLine.hasOption('y'); 
		
		try{
			i5Indices = new BarcodeMatcher(commandLine.getOptionValue("i5-indices"), maxHammingDistance);
			i7Indices = new BarcodeMatcher(commandLine.getOptionValue("i7-indices"), maxHammingDistance);
			barcodes = new BarcodeMatcher(commandLine.getOptionValue('b'), maxHammingDistance);
		} catch(IOException e){
			System.exit(1);
		}
		
		String[] remainingArgs = commandLine.getArgs();
		// If there are fixed indices, we use those and do not have index reads
		if(commandLine.hasOption("fixed-i5") && commandLine.hasOption("fixed-i7")) {
			String i5Label = commandLine.getOptionValue("fixed-i5");
			String i7Label = commandLine.getOptionValue("fixed-i7");
			if(i5Indices.getBarcodeLength(i5Label) == 0) {
				throw new RuntimeException("Bad index label: " + i5Label);
			}
			if(i7Indices.getBarcodeLength(i7Label) == 0) {
				throw new RuntimeException("Bad index label: " + i7Label);
			}

			try(
					FileInputStream r1File = new FileInputStream(remainingArgs[0]);
					FileInputStream r2File = new FileInputStream(remainingArgs[1]);

					FastqReader r1Reader = new FastqReader(new BufferedReader(new InputStreamReader(new GZIPInputStream(r1File))));
					FastqReader r2Reader = new FastqReader(new BufferedReader(new InputStreamReader(new GZIPInputStream(r2File))));
					){
				while(r1Reader.hasNext() && r2Reader.hasNext() ){
					Read r1 = new Read(r1Reader.next());
					Read r2 = new Read(r2Reader.next());

					IndexAndBarcodeKey key = MergedRead.findExperimentKey(r1, r2, i5Label, i7Label, barcodes, -1);
					updateCounters(key, sampleSetCounter);
				}
			} catch(IOException e){
				System.err.println(e);
				System.exit(1);
			}
		}
		else { // no fixed index labels, so use index reads 
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

				while(r1Reader.hasNext() && r2Reader.hasNext() && i1Reader.hasNext() && i2Reader.hasNext()){
					Read r1 = new Read(r1Reader.next());
					Read r2 = new Read(r2Reader.next());
					Read i1 = new Read(i1Reader.next());
					Read i2 = new Read(i2Reader.next());

					if(reverseComplementI5) {
						i2 = i2.reverseComplement();
					}

					IndexAndBarcodeKey key = MergedRead.findExperimentKey(r1, r2, i1, i2, 
							i5Indices, i7Indices, barcodes, -1);
					updateCounters(key, sampleSetCounter);
				}
			} catch(IOException e){
				System.err.println(e);
				System.exit(1);
			}
		}
		// output map statistics
		PrintStream statisticsOutput = System.out;
		statisticsOutput.println(sampleSetCounter.toString());
	}
	
	public static void updateCounters(IndexAndBarcodeKey key, SampleSetsCounter sampleSetCounter) {
		IndexAndBarcodeKey keyFlattened = null;
		sampleSetCounter.increment(); // statistics recording
		if(key != null){
			keyFlattened = key.flatten();
			// differentiate between keys with and without barcodes
			String i5 = keyFlattened.getI5Label();
			String i7 = keyFlattened.getI7Label();
			String p5 = keyFlattened.getP5Label();
			String p7 = keyFlattened.getP7Label();
			
			IndexAndBarcodeKey keyIndexOnly = new IndexAndBarcodeKey(i5, i7, null, null);
			if(p5 != null && p7 != null){
				sampleSetCounter.increment(keyIndexOnly.toString(), p5 + IndexAndBarcodeKey.FIELD_SEPARATOR + p7);
			} else{
				sampleSetCounter.increment(keyIndexOnly.toString(), WITHOUT_BARCODES);
			}
		}
	}
}
