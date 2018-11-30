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
		CommandLine commandLine	= parser.parse(options, args);
		
		BarcodeMatcher i5Indices = null, i7Indices = null;
		BarcodeMatcher barcodes = null;
		// We keep statistics for each 4-tuple of indices and barcodes
		SampleSetsCounter sampleSetCounter = new SampleSetsCounter();
		final int maxHammingDistance = Integer.valueOf(commandLine.getOptionValue('h', "1"));
		
		try{
			i5Indices = new BarcodeMatcher(commandLine.getOptionValue("i5-indices"), maxHammingDistance);
			i7Indices = new BarcodeMatcher(commandLine.getOptionValue("i7-indices"), maxHammingDistance);
			barcodes = new BarcodeMatcher(commandLine.getOptionValue('b'), maxHammingDistance);
		} catch(IOException e){
			System.exit(1);
		}

		String[] remainingArgs = commandLine.getArgs();
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
				
				IndexAndBarcodeKey key = MergedRead.findExperimentKey(r1, r2, i1, i2, 
						i5Indices, i7Indices, barcodes, -1);
				IndexAndBarcodeKey keyFlattened = null;
				sampleSetCounter.increment(); // statistics recording
				if(key != null){
					keyFlattened = key.flatten();
					// differentiate between keys with and without barcodes
					String i5 = keyFlattened.getI5Label();
					String i7 = keyFlattened.getI7Label();
					String p5 = keyFlattened.getP5Label();
					String p7 = keyFlattened.getP7Label();
					//int barcodeLength = barcodes.getBarcodeLength(p5);
					IndexAndBarcodeKey keyIndexOnly = new IndexAndBarcodeKey(i5, i7, null, null);
					if(p5 != null && p7 != null){
						sampleSetCounter.increment(keyIndexOnly.toString(), p5 + IndexAndBarcodeKey.FIELD_SEPARATOR + p7);
					} else{
						sampleSetCounter.increment(keyIndexOnly.toString(), WITHOUT_BARCODES);
					}
				}
			}
			// output map statistics
			PrintStream statisticsOutput = System.out;
			statisticsOutput.println(sampleSetCounter.toString());
		} catch(IOException e){
			System.err.println(e);
			System.exit(1);
		}
	}
}
