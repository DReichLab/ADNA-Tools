package adnascreen;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.Map;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import org.apache.commons.lang3.mutable.MutableInt;

import htsjdk.samtools.fastq.FastqReader;

public class IndexAndBarcodeScreener {

	public static void main(String []args) throws IOException{
		IndexMatcher i5Indices = null, i7Indices = null;
		BarcodeMatcher barcodes = null;
		Map<IndexAndBarcodeKey, MutableInt> SampleCounter = new HashMap<IndexAndBarcodeKey, MutableInt>();
		final int maxPenalty = 3;
		final int minOverlap = 10;
		final int minMergedLength = 30;
		final int numOutputFiles = 25;

		try{
			i5Indices = new IndexMatcher(args[0], 1);
			i7Indices = new IndexMatcher(args[1], 1);
			barcodes = new BarcodeMatcher(args[2], 1);
		} catch(IOException e){
			System.exit(1);
		}

		PrintWriter [] fileOutputs = new PrintWriter[numOutputFiles];;
		try(
				FileInputStream r1File = new FileInputStream(args[3]);
				FileInputStream r2File = new FileInputStream(args[4]);
				FileInputStream i1File = new FileInputStream(args[5]);
				FileInputStream i2File = new FileInputStream(args[6]);

				FastqReader r1Reader = new FastqReader(new BufferedReader(new InputStreamReader(new GZIPInputStream(r1File))));
				FastqReader r2Reader = new FastqReader(new BufferedReader(new InputStreamReader(new GZIPInputStream(r2File))));
				FastqReader i1Reader = new FastqReader(new BufferedReader(new InputStreamReader(new GZIPInputStream(i1File))));
				FastqReader i2Reader = new FastqReader(new BufferedReader(new InputStreamReader(new GZIPInputStream(i2File))));
				){
			// prepare output files for multiple parallel processing jobs downstream
			// for load balancing purposes, these are not demultiplexed
			String outputFilenameRoot = args[7];
			for(int i = 0; i < numOutputFiles; i++){
				// start counting from 1 for filenames
				String outputFilename = String.format("%s_%03d.fastq.gz", outputFilenameRoot, i + 1);
				fileOutputs[i] = new PrintWriter(new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outputFilename)))));
			}
			
			int pairedReadInputCount = 0;
			int pairedReadOutputCount = 0;
			while(r1Reader.hasNext() && r2Reader.hasNext() && i1Reader.hasNext() && i2Reader.hasNext()){
				pairedReadInputCount++;
				
				Read r1 = new Read(r1Reader.next());
				Read r2 = new Read(r2Reader.next());
				Read i1 = new Read(i1Reader.next());
				Read i2 = new Read(i2Reader.next());
				
				PairedRead merged = PairedRead.mergePairedSequences(r1, r2, i1, i2, 
						i5Indices, i7Indices, barcodes, maxPenalty, minOverlap, minMergedLength);
				if(merged != null){
					// TODO histogram of length distribution
					// separate into different files
					fileOutputs[pairedReadOutputCount % numOutputFiles].println(merged.toString());
					pairedReadOutputCount++;
				}
			}
			// output map statistics
			System.err.println("Number of paired reads: " + pairedReadInputCount);
			for(IndexAndBarcodeKey key : SampleCounter.keySet()){
				System.err.print(key);
				System.err.print('\t');
				System.err.println(SampleCounter.get(key).intValue());
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
