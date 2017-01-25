package adnascreen;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.util.HashMap;
import java.util.Map;
import java.util.zip.GZIPInputStream;

import org.apache.commons.lang3.mutable.MutableInt;

public class IndexAndBarcodeScreener {

	public static void main(String []args) throws IOException{
		IndexMatcher i5Indices = null, i7Indices = null;
		BarcodeMatcher barcodes = null;
		Map<IndexAndBarcodeKey, MutableInt> SampleCounter = new HashMap<IndexAndBarcodeKey, MutableInt>();

		try{
			i5Indices = new IndexMatcher(args[0], 1);
			i7Indices = new IndexMatcher(args[1], 1);
			barcodes = new BarcodeMatcher(args[2], 1);
		} catch(IOException e){
			System.exit(1);
		}

		try(			
				FileInputStream r1File = new FileInputStream(args[3]);
				FileInputStream r2File = new FileInputStream(args[4]);
				FileInputStream i1File = new FileInputStream(args[5]);
				FileInputStream i2File = new FileInputStream(args[6]);

				BufferedReader r1 = new BufferedReader(new InputStreamReader(new GZIPInputStream(r1File)));
				BufferedReader r2 = new BufferedReader(new InputStreamReader(new GZIPInputStream(r2File)));
				BufferedReader i1 = new BufferedReader(new InputStreamReader(new GZIPInputStream(i1File)));
				BufferedReader i2 = new BufferedReader(new InputStreamReader(new GZIPInputStream(i2File)));
				){
			
			int pairedReadCount = 0;
			String endTest;
			while(((endTest = r1.readLine()) != null) && endTest.length() > 0){
				pairedReadCount++;
				// FASTQ line 1 metadata
				FASTQHeader sequenceMetadataR1 = new FASTQHeader(endTest);
				FASTQHeader sequenceMetadataR2 = new FASTQHeader(r2.readLine());
				FASTQHeader sequenceMetadataI1 = new FASTQHeader(i1.readLine());
				FASTQHeader sequenceMetadataI2 = new FASTQHeader(i2.readLine());

				// check for metadata consistency
				if(!sequenceMetadataR1.equalsExceptRead(sequenceMetadataR2)
						|| !sequenceMetadataR1.equalsExceptRead(sequenceMetadataI1)
						|| !sequenceMetadataR1.equalsExceptRead(sequenceMetadataI2))
					throw new IllegalArgumentException("FASTQ metadata mismatch");

				// FASTQ line 2 base pair sequence
				// Index 1 is i7, Index 2 is i5
				String basePairR1 = r1.readLine().trim();
				DNASequence i7IndexRaw = new DNASequence(i1.readLine().trim());
				DNASequence i5IndexRaw = new DNASequence(i2.readLine().trim());
				String basePairR2 = r2.readLine().trim();

				DNASequence p5BarcodeRaw = new DNASequence(basePairR1.substring(0, 7));
				DNASequence p7BarcodeRaw = new DNASequence(basePairR2.substring(0, 7));

				// FASTQ line 3 quality metadata
				String qualityMetadataR1 = r1.readLine();
				String qualityMetadataR2 = r2.readLine();
				String qualityMetadataI1 = i1.readLine();
				String qualityMetadataI2 = i2.readLine();

				if(!qualityMetadataR1.equals(qualityMetadataR2)
						|| !qualityMetadataR1.equals(qualityMetadataI1)
						|| !qualityMetadataR1.equals(qualityMetadataI2))
					throw new IllegalArgumentException("FASTQ quality metadata mismatch");

				// FASTQ line 4
				String qualitySequenceR1 = r1.readLine();
				String qualitySequenceR2 = r2.readLine();
				String qualitySequenceI1 = i1.readLine();
				String qualitySequenceI2 = i2.readLine();

				// 
				DNASequence i5Index = i5Indices.find(i5IndexRaw);
				DNASequence i7Index = i7Indices.find(i7IndexRaw);
				String p5BarcodeSet = barcodes.find(p5BarcodeRaw);
				String p7BarcodeSet = barcodes.find(p7BarcodeRaw);

				if(i5Index != null 
						&& i7Index != null
						&& p5BarcodeSet != null
						&& p7BarcodeSet != null){
					// update counter
					IndexAndBarcodeKey key = new IndexAndBarcodeKey(i5Index, i7Index, p5BarcodeSet, p7BarcodeSet);
					MutableInt count = SampleCounter.get(key);
					if(count == null){
						count = new MutableInt(0);
						SampleCounter.put(key, count);
					}
					count.increment();

					// output interleaved FASTQ entries with augmented headers
					PrintStream out = System.out;
					//String keyString = key.toString();
					out.print(sequenceMetadataR1);
					//out.print(" ");
					//out.print(keyString);
					out.println();
					out.println(basePairR1);
					out.println(qualityMetadataR1);
					out.println(qualitySequenceR1);

					out.print(sequenceMetadataR2);
					//out.print(" ");
					//out.print(keyString);
					out.println();
					out.println(basePairR2);
					out.println(qualityMetadataR2);
					out.println(qualitySequenceR2);
				}
				else{
					//System.err.println(i5Index + " " + i7Index + " " + p5BarcodeSet + " " + p7BarcodeSet);
				}
			}
			// output map statistics
			System.err.println("Number of paired reads: " + pairedReadCount);
			for(IndexAndBarcodeKey key : SampleCounter.keySet()){
				System.err.print(key);
				System.err.print('\t');
				System.err.println(SampleCounter.get(key).intValue());
			}
		}
	}
}
