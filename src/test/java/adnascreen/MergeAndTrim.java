package adnascreen;

import static org.junit.Assert.*;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.concurrent.ExecutionException;
import java.util.zip.GZIPInputStream;

import org.apache.commons.cli.ParseException;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TemporaryFolder;

import htsjdk.samtools.fastq.FastqReader;

public class MergeAndTrim {
	static final int HAMMING_DISTANCE = 1;
	
	@Rule
    public TemporaryFolder tempFolder= new TemporaryFolder();
	
	private void compareFastq(String f1Filename, String f2Filename) {
		FileInputStream f1File;
		FileInputStream f2File;
		FastqReader f1Reader = null;
		FastqReader f2Reader = null;
		try {
			f1File = new FileInputStream(f1Filename);
			f2File = new FileInputStream(f2Filename);
			
			f1Reader = new FastqReader(new BufferedReader(new InputStreamReader(new GZIPInputStream(f1File))));
			f2Reader = new FastqReader(new BufferedReader(new InputStreamReader(new GZIPInputStream(f2File))));
			
			while(f1Reader.hasNext() || f2Reader.hasNext() ){
				Read r1 = new Read(f1Reader.next());
				Read r2 = new Read(f2Reader.next());
				
				assertEquals(r1, r2);
			}
			
		} catch (IOException e) {
			fail();
		}
		finally {
			if (f1Reader != null)
				f1Reader.close();
			if (f2Reader != null)
				f2Reader.close();
		}
	}
	
	private void compareTextFiles(String f1Filename, String f2Filename) {
		BufferedReader f1;
		BufferedReader f2;
		try {
			f1 = new BufferedReader(new FileReader(f1Filename));
			f2 = new BufferedReader(new FileReader(f2Filename));
			
			String expected;
			while((expected = f1.readLine()) != null) {
				String actual = f2.readLine();
				assertEquals(expected, actual);
			}
			assertNull(f2.readLine());
		} catch (IOException e) {
			fail();
		}
	}
	
	// test threading
	@Test
	public void fromFastq() throws IOException, ParseException, InterruptedException, ExecutionException {
		ClassLoader classLoader = getClass().getClassLoader();
		String r1Filename = classLoader.getResource("fastq/r1.fastq.gz").getPath();
		String r2Filename = classLoader.getResource("fastq/r2.fastq.gz").getPath();
		String i1Filename = classLoader.getResource("fastq/i1.fastq.gz").getPath();
		String i2Filename = classLoader.getResource("fastq/i2.fastq.gz").getPath();
		
		String i5Filename = classLoader.getResource("fastq/i5").getPath();
		String i7Filename = classLoader.getResource("fastq/i7").getPath();
		String barcodeFilename = classLoader.getResource("fastq/barcodes").getPath();
		
		String barcodeCountsFilename = classLoader.getResource("fastq/barcodeCounts").getPath();
		String expectedCountsFilename = classLoader.getResource("fastq/expected").getPath();
		
		try {
			BarcodeMatcher i5 = new BarcodeMatcher(i5Filename, HAMMING_DISTANCE);
			BarcodeMatcher i7 = new BarcodeMatcher(i7Filename, HAMMING_DISTANCE);
			BarcodeMatcher barcodes = new BarcodeMatcher(barcodeFilename, HAMMING_DISTANCE);
			
			IndexAndBarcodeScreener screener = new IndexAndBarcodeScreener(i5, i7, barcodes);
			screener.setNumOutputFiles(2);
			File barcodeCountStatisticsFile = new File(barcodeCountsFilename);
			screener.setBarcodeCountStatistics(new SampleSetsCounter(barcodeCountStatisticsFile));
			screener.setPrintStream(new FileOutputStream(tempFolder.getRoot() + "/counts"));
			
			String outputFileBase = tempFolder.getRoot() + "/test";
			screener.performScreeningMergeTrim(4, outputFileBase, r1Filename, r2Filename, i1Filename, i2Filename, null, null);
		} catch (IOException | ParseException | InterruptedException | ExecutionException e) {
			fail();
		}
		
		String outputBaseFilename1 = "test_001.fastq.gz";
		String outputBaseFilename2 = "test_002.fastq.gz";
		String expectedOutputFilename1 = classLoader.getResource("fastq/" + outputBaseFilename1).getPath();
		String expectedOutputFilename2 = classLoader.getResource("fastq/" + outputBaseFilename2).getPath();
		String outputFilename1 = tempFolder.getRoot() + "/" + outputBaseFilename1;
		String outputFilename2 = tempFolder.getRoot() + "/" + outputBaseFilename2;
		// compare results with non-threaded version
		compareFastq(expectedOutputFilename1, outputFilename1);
		compareFastq(expectedOutputFilename2, outputFilename2);
		
		String countFilename = tempFolder.getRoot() + "/counts";
		compareTextFiles(expectedCountsFilename, countFilename);
	}
}
