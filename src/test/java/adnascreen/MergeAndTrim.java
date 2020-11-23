package adnascreen;

import static org.junit.Assert.*;

import java.io.BufferedReader;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.util.LinkedList;
import java.util.List;
import java.util.Scanner;
import java.util.concurrent.ExecutionException;
import java.util.zip.GZIPInputStream;

import org.apache.commons.cli.ParseException;
import org.junit.After;
import org.junit.Before;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TemporaryFolder;

import htsjdk.samtools.fastq.FastqReader;

public class MergeAndTrim {
	static final int HAMMING_DISTANCE = 1;
	
	@Rule
    public TemporaryFolder tempFolder= new TemporaryFolder();
	
	private final ByteArrayOutputStream outContent = new ByteArrayOutputStream();
	private final ByteArrayOutputStream errContent = new ByteArrayOutputStream();
	private final PrintStream originalOut = System.out;
	private final PrintStream originalErr = System.err;

	@Before
	public void setUpStreams() {
	    System.setOut(new PrintStream(outContent));
	    System.setErr(new PrintStream(errContent));
	}

	@After
	public void restoreStreams() {
	    System.setOut(originalOut);
	    System.setErr(originalErr);
	}
	
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
	
	/*
	 * Test the code to process shotgun samples that come without index reads. 
	 */
	@Test
	public void no_indices() {
		ClassLoader classLoader = getClass().getClassLoader();
		String r1Filename = classLoader.getResource("no_index_fastq/1.fastq.gz").getPath();
		String r2Filename = classLoader.getResource("no_index_fastq/2.fastq.gz").getPath();
		
		String i5Filename = classLoader.getResource("no_index_fastq/i5").getPath();
		String i7Filename = classLoader.getResource("no_index_fastq/i7").getPath();
		String barcodeFilename = classLoader.getResource("fastq/barcodes").getPath();
		
		// This file contains the barcode counts generated with version 1.9.0
		String barcodeCountsFilename = classLoader.getResource("no_index_fastq/barcode_counts").getPath();
		
		String i5Label = "GCATGGC";
		String i7Label = "TGACGTC";
		
		try {
			BarcodeMatcher i5 = new BarcodeMatcher(i5Filename, HAMMING_DISTANCE);
			BarcodeMatcher i7 = new BarcodeMatcher(i7Filename, HAMMING_DISTANCE);
			BarcodeMatcher barcodes = new BarcodeMatcher(barcodeFilename, HAMMING_DISTANCE);
			
			String [] count_args = {"--i5-indices", i5Filename, "--i7-indices", i7Filename, "--barcodes", barcodeFilename, 
					"--fixed-i5", i5Label, "--fixed-i7", i7Label, r1Filename, r2Filename};
			BarcodeCount.main(count_args);
			// checkout barcode count output
			Scanner s = new Scanner(outContent.toString());
			File barcodeCountStatisticsFile = new File(barcodeCountsFilename);
			BufferedReader reader = new BufferedReader(new FileReader(barcodeCountStatisticsFile));
			String expectedLine;
			while((expectedLine = reader.readLine()) != null){
				String line = s.nextLine();
				assertEquals(expectedLine, line);
			}
			s.close();
			reader.close();
			
			
			IndexAndBarcodeScreener screener = new IndexAndBarcodeScreener(i5, i7, barcodes);
			screener.setNumOutputFiles(2);
			
			screener.setBarcodeCountStatistics(new SampleSetsCounter(barcodeCountStatisticsFile));
			String mergedCountsFilename = tempFolder.getRoot() + "/counts";
			screener.setPrintStream(new FileOutputStream(mergedCountsFilename));
			
			String outputFileBase = tempFolder.getRoot() + "/test";
			screener.performScreeningMergeTrim(4, outputFileBase, r1Filename, r2Filename, null, null, i5Label, i7Label);
			// check merge and trim output 
			// TODO check merges
			SampleSetsCounter mergedCounts = new SampleSetsCounter(new File(mergedCountsFilename));
			String indexBarcodeKeyString = "GCATGGC_TGACGTC_TACGTTC:ACGTAAG:CGTACCT:GTACGGA_CAGGAAT:GCTTCCA:TGAAGGC:ATCCTTG";
			assertEquals(5, mergedCounts.get(indexBarcodeKeyString, "raw"));
			assertEquals(5, mergedCounts.get(indexBarcodeKeyString, "merged"));
		} catch (IOException | ParseException | InterruptedException | ExecutionException e) {
			fail();
		}
	}
	
	// If flowcell 
	@Test(expected = ExecutionException.class)
	public void fastqFlowcellMismatch() throws IOException, ParseException, InterruptedException, ExecutionException {
		ClassLoader classLoader = getClass().getClassLoader();
		String barcodeFilename = classLoader.getResource("Barcodes_5-7bp").getPath();
		String i5Filename = classLoader.getResource("i5").getPath();
		String i7Filename = classLoader.getResource("i7").getPath();

		String r1Filename = classLoader.getResource("bad_fastq_read_group/1.fastq.gz").getPath();
		String r2Filename = classLoader.getResource("bad_fastq_read_group/2.fastq.gz").getPath();

		List<String> args = new LinkedList<String>();
		args.add("--i5-indices");
		args.add(i5Filename);
		args.add("--i7-indices");
		args.add(i7Filename);
		args.add("--barcodes");
		args.add(barcodeFilename);
		args.add("--fixed-i5");
		args.add("CGAATAGA");
		args.add("--fixed-i7");
		args.add("AGTCTGGA");

		args.add(r1Filename);
		args.add(r2Filename);

		String parentDirectory = tempFolder.getRoot().toString();
		String outputRoot = parentDirectory + "/test";
		args.add(outputRoot);

		String [] args_array = new String[args.size()];
		int n = 0;
		for (String arg : args){
			args_array[n++] = arg;
		}

		IndexAndBarcodeScreener.main(args_array);
	}
	
	// If fastq files for a sample are demultiplexed, then we may have more than one flowcell lane in the same fastq
	@Test
	public void multipleFlowcellLanes() throws Exception {
		ClassLoader classLoader = getClass().getClassLoader();
		String barcodeFilename = classLoader.getResource("Barcodes_5-7bp").getPath();
		String i5Filename = classLoader.getResource("i5").getPath();
		String i7Filename = classLoader.getResource("i7").getPath();

		String r1Filename = classLoader.getResource("fastq_multiple_lanes/1.fastq.gz").getPath();
		String r2Filename = classLoader.getResource("fastq_multiple_lanes/2.fastq.gz").getPath();

		List<String> args = new LinkedList<String>();
		args.add("--i5-indices");
		args.add(i5Filename);
		args.add("--i7-indices");
		args.add(i7Filename);
		args.add("--barcodes");
		args.add(barcodeFilename);
		args.add("--fixed-i5");
		args.add("CGAATAGA");
		args.add("--fixed-i7");
		args.add("AGTCTGGA");
		args.add("--read-group-file");
		String readGroupFilename = tempFolder.newFile("readGroup").getPath();
		args.add(readGroupFilename);
		args.add("--disable-read-group-lane-check");

		args.add(r1Filename);
		args.add(r2Filename);

		String parentDirectory = tempFolder.getRoot().toString();
		String outputRoot = parentDirectory + "/test";
		args.add(outputRoot);

		String [] args_array = new String[args.size()];
		int n = 0;
		for (String arg : args){
			args_array[n++] = arg;
		}

		IndexAndBarcodeScreener.main(args_array);
	}
	
	@Test(expected = IllegalStateException.class)
	public void multipleFlowcellLanesWithFlowcellLaneCheck() throws Exception {
		ClassLoader classLoader = getClass().getClassLoader();
		String barcodeFilename = classLoader.getResource("Barcodes_5-7bp").getPath();
		String i5Filename = classLoader.getResource("i5").getPath();
		String i7Filename = classLoader.getResource("i7").getPath();

		String r1Filename = classLoader.getResource("fastq_multiple_lanes/1.fastq.gz").getPath();
		String r2Filename = classLoader.getResource("fastq_multiple_lanes/2.fastq.gz").getPath();

		List<String> args = new LinkedList<String>();
		args.add("--i5-indices");
		args.add(i5Filename);
		args.add("--i7-indices");
		args.add(i7Filename);
		args.add("--barcodes");
		args.add(barcodeFilename);
		args.add("--fixed-i5");
		args.add("CGAATAGA");
		args.add("--fixed-i7");
		args.add("AGTCTGGA");

		args.add(r1Filename);
		args.add(r2Filename);

		String parentDirectory = tempFolder.getRoot().toString();
		String outputRoot = parentDirectory + "/test";
		args.add(outputRoot);

		String [] args_array = new String[args.size()];
		int n = 0;
		for (String arg : args){
			args_array[n++] = arg;
		}

		IndexAndBarcodeScreener.main(args_array);
	}
}
