package adnascreen;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.BufferedReader;
import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.zip.GZIPInputStream;

import org.junit.After;
import org.junit.Before;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TemporaryFolder;

import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;

public class FASTQFilterTests {
	protected int minimum_length = 30;
	
	private final ByteArrayOutputStream outContent = new ByteArrayOutputStream();
	private final PrintStream originalSystemOut = System.out;
	private final InputStream originalSystemIn = System.in;
	
	@Rule
	public TemporaryFolder testFolder = new TemporaryFolder();
	
	@Before
	public void setUpStreams() {
	    System.setOut(new PrintStream(outContent));
	}

	@After
	public void restoreStreams() {
	    System.setOut(originalSystemOut);
	    System.setIn(originalSystemIn);
	}
	
	@Test
	public void fileOutput() throws IOException{
		ClassLoader classLoader = getClass().getClassLoader();
		String fastqFilename = classLoader.getResource("fastq/ultima.fastq.gz").getPath();
		
		File filteredFile = testFolder.newFile("filtered.fastq.gz");
		String filteredFilename = filteredFile.getAbsolutePath();
		String [] args = {"-n", String.valueOf(minimum_length), "-i", fastqFilename, "-o", filteredFilename};
		
		FASTQFilter.main(args);
		
		// read back from disk
		FastqReader in = new FastqReader(new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(filteredFilename)))));
		checkFiltered(in);
		in.close();
	}
	
	@Test
	public void fileOutputAll() throws IOException{
		ClassLoader classLoader = getClass().getClassLoader();
		String fastqFilename = classLoader.getResource("fastq/ultima.fastq.gz").getPath();
		int mimimum_length_short = 10;
		
		File filteredFile = testFolder.newFile("filtered.fastq.gz");
		String filteredFilename = filteredFile.getAbsolutePath();
		String [] args = {"-n", String.valueOf(mimimum_length_short), "-i", fastqFilename, "-o", filteredFilename};
		
		FASTQFilter.main(args);
		
		// read back from disk
		FastqReader in = new FastqReader(new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(filteredFilename)))));
		
		ArrayList<FastqRecord> records = new ArrayList<FastqRecord>(); 
		while(in.hasNext()) {
			FastqRecord record = in.next();
			assertTrue(record.getReadLength() >= mimimum_length_short);
			records.add(record);
		}
		
		assertEquals(7, records.size());
		in.close();
	}
	
	@Test
	public void piped() throws IOException{
		ClassLoader classLoader = getClass().getClassLoader();
		String fastqFilename = classLoader.getResource("fastq/ultima.fastq.gz").getPath();
		
		String [] args = {"-n", String.valueOf(minimum_length)};
		
		System.setIn(new GZIPInputStream(new FileInputStream(fastqFilename)));
		FASTQFilter.main(args);
		
		FastqReader in = new FastqReader(new BufferedReader(new InputStreamReader(new ByteArrayInputStream(outContent.toByteArray()))));
		checkFiltered(in);
		in.close();
	}
	
	protected void checkFiltered(FastqReader in) {
		ArrayList<FastqRecord> records = new ArrayList<FastqRecord>(); 
		while(in.hasNext()) {
			FastqRecord record = in.next();
			assertTrue(record.getReadLength() >= minimum_length);
			records.add(record);
		}
		assertEquals("034278_1-Z0009-5345032633", records.get(0).getReadName());
		assertEquals("042417_2-Z0009-2743414696", records.get(1).getReadName());
		assertEquals("034275_2-Z0009-2499492645", records.get(2).getReadName());
		assertEquals("042417_2-Z0009-4720812868", records.get(3).getReadName());
		assertEquals(4, records.size());
	}
}
