package adnascreen;

import static org.junit.Assert.*;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;

import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TemporaryFolder;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class CacheTest {
	String a = "a";
	Integer A = 1;
	String b = "b";
	Integer B = 2;
	String c = "c";
	Integer C = 3;
	String d = "d";
	Integer D = 4;
	String e = "e";
	Integer E = 5;

	@Test
	public void singleEntry(){
		Cache<String, Integer> cache = new Cache<String, Integer>(3);
		Integer A1 = 1;
		Integer A2 = 2;
		cache.put(a, A1);
		assertTrue(cache.inCache(a));
		assertEquals(cache.get(a), A1);
		assertEquals(cache.size(), 1);
		cache.put(a, A1); // duplicate
		assertTrue(cache.inCache(a));
		assertEquals(cache.get(a), A1);
		assertEquals(cache.size(), 1);
		cache.put(a, A2); // overwrite
		assertTrue(cache.inCache(a));
		assertEquals(cache.get(a), A2);
		assertEquals(cache.size(), 1);
		cache.remove(a);
		assertFalse(cache.inCache(a));
		assertEquals(cache.get(a), null);
		assertEquals(cache.size(), 0);
	}

	@Test
	public void eviction(){
		Cache<String, Integer> cache = new Cache<String, Integer>(3);

		cache.put(a, A);
		assertEquals(cache.get(a), A);
		assertEquals(cache.size(), 1);
		cache.put(b, B);
		assertEquals(cache.get(a), A);
		assertEquals(cache.get(b), B);
		assertEquals(cache.size(), 2);
		cache.put(c, C);
		assertEquals(cache.get(a), A);
		assertEquals(cache.get(b), B);
		assertEquals(cache.get(c), C);
		assertEquals(cache.size(), 3);
		cache.put(d, D);
		assertFalse(cache.inCache(a));
		assertEquals(cache.get(a), null);
		assertEquals(cache.get(b), B);
		assertEquals(cache.get(c), C);
		assertEquals(cache.get(d), D);
		assertEquals(cache.size(), 3);
		cache.put(e, E);
		assertEquals(cache.get(c), C);
		assertEquals(cache.get(d), D);
		assertEquals(cache.get(e), E);
		assertEquals(cache.size(), 3);
	}

	@Test
	public void reorder(){
		Cache<String, Integer> cache = new Cache<String, Integer>(3);
		cache.put(a, A);
		cache.put(b, B);
		cache.put(c, C);
		assertEquals(cache.get(a), A); // b is last accessed
		cache.put(d, D);
		assertFalse(cache.inCache(b));
		assertEquals(cache.get(a), A);
		assertEquals(cache.get(b), null);
		assertEquals(cache.get(c), C);
		assertEquals(cache.get(d), D);
		assertEquals(cache.size(), 3);
	}

	@Test
	public void reorder2(){
		Cache<String, Integer> cache = new Cache<String, Integer>(3);
		cache.put(a, A);
		cache.put(b, B);
		cache.put(c, C);
		assertEquals(cache.get(a), A); // b is last accessed
		assertEquals(cache.get(b), B); // C is last accessed
		assertEquals(cache.size(), 3);
		cache.put(d, D);
		assertTrue(cache.inCache(a));
		assertTrue(cache.inCache(b));
		assertFalse(cache.inCache(c));
		assertEquals(cache.get(c), null);
		assertEquals(cache.get(a), A);
		assertEquals(cache.get(b), B);
		assertEquals(cache.get(d), D);
		assertEquals(cache.size(), 3);
	}

	@Rule
	public TemporaryFolder tempFolder = new TemporaryFolder();
	// Write sequence of numbers to files, and read numbers back again
	@Test
	public void fileTest() throws IOException{
		final int NUM_FILES = 4;
		final int NUM_WRITES = 10000;
		Cache<Integer, PrintWriter> cache = new Cache<Integer, PrintWriter>(NUM_FILES-1);
		String [] filePaths = new String[NUM_FILES];
		// create files in temporary folder
		for(int fileIndex = 0; fileIndex < NUM_FILES; fileIndex++){
			File f = tempFolder.newFile(String.valueOf(fileIndex));
			filePaths[fileIndex] = f.getAbsolutePath(); // save filenames
			assertTrue(f.exists());
		}
		// write integers to files
		// this is setup to require a new file opening for each write
		// there are no cache hits, only misses
		for (int n = 0; n < NUM_WRITES; n++){
			for(int fileIndex = 0; fileIndex < NUM_FILES; fileIndex++){
				PrintWriter output = cache.get(fileIndex);
				if(output == null){
					// reopen file
					output = new PrintWriter(new BufferedOutputStream(new FileOutputStream(filePaths[fileIndex], true)));
					cache.put(fileIndex, output);
				}
				output.println(n);
			}
		}
		cache.clear();
		assertEquals(0l, cache.getHits());
		assertEquals((long) (NUM_FILES * NUM_WRITES), cache.getMisses());
		assertEquals((long) (NUM_FILES * NUM_WRITES - NUM_FILES + 1), cache.getForcedCloses());

		// check file contents
		for(int fileIndex = 0; fileIndex < NUM_FILES; fileIndex++){
			try(BufferedReader in = new BufferedReader(new FileReader(filePaths[fileIndex]))){
				for(int n = 0; n < NUM_WRITES; n++){
					String line = in.readLine();
					int value = Integer.valueOf(line);
					assertEquals(n, value);
				}
			} catch(IOException e){
				fail(e.toString());
			}
		}
	}

	// This turns out to not be a useful test because SAM and BAM files cannot be reopened for writing. 
	@Test
	public void samMemoryLeak() throws IOException{
		ClassLoader classLoader = getClass().getClassLoader();
		String filename = classLoader.getResource("test.sam").getPath();
		SAMFileHeader header = null;
		SAMRecord record = null;
		SamInputResource bufferedSAMFile = SamInputResource.of(new BufferedInputStream(new FileInputStream(filename)));
		try(
				SamReader reader = SamReaderFactory.makeDefault().open(bufferedSAMFile);
				){
			header = reader.getFileHeader();
			SAMRecordIterator i = reader.iterator();
			record = i.next();
		}

		File f = tempFolder.newFile("samMemoryLeakOutput");
		final int NUM_FILES = 1000;
		SAMFileWriterFactory outputFileFactory = new SAMFileWriterFactory();
		for(int n = 0; n < NUM_FILES; n++){
			try(SAMFileWriter output = outputFileFactory.makeSAMWriter(header, true, f);){
				SAMRecord copy = record.deepCopy();
				output.addAlignment(copy);
			}
		}
	}
}
