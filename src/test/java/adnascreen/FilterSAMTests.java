package adnascreen;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;

import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TemporaryFolder;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class FilterSAMTests {
	@Rule
	public TemporaryFolder testFolder = new TemporaryFolder();
	
	@Test
	public void testHit() {
		ClassLoader classLoader = getClass().getClassLoader();
		String bedFilename = classLoader.getResource("filter_with_base_quality/hit.bed").getPath();
		String samFilename = classLoader.getResource("filter_with_base_quality/filter.sam").getPath();
		
		try {
			FilterSAM filter = new FilterSAM(bedFilename);
			
			SamInputResource bufferedSAMFile = SamInputResource.of(new BufferedInputStream(new FileInputStream(samFilename)));
			SamReader reader = SamReaderFactory.makeDefault().open(bufferedSAMFile);
			SAMRecordIterator i = reader.iterator();
			SAMRecord record = i.next();
			
			assertTrue(filter.filter(record, 0, 20));
			assertTrue(filter.filter(record, 0, 32));
			// required mapping quality too high
			assertFalse(filter.filter(record, 30, 20));
			// required base quality is too high
			assertFalse(filter.filter(record, 0, 33));
		}catch(Exception e) {
			fail();
		}
	}
	
	@Test
	public void testMiss() {
		ClassLoader classLoader = getClass().getClassLoader();
		String bedFilename = classLoader.getResource("filter_with_base_quality/miss.bed").getPath();
		String samFilename = classLoader.getResource("filter_with_base_quality/filter.sam").getPath();
		
		try {
			FilterSAM filter = new FilterSAM(bedFilename);
			
			SamInputResource bufferedSAMFile = SamInputResource.of(new BufferedInputStream(new FileInputStream(samFilename)));
			SamReader reader = SamReaderFactory.makeDefault().open(bufferedSAMFile);
			SAMRecordIterator i = reader.iterator();
			SAMRecord record = i.next();
			
			assertFalse(filter.filter(record, 0, 20));
		}catch(Exception e) {
			fail();
		}
	}
	
	@Test
	public void testHitAndMiss() {
		ClassLoader classLoader = getClass().getClassLoader();
		String bedFilename = classLoader.getResource("filter_with_base_quality/hit_and_miss.bed").getPath();
		String samFilename = classLoader.getResource("filter_with_base_quality/filter.sam").getPath();
		
		try {
			FilterSAM filter = new FilterSAM(bedFilename);
			
			SamInputResource bufferedSAMFile = SamInputResource.of(new BufferedInputStream(new FileInputStream(samFilename)));
			SamReader reader = SamReaderFactory.makeDefault().open(bufferedSAMFile);
			SAMRecordIterator i = reader.iterator();
			SAMRecord record = i.next();
			
			assertTrue(filter.filter(record, 0, 20));
			assertTrue(filter.filter(record, 0, 32));
			// required mapping quality too high
			assertFalse(filter.filter(record, 30, 20));
			// required base quality is too high
			assertFalse(filter.filter(record, 0, 33));
		}catch(Exception e) {
			fail();
		}
	}
	
	@Test
	public void testTwoHits() {
		ClassLoader classLoader = getClass().getClassLoader();
		String bedFilename = classLoader.getResource("filter_with_base_quality/two_hits.bed").getPath();
		String samFilename = classLoader.getResource("filter_with_base_quality/filter.sam").getPath();
		
		try {
			FilterSAM filter = new FilterSAM(bedFilename);
			
			SamInputResource bufferedSAMFile = SamInputResource.of(new BufferedInputStream(new FileInputStream(samFilename)));
			SamReader reader = SamReaderFactory.makeDefault().open(bufferedSAMFile);
			SAMRecordIterator i = reader.iterator();
			SAMRecord record = i.next();
			
			assertTrue(filter.filter(record, 0, 20));
			assertTrue(filter.filter(record, 0, 32));
			assertTrue(filter.filter(record, 0, 33)); // one hit above this quality threshold
			// required mapping quality too high
			assertFalse(filter.filter(record, 30, 20));
		}catch(Exception e) {
			fail();
		}
	}
	
	@Test
	public void testHitLeadingEdge() {
		ClassLoader classLoader = getClass().getClassLoader();
		String bedFilename = classLoader.getResource("filter_with_base_quality/edge.bed").getPath();
		String samFilename = classLoader.getResource("filter_with_base_quality/filter.sam").getPath();
		
		try {
			FilterSAM filter = new FilterSAM(bedFilename);
			
			SamInputResource bufferedSAMFile = SamInputResource.of(new BufferedInputStream(new FileInputStream(samFilename)));
			SamReader reader = SamReaderFactory.makeDefault().open(bufferedSAMFile);
			SAMRecordIterator i = reader.iterator();
			SAMRecord record = i.next();
			
			assertTrue(filter.filter(record, 0, 20));
		}catch(Exception e) {
			fail();
		}
	}
	
	@Test
	public void testHitTrailingEdgeWithDeletion() {
		ClassLoader classLoader = getClass().getClassLoader();
		String bedFilename = classLoader.getResource("filter_with_base_quality/deletion_end.bed").getPath();
		String samFilename = classLoader.getResource("filter_with_base_quality/deletion_end.sam").getPath();
		
		try {
			FilterSAM filter = new FilterSAM(bedFilename);
			
			SamInputResource bufferedSAMFile = SamInputResource.of(new BufferedInputStream(new FileInputStream(samFilename)));
			SamReader reader = SamReaderFactory.makeDefault().open(bufferedSAMFile);
			SAMRecordIterator i = reader.iterator();
			SAMRecord record = i.next();
			
			assertTrue(filter.filter(record, 30, 20));
		}catch(IOException e) {
			fail();
		}
	}
	
	// Deletion at the SNP
	@Test
	public void testDeletion() {
		ClassLoader classLoader = getClass().getClassLoader();
		String bedFilename = classLoader.getResource("filter_with_base_quality/deletion.bed").getPath();
		String samFilename = classLoader.getResource("filter_with_base_quality/deletion.sam").getPath();
		
		try {
			FilterSAM filter = new FilterSAM(bedFilename);
			
			SamInputResource bufferedSAMFile = SamInputResource.of(new BufferedInputStream(new FileInputStream(samFilename)));
			SamReader reader = SamReaderFactory.makeDefault().open(bufferedSAMFile);
			SAMRecordIterator i = reader.iterator();
			SAMRecord record = i.next();
			
			assertFalse(filter.filter(record, 0, 20));
		}catch(IOException e) {
			fail();
		}
	}
	
	@Test
	public void testMainNoClipping() {		
		ClassLoader classLoader = getClass().getClassLoader();
		String filename = classLoader.getResource("filter_with_base_quality/filter.sam").getPath();
		String bedFilename = classLoader.getResource("filter_with_base_quality/edge.bed").getPath();
		try {
			File outputFile = testFolder.newFile("clipped.sam");
			String outputFilename = outputFile.getAbsolutePath();
			String[] args = {"-i", filename, "-o", outputFilename, "-m", "0", "-q", "20", "-p", bedFilename};
			FilterSAM.main(args);
			
			// verify there is one read, and which contains target
			SamInputResource bufferedSAMFile = SamInputResource.of(new BufferedInputStream(new FileInputStream(outputFilename)));
			SamReader reader = SamReaderFactory.makeDefault().open(bufferedSAMFile);
			SAMRecordIterator i = reader.iterator();
			SAMRecord record = i.next();
			
			assertEquals(1, record.getReadPositionAtReferencePosition(31445));
		}
		catch(Exception e) {
			fail();
		}
	}

	@Test
	public void testMainWithClipping() {
		ClassLoader classLoader = getClass().getClassLoader();
		String filename = classLoader.getResource("filter_with_base_quality/filter.sam").getPath();
		String bedFilename = classLoader.getResource("filter_with_base_quality/edge.bed").getPath();
		try {
			File outputFile = testFolder.newFile("clipped.sam");
			String outputFilename = outputFile.getAbsolutePath();
			String[] args = {"-i", filename, "-o", outputFilename, "-m", "0", "-q", "20", "-p", bedFilename, "-n", "2"};
			FilterSAM.main(args);
			
			// verify there is are no reads
			SamInputResource bufferedSAMFile = SamInputResource.of(new BufferedInputStream(new FileInputStream(outputFilename)));
			SamReader reader = SamReaderFactory.makeDefault().open(bufferedSAMFile);
			SAMRecordIterator i = reader.iterator();
			assertFalse(i.hasNext());
		}
		catch(Exception e) {
			fail();
		}
	}
}
