package adnascreen;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.fail;

import org.junit.Test;

public class DuplicateHistogramTests {
	@Test
	public void testSingleRead() {
		ClassLoader classLoader = getClass().getClassLoader();
		String samFilename = classLoader.getResource("filter_with_base_quality/filter.sam").getPath();
		
		try {
			DuplicatesHistogram duplicatesHistogram = new DuplicatesHistogram(samFilename);
			int [] histogram = duplicatesHistogram.getHistogram();
			assertEquals(1, histogram.length);
			assertEquals(1, histogram[0]);
		}
		catch(Exception e) {
			fail();
		}
	}
	
	@Test
	public void testEmpty() {
		DuplicatesHistogram duplicatesHistogram = new DuplicatesHistogram();
		int[] histogram = duplicatesHistogram.getHistogram();
		assertEquals(0, histogram.length);
	}
	
	@Test
	public void testDoctoredMT() {
		ClassLoader classLoader = getClass().getClassLoader();
		String samFilename = classLoader.getResource("unique_read_histogram/doctored_histogram.sam").getPath();
		
		DuplicatesHistogram duplicatesHistogram;
		try {
			duplicatesHistogram = new DuplicatesHistogram(samFilename);
		}
		catch(Exception e) {
			fail();
			return;
		}
		int [] histogram = duplicatesHistogram.getHistogram();
		assertEquals(4, histogram.length);
		assertEquals(23, histogram[0]);
		assertEquals(0, histogram[1]);
		assertEquals(1, histogram[2]);
		assertEquals(1, histogram[3]);
	}
}
