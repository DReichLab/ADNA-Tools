package adnascreen;

import static org.junit.Assert.*;

import java.io.IOException;

import org.junit.Test;

public class TwoPassBarcodeTests {
	@Test
	public void priorPassMax(){
		SampleSetsCounter barcodeCounts = new SampleSetsCounter();
		
		ClassLoader classLoader = getClass().getClassLoader();
		String filename = classLoader.getResource("Barcodes_5-7bp").getPath();		
		try {
			BarcodeMatcher barcodeMatcher = new BarcodeMatcher(filename, 1);
			IndexAndBarcodeKey indexOnlyKey = new IndexAndBarcodeKey("1", "2", null, null);
			barcodeCounts.add(indexOnlyKey, "Q2_Q1", 1);
			barcodeCounts.add(indexOnlyKey, "Q1_Q2", 9);
			barcodeCounts.add(indexOnlyKey, "Q3_Q1", 5);
			
			int length;
			length = IndexAndBarcodeScreener.barcodeLengthFromPriorPassCounts(barcodeCounts, indexOnlyKey, barcodeMatcher);
			assertEquals(7, length);
			
			barcodeCounts.add(indexOnlyKey, "5-MM-1_5-MM-2", 100);
			length = IndexAndBarcodeScreener.barcodeLengthFromPriorPassCounts(barcodeCounts, indexOnlyKey, barcodeMatcher);
			assertEquals(5, length);
			
			barcodeCounts.add(indexOnlyKey, BarcodeCount.WITHOUT_BARCODES, 1000);
			length = IndexAndBarcodeScreener.barcodeLengthFromPriorPassCounts(barcodeCounts, indexOnlyKey, barcodeMatcher);
			assertEquals(0, length);
		} catch (IOException e) {
			fail();
		}
	}

}
