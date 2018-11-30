package adnascreen;

import static org.junit.Assert.*;

import java.io.IOException;

import org.junit.Test;

public class TwoPassBarcodeTests {
	@Test
	public void priorPassMax(){
		SampleSetsCounter barcodeCounts = new SampleSetsCounter();
		float threshold = 0.05f;
		
		ClassLoader classLoader = getClass().getClassLoader();
		String filename = classLoader.getResource("Barcodes_5-7bp").getPath();		
		try {
			BarcodeMatcher barcodeMatcher = new BarcodeMatcher(filename, 1);
			IndexAndBarcodeKey indexOnlyKey = new IndexAndBarcodeKey("1", "2", null, null);
			String indexOnlyKeyString = indexOnlyKey.toString();
			barcodeCounts.add(indexOnlyKeyString, "Q2_Q1", 1);
			barcodeCounts.add(indexOnlyKeyString, "Q1_Q2", 9);
			barcodeCounts.add(indexOnlyKeyString, "Q3_Q1", 5);
			
			int length;
			length = IndexAndBarcodeScreener.barcodeLengthFromPriorPassCounts(barcodeCounts, indexOnlyKey, barcodeMatcher, threshold);
			assertEquals(7, length);
			
			barcodeCounts.add(indexOnlyKeyString, "5-MM-1_5-MM-2", 100);
			length = IndexAndBarcodeScreener.barcodeLengthFromPriorPassCounts(barcodeCounts, indexOnlyKey, barcodeMatcher, threshold);
			assertEquals(5, length);
			
			barcodeCounts.add(indexOnlyKeyString, BarcodeCount.WITHOUT_BARCODES, 100000);
			length = IndexAndBarcodeScreener.barcodeLengthFromPriorPassCounts(barcodeCounts, indexOnlyKey, barcodeMatcher, threshold);
			assertEquals(0, length);
		} catch (IOException e) {
			fail();
		}
	}

	@Test
	public void lowBarcodePass(){
		SampleSetsCounter barcodeCounts = new SampleSetsCounter();
		float threshold = 0.05f;
		
		ClassLoader classLoader = getClass().getClassLoader();
		String filename = classLoader.getResource("Barcodes_5-7bp").getPath();		
		try {
			BarcodeMatcher barcodeMatcher = new BarcodeMatcher(filename, 1);
			IndexAndBarcodeKey indexOnlyKey = new IndexAndBarcodeKey("1", "2", null, null);
			int length;
			
			// barcode count is over the threshold
			barcodeCounts.add(indexOnlyKey.toString(), "5-MM-1_5-MM-2", 1000);
			barcodeCounts.add(indexOnlyKey.toString(), BarcodeCount.WITHOUT_BARCODES, 10000);
			length = IndexAndBarcodeScreener.barcodeLengthFromPriorPassCounts(barcodeCounts, indexOnlyKey, barcodeMatcher, threshold);
			assertEquals(5, length);
		} catch (IOException e) {
			fail();
		}
	}
	
	@Test
	public void noBarcodes(){
		SampleSetsCounter barcodeCounts = new SampleSetsCounter();
		float threshold = 0.05f;
		
		ClassLoader classLoader = getClass().getClassLoader();
		String filename = classLoader.getResource("Barcodes_5-7bp").getPath();		
		try {
			BarcodeMatcher barcodeMatcher = new BarcodeMatcher(filename, 1);
			IndexAndBarcodeKey indexOnlyKey = new IndexAndBarcodeKey("1", "2", null, null);
			int length;
			
			barcodeCounts.add(indexOnlyKey.toString(), BarcodeCount.WITHOUT_BARCODES, 10000);
			length = IndexAndBarcodeScreener.barcodeLengthFromPriorPassCounts(barcodeCounts, indexOnlyKey, barcodeMatcher, threshold);
			assertEquals(0, length);
		} catch (IOException e) {
			fail();
		}
	}
}
