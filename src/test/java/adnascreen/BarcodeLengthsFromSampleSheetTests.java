package adnascreen;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertThrows;
import static org.junit.Assert.fail;

import java.io.IOException;
import java.util.Map;

import org.apache.commons.cli.ParseException;
import org.junit.Test;

public class BarcodeLengthsFromSampleSheetTests {	
	private Map<IndexAndBarcodeKey, Integer> setupTest() throws IOException, ParseException {
		ClassLoader classLoader = getClass().getClassLoader();
		String sampleSheetFilename = classLoader.getResource("example.index_barcode_keys").getPath();
		
		BarcodeMatcher barcodeMatcher = new BarcodeMatcher();
		barcodeMatcher.setMaxHammingDistance(1);
		barcodeMatcher.addReferenceSet("ACAACC", "ACAACC");
		barcodeMatcher.addReferenceSet("CGCCATG:GTGGCAT:TATTGCA:ACAATGC", "CGCCATG:GTGGCAT:TATTGCA:ACAATGC");
		barcodeMatcher.addReferenceSet("CATAGGC:GCACTTG:TGCGAAT:ATGTCCA", "CATAGGC:GCACTTG:TGCGAAT:ATGTCCA");
		barcodeMatcher.addReferenceSet("GGTATCG:TTACAGT:AACGCTA:CCGTGAC", "GGTATCG:TTACAGT:AACGCTA:CCGTGAC"); // Q2
		
		Map<IndexAndBarcodeKey, Integer> lengths = IndexAndBarcodeScreener.barcodeLengthsByIndexPair(sampleSheetFilename, barcodeMatcher);
		return lengths;
	}
	
	@Test
	public void testLength6() {
		try{
			Map<IndexAndBarcodeKey, Integer> lengths = setupTest();
			int queryLength;
			
			IndexAndBarcodeKey query6 = new IndexAndBarcodeKey("CGCCGTC_AGCGCCA__");
			queryLength = lengths.get(query6);
			assertEquals(6, queryLength);
		}
		catch(Exception e){
			fail(e.toString());
		}
	}
	
	@Test
	public void testLength7() {
		try{
			Map<IndexAndBarcodeKey, Integer> lengths = setupTest();
			int queryLength;
			
			IndexAndBarcodeKey query7 = new IndexAndBarcodeKey("CGCCGTC_TTGGTCA__");
			queryLength = lengths.get(query7);
			assertEquals(7, queryLength);
		}
		catch(Exception e){
			fail(e.toString());
		}
	}
	
	// barcode is part of Q set, but not all four
	@Test
	public void testPartialQSet() {
		try{
			Map<IndexAndBarcodeKey, Integer> lengths = setupTest();
			int queryLength;
			
			IndexAndBarcodeKey query7 = new IndexAndBarcodeKey("GATCCAA_GCTTCAG__");
			queryLength = lengths.get(query7);
			assertEquals(7, queryLength);
		}
		catch(Exception e){
			fail(e.toString());
		}
	}
	
	@Test
	public void empty() {
		try{
			Map<IndexAndBarcodeKey, Integer> lengths = setupTest();
			
			IndexAndBarcodeKey queryNotFound = new IndexAndBarcodeKey("AAAAAAA_AAAAAAA__");
			assertNull(lengths.get(queryNotFound));
		}
		catch(Exception e){
			fail(e.toString());
		}
	}
	
	@Test
	public void illegalSampleSheet() { // one index pair contains multiple barcode lengths	
		ClassLoader classLoader = getClass().getClassLoader();
		String sampleSheetFilename = classLoader.getResource("illegal.index_barcode_keys").getPath();
		
		BarcodeMatcher barcodeMatcher = new BarcodeMatcher();
		barcodeMatcher.setMaxHammingDistance(1);
		barcodeMatcher.addReferenceSet("ACAACC", "ACAACC");
		barcodeMatcher.addReferenceSet("CGCCATG", "CGCCATG");
		barcodeMatcher.addReferenceSet("CATAGGC:GCACTTG:TGCGAAT:ATGTCCA", "CATAGGC:GCACTTG:TGCGAAT:ATGTCCA");
		assertThrows(IllegalStateException.class, () -> IndexAndBarcodeScreener.barcodeLengthsByIndexPair(sampleSheetFilename, barcodeMatcher));
	}
	
	@Test
	public void no_barcodes() {
		ClassLoader classLoader = getClass().getClassLoader();
		String sampleSheetFilename = classLoader.getResource("no_barcodes.index_barcode_keys").getPath();
		
		BarcodeMatcher barcodeMatcher = new BarcodeMatcher();
		int queryLength;
		
		try {
			Map<IndexAndBarcodeKey, Integer> lengths = IndexAndBarcodeScreener.barcodeLengthsByIndexPair(sampleSheetFilename, barcodeMatcher);
			IndexAndBarcodeKey query0 = new IndexAndBarcodeKey("CGAGATC_CCGTTGA__");
			queryLength = lengths.get(query0);
			assertEquals(0, queryLength);
		}
		catch(IOException | ParseException e) {
			fail("");
		}
	}
}
