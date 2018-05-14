package adnascreen;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.fail;

import java.io.IOException;
import java.util.Map;

import org.apache.commons.cli.ParseException;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.ExpectedException;

public class BarcodeLengthsFromSampleSheetTests {
	@Rule
	public ExpectedException thrown = ExpectedException.none();
	
	private Map<IndexAndBarcodeKey, Integer> setupTest() throws IOException, ParseException {
		ClassLoader classLoader = getClass().getClassLoader();
		String sampleSheetFilename = classLoader.getResource("example.index_barcode_keys").getPath();
		
		BarcodeMatcher barcodeMatcher = new BarcodeMatcher();
		barcodeMatcher.setMaxHammingDistance(1);
		barcodeMatcher.addReferenceSet("ACAACC", "ACAACC");
		barcodeMatcher.addReferenceSet("CGCCATG:GTGGCAT:TATTGCA:ACAATGC", "CGCCATG:GTGGCAT:TATTGCA:ACAATGC");
		barcodeMatcher.addReferenceSet("CATAGGC:GCACTTG:TGCGAAT:ATGTCCA", "CATAGGC:GCACTTG:TGCGAAT:ATGTCCA");
		
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
		try {
			thrown.expect(IllegalStateException.class);
			IndexAndBarcodeScreener.barcodeLengthsByIndexPair(sampleSheetFilename, barcodeMatcher);
			fail("Invalid barcodes for index pair were accepted");
		}
		catch(IOException | ParseException e) {
			fail("Wrong exception " + e.toString());
		}
	}
}
