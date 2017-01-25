package adnascreen;

import static org.junit.Assert.*;

import java.io.IOException;

import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.ExpectedException;

public class BarcodeMatcherTests {
	@Rule
	public ExpectedException thrown = ExpectedException.none();
	
	@Test
	public void exactMatch(){
		String barcodeSet = "ATCGATT:CAGTCAA:GCTAGCC:TGACTGG";

		DNASequence query = new DNASequence("ATCGATT");

		BarcodeMatcher barcodeMatcher = new BarcodeMatcher();
		barcodeMatcher.addReferenceSet(barcodeSet);
		String result = barcodeMatcher.find(query);
		assertEquals(result, barcodeSet);
	}
	
	@Test
	public void noMatch(){
		String barcodeSet = "ATCGATT:CAGTCAA:GCTAGCC:TGACTGG";

		DNASequence query = new DNASequence("AAAAAAA");

		BarcodeMatcher barcodeMatcher = new BarcodeMatcher();
		barcodeMatcher.addReferenceSet(barcodeSet);
		String result = barcodeMatcher.find(query);
		assertNull(result);
	}
	
	@Test
	public void matchDiff1(){
		String barcodeSet = "ATCGATT:CAGTCAA:GCTAGCC:TGACTGG";

		DNASequence query = new DNASequence("ATCGATG");

		BarcodeMatcher barcodeMatcher = new BarcodeMatcher();
		barcodeMatcher.setMaxHammingDistance(1);
		barcodeMatcher.addReferenceSet(barcodeSet);
		String result = barcodeMatcher.find(query);
		assertEquals(result, barcodeSet);
	}
	
	@Test
	public void invalidReferenceSet(){
		thrown.expect(IllegalArgumentException.class);
		
		String barcodeSet = "ATCGATT:CAGTCAA:GCTAGCC:TGACTGG:A";
		BarcodeMatcher barcodeMatcher = new BarcodeMatcher();
		barcodeMatcher.addReferenceSet(barcodeSet);
		fail("Invalid barcode set was accepted");
	}
	
	@Test
	public void file_barcodes(){
		try{
			ClassLoader classLoader = getClass().getClassLoader();
			String filename = classLoader.getResource("barcodes").getPath();
			
			BarcodeMatcher barcodeMatcher = new BarcodeMatcher(filename, 1);
			String result;
			// #25
			String expected = "CTTCCGA:GAAGGTC:TCCTTAG:AGGAACT";
			
			DNASequence query1 = new DNASequence("GAAGGTC");
			result = barcodeMatcher.find(query1);
			assertEquals(expected, result);
			
			DNASequence query2 = new DNASequence("GAAGGTT");
			result = barcodeMatcher.find(query2);
			assertEquals(expected, result);
			
			// repeat to use cache
			result = barcodeMatcher.find(query2);
			assertEquals(expected, result);
			
			DNASequence queryInteractive = new DNASequence("ACTGCNT");
			result = barcodeMatcher.find(queryInteractive);
			System.out.println(result);
		}
		catch(IOException e){
			fail(e.toString());
		}
	}
}
