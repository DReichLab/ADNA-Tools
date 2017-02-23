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
		String label = "Q22";
		DNASequence query = new DNASequence("ATCGATT");

		BarcodeMatcher barcodeMatcher = new BarcodeMatcher();
		barcodeMatcher.addReferenceSet(barcodeSet, label);
		String result = barcodeMatcher.find(query);
		assertEquals(label, result);
	}
	
	@Test
	public void noMatch(){
		String barcodeSet = "ATCGATT:CAGTCAA:GCTAGCC:TGACTGG";
		String label = "Q22";
		DNASequence query = new DNASequence("AAAAAAA");

		BarcodeMatcher barcodeMatcher = new BarcodeMatcher();
		barcodeMatcher.addReferenceSet(barcodeSet, label);
		String result = barcodeMatcher.find(query);
		assertNull(result);
	}
	
	@Test
	public void matchDiff1(){
		String barcodeSet = "ATCGATT:CAGTCAA:GCTAGCC:TGACTGG";
		String label = "Q22";
		DNASequence query = new DNASequence("ATCGATG");

		BarcodeMatcher barcodeMatcher = new BarcodeMatcher();
		barcodeMatcher.setMaxHammingDistance(1);
		
		barcodeMatcher.addReferenceSet(barcodeSet, label);
		String result = barcodeMatcher.find(query);
		assertEquals(label, result);
	}
	
	@Test
	public void invalidReferenceSet(){
		thrown.expect(IllegalArgumentException.class);
		
		String barcodeSet = "ATCGATT:CAGTCAA:GCTAGCC:TGACTGG:A";
		String label = "Q22";
		BarcodeMatcher barcodeMatcher = new BarcodeMatcher();
		barcodeMatcher.addReferenceSet(barcodeSet, label);
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
			String expected = "Q25"; //CTTCCGA:GAAGGTC:TCCTTAG:AGGAACT
			
			DNASequence query1 = new DNASequence("GAAGGTC");
			result = barcodeMatcher.find(query1);
			assertEquals(expected, result);
			
			DNASequence query2 = new DNASequence("GAAGGTT");
			result = barcodeMatcher.find(query2);
			assertEquals(expected, result);
			
			// repeat to use cache
			result = barcodeMatcher.find(query2);
			assertEquals(expected, result);
			
			/*
			DNASequence queryInteractive = new DNASequence("ACCTTGT");
			result = barcodeMatcher.find(queryInteractive);
			System.out.println(result);
			*/
		}
		catch(IOException e){
			fail(e.toString());
		}
	}
	
	@Test
	public void duplicateLabel(){
		thrown.expect(IllegalArgumentException.class);
		
		String barcodeSet = "ATCGATT:CAGTCAA:GCTAGCC:TGACTGG";
		String label = "Q22";
		BarcodeMatcher barcodeMatcher = new BarcodeMatcher();
		barcodeMatcher.addReferenceSet(barcodeSet, label);
		String barcodeSet2 = "ATCGACC";
		barcodeMatcher.addReferenceSet(barcodeSet2, label);
		
		fail("Label was allowed twice");
	}
}
