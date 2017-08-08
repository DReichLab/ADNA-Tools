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
		String label = "Q1";
		String expectedResult = "Q1" + BarcodeMatcher.INDEX_DELIMITER + "1";
		DNASequence query = new DNASequence("ATCGATT");

		BarcodeMatcher barcodeMatcher = new BarcodeMatcher();
		barcodeMatcher.addReferenceSet(barcodeSet, label);
		String result = barcodeMatcher.find(query);
		assertEquals(expectedResult, result);
	}
	
	@Test
	public void noMatch(){
		String barcodeSet = "ATCGATT:CAGTCAA:GCTAGCC:TGACTGG";
		String label = "Q1";
		DNASequence query = new DNASequence("AAAAAAA");

		BarcodeMatcher barcodeMatcher = new BarcodeMatcher();
		barcodeMatcher.addReferenceSet(barcodeSet, label);
		String result = barcodeMatcher.find(query);
		assertNull(result);
	}
	
	@Test
	public void matchDiff1(){
		String barcodeSet = "ATCGATT:CAGTCAA:GCTAGCC:TGACTGG";
		String label = "Q1";
		String expectedResult = "Q1" + BarcodeMatcher.INDEX_DELIMITER + "1";
		DNASequence query = new DNASequence("ATCGATG");

		BarcodeMatcher barcodeMatcher = new BarcodeMatcher();
		barcodeMatcher.setMaxHammingDistance(1);
		
		barcodeMatcher.addReferenceSet(barcodeSet, label);
		String result = barcodeMatcher.find(query);
		assertEquals(expectedResult, result);
	}
	
	@Test
	public void matchDiff1_fourthElement(){
		String barcodeSet = "ATCGATT:CAGTCAA:GCTAGCC:TGACTGG";
		String label = "Q1";
		String expectedResult = "Q1" + BarcodeMatcher.INDEX_DELIMITER + "4";
		DNASequence query = new DNASequence("TGACTTG");

		BarcodeMatcher barcodeMatcher = new BarcodeMatcher();
		barcodeMatcher.setMaxHammingDistance(1);
		
		barcodeMatcher.addReferenceSet(barcodeSet, label);
		String result = barcodeMatcher.find(query);
		assertEquals(expectedResult, result);
	}
	
	@Test
	public void invalidReferenceSet(){
		thrown.expect(IllegalArgumentException.class);
		
		String barcodeSet = "ATCGATT:CAGTCAA:GCTAGCC:TGACTGG:A";
		String label = "Q1";
		BarcodeMatcher barcodeMatcher = new BarcodeMatcher();
		barcodeMatcher.addReferenceSet(barcodeSet, label);
		fail("Invalid barcode set was accepted");
	}
	
	@Test
	public void file_barcodes(){
		try{
			ClassLoader classLoader = getClass().getClassLoader();
			String filename = classLoader.getResource("7bpBarcodes_Reich20170725").getPath();
			
			BarcodeMatcher barcodeMatcher = new BarcodeMatcher(filename, 1);
			String result;
			//CTTCCGA:GAAGGTC:TCCTTAG:AGGAACT
			String expected = "Q34" + BarcodeMatcher.INDEX_DELIMITER + "2"; 
			
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
		String label = "Q1";
		BarcodeMatcher barcodeMatcher = new BarcodeMatcher();
		barcodeMatcher.addReferenceSet(barcodeSet, label);
		String barcodeSet2 = "ATCGACC";
		barcodeMatcher.addReferenceSet(barcodeSet2, label);
		
		fail("Label was allowed twice");
	}
	
	@Test
	public void invalidLabelWithDelimiter(){
		thrown.expect(IllegalArgumentException.class);
		
		String barcodeSet = "ATCGATT:CAGTCAA:GCTAGCC:TGACTGG" + BarcodeMatcher.INDEX_DELIMITER;
		String label = "Q1";
		BarcodeMatcher barcodeMatcher = new BarcodeMatcher();
		barcodeMatcher.addReferenceSet(barcodeSet, label);
		fail("Illegal label was allowed");
	}
	
	@Test
	public void duplicateBarcode(){
		thrown.expect(IllegalArgumentException.class);
		
		String barcodeSet = "ATCGATT:CAGTCAA:GCTAGCC:TGACTGG";
		String label1 = "Q1";
		BarcodeMatcher barcodeMatcher = new BarcodeMatcher();
		barcodeMatcher.addReferenceSet(barcodeSet, label1);
		String barcodeSet2 = "ATCGATT"; // matches first barcode from above set
		String label2 = "duplicate";
		barcodeMatcher.addReferenceSet(barcodeSet2, label2);
		
		fail("Barcode was allowed twice");
	}
}
