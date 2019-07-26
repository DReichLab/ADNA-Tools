package adnascreen;

import static org.junit.Assert.*;

import java.io.IOException;
import java.util.List;

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
			String filename = classLoader.getResource("Barcodes_5-7bp").getPath();
			
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
	
	@Test
	public void basicLength(){
		String barcodeSet = "ATCGATT:CAGTCAA:GCTAGCC:TGACTGG";
		String label = "Q1";
		final int expectedLength = 7;

		BarcodeMatcher barcodeMatcher = new BarcodeMatcher();
		barcodeMatcher.setMaxHammingDistance(1);

		barcodeMatcher.addReferenceSet(barcodeSet, label);
		
		List<Integer> lengths = barcodeMatcher.getBarcodeLengths();
		assertEquals(1, lengths.size());
		assertEquals(expectedLength, (int) lengths.get(0));
		assertEquals(expectedLength, barcodeMatcher.getBarcodeLength(label));
		
		// check cache too
		lengths = barcodeMatcher.getBarcodeLengths();
		assertEquals(1, lengths.size());
		assertEquals(expectedLength, (int) lengths.get(0));
		assertEquals(expectedLength, barcodeMatcher.getBarcodeLength(label));
	}
	
	@Test
	public void lengthNotPresent(){
		String barcodeSet = "ATCGATT:CAGTCAA:GCTAGCC:TGACTGG";
		String label = "Q1";
		final int expectedLength = 0;
		
		BarcodeMatcher barcodeMatcher = new BarcodeMatcher();
		barcodeMatcher.setMaxHammingDistance(1);

		barcodeMatcher.addReferenceSet(barcodeSet, label);
		int length = barcodeMatcher.getBarcodeLength("not_present");
		assertEquals(expectedLength, length);
		
		int nullLength = barcodeMatcher.getBarcodeLength(null);
		assertEquals(expectedLength, nullLength);
	}
	
	@Test
	public void lengthBarcodePairLabel(){
		try{
			ClassLoader classLoader = getClass().getClassLoader();
			String filename = classLoader.getResource("Barcodes_5-7bp").getPath();

			BarcodeMatcher barcodeMatcher = new BarcodeMatcher(filename, 1);
			
			int length;
			length = barcodeMatcher.getBarcodePairLength("Q1" + IndexAndBarcodeKey.FIELD_SEPARATOR + "Q2");
			assertEquals(7, length);
			
			length = barcodeMatcher.getBarcodePairLength("6-1" + IndexAndBarcodeKey.FIELD_SEPARATOR + "6-1");
			assertEquals(6, length);
			
			length = barcodeMatcher.getBarcodePairLength("5-4" + IndexAndBarcodeKey.FIELD_SEPARATOR + "5-6");
			assertEquals(5, length);
		}
		catch(Exception e){
			fail();
		}
	}
	
	@Test
	public void illegalBarcodePairLabel(){
		thrown.expect(IllegalArgumentException.class);
		try{
			ClassLoader classLoader = getClass().getClassLoader();
			String filename = classLoader.getResource("Barcodes_5-7bp").getPath();

			BarcodeMatcher barcodeMatcher = new BarcodeMatcher(filename, 1);
			
			barcodeMatcher.getBarcodePairLength("Q1" + IndexAndBarcodeKey.FIELD_SEPARATOR + "6-1");
			fail();
		}
		catch(IOException e){
			fail();
		}
	}
	
	@Test
	public void caseChange() {
		String barcodeSet = "atcgatt:CAGTCAA";
		String label = "A1";
		String expectedResult = "A1" + BarcodeMatcher.INDEX_DELIMITER + "1";
		DNASequence query = new DNASequence("ATCGATT");

		BarcodeMatcher barcodeMatcher = new BarcodeMatcher();
		barcodeMatcher.addReferenceSet(barcodeSet, label);
		String result = barcodeMatcher.find(query);
		assertEquals(expectedResult, result);
	}
	
	@Test
	public void labelToBarcodeNullQuery(){
		String barcodeSet = "ATCGATT:CAGTCAA:GCTAGCC:TGACTGG";
		String label = "Q1";
		String query = null;
		DNASequence expectedResult = null;

		BarcodeMatcher barcodeMatcher = new BarcodeMatcher();
		barcodeMatcher.addReferenceSet(barcodeSet, label);
		DNASequence result = barcodeMatcher.getBarcode(query);
		assertEquals(expectedResult, result);
	}
	
	@Test
	public void labelToBarcodeMatch1(){
		String barcodeSet = "ATCGATT:CAGTCAA:GCTAGCC:TGACTGG";
		String label = "Q1";
		String query = "Q1" + BarcodeMatcher.INDEX_DELIMITER + "1";
		DNASequence expectedResult = new DNASequence("ATCGATT");

		BarcodeMatcher barcodeMatcher = new BarcodeMatcher();
		barcodeMatcher.addReferenceSet(barcodeSet, label);
		DNASequence result = barcodeMatcher.getBarcode(query);
		assertEquals(expectedResult, result);
	}
	
	@Test
	public void labelToBarcodeMatch2(){
		String barcodeSet = "ATCGATT:CAGTCAA:GCTAGCC:TGACTGG";
		String label = "Q1";
		String query = "Q1" + BarcodeMatcher.INDEX_DELIMITER + "2";
		DNASequence expectedResult = new DNASequence("CAGTCAA");

		BarcodeMatcher barcodeMatcher = new BarcodeMatcher();
		barcodeMatcher.addReferenceSet(barcodeSet, label);
		DNASequence result = barcodeMatcher.getBarcode(query);
		assertEquals(expectedResult, result);
	}
	
	@Test
	public void labelToBarcodeMatchNonset() {
		String barcodeSet = "ATCGATT";
		String label = "test";
		DNASequence expectedResult = new DNASequence(barcodeSet);

		BarcodeMatcher barcodeMatcher = new BarcodeMatcher();
		barcodeMatcher.addReferenceSet(barcodeSet, label);
		DNASequence result = barcodeMatcher.getBarcode(label);
		assertEquals(expectedResult, result);
	}
	
	@Test
	public void labelToBarcodeNoMatch() {
		String barcodeSet = "ATCGATT:CAGTCAA:GCTAGCC:TGACTGG";
		String label = "Q1";
		String query = "Q1" + BarcodeMatcher.INDEX_DELIMITER + "5";

		BarcodeMatcher barcodeMatcher = new BarcodeMatcher();
		barcodeMatcher.addReferenceSet(barcodeSet, label);
		
		DNASequence result = barcodeMatcher.getBarcode(query);
		assertNull(result);
		
		assertNull(barcodeMatcher.getBarcode("missing_label"));
	}
	
	@Test
	public void badLabel() {
		thrown.expect(IllegalArgumentException.class);
		String barcodeSet = "ATCGATT:CAGTCAA:GCTAGCC:TGACTGG";
		String label = "Q1" + BarcodeMatcher.INDEX_DELIMITER;
		
		BarcodeMatcher barcodeMatcher = new BarcodeMatcher();
		barcodeMatcher.addReferenceSet(barcodeSet, label);
	}
	
	@Test
	public void changingHammingDistance() {
		String barcodeSet = "ATCGATT:CAGTCAA:GCTAGCC:TGACTGG";
		String label = "Q1";
		String expectedResult = "Q1" + BarcodeMatcher.INDEX_DELIMITER + "1";
		DNASequence query = new DNASequence("ATCGACC"); // off by 2

		BarcodeMatcher barcodeMatcher = new BarcodeMatcher();
		barcodeMatcher.setMaxHammingDistance(2);
		barcodeMatcher.addReferenceSet(barcodeSet, label);
		String result = barcodeMatcher.find(query);
		assertEquals(expectedResult, result);
		
		barcodeMatcher.setMaxHammingDistance(1);
		assertNull(barcodeMatcher.find(query));
	}
	
	@Test
	public void badHammingDistance() {
		thrown.expect(IllegalArgumentException.class);
		BarcodeMatcher barcodeMatcher = new BarcodeMatcher();
		barcodeMatcher.setMaxHammingDistance(-1);
	}
	
	@Test
	public void emptyFile() {
		ClassLoader classLoader = getClass().getClassLoader();
		String filename = classLoader.getResource("empty").getPath();
		
		BarcodeMatcher barcodeMatcher = new BarcodeMatcher();
		try {
			barcodeMatcher.loadFile(filename);
			assertEquals(0, barcodeMatcher.getBarcodeLengths().size());
		}
		catch(Exception e) {
			fail();
		}
	}
	
	private DNASequence dnaStringFromInt(int n, int length){
		final char[] BASES = {'A', 'C', 'G', 'T'}; 
		StringBuilder b = new StringBuilder();
		for(int i = 0; i < length; i++){
			char c = BASES[n % BASES.length];
			n /= BASES.length;
			b.append(c);
		}
		String dnaString = b.reverse().toString();
		return new DNASequence(dnaString);
	}
	
	//@Test
	public void barcodeDensity(){
		try{
			ClassLoader classLoader = getClass().getClassLoader();
			String filename = classLoader.getResource("Barcodes_5-7bp").getPath();
			
			SampleCounter whichBarcodesMatch = new SampleCounter();
			
			BarcodeMatcher barcodeMatcher = new BarcodeMatcher(filename, 1);
			for(int n = 0; n < 16384; n++){
				DNASequence query_7bp = dnaStringFromInt(n, 7);
				DNASequence query_6bp = query_7bp.subsequence(0, 6);
				assertEquals(6, query_6bp.length());
				DNASequence query_5bp = query_7bp.subsequence(0, 5);
				assertEquals(5, query_5bp.length());
				
				String result_7bp = barcodeMatcher.find(query_7bp);
				String result_6bp = barcodeMatcher.find(query_6bp);
				String result_5bp = barcodeMatcher.find(query_5bp);
				
				String result = (result_5bp == null ? "0" : "1")
						+ 		(result_6bp == null ? "0" : "1")
						+ 		(result_7bp == null ? "0" : "1");
				whichBarcodesMatch.increment(result);
				if(result.equals("101")){
					System.out.println(n + "\t" + query_7bp + "\t" + result_5bp + "\t" + result_6bp + "\t" + result_7bp);
				}
			}
			System.out.println(whichBarcodesMatch.toString());
		}
		catch(IOException e){
			fail(e.toString());
		}
	}
}
