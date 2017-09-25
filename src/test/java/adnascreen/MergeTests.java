package adnascreen;

import static org.junit.Assert.*;

import java.io.IOException;
import java.util.List;

import org.junit.Test;

public class MergeTests {
	static final int maxPenalty = 3;
	static final int minOverlapLength = 10;
	static final int maxPositions = 4;
	static final int minResultLength = 30;
	
	static final int mismatchPenaltyHigh = 3;
	static final int mismatchPenaltyLow = 1;
	static final int mismatchBaseQualityThreshold = 20;
	
	
	@Test
	public void penaltyAndSmallOverlap(){
		// overlap is 10 A bases
		Read r1 = new Read("", "ACCTGATGCGTCAAAAAAAAAA", "EEEEEEEEEEEEEEEEEEEEEE");
		Read r2 = new Read("", "AAAAAAAAAATGTTGGTACCAT", "EEEEEEEEEEEEEEEEEEEEEE");
		
		// overlap starting from A's should be full overlap
		boolean expected = true;
		boolean result = Read.alignmentAssessment(r1, r2, 12, 0, 5, 10, 
				mismatchPenaltyHigh, mismatchPenaltyLow, mismatchBaseQualityThreshold);
		assertEquals(expected, result);
		result = Read.alignmentAssessment(r1, r2, 12, 0, 10, 3,
				mismatchPenaltyHigh, mismatchPenaltyLow, mismatchBaseQualityThreshold);
		assertEquals(expected, result);
		
		expected = false;
		result = Read.alignmentAssessment(r1, r2, 4, 0, 5, 10,
				mismatchPenaltyHigh, mismatchPenaltyLow, mismatchBaseQualityThreshold);
		assertEquals(expected, result);
		
		int thisShortMinResultLength = 15; // this test's reads have length 22, which is too short for standard 30 length
		List<Integer> alignments = Read.findBestAlignment(r1, r1, maxPenalty, minOverlapLength, thisShortMinResultLength, maxPositions,
				mismatchPenaltyHigh, mismatchPenaltyLow, mismatchBaseQualityThreshold);
		assertEquals(1, alignments.size());
		int offset = alignments.get(0);
		assertEquals(0, offset);
	}
	
	@Test
	public void duplicateAlignment(){
		// reads minus barcodes
		Read r1 = new Read("@NS500217:348:HTW2FBGXY:1:11101:22352:1064 1:N:0:0",
				"CTAGCATTACTTATATGATATGTCTCCATACCAATTACAATCTCCAAGTGAACGAGATCGGAAGAGCAC",
				"EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE");
		
		List<Integer> alignments = Read.findBestAlignment(r1, r1, maxPenalty, minOverlapLength, minResultLength, maxPositions,
				mismatchPenaltyHigh, mismatchPenaltyLow, mismatchBaseQualityThreshold);
		assertEquals(1, alignments.size());
		int offset = alignments.get(0);
		assertEquals(0, offset);
	}
	
	@Test
	public void testTrimTrailingN(){
		Read r1 = new Read("@NS500217:348:HTW2FBGXY:1:11101:22352:1064 1:N:0:0",
				"CNAGCATTANNnN",
				"EEEEEEEEEEEEE");
		Read trimmedR1 = r1.trimTrailingUnknownBases();
		DNASequence expected = new DNASequence("CNAGCATTA");
		assertEquals(expected, trimmedR1.getDNASequence());
		for(int i = 0; i < trimmedR1.length(); i++){
			assertEquals(r1.getQualitySequence().getQuality(i), trimmedR1.getQualitySequence().getQuality(i));
		}
		Read doubleTrimmedR1 = trimmedR1.trimTrailingUnknownBases();
		assertEquals(expected, doubleTrimmedR1.getDNASequence());
	}
	
	@Test
	public void trimEmpty(){
		Read r1 = new Read("@NS500217:348:HTW2FBGXY:1:11101:22352:1064 1:N:0:0", "", "");
		Read trimmedR1 = r1.trimTrailingUnknownBases();
		assertEquals(0, r1.length());
		assertEquals(0, trimmedR1.length());
	}
	
	@Test
	public void MergedReadToString(){
		Read r1 = new Read("@NS500217:348:HTW2FBGXY:1:11101:22352:1064 1:N:0:0",
				"CTAGCATTACTTATATGATATGTCTCCATACCAATTACAATCTCCAAGTGAACGAGATCGGAAGAGCAC",
				"EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE");
		IndexAndBarcodeKey key = new IndexAndBarcodeKey("1", "2",
				"Q1", "Q2");
		MergedRead p = new MergedRead(r1, key);
		String s = p.toString();
		//System.out.println(s);;
		// TODO
	}
	
	@Test
	public void mergeMatchingBases(){
		char b1 = 'A';
		char b2 = 'A';
		int quality1 = 30;
		int quality2 = 40;
		MergedRead.BaseWithQuality result = MergedRead.mergeBases(50, b1, quality1, b2, quality2);
		assertEquals(b2, result.base);
		assertEquals(quality2, result.quality);
	}
	
	@Test
	public void mergeNonmatchingBases(){
		char b1 = 'A';
		char b2 = 'C';
		int quality1 = 30;
		int quality2 = 40;
		MergedRead.BaseWithQuality result = MergedRead.mergeBases(50, b1, quality1, b2, quality2);
		assertEquals(b2, result.base);
		assertEquals(quality2 - quality1, result.quality);
	}
	
	@Test
	public void MergeSequencesPartialOverlap(){
		try{
			ClassLoader classLoader = getClass().getClassLoader();
			String filename = classLoader.getResource("Barcodes_5-7bp").getPath();
			BarcodeMatcher barcodeMatcher = new BarcodeMatcher(filename, 1);

			BarcodeMatcher indexMatcher = new BarcodeMatcher();
			DNASequence indexReference = new DNASequence("AAAAAAA");
			indexMatcher.addReferenceSet(indexReference.toString(), "1");
			Read indexRead5 = new Read("@NS500217:348:HTW2FBGXY:1:11101:13418:1065 1:N:0:0", "AAAAAAA", "EEEEEEE");
			Read indexRead7 = new Read("@NS500217:348:HTW2FBGXY:1:11101:13418:1065 2:N:0:0", "AAAAAAA", "EEEEEEE");
			/* splitting this single read into two
		@NS500217:348:HTW2FBGXY:1:11101:13418:1065 1:N:0:0
		CAGCTACTCAGGAGGCT GAGGCAGGAGAATTAACTT GAGATCGGAAGAGCACACGTCTGAACTCCAGTC
		A only           | overlap A and B   | B only
		+
		EEEAEEEAAEEEAEEEE EEE/EEEEEEEEAEEEEEE EEEEEAAEEEEEEEEEEEEEEEEEEEAEEEEEE
			 */
			DNASequence expectedDNA = new DNASequence("CAGCTACTCAGGAGGCTGAGGCAGGAGAATTAACTTGAGATCGGAAGAGCACACGTCTGAACTCCAGTC");
			QualitySequence expectedQuality = new QualitySequence("EEEAEEEAAEEEAEEEEEEE/EEEEEEEEAEEEEEEEEEEEAAEEEEEEEEEEEEEEEEEEEAEEEEEE");
			// ACAATNC barcode prepended
			Read r1 = new Read("@NS500217:348:HTW2FBGXY:1:11101:13418:1065 1:N:0:0",
					"ACAATNCCAGCTACTCAGGAGGCTGAGGCAGGAGAATTAACTT",
					"AAAAA#EEEEAEEEAAEEEAEEEEEEE/EEEEEEEEAE!EEE!"); // ! values are tweaked
			/* non-reversed 
			 * no barcode 
			Read r2 = new Read("@NS500217:348:HTW2FBGXY:1:11101:13418:1065 2:N:0:0",
					"GAGGCAGGAGAATTAACTTGAGATCGGAAGAGCACACGTCTGAACTCCAGTC",
					"EEE/EEEEEEEEAEEEEEEEEEEEAAEEEEEEEEEEEEEEEEEEEAEEEEEE");
			*/
			// reversed, barcode GTCTCAA prepended
			Read r2 = new Read("@NS500217:348:HTW2FBGXY:1:11101:13418:1065 2:N:0:0",
					"GTCTCAAGACTGGAGTTCAGACGTGTGCTCTTCCGATCTCAAGTTAATTCTCCTGCCTC",
					"BBBBBBBEEEEEEAEEEEEEEEEEEEEEEEEEEAAEEEEEEEEEEEAEEEEEEEE/EEE");
			
			// Test search with unrestricted barcode length
			IndexAndBarcodeKey key = MergedRead.findExperimentKey(r1, r2, indexRead5, indexRead7, 
					indexMatcher, indexMatcher, barcodeMatcher, -1);
			// Test that restricting to correct barcode length works too
			IndexAndBarcodeKey key7 = MergedRead.findExperimentKey(r1, r2, indexRead5, indexRead7, 
					indexMatcher, indexMatcher, barcodeMatcher, 7);
			assertEquals(key, key7);
			
			int r1BarcodeLength = barcodeMatcher.getBarcodeLength(key.getP5Label());
			int r2BarcodeLength = barcodeMatcher.getBarcodeLength(key.getP7Label());
			assertEquals(r1BarcodeLength, 7);
			assertEquals(r2BarcodeLength, 7);
			MergedRead merged = MergedRead.mergePairedSequences(r1, r2, key, 
					r1BarcodeLength, r2BarcodeLength, maxPenalty, minOverlapLength, minResultLength,
					mismatchPenaltyHigh, mismatchPenaltyLow, mismatchBaseQualityThreshold);
			
			assertEquals(expectedDNA, merged.getDNASequence());
			assertEquals(expectedQuality, merged.getQualitySequence());
		} catch(IOException e){
			fail();
		}
	}
	
	// This test does not use adapters and barcodes where the simulated data overruns
	@Test
	public void overlappingMerge(){
		try{
			ClassLoader classLoader = getClass().getClassLoader();
			String filename = classLoader.getResource("Barcodes_5-7bp").getPath();
			BarcodeMatcher barcodeMatcher = new BarcodeMatcher(filename, 1);

			BarcodeMatcher indexMatcher = new BarcodeMatcher();
			DNASequence indexReference = new DNASequence("AAAAAAA");
			indexMatcher.addReferenceSet(indexReference.toString(), "1");
			Read indexRead5 = new Read("@NS500217:348:HTW2FBGXY:1:11101:13418:1065 1:N:0:0", "AAAAAAA", "EEEEEEE");
			Read indexRead7 = new Read("@NS500217:348:HTW2FBGXY:1:11101:13418:1065 2:N:0:0", "AAAAAAA", "EEEEEEE");
			/* splitting this single read into two
		@NS500217:348:HTW2FBGXY:1:11101:13418:1065 1:N:0:0
		CAGCTACTC AGGAGGCTGAGGCAGGAGAATTAACTTGAGATCGGAAGAGCACAC GTCTGAACTCCAGTC
		B only   |         overlap A and B                     | A only
		+
		EEEAEEEAA EEEAEEEEEEE/EEEEEEEEAEEEEEEEEEEEAAEEEEEEEEEEE EEEEEEEEAEEEEEE
			 */
			DNASequence expectedDNA = new DNASequence("AGGAGGCTGAGGCAGGAGAATTAACTTGAGATCGGAAGAGCACAC");
			QualitySequence expectedQuality = new QualitySequence("EEEAEEEEEEE/EEEEEEEEAEEEEEEEEEEEAAEEEEEEEEEEE");
			// ACAATNC barcode prepended
			Read r1 = new Read("@NS500217:348:HTW2FBGXY:1:11101:13418:1065 1:N:0:0",
					"ACAATNCAGGAGGCTGAGGCAGGAGAATTAACTTGAGATCGGAAGAGCACACGTCTGAACTCCAGTC",
					"AAAAA#EEEEAEEEEEEE/EEEEEEEEAEEEEEEEEEEEAAEEEEEEEEEEEEEEEEEEEAEEEEEE");
			/* non-reversed 
			 * no barcode
			Read r2 = new Read("@NS500217:348:HTW2FBGXY:1:11101:13418:1065 2:N:0:0",
					"CAGCTACTCAGGAGGCTGAGGCAGGAGAATTAACTTGAGATCGGAAGAGCACAC",
					"EEEAEEEAAEEEAEEEEEEE/EEEEEEEEAEEEEEEEEEEEAAEEEEEEEEEEE");
			*/
			// reversed, barcode GTCTCAA prepended
			Read r2 = new Read("@NS500217:348:HTW2FBGXY:1:11101:13418:1065 2:N:0:0",
					"GTCTCAAGTGTGCTCTTCCGATCTCAAGTTAATTCTCCTGCCTCAGCCTCCTGAGTAGCTG",
					"BBBBBBBEEEEEEEEEEEAAEEEEEEEEEEEAEEEEEEEE/EEEEEEEAEEEAAEEEAEEE");
			
			IndexAndBarcodeKey key = MergedRead.findExperimentKey(r1, r2, indexRead7, indexRead5, 
					indexMatcher, indexMatcher, barcodeMatcher, -1);
			int r1BarcodeLength = barcodeMatcher.getBarcodeLength(key.getP5Label());
			int r2BarcodeLength = barcodeMatcher.getBarcodeLength(key.getP7Label());
			assertEquals(r1BarcodeLength, 7);
			assertEquals(r2BarcodeLength, 7);
			MergedRead merged = MergedRead.mergePairedSequences(r1, r2, key, 
					r1BarcodeLength, r2BarcodeLength, maxPenalty, minOverlapLength, minResultLength,
					mismatchPenaltyHigh, mismatchPenaltyLow, mismatchBaseQualityThreshold);
			
			assertEquals(expectedDNA, merged.getDNASequence());
			assertEquals(expectedQuality, merged.getQualitySequence());
		} catch(IOException e){
			fail();
		}
	}
	
	@Test
	public void ambiguousMerge(){
		StringBuilder qualityBuilder = new StringBuilder();
		for(int i = 0; i < 60; i++){
			qualityBuilder.append('E');
		}
		String qualityString = qualityBuilder.toString();
		
		Read r1 = new Read("", "ACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTG", qualityString);
		Read r2 = new Read("", "CTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGACTGA", qualityString);
		
		List<Integer> alignments = Read.findBestAlignment(r1, r2, maxPenalty, minOverlapLength, minResultLength, maxPositions,
				mismatchPenaltyHigh, mismatchPenaltyLow, mismatchBaseQualityThreshold);
		assertTrue(alignments.size() > 1);
	}
	
	@Test
	public void noBarcodeMerge(){
		try{
			ClassLoader classLoader = getClass().getClassLoader();
			String barcodeFilename = classLoader.getResource("Barcodes_5-7bp").getPath();
			BarcodeMatcher barcodeMatcher = new BarcodeMatcher(barcodeFilename, 1);
			
			String i5Filename = classLoader.getResource("P5indices_Reich20170725").getPath();
			BarcodeMatcher i5Matcher = new BarcodeMatcher(i5Filename, 1);
			String i7Filename = classLoader.getResource("P7indices_Reich20170725").getPath();
			BarcodeMatcher i7Matcher = new BarcodeMatcher(i7Filename, 1);
			Read i1Read = new Read("@NS500217:33:H0NW5AGXX:2:11108:10995:8106 1:N:0:0", "GCTCCGT", "AAAAAFF");
			Read i2Read = new Read("@NS500217:33:H0NW5AGXX:2:11108:10995:8106 2:N:0:0", "ACGGAGC", "<AAAAFF");
			
			Read r1 = new Read("@NS500217:33:H0NW5AGXX:2:11108:10995:8106 1:N:0:0",
					"ATAGAAGCAGAGAATAGTATGATGGTTACTAGAGATCGGAAGAGCACACGTCTGAACTCCAGTCACGCTCCGTATC",
					"AA<AAFFFFFFFFFFAFFFFFFFFFFFFFFFFFFFFFFFFFFFFFAFAFFFFFFAFFFFFFFFF7FFFFFFFFAFF");
			Read r2 = new Read("@NS500217:33:H0NW5AGXX:2:11108:10995:8106 2:N:0:0",
					"CTAGTAACCATCATACTATTCTCTGCTTCTATAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTACGGAGCGTGT",
					"AAAAAF<FFFFFFFFFFFFFFFFFFFFFFFAF<FFFFFFFFFFFFFFFFFF.)FFF.AFFFFFFFF<FF.FAA.FA");
/*
                                            ATAGAAGCAGAGAATAGTATGATGGTTACTAGAGATCGGAAGAGCACACGTCTGAACTCCAGTCACGCTCCGTATC
ACACGCTCCGTACACTCTTTCCCTACACGACGCTCTTCCGATCTATAGAAGCAGAGAATAGTATGATGGTTACTAG
                                            AA<AAFFFFFFFFFFAFFFFFFFFFFFFFFFFFFFFFFFFFFFFFAFAFFFFFFAFFFFFFFFF7FFFFFFFFAFF
AF.AAF.FF<FFFFFFFFA.FFF).FFFFFFFFFFFFFFFFFF<FAFFFFFFFFFFFFFFFFFFFFFFF<FAAAAA
                                            FAFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
*/
			// Test search 0 barcode length
			IndexAndBarcodeKey key = MergedRead.findExperimentKey(r1, r2, i1Read, i2Read, 
					i5Matcher, i7Matcher, barcodeMatcher, 0);
			assertEquals(key.getI5Label(), "CR-P5-i01");
			assertEquals(key.getI7Label(), "CR-P7-i01");
			
			MergedRead merged = MergedRead.mergePairedSequences(r1, r2, key, 
					0, 0, maxPenalty, minOverlapLength, minResultLength,
					mismatchPenaltyHigh, mismatchPenaltyLow, mismatchBaseQualityThreshold);
			assertEquals("ATAGAAGCAGAGAATAGTATGATGGTTACTAG", merged.getDNASequence().toString());
			assertEquals("FAFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF", merged.getQualitySequence().toString());
		}
		catch(Exception e){
			fail();
		}
	}
}
