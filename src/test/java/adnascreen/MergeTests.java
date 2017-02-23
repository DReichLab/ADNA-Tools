package adnascreen;

import static org.junit.Assert.*;

import java.io.IOException;
import java.util.List;

import org.junit.Test;

public class MergeTests {
	static int maxPenalty = 3;
	static int minOverlapLength = 10;
	static int maxPositions = 4;
	static int minResultLength = 30;
	
	@Test
	public void penaltyAndSmallOverlap(){
		// overlap is 10 A bases
		Read r1 = new Read("", "ACCTGATGCGTCAAAAAAAAAA", "EEEEEEEEEEEEEEEEEEEEEE");
		Read r2 = new Read("", "AAAAAAAAAATGTTGGTACCAT", "EEEEEEEEEEEEEEEEEEEEEE");
		
		// overlap starting from A's should be full overlap
		boolean expected = true;
		boolean result = Read.alignmentAssessment(r1, r2, 12, 0, 5, 10);
		assertEquals(expected, result);
		result = Read.alignmentAssessment(r1, r2, 12, 0, 10, 3);
		assertEquals(expected, result);
		
		expected = false;
		result = Read.alignmentAssessment(r1, r2, 4, 0, 5, 10);
		assertEquals(expected, result);
		
		int thisShortMinResultLength = 15; // this test's reads have length 22, which is too short for standard 30 length
		List<Integer> alignments = Read.findBestAlignment(r1, r1, maxPenalty, minOverlapLength, thisShortMinResultLength, maxPositions);
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
		
		List<Integer> alignments = Read.findBestAlignment(r1, r1, maxPenalty, minOverlapLength, minResultLength, maxPositions);
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
	public void MergedReadToString(){
		Read r1 = new Read("@NS500217:348:HTW2FBGXY:1:11101:22352:1064 1:N:0:0",
				"CTAGCATTACTTATATGATATGTCTCCATACCAATTACAATCTCCAAGTGAACGAGATCGGAAGAGCAC",
				"EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE");
		IndexAndBarcodeKey key = new IndexAndBarcodeKey(new DNASequence("AAAAAAA"), new DNASequence("CCCCCCC"),
				"TGACGCA:ATCGTGC:CAGTATG:GCTACAT", "GTCTCAA:TAGAGCC:ACTCTGG:CGAGATT");
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
			String filename = classLoader.getResource("barcodes").getPath();
			BarcodeMatcher barcodeMatcher = new BarcodeMatcher(filename, 1);

			IndexMatcher indexMatcher = new IndexMatcher();
			DNASequence indexReference = new DNASequence("AAAAAAA");
			indexMatcher.addReference(indexReference);
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
			
			MergedRead merged = MergedRead.mergePairedSequences(r1, r2, indexRead5, indexRead7, 
					indexMatcher, indexMatcher, barcodeMatcher, maxPenalty, minOverlapLength, minResultLength);
			
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
			String filename = classLoader.getResource("barcodes").getPath();
			BarcodeMatcher barcodeMatcher = new BarcodeMatcher(filename, 1);

			IndexMatcher indexMatcher = new IndexMatcher();
			DNASequence indexReference = new DNASequence("AAAAAAA");
			indexMatcher.addReference(indexReference);
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
			
			MergedRead merged = MergedRead.mergePairedSequences(r1, r2, indexRead5, indexRead7, 
					indexMatcher, indexMatcher, barcodeMatcher, maxPenalty, minOverlapLength, minResultLength);
			
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
		
		List<Integer> alignments = Read.findBestAlignment(r1, r2, maxPenalty, minOverlapLength, minResultLength, maxPositions);
		assertTrue(alignments.size() > 1);
	}
}
