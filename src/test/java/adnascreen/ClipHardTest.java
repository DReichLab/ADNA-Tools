package adnascreen;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TemporaryFolder;

import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

public class ClipHardTest {
	@Rule
	public TemporaryFolder testFolder = new TemporaryFolder();
	
	// Test whether cigar bases are consumed as expected in the reference
	@Test
	public void testCigarConsumes() {
		CigarOperator x;
		char [] queryConsuming = {'M', 'I', 'S', '=', 'X'};
		for (char c : queryConsuming) {
			x = CigarOperator.characterToEnum(c);
			assertTrue(x.consumesReadBases());
		}
		
		char [] queryNonConsuming = {'D', 'N', 'H', 'P'};
		for (char c : queryNonConsuming) {
			x = CigarOperator.characterToEnum(c);
			assertFalse(x.consumesReadBases());
		}
		
		char [] referenceConsuming = {'M', 'D', 'N', '=', 'X'};
		for (char c : referenceConsuming) {
			x = CigarOperator.characterToEnum(c);
			assertTrue(x.consumesReferenceBases());
		}
		
		char [] referenceNonConsuming = {'I', 'S', 'H', 'P'};
		for (char c : referenceNonConsuming) {
			x = CigarOperator.characterToEnum(c);
			assertFalse(x.consumesReferenceBases());
		}
	}
	
	@Test
	public void clipRecord2(){
		int basesToClip = 2;
		/*
		 * NS500217:348:HTW2FBGXY:1:12201:4978:19998	16	MT	12817	37	46M1I12M1D5M	*	0	0	
		 * CGAGCAGATGCCAACACAGCAGCCATTCAAGCAATCCTATAAAACCAATATCGCCGATACGGTT	
		 * EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEEEEEEEEEEEEEEEEEEEEEEEEEEE	
		 * X0:i:1	X1:i:0	MD:Z:41C4G5G5^T5	XD:Z:3-14-Q17.2-Q38.4_64	PG:Z:MarkDuplicates	XG:i:2	NM:i:5	XM:i:3	XO:i:2	XT:A:U
		 */
		String readName = "test";
		String referenceName = "MT";
		int flags = 16;
		int alignmentStart = 12817;
		String readString = "CGAGCAGATGCCAACACAGCAGCCATTCAAGCAATCCTATAAAACCAATATCGCCGATACGGTT";
		String qualityString = "EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEEEEEEEEEEEEEEEEEEEEEEEEEEE";
		
		SAMRecord record = new SAMRecord(null);
		record.setReadName(readName);
		record.setReferenceName(referenceName);
		record.setFlags(flags);
		record.setAlignmentStart(alignmentStart);
		record.setReadString(readString);
		record.setBaseQualityString(qualityString);
		record.setMappingQuality(37);
		record.setCigarString("46M1I12M1D5M");
		record.setAttribute("MD", "41C4G5G5^T5");
		record.setAttribute("NM", 5);
		
		Clipping.hardClipBothEndsOfRead(record, 2);
		String shortenedReadString = readString.substring(2, readString.length() - 2);
		String shortenedQualityString = qualityString.substring(2, qualityString.length() - 2);
		assertEquals(readName, record.getReadName());
		assertEquals(referenceName, record.getReferenceName());
		assertEquals(flags, record.getFlags());
		assertEquals(alignmentStart + basesToClip, record.getAlignmentStart());
		assertEquals(shortenedReadString, record.getReadString());
		assertEquals(shortenedQualityString, record.getBaseQualityString());
		assertEquals("2H44M1I12M1D3M2H", record.getCigarString());
		assertEquals("39C4G5G5^T3", record.getAttribute("MD"));
		assertEquals(5, record.getAttribute("NM"));
		
		record.validateCigar(-1);
		assertNull(record.isValid());
	}
	
	@Test
	public void deletion() {
		/*
		NS500217:488:HWCHLBGX3:3:22409:9911:4321        16      19      38504979        37      1M2D2M4I77M     *       0       0       
		CTGATAAAAGCGTGCAGCACCACCAAATCTGCCACTCATGAGACCAGCACAGAGGAATTACAGAAATGTTGAAGTCAGGACATC        
		AAAAAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAAAAA    
		X0:i:1      X1:i:0  MD:Z:1^CG79     XD:Z:__84       XG:i:2  NM:i:6  XM:i:3  XO:i:1  XT:A:U
		*/
		String readName = "test";
		String referenceName = "19";
		int flags = 16;
		int alignmentStart = 38504979;
		String readString = "CTGATAAAAGCGTGCAGCACCACCAAATCTGCCACTCATGAGACCAGCACAGAGGAATTACAGAAATGTTGAAGTCAGGACATC";
		String qualityString = "AAAAAEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAAAAA";
		
		int readLength = 84;
		
		SAMRecord record = new SAMRecord(null);
		record.setReadName(readName);
		record.setReferenceName(referenceName);
		record.setFlags(flags);
		record.setAlignmentStart(alignmentStart);
		record.setReadString(readString);
		record.setBaseQualityString(qualityString);
		record.setMappingQuality(37);
		record.setCigarString("1M2D2M4I77M");
		record.setAttribute("MD", "1^CG79");
		record.setAttribute("NM", 6);
		
		SAMRecord copy = new SAMRecord(null);
		copy.setReadName(readName);
		copy.setReferenceName(referenceName);
		copy.setFlags(flags);
		copy.setAlignmentStart(alignmentStart);
		copy.setReadString(readString);
		copy.setBaseQualityString(qualityString);
		copy.setMappingQuality(37);
		copy.setCigarString("1M2D2M4I77M");
		copy.setAttribute("MD", "1^CG79");
		copy.setAttribute("NM", 6);
		
		int basesToClip = 2;
		Clipping.hardClipBothEndsOfRead(record, basesToClip);
		assertEquals(readName, record.getReadName());
		assertEquals(referenceName, record.getReferenceName());
		assertEquals(flags, record.getFlags());
		assertEquals(alignmentStart + basesToClip + 2, record.getAlignmentStart()); // two bases clipped plus two deletions
		assertEquals(readString.substring(2, readString.length()-2), record.getReadString());
		assertEquals(qualityString.substring(2, qualityString.length()-2), record.getBaseQualityString());
		assertEquals("4H1M4I75M2H", record.getCigarString());
		assertEquals("76", record.getAttribute("MD"));
		assertEquals(4, record.getAttribute("NM"));
		
		record.validateCigar(-1);
		assertNull(record.isValid());
		
		// check that bases and qualities are still the same at the same reference position
		int referencePosition = record.getAlignmentStart();
		int copyPosition = copy.getReadPositionAtReferencePosition(referencePosition);
		for (int i = 0; i < record.getReadLength(); i++) {
			assertEquals(record.getReadBases()[i], copy.getReadBases()[copyPosition-1+i]);
			assertEquals(record.getBaseQualities()[i], copy.getBaseQualities()[copyPosition-1+i]);
		}
	}
}
