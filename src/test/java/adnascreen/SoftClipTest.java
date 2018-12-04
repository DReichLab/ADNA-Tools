package adnascreen;

import static org.junit.Assert.*;

import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.commons.cli.ParseException;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TemporaryFolder;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.TextCigarCodec;
import htsjdk.samtools.util.CigarUtil;

public class SoftClipTest {
	@Rule
	public TemporaryFolder testFolder = new TemporaryFolder();
	
	@Test
	public void fileTest(){
		ClassLoader classLoader = getClass().getClassLoader();
		String filename = classLoader.getResource("target-test.sam").getPath();
		int numSoftClipBases = 2;
		try {
			SamInputResource bufferedSAMFile = SamInputResource.of(new BufferedInputStream(new FileInputStream(filename)));
			SamReader reader = SamReaderFactory.makeDefault().open(bufferedSAMFile);
			SAMRecordIterator i = reader.iterator();
			while(i.hasNext()){
				SAMRecord record = i.next();
				int alignmentStart = record.getAlignmentStart();

				Cigar unmodified = record.getCigar();
				assertNull(unmodified.isValid(null, -1));
				assertFalse(record.getCigar().isClipped());
				assertFalse(record.getCigar().isRightClipped());
				assertFalse(record.getCigar().isLeftClipped());

				SoftClip.softClipBothEndsOfRead(record, numSoftClipBases);

				//System.out.println(unmodified);
				//System.out.println(record.getCigar());
				assertNull(record.getCigar().isValid(null, -1));
				assertTrue(record.getCigar().isClipped());
				assertTrue(record.getCigar().isLeftClipped());
				assertTrue(record.getCigar().isRightClipped());
				
				CigarElement firstCigarElement = record.getCigar().getFirstCigarElement();
				assertEquals(CigarOperator.SOFT_CLIP, firstCigarElement.getOperator());
				assertTrue(firstCigarElement.getLength() >= numSoftClipBases);
				CigarElement lastCigarElement = record.getCigar().getFirstCigarElement();
				assertEquals(CigarOperator.SOFT_CLIP, lastCigarElement.getOperator());
				assertTrue(lastCigarElement.getLength() >= numSoftClipBases);

				assertEquals(unmodified.getReadLength(), record.getCigar().getReadLength());
				// all of the examples in this file have 2 cigar operators that consume reference bases
				// at the start, so they will adjust the alignment start
				assertEquals(alignmentStart + numSoftClipBases, record.getAlignmentStart());		
			}
		} catch (Exception e) {
			fail();
		}
	}

	@Test
	public void cigarOperatorBehavior(){
		CigarOperator [] cigarOperators = CigarOperator.values();
		for(CigarOperator cigarOperator : cigarOperators){
			char c = (char) CigarOperator.enumToCharacter(cigarOperator);
			CigarOperator current = CigarOperator.characterToEnum(c);
			assertEquals(current, cigarOperator);
		}
	}
	
	// make sure the regex matches valid sequences
	@Test
	public void testMD_validation(){
		Pattern md_regex = SAM_MD.VALIDATION_REGEX;
		Matcher valid_md;
		
		valid_md = md_regex.matcher("70");
		assertTrue(valid_md.matches());
		valid_md = md_regex.matcher("6^G106");
		assertTrue(valid_md.matches());
		valid_md = md_regex.matcher("0C0G84");
		assertTrue(valid_md.matches());
		valid_md = md_regex.matcher("0G108");
		assertTrue(valid_md.matches());
		valid_md = md_regex.matcher("11G10G3G5G36G2G42");
		assertTrue(valid_md.matches());
		valid_md = md_regex.matcher("48G15G29A5G3G0");
		assertTrue(valid_md.matches());
		valid_md = md_regex.matcher("3C55G0");
		assertTrue(valid_md.matches());
	}
	
	@Test
	public void testMD_70(){
		String md = "70";
		SAM_MD test = new SAM_MD(md);
		assertEquals(md, test.toString());
		assertEquals(0, test.editDistance());
	}
	@Test
	public void testMD_6_G106(){
		String md = "6^G106";
		SAM_MD test = new SAM_MD(md);
		assertEquals(md, test.toString());
		assertEquals(1, test.editDistance());
	}
	@Test
	public void testMD_0C0G84(){
		String md = "0C0G84";
		SAM_MD test = new SAM_MD(md);
		assertEquals(md, test.toString());
		assertEquals(2, test.editDistance());
	}
	@Test
	public void testMD_0G108(){
		String md = "0G108";
		SAM_MD test = new SAM_MD(md);
		assertEquals(md, test.toString());
		assertEquals(1, test.editDistance());
	}
	@Test
	public void testMD_11G10G3G5G36G2G42(){
		String md = "11G10G3G5G36G2G42";
		SAM_MD test = new SAM_MD(md);
		assertEquals(md, test.toString());
		assertEquals(6, test.editDistance());
	}
	@Test
	public void testMD_48G15G29A5G3G0(){
		String md = "48G15G29A5G3G0";
		SAM_MD test = new SAM_MD(md);
		assertEquals(md, test.toString());
		assertEquals(5, test.editDistance());
	}
	
	@Test
	public void testMD_3C55G0(){
		String md = "3C55G0";
		SAM_MD test = new SAM_MD(md);
		assertEquals(md, test.toString());
		assertEquals(2, test.editDistance());
	}
	
	@Test
	public void testMD_42_C26(){
		String md = "42^C26";
		SAM_MD test = new SAM_MD(md);
		assertEquals(md, test.toString());
		assertEquals(1, test.editDistance());
	}
	
	@Test
	public void testMD_19A15G6_T35(){
		String md = "19A15G6^T35";
		SAM_MD test = new SAM_MD(md);
		assertEquals(md, test.toString());
		assertEquals(3, test.editDistance());
	}
	
	@Test
	public void testMD_0C53_CG55(){
		String md = "0C53^CG55";
		SAM_MD test = new SAM_MD(md);
		assertEquals(md, test.toString());
		assertEquals(3, test.editDistance());
	}
	
	@Test
	public void testMD_50_clip(){
		String md = "50";
		String expected = "46";
		SAM_MD test = new SAM_MD(md);
		test.clip(2, 2);
		assertEquals(expected, test.toString());
	}
	
	@Test
	public void testMD_3C55G0_clip2(){
		String md = "3C55G0";
		String expected = "1C54";
		SAM_MD test = new SAM_MD(md);
		test.clip(2, 2);
		assertEquals(expected, test.toString());
	}
	
	@Test
	public void testMD_3C55G0_clip3(){
		String md = "3C55G0";
		String expected = "0C53";
		SAM_MD test = new SAM_MD(md);
		test.clip(3, 3);
		assertEquals(expected, test.toString());
	}
	
	@Test
	public void testMD_3C55G0_clip4(){
		String md = "3C55G0";
		String expected = "52";
		SAM_MD test = new SAM_MD(md);
		test.clip(4, 4);
		assertEquals(expected, test.toString());
	}
	
	@Test
	public void cigarInsertLength(){
		Cigar c = TextCigarCodec.decode("52M2I5M");
		assertEquals(2, SoftClip.editDistance(c));
	}
	
	@Test
	public void cigarInsertLength2(){
		Cigar c = TextCigarCodec.decode("52M1D5M");
		assertEquals(0, SoftClip.editDistance(c));
	}
	
	@Test
	public void clipRecord(){
		int basesToClip = 2;
		/*NS500217:348:HTW2FBGXY:2:23107:4108:10621	16	MT	5510	37	20M1D31M	*	0	0	
		 * CTAGAAATTTAGGTTAAACAAGACCAAGAGCCTTCAAAGCCCTTAGTAAGT	
		 * EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE	
		 * X0:i:1	X1:i:0	MD:Z:0A17T1^C23C7	XD:Z:3-14-Q17.4-Q38.2_51	PG:Z:MarkDuplicates	XG:i:1	NM:i:4	XM:i:3	XO:i:1	XT:A:U
		 */
		String readName = "test";
		String referenceName = "MT";
		int flags = 16;
		int alignmentStart = 5510;
		String readString = "CTAGAAATTTAGGTTAAACAAGACCAAGAGCCTTCAAAGCCCTTAGTAAGT";
		String qualityString = "EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE";
		
		SAMRecord record = new SAMRecord(null);
		record.setReadName(readName);
		record.setReferenceName(referenceName);
		record.setFlags(flags);
		record.setAlignmentStart(alignmentStart);
		record.setReadString(readString);
		record.setBaseQualityString(qualityString);
		record.setMappingQuality(37);
		record.setCigarString("20M1D31M");
		record.setAttribute("MD", "0A17T1^C23C7");
		record.setAttribute("NM", 4);
		
		SoftClip.softClipBothEndsOfRead(record, 2);
		assertEquals(readName, record.getReadName());
		assertEquals(referenceName, record.getReferenceName());
		assertEquals(flags, record.getFlags());
		assertEquals(alignmentStart + basesToClip, record.getAlignmentStart());
		assertEquals(readString, record.getReadString());
		assertEquals(qualityString, record.getBaseQualityString());
		assertEquals("2S18M1D29M2S", record.getCigarString());
		assertEquals("16T1^C23C5", record.getAttribute("MD"));
		assertEquals(3, record.getAttribute("NM"));
		
		record.validateCigar(-1);
		assertNull(record.isValid());
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
		
		SoftClip.softClipBothEndsOfRead(record, 2);
		assertEquals(readName, record.getReadName());
		assertEquals(referenceName, record.getReferenceName());
		assertEquals(flags, record.getFlags());
		assertEquals(alignmentStart + basesToClip, record.getAlignmentStart());
		assertEquals(readString, record.getReadString());
		assertEquals(qualityString, record.getBaseQualityString());
		assertEquals("2S44M1I12M1D3M2S", record.getCigarString());
		assertEquals("39C4G5G5^T3", record.getAttribute("MD"));
		assertEquals(5, record.getAttribute("NM"));
		
		record.validateCigar(-1);
		assertNull(record.isValid());
	}
	
	
	@Test
	public void clipComparisonWithOriginalReadByPosition(){
		int basesToClip = 2;
		/*NS500217:348:HTW2FBGXY:2:23107:4108:10621	16	MT	5510	37	20M1D31M	*	0	0	
		 * CTAGAAATTTAGGTTAAACAAGACCAAGAGCCTTCAAAGCCCTTAGTAAGT	
		 * EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE	
		 * X0:i:1	X1:i:0	MD:Z:0A17T1^C23C7	XD:Z:3-14-Q17.4-Q38.2_51	PG:Z:MarkDuplicates	XG:i:1	NM:i:4	XM:i:3	XO:i:1	XT:A:U
		 */
		String readName = "test";
		String referenceName = "MT";
		int flags = 16;
		int alignmentStart = 5510;
		String readString = "CTAGAAATTTAGGTTAAACAAGACCAAGAGCCTTCAAAGCCCTTAGTAAGT";
		String qualityString = "EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE";
		
		int readLength = 51;
		
		SAMRecord record = new SAMRecord(null);
		record.setReadName(readName);
		record.setReferenceName(referenceName);
		record.setFlags(flags);
		record.setAlignmentStart(alignmentStart);
		record.setReadString(readString);
		record.setBaseQualityString(qualityString);
		record.setMappingQuality(37);
		record.setCigarString("20M1D31M");
		record.setAttribute("MD", "0A17T1^C23C7");
		record.setAttribute("NM", 4);
		
		assertEquals(readLength, record.getReadLength());
		
		SAMRecord original = new SAMRecord(null);
		original.setReadName(readName);
		original.setReferenceName(referenceName);
		original.setFlags(flags);
		original.setAlignmentStart(alignmentStart);
		original.setReadString(readString);
		original.setBaseQualityString(qualityString);
		original.setMappingQuality(37);
		original.setCigarString("20M1D31M");
		original.setAttribute("MD", "0A17T1^C23C7");
		original.setAttribute("NM", 4);
		
		SoftClip.softClipBothEndsOfRead(record, basesToClip);
		assertEquals(readName, record.getReadName());
		assertEquals(referenceName, record.getReferenceName());
		assertEquals(flags, record.getFlags());
		assertEquals(alignmentStart + basesToClip, record.getAlignmentStart());
		assertEquals(readString, record.getReadString());
		assertEquals(qualityString, record.getBaseQualityString());
		assertEquals("2S18M1D29M2S", record.getCigarString());
		assertEquals("16T1^C23C5", record.getAttribute("MD"));
		assertEquals(3, record.getAttribute("NM"));
		//assertEquals(readLength - 2 * basesToClip, record.getReadLength());
		//System.out.println(record.getAlignmentEnd());
		
		record.validateCigar(-1);
		assertNull(record.isValid());
		
		//System.out.println(record.getSAMString());
		
		for (int i = record.getAlignmentStart(); i < record.getAlignmentEnd(); i++){
			assertEquals(original.getReadPositionAtReferencePosition(i), record.getReadPositionAtReferencePosition(i));
		}
	}
	
	@Test
	public void sample(){
		/* NS500217:348:HTW2FBGXY:1:21302:14243:3111	4	GL000248.1	39786	0	7M4I42M	*	0	0	
		 * CGATCTGTATCGGCACCTTGGCCTCCCAAAGTGCTGGGATTACAGGCATGAGC	
		 * EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE	
		 * X0:i:2	X1:i:0	XA:Z:16,+75842244,8M4I41M,4;	MD:Z:7C41	XD:Z:48-91-Q11.2-Q43.4_53	XG:i:4	NM:i:5	XM:i:0	XO:i:2	XT:A:R
		 * */
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
		SoftClip.softClipBothEndsOfRead(record, basesToClip);
		assertEquals(readName, record.getReadName());
		assertEquals(referenceName, record.getReferenceName());
		assertEquals(flags, record.getFlags());
		assertEquals(alignmentStart + basesToClip + 2, record.getAlignmentStart()); // two bases clipped plus two deletions
		assertEquals(readString, record.getReadString());
		assertEquals(qualityString, record.getBaseQualityString());
		//assertEquals("1S2H1S1M4I75M2S", record.getCigarString()); // expected, but illegal because hard clip after soft clip
		assertEquals("2H2S1M4I75M2S", record.getCigarString());
		assertEquals("76", record.getAttribute("MD"));
		assertEquals(4, record.getAttribute("NM"));
		
		record.validateCigar(-1);
		assertNull(record.isValid());
		
		for (int i = record.getAlignmentStart(); i < record.getAlignmentEnd(); i++){
			assertEquals(copy.getReadPositionAtReferencePosition(i), record.getReadPositionAtReferencePosition(i));
		}
	}
	
	// For very short reads, clipping may result in read disappearing. Handle this edge case.
	@Test
	public void shortReadEntirelyClipped() {
		ClassLoader classLoader = getClass().getClassLoader();
		String samToClip = classLoader.getResource("shortReadForClipping.sam").getPath();
		final int basesToClip = 20; // Clip more bases than are in a read
		try {
			// clip test file
			File clippledSam = testFolder.newFile("clipped.sam");
			String clippedSamFilename = clippledSam.getAbsolutePath();
			String [] args = {"-n", String.valueOf(basesToClip), "-i", samToClip, "-o", clippedSamFilename};
			SoftClip.main(args);
			
			SamInputResource bufferedSAMFile;
			SamReader reader;
			SAMRecordIterator i;
			// count reads in input file
			int inputCount = 0;
			bufferedSAMFile = SamInputResource.of(new BufferedInputStream(new FileInputStream(samToClip)));
			reader = SamReaderFactory.makeDefault().open(bufferedSAMFile);
			i = reader.iterator();
			while(i.hasNext()){
				SAMRecord record = i.next();
				assertNull(record.isValid()); // check that input reads
				Cigar cigar = record.getCigar();
				assertTrue(cigar.containsOperator(CigarOperator.MATCH_OR_MISMATCH)); // this is our definition of nonempty
				inputCount++;
			}
			
			// verify output has no empty reads
			int outputCount = 0;
			bufferedSAMFile = SamInputResource.of(new BufferedInputStream(new FileInputStream(clippedSamFilename)));
			reader = SamReaderFactory.makeDefault().open(bufferedSAMFile);
			i = reader.iterator();
			assertTrue(i.hasNext());
			while(i.hasNext()){
				SAMRecord record = i.next();
				assertNull(record.isValid());
				Cigar cigar = record.getCigar();
				assertTrue(cigar.containsOperator(CigarOperator.MATCH_OR_MISMATCH)); // this is our definition of nonempty
				outputCount++;
			}
			
			// empty reads are dropped; there should be fewer output reads than input reads
			assertTrue(outputCount < inputCount);
		} catch (IOException | ParseException e) {
			fail();
		}
	}
	
	@Test
	public void multiLengthClipping() {
		ClassLoader classLoader = getClass().getClassLoader();
		String filename = classLoader.getResource("multi_lib.sam").getPath();
		try {
			File clippedSam = testFolder.newFile("clipped.sam");
			String clippedSamFilename = clippedSam.getAbsolutePath();
			String [] args = {"-n", "2", "-i", filename, "-o", clippedSamFilename, "-x", "10", "-s", "Lib2", "-s", "Lib3", "-y", "-0", "-t", "Lib4"};
			SoftClip.main(args);
			
			// check clipping in new file
			SamInputResource bufferedSAMFile = SamInputResource.of(new BufferedInputStream(new FileInputStream(clippedSamFilename)));
			SamReader reader = SamReaderFactory.makeDefault().open(bufferedSAMFile);
			SAMRecordIterator i = reader.iterator();
			assertTrue(i.hasNext());
			while(i.hasNext()){
				SAMRecord record = i.next();
				SAMReadGroupRecord readGroup = record.getReadGroup();
				assertNull(record.isValid());
				
				int expectedClippedBases = -1;
				String library = readGroup.getLibrary();
				switch(library) {
				case "Lib1":
					expectedClippedBases = 2;
					break;
				case "Lib2":
				case "Lib3":
					expectedClippedBases = 10;
					break;
				case "Lib4":
					expectedClippedBases = 0;
					break;
				default:
					fail(); // all reads should have a read group and library
				}
				Cigar cigar = record.getCigar();
				char [] cigarUnrolled = CigarUtil.cigarArrayFromElements(cigar.getCigarElements() );
				// ends of reads should be soft clipped
				for(int n = 0; n < expectedClippedBases; n++) {
					assertEquals(cigarUnrolled[n], 'S');
					assertEquals(cigarUnrolled[cigarUnrolled.length - (1 + n)], 'S');
				}
				// further in, read should not be soft clipped
				assertNotEquals(cigarUnrolled[expectedClippedBases], 'S');
				assertNotEquals(cigarUnrolled[cigarUnrolled.length - (1 + expectedClippedBases)], 'S');
			}
		} catch (ParseException e) {
			fail();
		} catch (IOException e) {
			fail();
		}
	}
}
