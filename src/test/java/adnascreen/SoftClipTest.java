package adnascreen;

import static org.junit.Assert.*;

import java.io.BufferedInputStream;
import java.io.FileInputStream;

import org.junit.Test;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CigarUtil;

public class SoftClipTest {
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
}
