package adnascreen;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import java.io.BufferedInputStream;
import java.io.FileInputStream;
import java.io.IOException;

import org.apache.commons.cli.ParseException;
import org.junit.Test;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMRecordSetBuilder;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class AlignmentComparisonTests {
	String [] samFilenames = 
			{"end_clipping.sam", "heng_shotgun.sam", "multi_lib.sam", 
			 "shortReadForClipping.sam", "test.sam", "target-test.sam"};
	
	public SAMRecord getRecord(String filename, int zeroBasedIndex) throws IOException {
		ClassLoader classLoader = getClass().getClassLoader();
		String filepath = classLoader.getResource(filename).getPath();
		SamReaderFactory samReaderFactory = SamReaderFactory.makeDefault();
		try(
		SamReader reader = samReaderFactory.open(SamInputResource.of(new BufferedInputStream(new FileInputStream(filepath))));
			){
			SAMRecordIterator outputIterator = reader.iterator();
			for (int i = 0; i < zeroBasedIndex; i++) {
				outputIterator.next();
			}
			return outputIterator.next();
		}
	}
	
	@Test
	public void testRecordReflexive() throws IOException {
		SAMRecord x = getRecord("multi_lib.sam", 0);
		SAMRecord y = x.deepCopy();
		String tags[] = {"MD", "XD", "RG", "NM"};
		assertTrue(AlignmentComparison.compareFullAlignments(x, y, tags));
	}
	
	@Test
	public void testRecordDuplicateFlag() throws IOException {
		SAMRecord x = getRecord("multi_lib.sam", 0);
		SAMRecord y = x.deepCopy();
		String tags[] = {"MD", "XD", "RG", "NM"};
		x.setDuplicateReadFlag(true);
		assertTrue(AlignmentComparison.compareFullAlignments(x, y, tags));
	}
	
	@Test
	public void testRecordVendorQualityFlagFail() throws IOException {
		SAMRecord x = getRecord("multi_lib.sam", 0);
		SAMRecord y = x.deepCopy();
		String tags[] = {"MD", "XD", "RG", "NM"};
		x.setReadFailsVendorQualityCheckFlag(true);
		assertFalse(AlignmentComparison.compareFullAlignments(x, y, tags));
	}
	
	@Test public void testRecordBaseFail() throws IOException {
		SAMRecord x = getRecord("multi_lib.sam", 0);
		SAMRecord y = x.deepCopy();
		String tags[] = {"MD", "XD", "RG", "NM"};
		byte[] bases = x.getReadBases();
		bases[1] = 'T';
		x.setReadBases(bases);
		assertFalse(AlignmentComparison.compareFullAlignments(x, y, tags));
	}
	
	@Test public void testRecordBaseQualityFail() throws IOException {
		SAMRecord x = getRecord("multi_lib.sam", 0);
		SAMRecord y = x.deepCopy();
		String tags[] = {"MD", "XD", "RG", "NM"};
		byte[] qualities = x.getBaseQualities();
		qualities[2] = 'A';
		x.setBaseQualities(qualities);
		assertFalse(AlignmentComparison.compareFullAlignments(x, y, tags));
	}
	
	@Test public void testRecordXDFail() throws IOException {
		SAMRecord x = getRecord("multi_lib.sam", 0);
		SAMRecord y = x.deepCopy();
		String tags[] = {"MD", "XD", "RG", "NM"};
		x.setAttribute("XD", "ATTGGCA_AAGCATC_40");
		assertFalse(AlignmentComparison.compareFullAlignments(x, y, tags));
	}
	
	@Test
	public void testRecordMappingQualityFlag() throws IOException {
		SAMRecord x = getRecord("multi_lib.sam", 0);
		SAMRecord y = x.deepCopy();
		String tags[] = {"MD", "XD", "RG", "NM"};
		x.setMappingQuality(1);
		assertFalse(AlignmentComparison.compareFullAlignments(x, y, tags));
	}
	
	@Test
	public void testRecordSetReflexive() {
		SAMRecordSetBuilder testSet = new SAMRecordSetBuilder();
		SAMRecord test = testSet.addFrag("abc", 1, 10000, false);
		String tags[] = {"MD"};
		assertTrue(AlignmentComparison.compareFullAlignments(test, test, tags));
	}
	
	private void SAMFileReflexive(String filename, String[] tags) {
		ClassLoader classLoader = getClass().getClassLoader();
		String samFilename = classLoader.getResource(filename).getPath();
		String[] inputs = new String[1];
		inputs[0] = samFilename;
		try {
			assertTrue(AlignmentComparison.compareAlignmentFiles(samFilename, inputs, tags, null));
		} catch (IOException e) {
			fail(e.getMessage());
		}
	}
	
	@Test
	public void testSAMFileReflexiveSingle() {
		String[] tags = {"MD"};
		for (String samFilename : samFilenames) {
			SAMFileReflexive(samFilename, tags);
		}
	}
	
	private void differentSAM(String filename1, String filename2, String[] tags) {
		ClassLoader classLoader = getClass().getClassLoader();
		String samFilename1 = classLoader.getResource(filename1).getPath();
		String samFilename2 = classLoader.getResource(filename2).getPath();
		String[] inputs = new String[1];
		try {
			inputs[0] = samFilename2;
			assertFalse(AlignmentComparison.compareAlignmentFiles(samFilename1, inputs, tags, null));
			inputs[0] = samFilename1;
			assertFalse(AlignmentComparison.compareAlignmentFiles(samFilename2, inputs, tags, null));
		} catch (IOException e) {
			fail(e.getMessage());
		}
	}
	
	@Test
	public void differentSAMFiles() {
		String[] tags = {"MD"};
		for (int x = 0; x < samFilenames.length; x++) {
			for (int y = x+1; y < samFilenames.length; y++) {
				differentSAM(samFilenames[x], samFilenames[y], tags);
			}
		}
	}
	
	private void differentSAMCommandLine(String filename1, String filename2) throws IOException, ParseException {
		ClassLoader classLoader = getClass().getClassLoader();
		String samFilename1 = classLoader.getResource(filename1).getPath();
		String samFilename2 = classLoader.getResource(filename2).getPath();
		
		String[] commandArray = new String[5];
		commandArray[0] = "AlignmentComparison";
		commandArray[1] = "-c";
		commandArray[3] = "-i";
		commandArray[2] = samFilename1;
		commandArray[4] = samFilename2;
		assertEquals(AlignmentComparison.NO_MATCH, AlignmentComparison.main(commandArray));
		commandArray[2] = samFilename2;
		commandArray[4] = samFilename1;
		assertEquals(AlignmentComparison.NO_MATCH, AlignmentComparison.main(commandArray));
	}
	
	@Test
	public void differentSAMFilesCommandLine() throws IOException, ParseException {
		for (int x = 0; x < samFilenames.length; x++) {
			for (int y = x+1; y < samFilenames.length; y++) {
				differentSAMCommandLine(samFilenames[x], samFilenames[y]);
			}
		}
	}
	
	@Test
	public void testSelfPlusSecond() throws IOException {
		ClassLoader classLoader = getClass().getClassLoader();
		String samFilename1 = classLoader.getResource("test.sam").getPath();
		String samFilename2 = classLoader.getResource("target-test.sam").getPath();
		String[] inputs = new String[2];
		String[] tags = new String[0];
		// try both orders
		inputs[0] = samFilename1;
		inputs[1] = samFilename2;
		assertFalse(AlignmentComparison.compareAlignmentFiles(samFilename1, inputs, tags, null));
		assertFalse(AlignmentComparison.compareAlignmentFiles(samFilename2, inputs, tags, null));
		inputs[0] = samFilename2;
		inputs[1] = samFilename1;
		assertFalse(AlignmentComparison.compareAlignmentFiles(samFilename1, inputs, tags, null));
		assertFalse(AlignmentComparison.compareAlignmentFiles(samFilename2, inputs, tags, null));
	}
	
	@Test
	public void testMerge() throws IOException {
		ClassLoader classLoader = getClass().getClassLoader();
		String samFilename1 = classLoader.getResource("alignment_comparison/multi_lib123.sam").getPath();
		String samFilename2 = classLoader.getResource("alignment_comparison/multi_lib45.sam").getPath();
		String merged = classLoader.getResource("multi_lib.sam").getPath();
		String[] inputs = {samFilename1, samFilename2};
		String[] tags = {"RG", "XD", "MD", "NM"};
		
		assertTrue(AlignmentComparison.compareAlignmentFiles(merged, inputs, tags, null));
	}
}
