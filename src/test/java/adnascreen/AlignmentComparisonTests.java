package adnascreen;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertThrows;
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
	String tags[] = {"MD", "XD", "RG", "NM"};
	
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
		assertTrue(AlignmentComparison.compareFullAlignments(x, y, tags));
	}
	
	@Test
	public void testRecordDuplicateFlag() throws IOException {
		SAMRecord x = getRecord("multi_lib.sam", 0);
		SAMRecord y = x.deepCopy();
		x.setDuplicateReadFlag(true);
		assertTrue(AlignmentComparison.compareFullAlignments(x, y, tags));
	}
	
	@Test
	public void testRecordVendorQualityFlagFail() throws IOException {
		SAMRecord x = getRecord("multi_lib.sam", 0);
		SAMRecord y = x.deepCopy();
		x.setReadFailsVendorQualityCheckFlag(true);
		assertFalse(AlignmentComparison.compareFullAlignments(x, y, tags));
	}
	
	@Test public void testRecordBaseFail() throws IOException {
		SAMRecord x = getRecord("multi_lib.sam", 0);
		SAMRecord y = x.deepCopy();
		byte[] bases = x.getReadBases();
		bases[1] = 'T';
		x.setReadBases(bases);
		assertFalse(AlignmentComparison.compareFullAlignments(x, y, tags));
	}
	
	@Test public void testRecordBaseQualityFail() throws IOException {
		SAMRecord x = getRecord("multi_lib.sam", 0);
		SAMRecord y = x.deepCopy();
		byte[] qualities = x.getBaseQualities();
		qualities[2] = 'A';
		x.setBaseQualities(qualities);
		assertFalse(AlignmentComparison.compareFullAlignments(x, y, tags));
	}
	
	@Test public void testRecordXDFail() throws IOException {
		SAMRecord x = getRecord("multi_lib.sam", 0);
		SAMRecord y = x.deepCopy();
		x.setAttribute("XD", "ATTGGCA_AAGCATC_40");
		assertFalse(AlignmentComparison.compareFullAlignments(x, y, tags));
	}
	
	@Test
	public void testRecordMappingQualityFlag() throws IOException {
		SAMRecord x = getRecord("multi_lib.sam", 0);
		SAMRecord y = x.deepCopy();
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
		
		String[] commandArray = new String[4];
		commandArray[0] = "-c";
		commandArray[2] = "-i";
		commandArray[1] = samFilename1;
		commandArray[3] = samFilename2;
		assertEquals(AlignmentComparison.NO_MATCH, AlignmentComparison.main(commandArray));
		commandArray[1] = samFilename2;
		commandArray[3] = samFilename1;
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
		
		assertTrue(AlignmentComparison.compareAlignmentFiles(merged, inputs, tags, null));
	}
	
	@Test
	public void testMergeCommandLine() throws IOException, ParseException {
		ClassLoader classLoader = getClass().getClassLoader();
		String samFilename1 = classLoader.getResource("alignment_comparison/multi_lib123.sam").getPath();
		String samFilename2 = classLoader.getResource("alignment_comparison/multi_lib45.sam").getPath();
		String merged = classLoader.getResource("multi_lib.sam").getPath();
		
		String[] commandArray = new String[9];
		commandArray[0] = "-i";
		commandArray[1] = samFilename1;
		commandArray[2] = samFilename2;
		commandArray[3] = "-c";
		commandArray[4] = merged;
		for (int i = 0; i  < tags.length; i++)
			commandArray[i+5] = tags[i];
		int result = AlignmentComparison.main(commandArray);
		assertEquals(AlignmentComparison.MATCH, result);
	}
	
	@Test
	public void testNoInputsFail() {
		String[] commandArray = new String[0];
		assertThrows(org.apache.commons.cli.MissingOptionException.class, () -> AlignmentComparison.main(commandArray));
	}
	
	@Test
	public void testCheckOnlyNoInputsFail() {
		ClassLoader classLoader = getClass().getClassLoader();
		String merged = classLoader.getResource("multi_lib.sam").getPath();
		
		String[] commandArray = new String[2];
		commandArray[0] = "-c";
		commandArray[1] = merged;
	
		assertThrows(org.apache.commons.cli.MissingOptionException.class, () -> AlignmentComparison.main(commandArray));
	}
	
	@Test
	public void testNoCheckFail() {
		ClassLoader classLoader = getClass().getClassLoader();
		String samFilename1 = classLoader.getResource("alignment_comparison/multi_lib123.sam").getPath();
		String samFilename2 = classLoader.getResource("alignment_comparison/multi_lib45.sam").getPath();
		
		String[] commandArray = new String[3];
		commandArray[0] = "-i";
		commandArray[1] = samFilename1;
		commandArray[2] = samFilename2;
	
		assertThrows(org.apache.commons.cli.MissingOptionException.class, () -> AlignmentComparison.main(commandArray));
	}
	
	@Test
	public void testDifferentReferenceDictionary() {
		ClassLoader classLoader = getClass().getClassLoader();
		String samFilename1 = classLoader.getResource("alignment_comparison/multi_lib45.sam").getPath();
		String samFilename2 = classLoader.getResource("alignment_comparison/multi_lib45_extraSeqDict.sam").getPath();
		
		String[] inputs = new String[1];
		inputs[0] = samFilename2;
		try {
			assertFalse(AlignmentComparison.compareAlignmentFiles(samFilename1, inputs, tags, null));
		} catch (IOException e) {
			fail(e.getMessage());
		}
	}
	
	@Test
	public void testDifferentReferenceDictionaryCommandLine() throws IOException, ParseException {
		String samFilename1 = "alignment_comparison/multi_lib45.sam";
		String samFilename2 = "alignment_comparison/multi_lib45_extraSeqDict.sam";
		differentSAMCommandLine(samFilename1, samFilename2);
	}
}
