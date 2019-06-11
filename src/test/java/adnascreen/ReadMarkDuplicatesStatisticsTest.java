package adnascreen;

import static org.junit.Assert.*;

import org.junit.Test;

public class ReadMarkDuplicatesStatisticsTest {
	@Test
	public void readNumberDuplicatesTest(){
		int numDuplicates = -2;
		ClassLoader classLoader = getClass().getClassLoader();
		String filename = classLoader.getResource("1-2-Q20-Q41.sam.dedup_stats").getPath();
		try{
			numDuplicates = ReadMarkDuplicatesStatistics.readMarkDuplicateStatistics(filename);
			assertEquals(449, numDuplicates);
		}
		catch(Exception e){
			fail();
		}
	}
	
	@Test
	public void keyParseFromFilename1() {
		String filename = "/test/abc/GTCGCAG_TTGGATC_CGCTGAG-GTGATCT-TATCAGA-ACAGCTC_AGCATCA-CTGCAGC-GATGCTG-TCATGAT.bam.stats";
		String expectedKey = "GTCGCAG_TTGGATC_CGCTGAG-GTGATCT-TATCAGA-ACAGCTC_AGCATCA-CTGCAGC-GATGCTG-TCATGAT";
		String result = ReadMarkDuplicatesStatistics.keyFromFilename(filename);
		assertEquals(expectedKey, result);
	}
	
	@Test
	public void keyParseFromFilename2() {
		String filename = "/test/dir/S6913.E1.L6.bam.stats";
		String expectedKey = "S6913.E1.L6";
		String result = ReadMarkDuplicatesStatistics.keyFromFilename(filename);
		assertEquals(expectedKey, result);
	}
}
