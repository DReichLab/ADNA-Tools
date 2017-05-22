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
}
