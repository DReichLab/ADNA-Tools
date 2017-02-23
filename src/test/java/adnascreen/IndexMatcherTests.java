package adnascreen;

import static org.junit.Assert.*;

import java.io.IOException;

import org.junit.Test;

public class IndexMatcherTests {
	@Test
	public void noMatch(){
		BarcodeMatcher indexMatcher = new BarcodeMatcher();
		DNASequence reference = new DNASequence("AAAAAAA");
		indexMatcher.addReferenceSet(reference.toString(), "1");
		
		String result;
		DNASequence query1 = new DNASequence("TTTTTTT");
		result = indexMatcher.find(query1);
		assertNull(result);
		// run a second time for cache effects
		result = indexMatcher.find(query1);
		assertNull(result);
		
		// default Hamming distance is 0, so off by one query should return empty
		DNASequence query2 = new DNASequence("AAAAAAT");
		result = indexMatcher.find(query2);
		assertNull(result);
	}
	
	@Test
	public void matchDiff1(){
		BarcodeMatcher indexMatcher = new BarcodeMatcher();
		DNASequence reference = new DNASequence("AAAAAAA");
		String label = "1";
		indexMatcher.addReferenceSet(reference.toString(), label);
		indexMatcher.setMaxHammingDistance(1);
		
		DNASequence query = new DNASequence("AAAAAAT");
		String result = indexMatcher.find(query);
		assertEquals(label, result);
	}
	
	@Test
	public void i7_search(){
		ClassLoader classLoader = getClass().getClassLoader();
		String filename = classLoader.getResource("i7").getPath();
		try{
			BarcodeMatcher indexMatcher = new BarcodeMatcher(filename, 1);
			DNASequence queryExact = new DNASequence("TCGCAGG");
			String expectedLabel = "1";
			String result = indexMatcher.find(queryExact);
			assertEquals(expectedLabel, result);
			
			DNASequence queryOffByOne = new DNASequence("TCGCAGT");
			result = indexMatcher.find(queryOffByOne);
			assertEquals(expectedLabel, result);
			
			// repeat for cache
			result = indexMatcher.find(queryOffByOne);
			assertEquals(expectedLabel, result);
		}
		catch(IOException e){
			fail();
		}
	}
}
