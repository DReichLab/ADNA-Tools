package adnascreen;

import static org.junit.Assert.*;

import java.io.IOException;

import org.junit.Test;

public class IndexMatcherTests {
	@Test
	public void noMatch(){
		IndexMatcher indexMatcher = new IndexMatcher();
		DNASequence reference = new DNASequence("AAAAAAA");
		indexMatcher.addReference(reference);
		
		DNASequence result;
		DNASequence query1 = new DNASequence("TTTTTTT");
		result = indexMatcher.find(query1);
		assertNull(result);
		
		DNASequence query2 = new DNASequence("AAAAAAT");
		result = indexMatcher.find(query2);
		assertNull(result);
	}
	
	@Test
	public void matchDiff1(){
		IndexMatcher indexMatcher = new IndexMatcher();
		DNASequence reference = new DNASequence("AAAAAAA");
		indexMatcher.addReference(reference);
		indexMatcher.setMaxHammingDistance(1);
		
		DNASequence query = new DNASequence("AAAAAAT");
		DNASequence result = indexMatcher.find(query);
		assertEquals(result, reference);
	}
	
	@Test
	public void i7_search(){
		ClassLoader classLoader = getClass().getClassLoader();
		String filename = classLoader.getResource("i7").getPath();
		try{
			IndexMatcher indexMatcher = new IndexMatcher(filename, 1);
			DNASequence queryExact = new DNASequence("TCGCAGG");
			DNASequence result = indexMatcher.find(queryExact);
			assertEquals(queryExact, result);
			
			DNASequence queryOffByOne = new DNASequence("TCGCAGT");
			result = indexMatcher.find(queryOffByOne);
			assertEquals(queryExact, result);
			
			// repeat for cache
			result = indexMatcher.find(queryOffByOne);
			assertEquals(queryExact, result);
		}
		catch(IOException e){
			fail();
		}
	}
}
