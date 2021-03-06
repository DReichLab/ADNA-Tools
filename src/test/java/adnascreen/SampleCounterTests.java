package adnascreen;

import org.junit.Test;

import static org.junit.Assert.*;

import java.util.List;

public class SampleCounterTests {
	@Test
	public void basic(){
		SampleCounter sampleCounter = new SampleCounter();
		
		String test = "test";
		String dne = "empty";
		sampleCounter.increment(test);
		assertEquals(1, sampleCounter.get(test));
		assertEquals(0, sampleCounter.get(dne));
		
		List<String> labels = sampleCounter.getLabelList();
		assertEquals(1, labels.size());
		assertEquals(test, labels.get(0));
	}
	
	@Test
	public void equals(){
		SampleCounter a = new SampleCounter();
		SampleCounter b = new SampleCounter();
		assertEquals(a, b);
		
		String test = "test";
		a.increment(test);
		assertNotEquals(a, b);
		b.increment(test);
		assertEquals(a, b);
		
		a.increment("a_only");
		assertNotEquals(a, b);
	}
	
	@Test
	public void combine(){
		SampleCounter a = new SampleCounter();
		SampleCounter b = new SampleCounter();
		
		a.increment("both");
		b.increment("both");
		b.increment("both");
		
		a.increment("a_only");
		b.increment("b_only");
		
		SampleCounter c = new SampleCounter(a);
		assertEquals(a, c);
		assertEquals(c, a);
		
		c.combine(b);
		assertEquals(3, c.get("both"));
		assertEquals(1, c.get("a_only"));
		assertEquals(1, c.get("b_only"));
	}
	
	@Test
	public void counterToString(){
		SampleCounter s = new SampleCounter();
		s.add("raw", 10);
		s.add("merged", 8);
		s.add("aligned", 7);
		assertEquals("raw\t10\tmerged\t8\taligned\t7", s.toString());
		assertEquals(10, s.get("raw"));
		assertEquals(8, s.get("merged"));
		assertEquals(7, s.get("aligned"));
		
		// also test constructor from string
		SampleCounter fromString = new SampleCounter(s.toString());
		assertEquals(s, fromString);
		
		List<String> labels = s.getLabelList();
		assertEquals(3, labels.size());
		// order should be the same as the added order
		assertEquals("raw", labels.get(0));
		assertEquals("merged", labels.get(1));
		assertEquals("aligned", labels.get(2));
	}
	
	@Test
	// Test that we can handle long-sized values
	public void longCount() {
		SampleCounter s = new SampleCounter();
		long expected = 1L + (long) Integer.MAX_VALUE;
		
		String TEST = "long_test";
		s.add(TEST, Integer.MAX_VALUE);
		s.increment(TEST);
		
		assertEquals(expected, s.get(TEST));
	}
	
	@Test
	// test copy constructor with long-sized value
	public void copy() {
		SampleCounter s = new SampleCounter();
		long expected = 1L + (long) Integer.MAX_VALUE;
		String TEST = "long_test";
		s.add(TEST, expected);
		
		SampleCounter s2 = new SampleCounter(s);
		assertEquals(s.get(TEST), s2.get(TEST));
	}
	
	private final String RAW = "raw";
	private final String MERGED = "merged";
	private final String ALIGNED = "aligned";
	@Test
	public void basicSets(){
		String key1 = new IndexAndBarcodeKey("1_2_Q1_Q2").toString();
		String key2 = new IndexAndBarcodeKey("3_4_Q3_Q4").toString();
		
		SampleSetsCounter set1 = new SampleSetsCounter();
		SampleSetsCounter set2 = new SampleSetsCounter();
		
		set1.increment();
		set2.increment();
		set1.increment(key1, RAW);
		set2.increment(key2, RAW);
		assertNotEquals(set1, set2);
		
		set1.increment();
		set2.increment();
		set1.increment(key2, RAW);
		set2.increment(key1, RAW);
		assertEquals(set1, set2);
		
		assertEquals(1, set1.get(key1, RAW));
		assertEquals(1, set1.get(key2, RAW));
	}
	
	@Test
	public void combineSets(){
		String key1 = new IndexAndBarcodeKey("1_2_Q1_Q2").toString();
		String key2 = new IndexAndBarcodeKey("3_4_Q3_Q4").toString();
		
		SampleSetsCounter set1 = new SampleSetsCounter();
		SampleSetsCounter set2 = new SampleSetsCounter();
		SampleSetsCounter set3 = new SampleSetsCounter();
		
		set1.increment();
		set1.increment(key1, RAW);
		set1.increment();
		set1.increment(key1, RAW);
		set1.increment();
		set1.increment(key2, RAW);
		
		set2.increment();
		set2.increment(key1, RAW);
		set2.increment();
		set2.increment(key2, RAW);
		
		set3.increment();
		set3.increment(key1, RAW);
		assertNotEquals(set1, set2);
		assertNotEquals(set2, set3);
		assertNotEquals(set1, set3);
		
		SampleSetsCounter set4 = new SampleSetsCounter(set2);
		assertEquals(set2, set4);
		set4.combine(set3);
		assertNotEquals(set2, set4);
		assertEquals(set1, set4);
		assertEquals(set1.toString(), set4.toString());
	}
	
	@Test
	public void setsToString(){
		String key1 = new IndexAndBarcodeKey("1_2_Q1_Q2").toString();
		String key2 = new IndexAndBarcodeKey("3_4_Q3_Q4").toString();
		
		SampleSetsCounter set1 = new SampleSetsCounter();
		
		set1.increment();
		set1.increment(key1, RAW);
		set1.increment();
		set1.increment(key1, RAW);
		set1.increment(key1, MERGED);
		set1.increment(key1, ALIGNED);
		set1.increment();
		set1.increment(key2, RAW);
		
		String s = set1.toString();
		SampleSetsCounter fromString = new SampleSetsCounter(s);
		assertEquals(set1, fromString);
	}
	
	@Test
	public void sortedToString(){
			String key1 = new IndexAndBarcodeKey("1_4_Q1_Q4").toString();
			String key2 = new IndexAndBarcodeKey("2_5_Q2_Q5").toString();
			String key3 = new IndexAndBarcodeKey("3_6_Q3_Q6").toString();
			
			SampleSetsCounter set1 = new SampleSetsCounter();
			
			set1.add(9);
			set1.add(key1, RAW, 9);
			set1.add(key1, MERGED, 4);
			set1.add(key1, ALIGNED, 3);
			
			set1.add(10);	
			set1.add(key2, RAW, 10);
			set1.add(key2, MERGED, 8);
			set1.add(key2, ALIGNED, 7);
			
			set1.add(10);	
			set1.add(key3, RAW, 10);
			set1.add(key3, MERGED, 5);
			set1.add(key3, ALIGNED, 2);
			
			String s = set1.toString();
			SampleSetsCounter fromString = new SampleSetsCounter(s);
			assertEquals(set1, fromString);
			
			// sorting section
			StringBuilder expectedBuilder = new StringBuilder();
			expectedBuilder.append("29\n");
			
			expectedBuilder.append(key2.toString());
			expectedBuilder.append("\t");
			expectedBuilder.append(RAW);
			expectedBuilder.append("\t");
			expectedBuilder.append(10);
			expectedBuilder.append("\t");
			expectedBuilder.append(MERGED);
			expectedBuilder.append("\t");
			expectedBuilder.append(8);
			expectedBuilder.append("\t");
			expectedBuilder.append(ALIGNED);
			expectedBuilder.append("\t");
			expectedBuilder.append(7);
			expectedBuilder.append("\n");
			
			expectedBuilder.append(key1.toString());
			expectedBuilder.append("\t");
			expectedBuilder.append(RAW);
			expectedBuilder.append("\t");
			expectedBuilder.append(9);
			expectedBuilder.append("\t");
			expectedBuilder.append(MERGED);
			expectedBuilder.append("\t");
			expectedBuilder.append(4);
			expectedBuilder.append("\t");
			expectedBuilder.append(ALIGNED);
			expectedBuilder.append("\t");
			expectedBuilder.append(3);
			expectedBuilder.append("\n");
			
			expectedBuilder.append(key3.toString());
			expectedBuilder.append("\t");
			expectedBuilder.append(RAW);
			expectedBuilder.append("\t");
			expectedBuilder.append(10);
			expectedBuilder.append("\t");
			expectedBuilder.append(MERGED);
			expectedBuilder.append("\t");
			expectedBuilder.append(5);
			expectedBuilder.append("\t");
			expectedBuilder.append(ALIGNED);
			expectedBuilder.append("\t");
			expectedBuilder.append(2);
			
			assertEquals(expectedBuilder.toString(), set1.toStringSorted(ALIGNED));
	}
}
