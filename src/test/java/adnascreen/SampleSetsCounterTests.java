package adnascreen;

import org.junit.Test;
import static org.junit.Assert.*;

public class SampleSetsCounterTests {
	@Test
	public void rawLongCount() {
		SampleSetsCounter s = new SampleSetsCounter();
		long expected = 1L + (long) Integer.MAX_VALUE;
		long result = s.add(expected);
		
		assertEquals(expected, result);
	}
	
	@Test
	public void retrieval() {
		SampleSetsCounter s = new SampleSetsCounter();
		
		final String key = "TEST";
		final String label = "raw";
		long value = 1L + (long) Integer.MAX_VALUE;
		
		s.add(key, label, value);
		assertEquals(value, s.get(key, label));
		
		// increment using add and retest
		s.add(key, label, 1);
		assertEquals(value + 1, s.get(key, label));
		
		// increment and retest
		s.increment(key, label);
		assertEquals(value + 2, s.get(key, label));
	}
	
	@Test
	public void equalsTest() {
		SampleSetsCounter s = new SampleSetsCounter();
		
		final String key = "TEST";
		final String label = "raw";
		long value = 1L + (long) Integer.MAX_VALUE;
		
		s.add(key, label, value);
		
		SampleSetsCounter s2 = new SampleSetsCounter(s);
		
		assertEquals(s, s2);
		
		// now alter one and retest
		s.increment(key, label);
		assertNotEquals(s, s2);
	}
	
	@Test
	public void equalsTestRawNotEqual() {
		SampleSetsCounter s = new SampleSetsCounter();
		
		final String key = "TEST";
		final String label = "raw";
		long value = 1L + (long) Integer.MAX_VALUE;
		
		s.add(key, label, value);
		
		SampleSetsCounter s2 = new SampleSetsCounter(s);
		
		assertEquals(s, s2);
		
		// now alter raw and retest
		s.increment();
		assertNotEquals(s, s2);
	}

	
	@Test
	public void stringAndBack() {
		SampleSetsCounter s = new SampleSetsCounter();
		final String key = "TEST";
		final String label = "raw";
		long value = 1L + (long) Integer.MAX_VALUE;
		
		s.add(key, label, value);
		
		SampleSetsCounter s2 = new SampleSetsCounter(s.toString());
		
		assertEquals(s, s2);
	}
	
	@Test
	public void combine() {
		SampleSetsCounter s1 = new SampleSetsCounter();
		final String key1 = "TEST1";
		final String label1A = "1A";
		final String label1B = "1B";
		final long A = 100;
		final long B = 1000;
		s1.add(A);
		s1.add(B);
		s1.add(key1, label1A, A);
		s1.add(key1, label1B, B);
		
		SampleSetsCounter s2 = new SampleSetsCounter();
		final String key2 = "TEST2";
		final long A2 = 50000;
		s2.add(A2);
		s1.add(key1, label1A, 1);
		s2.add(key2, label1A, A2);
		
		s1.combine(s2);
		
		assertEquals(s1.get(key1, label1A), A+1);
		assertEquals(s1.get(key1, label1B), B);
		assertEquals(s1.get(key2, label1A), A2);
		
		// test from string and back with this more complicated examples
		SampleSetsCounter fromString = new SampleSetsCounter(s1.toString());
		assertEquals(s1, fromString);
	}
}
