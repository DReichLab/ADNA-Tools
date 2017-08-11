package adnascreen;

import static org.junit.Assert.*;
import org.junit.Test;

public class KeyTests {
	private static final String i5 = "5";
	private static final String i7 = "7";
	private static final String p5 = "Q5";
	private static final String p7 = "Q7";
	
	@Test
	public void basic(){		
		IndexAndBarcodeKey key = new IndexAndBarcodeKey(i5, i7, p5, p7);
		assertEquals(i5, key.getI5Label());
		assertEquals(i7, key.getI7Label());
		assertEquals(p5, key.getP5Label());
		assertEquals(p7, key.getP7Label());
	}
	
	@Test
	public void serialize(){
		IndexAndBarcodeKey key = new IndexAndBarcodeKey(i5, i7, p5, p7);
		String serialized = key.toString();
		IndexAndBarcodeKey fromString = new IndexAndBarcodeKey(serialized);
		assertEquals(key, fromString);
	}
	
	@Test
	public void flatten(){
		char separator = BarcodeMatcher.INDEX_DELIMITER;
		IndexAndBarcodeKey key1 = new IndexAndBarcodeKey("1", "2", "Q3" + separator + "3", "Q4" + separator + "4");
		IndexAndBarcodeKey key2 = new IndexAndBarcodeKey("1", "2", "Q3" + separator + "2", "Q4" + separator + "1");
		IndexAndBarcodeKey expectedFlattened = new IndexAndBarcodeKey("1", "2", "Q3", "Q4");
		
		assertNotEquals(key1, key2);
		
		IndexAndBarcodeKey key1Flattened = key1.flatten();
		IndexAndBarcodeKey key2Flattened = key2.flatten();
		assertNotEquals(key1, key1Flattened);
		assertNotEquals(key2, key2Flattened);
		assertEquals(key1Flattened, key2Flattened);
		assertEquals(expectedFlattened, key1Flattened);
		assertEquals(expectedFlattened, key2Flattened);
	}
	
	@Test
	public void nullValues(){
		IndexAndBarcodeKey key1 = new IndexAndBarcodeKey(i5, i7, p5, null);
		assertNotNull(key1.toString());
		IndexAndBarcodeKey flattened1 = key1.flatten();
		assertNotNull(flattened1.toString());
		assertEquals(key1, flattened1);
		
		IndexAndBarcodeKey key2 = new IndexAndBarcodeKey(i5, i7, null, p7);
		assertNotNull(key2.toString());
		IndexAndBarcodeKey flattened2 = key2.flatten();
		assertNotNull(flattened2.toString());
		assertEquals(key2, flattened2);
		
		assertNotEquals(key1, key2);
	}
}
