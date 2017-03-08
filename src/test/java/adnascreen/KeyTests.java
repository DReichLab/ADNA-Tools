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
}
