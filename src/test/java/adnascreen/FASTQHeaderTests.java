package adnascreen;

import static org.junit.Assert.*;

import org.junit.Test;

public class FASTQHeaderTests {
	@Test
	public void identical(){
		String header = "@NS500217:348:HTW2FBGXY:1:11101:22352:1064 1:N:0:0";
		FASTQHeader h1 = new FASTQHeader(header);
		FASTQHeader h2 = new FASTQHeader(header);
		assertTrue(h1.equalsExceptRead(h2));
		assertEquals(h1, h2);
	}
	
	@Test
	public void fields(){
		String header = "@NS500217:348:HTW2FBGXY:1:11101:22352:1064 1:N:0:0";
		FASTQHeader h1 = new FASTQHeader(header);
		assertEquals("NS500217", h1.getInstrument());
		assertEquals(348, h1.getRunNumber());
		assertEquals("HTW2FBGXY", h1.getFlowcellID());
		assertEquals(1, h1.getLane());
		assertEquals(11101, h1.getTile());
		assertEquals(22352, h1.getX());
		assertEquals(1064, h1.getY());
		assertEquals(1, h1.getRead());
		assertEquals(false, h1.isFiltered());
		assertEquals(0, h1.getControlNumber());
		assertEquals("0", h1.getIndex());
	}
	
	@Test
	public void to_string(){
		String header = "@NS500217:348:HTW2FBGXY:1:11101:22352:1064 1:N:0:0";
		FASTQHeader h1 = new FASTQHeader(header);
		String s = h1.toString();
		assertEquals(s, header);
	}
	
	@Test
	public void different_read(){
		FASTQHeader	h1 = new FASTQHeader("@NS500217:348:HTW2FBGXY:1:11101:22352:1064 1:N:0:0");
		FASTQHeader h2 = new FASTQHeader("@NS500217:348:HTW2FBGXY:1:11101:22352:1064 2:N:0:0");
		assertTrue(h1.equalsExceptRead(h2));
		assertNotEquals(h1, h2);
	}
	
	@Test
	public void different_instrument(){
		FASTQHeader h1 = new FASTQHeader("@NS500217:348:HTW2FBGXY:1:11101:22352:1064 1:N:0:0");
		FASTQHeader h2 = new FASTQHeader("@NS500218:348:HTW2FBGXY:1:11101:22352:1064 1:N:0:0");
		assertFalse(h1.equalsExceptRead(h2));
		assertNotEquals(h1, h2);
	}
	
	@Test
	public void bad_filteredField(){
		// Filtered field has invalid character
		assertThrows(IllegalArgumentException.class, () -> new FASTQHeader("@NS500217:348:HTW2FBGXY:1:11101:22352:1064 1:Z:0:0"));
	}
	
	@Test
	public void readGroupElements(){
		String header = "@NS500217:348:HTW2FBGXY:1:11101:22352:1064 1:N:0:0";
		FASTQHeader h1 = new FASTQHeader(header);
		String expected = "PM:NS500217\tPU:HTW2FBGXY.348.1";
		String s = h1.getReadGroupElements();
		assertEquals(expected, s);
	}
	
	@Test
	public void hasIndex(){
		String header = "@NS500217:33:H0NW5AGXX:1:11101:10568:3456 1:N:0:TAATTCG+GAAATAA";
		FASTQHeader h1 = new FASTQHeader(header);
		String expected = "TAATTCG+GAAATAA";
		String s = h1.getIndex();
		assertEquals(expected, s);
	}
	
	// Other tests use data from NextSeq runs. Following tests use data from a Broad HiSeq X10.
	@Test
	public void badBroadRunNumber() {
		String header = "@E00151:HY352CCXY190421:HY352CCXY:1:1101:10003:10029 1:N:0:";
		FASTQHeader h1 = new FASTQHeader(header);
		
		assertEquals("E00151", h1.getInstrument());
		assertEquals(0, h1.getRunNumber()); // bad run numbers are converted to 0
		assertEquals("HY352CCXY", h1.getFlowcellID());
		assertEquals(1, h1.getLane());
		assertEquals(1101, h1.getTile());
		assertEquals(10003, h1.getX());
		assertEquals(10029, h1.getY());
		assertEquals(1, h1.getRead());
		assertEquals(false, h1.isFiltered());
		assertEquals(0, h1.getControlNumber());
		assertEquals("", h1.getIndex());
	}
	
	@Test
	public void badBroadRunNumber2() {
		String header = "@E00151:HY352CCXY190421:HY352CCXY:1:1101:10003:10029 2:N:0:";
		FASTQHeader h1 = new FASTQHeader(header);
		
		assertEquals("E00151", h1.getInstrument());
		assertEquals(0, h1.getRunNumber()); // bad run numbers are converted to 0
		assertEquals("HY352CCXY", h1.getFlowcellID());
		assertEquals(1, h1.getLane());
		assertEquals(1101, h1.getTile());
		assertEquals(10003, h1.getX());
		assertEquals(10029, h1.getY());
		assertEquals(2, h1.getRead());
		assertEquals(false, h1.isFiltered());
		assertEquals(0, h1.getControlNumber());
		assertEquals("", h1.getIndex());
	}
	
	@Test
	public void badBroadReadNumber() {
		String header = "@E00151:HY352CCXY190421:HY352CCXY:1:1101:10003:10029 :N:0:";
		FASTQHeader h1 = new FASTQHeader(header);
		
		assertEquals("E00151", h1.getInstrument());
		assertEquals(0, h1.getRunNumber()); // bad run numbers are converted to 0
		assertEquals("HY352CCXY", h1.getFlowcellID());
		assertEquals(1, h1.getLane());
		assertEquals(1101, h1.getTile());
		assertEquals(10003, h1.getX());
		assertEquals(10029, h1.getY());
		assertEquals(0, h1.getRead()); // bad read numbers are converted to 0
		assertEquals(false, h1.isFiltered());
		assertEquals(0, h1.getControlNumber());
		assertEquals("", h1.getIndex());
	}
	
	@Test
	// The Picard SamToFastq program also produces non-standard fastq headers
	public void broad_samtofastq() {
		String header = "@H2NMFCCXY170825:1:1101:10003:10029/2";
		FASTQHeader h1 = new FASTQHeader(header);
		
		assertEquals("", h1.getInstrument());
		assertEquals(0, h1.getRunNumber());
		assertEquals("H2NMFCCXY", h1.getFlowcellID());
		assertEquals(1, h1.getLane());
		assertEquals(1101, h1.getTile());
		assertEquals(10003, h1.getX());
		assertEquals(10029, h1.getY());
		
		assertEquals(2, h1.getRead());
		assertEquals(false, h1.isFiltered());
		assertEquals(0, h1.getControlNumber());
		assertEquals("", h1.getIndex());
	}
	
	@Test
	public void fastqWithIndexAndBarcodeKey() {
		String header = "@NS500217:706:HGL7NBGXB:1:11101:5330:1044;CCGGATG_GGCGGTC_TGATCTC:ATCAGAG:CAGCTCT:GCTGAGA.4_CTACTCG:GACGAGT:TCGTCTA:AGTAGAC.2 1:N:0:0";
		
		FASTQHeader h1 = new FASTQHeader(header);
		assertEquals(header, h1.toString());
		
		assertEquals(5330, h1.getX());
		assertEquals(1044, h1.getY());
		
		assertEquals(1, h1.getRead());
		assertEquals(false, h1.isFiltered());
		assertEquals(0, h1.getControlNumber());
		
		IndexAndBarcodeKey key = h1.getKey();
		assertEquals("CCGGATG", key.getI5Label());
		assertEquals("GGCGGTC", key.getI7Label());
		assertEquals("TGATCTC:ATCAGAG:CAGCTCT:GCTGAGA.4", key.getP5Label());
		assertEquals("CTACTCG:GACGAGT:TCGTCTA:AGTAGAC.2", key.getP7Label());
	}
}
