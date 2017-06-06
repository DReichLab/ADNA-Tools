package adnascreen;

import static org.junit.Assert.*;

import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.ExpectedException;

public class FASTQHeaderTests {
	@Rule
	public ExpectedException thrown = ExpectedException.none();
	
	@Test
	public void identical(){
		String header = "@NS500217:348:HTW2FBGXY:1:11101:22352:1064 1:N:0:0";
		FASTQHeader h1 = new FASTQHeader(header);
		FASTQHeader h2 = new FASTQHeader(header);
		assertTrue(h1.equalsExceptRead(h2));
		assertEquals(h1, h2);
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
		thrown.expect(IllegalArgumentException.class);
		new FASTQHeader("@NS500217:348:HTW2FBGXY:1:11101:22352:1064 1:Z:0:0");
		fail("Filtered field has invalid character but was accepted"); // should not reach here
	}
	
	@Test
	public void readGroup(){
		String header = "@NS500217:348:HTW2FBGXY:1:11101:22352:1064 1:N:0:0";
		FASTQHeader h1 = new FASTQHeader(header);
		String readGroup = h1.getReadGroup();
		String readGroupHeaderSection = "NS500217:348:HTW2FBGXY:1";
		String expected = readGroupHeaderSection.replaceAll(":", String.valueOf(FASTQHeader.READ_GROUP_FIELD_DELIMITER));
		assertEquals(expected, readGroup);
	}
}
