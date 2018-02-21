package adnascreen;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

public class DriverTests {
	@Test
	public void testBAM() {
		assertTrue(Driver.isBAMFilename("test.bam"));
	}
	
	@Test
	public void testLongBAM() {
		assertTrue(Driver.isBAMFilename("AGTTGGT_AAGCTAA_GCCTAAC-TGGACCG-ATTCGGT-CAAGTTA_GGAGTAC-TTCTACG-AAGACGT-CCTCGTA.prededup.clipped.sorted.bam"));
	}
	
	@Test
	public void testSAM() {
		assertFalse(Driver.isBAMFilename("test.sam"));
	}
	
	@Test
	public void testIndexBarcodeBAM() {
		assertTrue(Driver.isBAMFilename("AGTTGGT_AAGCTAA_GCCTAAC-TGGACCG-ATTCGGT-CAAGTTA_GGAGTAC-TTCTACG-AAGACGT-CCTCGTA.bam"));
	}
}
