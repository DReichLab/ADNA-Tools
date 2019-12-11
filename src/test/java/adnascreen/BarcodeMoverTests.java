package adnascreen;

import static org.junit.Assert.*;

import java.io.BufferedInputStream;
import java.io.FileInputStream;
import java.io.IOException;

import org.junit.Test;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class BarcodeMoverTests {
	
	@Test
	public void file_barcodes(){
		try{
			ClassLoader classLoader = getClass().getClassLoader();
			String samFilename = classLoader.getResource("heng_shotgun.sam").getPath();
			
			SamInputResource bufferedSAMFile = SamInputResource.of(new BufferedInputStream(new FileInputStream(samFilename)));
			SamReader reader = SamReaderFactory.makeDefault().open(bufferedSAMFile);
			SAMRecordIterator i = reader.iterator();
			SAMRecord record = i.next();
			SAMRecord modified = BarcodeMover.barcodeToTag(record, DemultiplexSAM.duplicatesSAMTag);
			
			assertEquals("GTATCGG_TTACAGT_44", modified.getAttribute(DemultiplexSAM.duplicatesSAMTag));
			assertEquals("HL32WALXX170627:4:2119:7050:54102_21:GTATCGGTTACAGT", record.getReadName());
			assertEquals("HL32WALXX170627:4:2119:7050:54102_21", modified.getReadName());
			assertEquals(record.getBaseQualityString(), modified.getBaseQualityString());
			assertEquals(record.getCigarString(), modified.getCigarString());
			assertArrayEquals(record.getReadBases(), modified.getReadBases());
			
			record = i.next();
			modified = BarcodeMover.barcodeToTag(record, DemultiplexSAM.duplicatesSAMTag);
			
			assertEquals("CGTGACC_CCGTGAC_84", modified.getAttribute(DemultiplexSAM.duplicatesSAMTag));
			assertEquals("HL32WALXX170627:4:1211:3589:49795_21:CGTGACCCCGTGAC", record.getReadName());
			assertEquals("HL32WALXX170627:4:1211:3589:49795_21", modified.getReadName());
			assertEquals(record.getBaseQualityString(), modified.getBaseQualityString());
			assertEquals(record.getCigarString(), modified.getCigarString());
			assertArrayEquals(record.getReadBases(), modified.getReadBases());
			
			record = i.next();
			modified = BarcodeMover.barcodeToTag(record, DemultiplexSAM.duplicatesSAMTag);
			
			assertEquals("TACAGTT_AACGCTA_82", modified.getAttribute(DemultiplexSAM.duplicatesSAMTag));
			assertEquals("HL32WALXX170627:4:2203:7659:60114_21:TACAGTTAACGCTA", record.getReadName());
			assertEquals("HL32WALXX170627:4:2203:7659:60114_21", modified.getReadName());
			assertEquals(record.getBaseQualityString(), modified.getBaseQualityString());
			assertEquals(record.getCigarString(), modified.getCigarString());
			assertArrayEquals(record.getReadBases(), modified.getReadBases());
			
			record = i.next();
			modified = BarcodeMover.barcodeToTag(record, DemultiplexSAM.duplicatesSAMTag);
			
			assertEquals("TACAGTT_CCGTGAC_85", modified.getAttribute(DemultiplexSAM.duplicatesSAMTag));
			assertEquals("HL32WALXX170627:4:1216:3041:41602_21:TACAGTTCCGTGAC", record.getReadName());
			assertEquals("HL32WALXX170627:4:1216:3041:41602_21", modified.getReadName());
			assertEquals(record.getBaseQualityString(), modified.getBaseQualityString());
			assertEquals(record.getCigarString(), modified.getCigarString());
			assertArrayEquals(record.getReadBases(), modified.getReadBases());
			
			record = i.next();
			modified = BarcodeMover.barcodeToTag(record, DemultiplexSAM.duplicatesSAMTag);
			
			assertEquals("ACGCTAA_CCGTGAC_67", modified.getAttribute(DemultiplexSAM.duplicatesSAMTag));
			assertEquals("HL32WALXX170627:4:2106:18710:25534_21:ACGCTAACCGTGAC", record.getReadName());
			assertEquals("HL32WALXX170627:4:2106:18710:25534_21", modified.getReadName());
			assertEquals(record.getBaseQualityString(), modified.getBaseQualityString());
			assertEquals(record.getCigarString(), modified.getCigarString());
			assertArrayEquals(record.getReadBases(), modified.getReadBases());
		}
		catch(IOException e){
			fail(e.toString());
		}
	}
}
