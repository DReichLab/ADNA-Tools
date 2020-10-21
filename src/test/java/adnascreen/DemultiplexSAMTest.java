package adnascreen;

import static org.junit.Assert.*;

import java.io.BufferedInputStream;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Queue;
import java.util.Set;

import org.apache.commons.cli.ParseException;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TemporaryFolder;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class DemultiplexSAMTest {
	@Rule
	public TemporaryFolder testFolder = new TemporaryFolder();

	@Test
	public void testSelectTopSample() {
		ClassLoader classLoader = getClass().getClassLoader();
		String filename = classLoader.getResource("fastq/expected").getPath();
		
		try {
			int numTopSamples = 1;
			int minimumReads = 1;
			int thresholdReads = -1;
			Queue<IndexAndBarcodeKey> topQueue = DemultiplexSAM.selectTopSamples(filename, numTopSamples, minimumReads, thresholdReads);
			assertEquals(1, topQueue.size());
			IndexAndBarcodeKey key = topQueue.remove();
			assertEquals("CAGGTCG_GAATCTC_CGCTGAG:GTGATCT:TATCAGA:ACAGCTC_TCGCATT:AGTGCAA:CTATGCC:GACATGG", key.toString() );
		} catch (IOException e) {
			fail();
		}
	}
	
	@Test
	public void testSelectTopSamplesWithThreshold() {
		ClassLoader classLoader = getClass().getClassLoader();
		String filename = classLoader.getResource("fastq/expected").getPath();
		
		try {
			int numTopSamples = 2;
			int minimumReads = 1;
			int thresholdReads = 4;
			Queue<IndexAndBarcodeKey> topQueue = DemultiplexSAM.selectTopSamples(filename, numTopSamples, minimumReads, thresholdReads);
			assertEquals(6, topQueue.size());
			IndexAndBarcodeKey key;
			
			key = topQueue.remove();
			assertEquals("CAGGTCG_GAATCTC_CGCTGAG:GTGATCT:TATCAGA:ACAGCTC_TCGCATT:AGTGCAA:CTATGCC:GACATGG", key.toString() );
			key = topQueue.remove();
			assertEquals("CAGGTCG_ATACTGA_CTAGACA:GACTCGC:TCGAGTG:AGTCTAT_CTGGCTA:GATTGAC:TCAATCG:AGCCAGT", key.toString() );
			key = topQueue.remove();
			assertEquals("CAGGTCG_TACTTAG_ACTACTG:CGACGAT:GTCGTCA:TAGTAGC_CATAGGC:GCACTTG:TGCGAAT:ATGTCCA", key.toString() );
			key = topQueue.remove();
			assertEquals("CAGGTCG_CCGATTG_CGTGACC:GTATCGG:TACAGTT:ACGCTAA_AGTCACG:CTAGCGT:GACTGTA:TCGATAC", key.toString() );
			key = topQueue.remove();
			assertEquals("CAGGTCG_ATGGAGA_GTCTCAA:TAGAGCC:ACTCTGG:CGAGATT_GGAGTAC:TTCTACG:AAGACGT:CCTCGTA", key.toString() );
			key = topQueue.remove();
			assertEquals("CCGGATG_AATAAGC_CCTTACG:GGAACGT:TTCCGTA:AAGGTAC_TTGTCAG:AATAGCT:CCACTGA:GGCGATC", key.toString() );
		} catch (IOException e) {
			fail();
		}
	}
	
	/**
	 * Verify that the given bam exists and has the expected number of reads
	 * @param bamFilename
	 * @param expectedNumReads
	 * @return
	 * @throws FileNotFoundException 
	 */
	public void bamChecks(String parentDirectory, String filenameOnly, int expectedNumReads) throws FileNotFoundException {
		String bam_path = parentDirectory + "/" + filenameOnly;
		String [] filenameParts = filenameOnly.split("\\.");
		assertEquals(2, filenameParts.length);
		assertTrue(filenameParts[1].equals("bam") || filenameParts[1].equals("sam"));
		String [] indicesAndBarcodes = filenameParts[0].split(String.valueOf(IndexAndBarcodeKey.FIELD_SEPARATOR));
		String p5BarcodeString = indicesAndBarcodes[2];
		String p7BarcodeString = indicesAndBarcodes[3];
		Set<String> p5Barcodes = new HashSet<String>(Arrays.asList(p5BarcodeString.split(":|-")));
		assertTrue(p5Barcodes.size() == 4 || p5Barcodes.size() == 0);
		Set<String> p7Barcodes = new HashSet<String>(Arrays.asList(p7BarcodeString.split(":|-")));
		assertTrue(p7Barcodes.size() == 4 || p7Barcodes.size() == 0);
		
		SamInputResource bufferedSAMFile = SamInputResource.of(new BufferedInputStream(new FileInputStream(bam_path)));
		SamReader reader = SamReaderFactory.makeDefault().open(bufferedSAMFile);
		SAMRecordIterator i = reader.iterator();

		int readCount = 0;
		SAMRecord record;
		while (i.hasNext()){
			record = i.next();
			readCount++;
			// check that barcodes in XD tag match filename
			String xd = (String) record.getAttribute(DemultiplexSAM.duplicatesSAMTag);
			String [] fields = xd.split(String.valueOf(IndexAndBarcodeKey.FIELD_SEPARATOR));
			String p5Barcode = fields[0];
			String p7Barcode = fields[1];
			int length = Integer.valueOf(fields[2]);
			assertTrue(p5Barcodes.contains(p5Barcode));
			assertTrue(p7Barcodes.contains(p7Barcode));
			
			assertEquals(record.getBaseQualities().length, length);
		}
		assertEquals(expectedNumReads, readCount);
	}
	
/*
CAGGTCG_GAATCTC_CGCTGAG:GTGATCT:TATCAGA:ACAGCTC_TCGCATT:AGTGCAA:CTATGCC:GACATGG	raw	7	merged	5
CAGGTCG_ATACTGA_CTAGACA:GACTCGC:TCGAGTG:AGTCTAT_CTGGCTA:GATTGAC:TCAATCG:AGCCAGT	raw	5	merged	4
CAGGTCG_TACTTAG_ACTACTG:CGACGAT:GTCGTCA:TAGTAGC_CATAGGC:GCACTTG:TGCGAAT:ATGTCCA	raw	4	merged	4
CAGGTCG_CCGATTG_CGTGACC:GTATCGG:TACAGTT:ACGCTAA_AGTCACG:CTAGCGT:GACTGTA:TCGATAC	raw	4	merged	2
CAGGTCG_ATGGAGA_GTCTCAA:TAGAGCC:ACTCTGG:CGAGATT_GGAGTAC:TTCTACG:AAGACGT:CCTCGTA	raw	4	merged	3
CCGGATG_AATAAGC_CCTTACG:GGAACGT:TTCCGTA:AAGGTAC_TTGTCAG:AATAGCT:CCACTGA:GGCGATC	raw	4	merged	3
 */
	
	@Test
	public void testDemultiplexSam() {
		String parentDirectory = testFolder.getRoot().toString();
		testDemultiplexCommon(parentDirectory, 2, 1, 0, false);
		// verify output
		try {
			bamChecks(parentDirectory, "CAGGTCG_GAATCTC_CGCTGAG-GTGATCT-TATCAGA-ACAGCTC_TCGCATT-AGTGCAA-CTATGCC-GACATGG.sam", 5);
			bamChecks(parentDirectory, "CAGGTCG_ATACTGA_CTAGACA-GACTCGC-TCGAGTG-AGTCTAT_CTGGCTA-GATTGAC-TCAATCG-AGCCAGT.sam", 4);
		} catch (Exception e){
			fail();
		}
	}
	
	// Same as sam above, but with bam format
	@Test
	public void testDemultiplexBam() {
		String parentDirectory = testFolder.getRoot().toString();
		testDemultiplexCommon(parentDirectory, 2, 1, 0, true);
		// verify output
		try {
			bamChecks(parentDirectory, "CAGGTCG_GAATCTC_CGCTGAG-GTGATCT-TATCAGA-ACAGCTC_TCGCATT-AGTGCAA-CTATGCC-GACATGG.bam", 5);
			bamChecks(parentDirectory, "CAGGTCG_ATACTGA_CTAGACA-GACTCGC-TCGAGTG-AGTCTAT_CTGGCTA-GATTGAC-TCAATCG-AGCCAGT.bam", 4);
		} catch (Exception e){
			fail();
		}
		
		// This is an output bam that should be present with threshold reads = 4, but not with the top 2 most common combinations only
		boolean missing = false;
		try {
			bamChecks(parentDirectory, "CAGGTCG_TACTTAG_ACTACTG-CGACGAT-GTCGTCA-TAGTAGC_CATAGGC-GCACTTG-TGCGAAT-ATGTCCA.bam", 4);
		} catch(FileNotFoundException e) {
			missing = true; 
		} catch(Exception e) {
			fail();
		}
		assertTrue(missing);
	}
	
	
	// Use threshold reads option. This is for low-level shotgun sequencing
	@Test
	public void testDemultiplexThreshold() {
		String parentDirectory = testFolder.getRoot().toString();
		testDemultiplexCommon(parentDirectory, 2, 1, 4, true);
		// verify output
		try {
			bamChecks(parentDirectory, "CAGGTCG_GAATCTC_CGCTGAG-GTGATCT-TATCAGA-ACAGCTC_TCGCATT-AGTGCAA-CTATGCC-GACATGG.bam", 5);
			bamChecks(parentDirectory, "CAGGTCG_ATACTGA_CTAGACA-GACTCGC-TCGAGTG-AGTCTAT_CTGGCTA-GATTGAC-TCAATCG-AGCCAGT.bam", 4);
			
			bamChecks(parentDirectory, "CAGGTCG_TACTTAG_ACTACTG-CGACGAT-GTCGTCA-TAGTAGC_CATAGGC-GCACTTG-TGCGAAT-ATGTCCA.bam", 4);
			bamChecks(parentDirectory, "CAGGTCG_CCGATTG_CGTGACC-GTATCGG-TACAGTT-ACGCTAA_AGTCACG-CTAGCGT-GACTGTA-TCGATAC.bam", 2);
			bamChecks(parentDirectory, "CAGGTCG_ATGGAGA_GTCTCAA-TAGAGCC-ACTCTGG-CGAGATT_GGAGTAC-TTCTACG-AAGACGT-CCTCGTA.bam", 3);
			bamChecks(parentDirectory, "CCGGATG_AATAAGC_CCTTACG-GGAACGT-TTCCGTA-AAGGTAC_TTGTCAG-AATAGCT-CCACTGA-GGCGATC.bam", 3);
		} catch (Exception e){
			fail();
		}
	}

	protected void testDemultiplexCommon(String parentDirectory, int numSamplesToOutput, int maximumConcurrentOpenFiles, int thresholdReads, boolean useBAM) {
		ClassLoader classLoader = getClass().getClassLoader();
		String filename_bam1 = classLoader.getResource("fastq/aligned_001.bam").getPath();
		String filename_bam2 = classLoader.getResource("fastq/aligned_002.bam").getPath();
		String statistics = classLoader.getResource("fastq/expected").getPath();
		
		String barcodeFilename = classLoader.getResource("fastq/barcodes").getPath();
		
		try {
			List<String> args = new LinkedList<String>();
			args.add("-s");
			args.add(statistics);
			args.add("-f");
			args.add(barcodeFilename);
			args.add("-n");
			args.add((new Integer(numSamplesToOutput)).toString());
			args.add("-m");
			args.add((new Integer(maximumConcurrentOpenFiles)).toString());
			if (thresholdReads > 0) {
				args.add("--thresholdReads");
				args.add((new Integer(thresholdReads)).toString());
			}
			if(useBAM) {
				args.add("--BAM");
			}
			args.add("--outputDirectory");
			args.add(parentDirectory);
			args.add(filename_bam1);
			args.add(filename_bam2);
			
			String [] args_array = new String[args.size()];
			int n = 0;
			for (String arg : args){
				args_array[n++] = arg;
			}
			DemultiplexSAM.main(args_array);
		} catch (IOException | ParseException e) {
			fail();
		}
	}
}
