package adnascreen;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class DuplicatesTagRewrite {
	/**
	 * Take a bam, and rewrite each alignment's XD tag to the read length
	 * This is to combine and deduplicate bam for a single library
	 * If a bam has no barcode information, then we cannot use barcodes in deduplicating the library
	 * We treat this conservatively by removing barcodes and deduplicating using only the position and length
	 * @param args
	 * @throws ParseException 
	 * @throws IOException 
	 */
	public static void main(String [] args) throws ParseException, IOException {
		CommandLineParser parser = new DefaultParser();
		Options options = new Options();
		options.addRequiredOption("i", "input-filename", true, "input SAM/BAM filename");
		options.addRequiredOption("o", "output-filename", true, "output SAM/BAM filename");
		// Ultima wrote barcodes to separate tags, so adding an option to write
		// existing tags into deduplication tag
		options.addOption("l", "p5-barcode-tag", true, "Tag with p5 barcode to copy into deduplicating tag");
		options.addOption("k", "p7-barcode-tag", true, "Tag with p7 barcode to copy into deduplicating tag");
		options.addOption("b", "BAM", false, "Output bam file");
		
		CommandLine commandLine	= parser.parse( options, args );
		
		String inputFilename = commandLine.getOptionValue("input-filename");
		String outputFilename = commandLine.getOptionValue("output-filename");
		boolean useBAM = commandLine.hasOption("BAM") || Driver.isBAMFilename(outputFilename);
		
		SAMFileWriter output = null;
		SamInputResource bufferedSAMFile = SamInputResource.of(new BufferedInputStream(new FileInputStream(inputFilename)));
		try(
				SamReader reader = SamReaderFactory.makeDefault().open(bufferedSAMFile);
				){
			SAMFileHeader header = reader.getFileHeader();
			
			SAMFileWriterFactory outputFileFactory = new SAMFileWriterFactory();
			BufferedOutputStream outputFile = new BufferedOutputStream(new FileOutputStream(outputFilename));
			if(useBAM){
				output = outputFileFactory.makeBAMWriter(header, true, outputFile);
			} else {
				output = outputFileFactory.makeSAMWriter(header, true, outputFile);
			}
			
			String p5_barcode_tag = commandLine.getOptionValue("p5-barcode-tag");
			String p7_barcode_tag = commandLine.getOptionValue("p7-barcode-tag");

			SAMRecordIterator i = reader.iterator();
			while(i.hasNext()){
				// iterate through alignments
				SAMRecord record = i.next();
				StringBuilder builder = new StringBuilder();
				if (p5_barcode_tag != null) {
					String p5_barcode = record.getStringAttribute(p5_barcode_tag);
					builder.append(p5_barcode);
				}
				builder.append(IndexAndBarcodeKey.FIELD_SEPARATOR);
				if (p7_barcode_tag != null) {
					String p7_barcode = record.getStringAttribute(p7_barcode_tag);
					builder.append(p7_barcode);
				}
				builder.append(IndexAndBarcodeKey.FIELD_SEPARATOR);
				builder.append(record.getReadLength());
				String deduplicationCriterion = builder.toString();
				record.setAttribute(DemultiplexSAM.duplicatesSAMTag, deduplicationCriterion);
				output.addAlignment(record);
			}
			
		}
		finally {
			if (output != null)
				output.close();
		}
	}
}
