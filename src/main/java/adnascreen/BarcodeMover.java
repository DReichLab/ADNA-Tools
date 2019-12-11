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

import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMFormatException;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class BarcodeMover {
	public static void main(String [] args) throws ParseException, IOException {
		CommandLineParser parser = new DefaultParser();
		Options options = new Options();
		options.addRequiredOption("i", "input", true, "Input file with barcodes at end of read name");
		options.addRequiredOption("o", "output", true, "Output file with barcodes and length in separate tag");
		options.addOption("t", "tag", true, "SAM tag label to use for deduplicating. Barcodes and length go in this tag.");
		
		options.addOption("b", "BAM", false, "Use bam files for output");
		options.addOption(null, "bufferSize", true, "Output file buffer size for performance");
		CommandLine commandLine	= parser.parse(options, args);
		
		String tag = commandLine.getOptionValue("tag", DemultiplexSAM.duplicatesSAMTag);
		int bufferSize = Integer.valueOf(commandLine.getOptionValue("bufferSize", "1048576")); // 1 MB
		
		String inputFilename = commandLine.getOptionValue("input");
		String outputFilename = commandLine.getOptionValue("output");
		
		SAMFileWriterFactory outputFileFactory = new SAMFileWriterFactory();
		SamInputResource bufferedSAMFile = SamInputResource.of(new BufferedInputStream(new FileInputStream(inputFilename), bufferSize));
		try(
				SamReader reader = SamReaderFactory.makeDefault().open(bufferedSAMFile);
				SAMFileWriter output = outputFileFactory.makeBAMWriter(reader.getFileHeader(), true, 
						new BufferedOutputStream(new FileOutputStream(outputFilename), bufferSize));
				){
			SAMRecordIterator i = reader.iterator();
			while(i.hasNext()){
				// iterate through alignments
				try{
					SAMRecord record = i.next();
					SAMRecord modified = barcodeToTag(record, tag);
					
					output.addAlignment(modified);
				} catch (SAMFormatException e){
					System.err.println(e);
					// ignore this record and continue to the next
				}
			}
		}
	}
	
	/**
	 * Heng's processing of shotgun data puts the barcodes into the query name
	 * Move those barcodes into a SAM/BAM tag
	 * @param hengRead
	 * @return
	 */
	public static SAMRecord barcodeToTag(SAMRecord hengRead, String tag) {
		SAMRecord modified = hengRead.deepCopy();
		
		String readName = hengRead.getReadName();
		String [] fields = readName.split(":");
		
		// Last field should be a DNA sequence with an even number of bases (P5 and P7 barcodes concatenated)
		DNASequence barcodesConcatenated = new DNASequence(fields[fields.length-1]);
		if (barcodesConcatenated.length() % 2 != 0) {
			throw new IllegalArgumentException("DNA barcode length: " + barcodesConcatenated.length());
		}
		String barcodesString = barcodesConcatenated.toString();
		int singleBarcodeLength = barcodesConcatenated.length() / 2;
		String barcode1 = barcodesString.substring(0, singleBarcodeLength);
		String barcode2 = barcodesString.substring(singleBarcodeLength);
		
		// write duplicates tag with barcodes and length
		StringBuilder builder = new StringBuilder();
		builder.append(barcode1);
		builder.append(IndexAndBarcodeKey.FIELD_SEPARATOR);
		builder.append(barcode2);
		builder.append(IndexAndBarcodeKey.FIELD_SEPARATOR);
		builder.append(hengRead.getReadLength());
		String duplicatesTagContent = builder.toString();
		
		modified.setAttribute(DemultiplexSAM.duplicatesSAMTag, duplicatesTagContent);
		
		// remove barcode from read name for space saving
		int toRemoveIndex = readName.indexOf(":" + barcodesString);
		String trimmedReadName = readName.substring(0, toRemoveIndex);
		modified.setReadName(trimmedReadName);
		
		return modified;
	}
}
