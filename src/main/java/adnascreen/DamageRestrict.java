package adnascreen;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
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
import htsjdk.samtools.cram.ref.CRAMReferenceSource;
import htsjdk.samtools.cram.ref.ReferenceSource;

public class DamageRestrict {
	public static int main(String[] args) throws IOException{
		CommandLineParser parser = new DefaultParser();
		
		Options options = new Options();
		options.addRequiredOption("i", "input", true, "SAM/BAM/CRAM to be filtered for damage score, applied by Nick Patterson's damagescore program");
		options.addRequiredOption("o", "output", true, "SAM/BAM/CRAM filtered for damage score.");
		options.addRequiredOption("d", "damage", true, "Damage threshold > for retaining read");
		options.addOption("r", "reference", true, "CRAM reference: if any cram files are used, the reference must be specified and the same for all CRAM files.");
		options.addOption("t", "tag", true, "SAM tag to use, default: 'ds'");
		options.addOption("c", "compression", true, "HTSJDK compression parameter for BAM/CRAM: 0=none, 9=max, default 5");
		
		CommandLine commandLine;
		try {
			commandLine	= parser.parse(options, args);
		
			String inputFilename = commandLine.getOptionValue("input");
			String outputFilename = commandLine.getOptionValue("output");
			float damageThreshold = Float.valueOf(commandLine.getOptionValue("damage"));
			String damageTag = "ds";
			if (commandLine.hasOption("tag")){
				damageTag = commandLine.getOptionValue("tag");
			}
			int compression = Integer.valueOf(commandLine.getOptionValue("compression", "5"));
			
			// if there is a CRAM file, we need to make sure there is a reference
			boolean isCramPresent = AlignmentComparison.isCramFilename(inputFilename) ||  AlignmentComparison.isCramFilename(outputFilename);
			if (isCramPresent && !commandLine.hasOption("reference")) {
				throw new RuntimeException("Missing CRAM reference file");
			}
			String cramReference = null;
			if (commandLine.hasOption("reference")) {
				cramReference = commandLine.getOptionValue("reference");
			}
			filterDamageScore(inputFilename, outputFilename, damageTag, damageThreshold, cramReference, compression);
			return 0;
		} catch(ParseException e) {
			System.err.println(e.getMessage());
			HelpFormatter formatter = new HelpFormatter();
			formatter.printHelp("DamageRestrict", options);
			return -1;
		}
	}
	
	public static void filterDamageScore(String inputFilename, String outputFilename, String damageTag, float damageThreshold, String reference, int compression) throws IOException {
		SamReader inputReader = null;
		SamReaderFactory samReaderFactory = SamReaderFactory.makeDefault();
		CRAMReferenceSource cramReference = null;
		if (reference != null) {
			Path referencePath = Paths.get(reference);
			cramReference = new ReferenceSource(referencePath);
			samReaderFactory.referenceSource(cramReference);
		}
		if(AlignmentComparison.isCramFilename(inputFilename)) {
			inputReader = SamReaderFactory.makeDefault().referenceSource(cramReference).open(SamInputResource.of(new BufferedInputStream(new FileInputStream(inputFilename))));
		} else {
			inputReader = samReaderFactory.open(SamInputResource.of(new BufferedInputStream(new FileInputStream(inputFilename))));
		}
		SAMRecordIterator inputIterator = inputReader.iterator();
		SAMFileHeader header = inputReader.getFileHeader();
		
		SAMFileWriter output;
		SAMFileWriterFactory outputFileFactory = new SAMFileWriterFactory();
		outputFileFactory.setCompressionLevel(compression);
		BufferedOutputStream outputFile = new BufferedOutputStream(new FileOutputStream(outputFilename));
		if(outputFilename.toLowerCase().endsWith(".cram")){
			output = outputFileFactory.makeCRAMWriter(header, outputFile, new java.io.File(reference));
		} else if(outputFilename.toLowerCase().endsWith(".bam")) {
			output = outputFileFactory.makeBAMWriter(header, true, outputFile);
		} else {
			output = outputFileFactory.makeSAMWriter(header, true, outputFile);
		}
		
		while(inputIterator.hasNext()){
			// iterate through alignments
			SAMRecord record = null;
			try{
				record = inputIterator.next();
				if (record.hasAttribute(damageTag) && record.getFloatAttribute(damageTag) > damageThreshold) {
					output.addAlignment(record);
				}
			}
			catch(Exception e){
				System.err.println(e.toString());
				if(record != null) {
					System.err.println(record.toString());
				}
			}
		}
		output.close();
	}
}
