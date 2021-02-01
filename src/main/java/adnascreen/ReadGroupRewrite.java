package adnascreen;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.List;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMProgramRecord;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class ReadGroupRewrite {
	/**
	 * Rewrite sections of read groups for an existing bam into a new bam
	 * Usage examples:
	 * 1. cases of platforms that violate the SAM specification
	 * 2. merging bams for a sample requiring rewriting sample fields to agree
	 * @param args
	 * @throws ParseException
	 * @throws IOException
	 * @throws java.text.ParseException
	 */
	public static void main(String[] args) throws ParseException, IOException, java.text.ParseException{
		CommandLineParser parser = new DefaultParser();
		Options options = new Options();
		options.addRequiredOption("i", "input-filename", true, "input SAM/BAM filename");
		options.addRequiredOption("o", "output-filename", true, "output SAM/BAM filename");
		options.addRequiredOption("s", "sample", true, "sample ID");
		options.addOption("b", "BAM", false, "Output bam file");
		options.addOption("l", "library", true, "Library");
		options.addOption("q", "sequencing-center", true, "sequencing center producing read");
		options.addOption("p", "sequencing-platform", true, "sequencing platform producing read");

		CommandLine commandLine	= parser.parse( options, args );
		
		String filename = commandLine.getOptionValue("input-filename");
		String outputFilename = commandLine.getOptionValue("output-filename");
		String sampleID = commandLine.getOptionValue("sample");
		boolean useBAM = commandLine.hasOption("BAM") || Driver.isBAMFilename(outputFilename);
		String library = commandLine.getOptionValue("library");
		String sequencingCenter = commandLine.getOptionValue("sequencing-center");
		String sequencingPlatform = commandLine.getOptionValue("sequencing-platform");
		
		SAMFileWriterFactory outputFileFactory = new SAMFileWriterFactory();
		
		SamInputResource bufferedSAMFile = SamInputResource.of(new BufferedInputStream(new FileInputStream(filename)));
		SAMFileHeader header;
		try(
				SamReader reader = SamReaderFactory.makeDefault().open(bufferedSAMFile);
		){
			header = reader.getFileHeader();
			List<SAMReadGroupRecord> readGroups = header.getReadGroups();
			// alter read groups according to command line parameters
			for (SAMReadGroupRecord readGroup : readGroups) {
				// silently remove platform values that violate specification
				try {
					SAMReadGroupRecord.PlatformValue.valueOf(readGroup.getPlatform());
				}
				catch(IllegalArgumentException e) {
					readGroup.setPlatform(null);
				}
				
				// replace any specified options
				if(sampleID != null) {
					readGroup.setSample(sampleID);
				}
				if(library != null) {
					readGroup.setLibrary(library);
				}
				if(sequencingCenter != null) {
					readGroup.setSequencingCenter(sequencingCenter);
				}
				if(sequencingPlatform != null) {
					// check that any platform entered is valid
					SAMReadGroupRecord.PlatformValue.valueOf(sequencingPlatform);
					readGroup.setPlatform(sequencingPlatform);
				}
			}
			// add program group header
			List<SAMProgramRecord> programEntries = header.getProgramRecords();
			// ensure unique ID
			int count = 0;
			boolean found;
			String candidateID;
			do {
				count++;
				found = false;
				candidateID = Driver.PROGRAM_NAME + '.' + count;
				for(SAMProgramRecord entry : programEntries) {
					if(entry.getId().equals(candidateID)) {
						found = true;
						break;
					}
				}
			} while(found);
			
			SAMProgramRecord programEntry = new SAMProgramRecord(candidateID);
			programEntry.setProgramVersion(Driver.versionString());
			programEntry.setCommandLine("ReadGroupRewrite" + String.join(" ", args));
			programEntry.setProgramName("ReadGroupRewrite");
			// chain program entry to last entry
			if(programEntries.size() > 0)
				programEntry.setPreviousProgramGroupId(programEntries.get(programEntries.size()-1).getId());
			header.addProgramRecord(programEntry);
			
			// rewrite alignment file with new header
			SAMFileWriter output;
			BufferedOutputStream outputFile = new BufferedOutputStream(new FileOutputStream(outputFilename));
			if(useBAM){
				output = outputFileFactory.makeBAMWriter(header, false, outputFile);
			} else {
				output = outputFileFactory.makeSAMWriter(header, false, outputFile);
			}
			
			SAMRecordIterator i = reader.iterator();
			while(i.hasNext()){
				SAMRecord record = i.next();
				output.addAlignment(record);
			}
			output.close();
		}
	}

}
