package adnascreen;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.HashMap;
import java.util.Map;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

/**
 * Take a SAM/BAM file produced with no read group information for a single Illumina run
 * Add one read group for each lane to assign each read a read group
 * Sample ID is required
 * 
 * @author Matthew Mah
 *
 */
public class AssignReadGroups {
	private static final SimpleDateFormat dateFormat = new SimpleDateFormat("yyyyMMdd");
	
	public static void main(String [] args) throws ParseException, IOException, java.text.ParseException{
		CommandLineParser parser = new DefaultParser();
		Options options = new Options();
		options.addRequiredOption("i", "input-filename", true, "input SAM/BAM filename");
		options.addRequiredOption("o", "output-filename", true, "input SAM/BAM filename");
		options.addRequiredOption("s", "sample", true, "sample ID");
		options.addOption("b", "BAM", false, "Output bam file");
		options.addOption("x", "label", true, "");
		options.addRequiredOption("d", "date", true, "sequencing date in YYYYMMDD format");
		options.addOption("l", "library", true, "Library");
		options.addOption("q", "sequencing-center", true, "sequencing center producing read");
		options.addOption("p", "sequencing-platform", true, "sequencing platform producing read");
		
		CommandLine commandLine	= parser.parse( options, args );
		
		String filename = commandLine.getOptionValue("input-filename");
		String outputFilename = commandLine.getOptionValue("output-filename");
		String sampleID = commandLine.getOptionValue("sample");
		boolean useBAM = commandLine.hasOption("BAM") || Driver.isBAMFilename(outputFilename);
		String label = commandLine.getOptionValue("label", "");
		String dateString = commandLine.getOptionValue("date");
		Date date = dateFormat.parse(dateString);
		String library = commandLine.getOptionValue("library", null);
		String sequencingCenter = commandLine.getOptionValue("sequencing-center", null);
		String sequencingPlatform = commandLine.getOptionValue("sequencing-platform", "ILLUMINA");
		
		SAMFileWriterFactory outputFileFactory = new SAMFileWriterFactory();
		// key format is platform unit: flowcellID.lane
		// this will be unique because flowcell IDs are unique
		Map<String, SAMReadGroupRecord> readGroups = new HashMap<String, SAMReadGroupRecord>();
		
		// make two passes through the file
		// on first pass, discover which read groups we need
		SamInputResource bufferedSAMFile = SamInputResource.of(new BufferedInputStream(new FileInputStream(filename)));
		SAMFileHeader header;
		try(
				SamReader reader = SamReaderFactory.makeDefault().open(bufferedSAMFile);
		){
			header = reader.getFileHeader();
			SAMRecordIterator i = reader.iterator();
			while(i.hasNext()){
				SAMRecord record = i.next();
				String readName = record.getReadName();
				// parse read group information from read name
				String[] fields = readName.split(":");
				
				String instrument = fields[0];
				int runNumber = Integer.valueOf(fields[1]);
				String flowcellID = fields[2];
				int lane = Integer.valueOf(fields[3]);

				// generate a read group, if there is not one present already
				// use shortened platform unit
				String platformUnit = String.format("%s.%d", flowcellID, lane);
				if(!readGroups.containsKey(platformUnit)){
					// combine with label and data
					String readGroupID = assembleReadGroupID(label, date, flowcellID, lane);
					
					SAMReadGroupRecord group = new SAMReadGroupRecord(readGroupID);
					if(sequencingCenter != null)
						group.setSequencingCenter(sequencingCenter);
					group.setRunDate(date);
					if(library != null)
						group.setLibrary(library);
					group.setPlatform(sequencingPlatform);
					group.setPlatformModel(String.format("%s %d", instrument, runNumber));
					group.setPlatformUnit(platformUnit);
					group.setSample(sampleID);
					
					readGroups.put(platformUnit, group);
				}				
			}
		}
		// on second pass, write to file with read groups
		// first, add read groups to header
		for(SAMReadGroupRecord readGroup : readGroups.values() ){
			header.addReadGroup(readGroup);
		}
		bufferedSAMFile = SamInputResource.of(new BufferedInputStream(new FileInputStream(filename)));
		try(
				SamReader reader = SamReaderFactory.makeDefault().open(bufferedSAMFile);
		){
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
				// shorten read names by removing info common to read group (instrument, runNumber, flowcellID, lane)
				String readName = record.getReadName();
				// parse read group information from read name
				String[] fields = readName.split(":");
				
				String flowcellID = fields[2];
				int lane = Integer.valueOf(fields[3]);
				int tile = Integer.valueOf(fields[4]);
				int x = Integer.valueOf(fields[5]);
				int y = Integer.valueOf(fields[6]);
				
				String shortenedReadName = String.format("%d:%d:%d", tile, x, y);
				record.setReadName(shortenedReadName);
				
				// set read group
				String readGroupID = assembleReadGroupID(label, date, flowcellID, lane);
				record.setAttribute(SAMTag.RG.toString(), readGroupID);
				// write to file with read group
				output.addAlignment(record);
			}
			output.close();
		}
	}
	
	public static String assembleReadGroupID(String label, Date date, String flowcellID, int lane) {
		String readGroupID = String.format("%s_%s_%s_%d", label, dateFormat.format(date), flowcellID, lane);
		return readGroupID;
	}
}
