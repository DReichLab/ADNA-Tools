package adnascreen;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMFormatException;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class DemultiplexSAM {
	public static void main(String [] args) throws IOException {
		Cache<IndexAndBarcodeKey, SAMFileWriter> outputFiles = new Cache<IndexAndBarcodeKey, SAMFileWriter>(2000);
		
		SAMFileHeader masterHeader = null;
		SAMFileWriterFactory outputFileFactory = new SAMFileWriterFactory();

		int openedOutputFiles = 0;
		long alignmentsProcessed = 0;
		// iterate through input files
		for(String filename : args){
			SamInputResource bufferedSAMFile = SamInputResource.of(new BufferedInputStream(new FileInputStream(filename)));
			try(
					SamReader reader = SamReaderFactory.makeDefault().open(bufferedSAMFile);
			){
				SAMFileHeader header = reader.getFileHeader();
				if(masterHeader == null){
					masterHeader = header;
				} else if(!masterHeader.equals(header)){
					throw new IllegalArgumentException("SAM headers do not match");
				}

				SAMRecordIterator i = reader.iterator();
				while(i.hasNext()){
					// iterate through alignments
					try{
						SAMRecord record = i.next();
						// parse out key, which is a 4-tuple of indices and barcodes
						String readname = record.getReadName();
						String [] readnameParts = readname.split(String.valueOf(MergedRead.KEY_SEPARATOR));
						IndexAndBarcodeKey key = new IndexAndBarcodeKey(readnameParts[1]);
						// remove indices and barcodes from read name
						//record.setReadName(readnameParts[0]);
						
						// find file corresponding to this key
						SAMFileWriter output = outputFiles.get(key);
						if(output == null){ // open new file, if none exists for this key
							String outputFilename = key.toString() + ".sam";
							//BufferedOutputStream outputFile = new BufferedOutputStream(new FileOutputStream(outputFilename, true));
							File outputFile = new File(outputFilename);
							openedOutputFiles++;
							output = outputFileFactory.makeSAMWriter(header, true, outputFile);
							//output = outputFileFactory.makeBAMWriter(header, true, outputFile);
							outputFiles.put(key, output); // 
						}
						// write out to alignment to file
						output.addAlignment(record);
					} catch (SAMFormatException e){
						System.err.println(e);
						// ignore this record and continue to the next
					} catch (Exception e){
						System.err.println(e);
						System.err.println("Opened output files: " + openedOutputFiles);
					}
					alignmentsProcessed++;
					/*
					if(alignmentsProcessed % 1000 == 0){
						System.err.print(alignmentsProcessed);
						System.err.print('\t');
						System.err.print(outputFiles.getHits());
						System.err.print('\t');
						System.err.print(outputFiles.getMisses());
						System.err.print('\t');
						System.err.print(outputFiles.getForcedCloses());
						System.err.print('\n');
						
					}
					*/
				}
			}
		}
		// cleanup, close all files
		outputFiles.close();
	}
}
