This project is a collection of tools for bioinformatic processing of ancient DNA data at the [Reich Lab at Harvard Medical School](https://reich.hms.harvard.edu). 

Originally written for the screening process, this project is now used for general processing of ancient DNA data from Illumina NextSeq. It performs tasks including:

- Count and analyze barcodes and indices for determining contents of a sequencing run. Given an input set of barcodes and indices, these counts are used to infer samples omitted from sample sheets present in sequencing runs.
- With the expectation of short ancient DNA molecules, merge paired forward and reverse reads into a single sequence. This step also implicitly trims adapters. 
- Write single or double barcodes into fastq headers that are preserved by bwa in header
- Demultiplex aligned SAM/BAM files based on barcodes in headers. Demultiplexing after aligning allows more balanced load sharing. The insert size is not an issue for ancient DNA because the DNA fragments are small and paired forward and reverse reads are expected to overlap. 
- Soft clip SAM/BAM files to remove deamination damage. 
- Assign read groups

---

SOFTWARE COPYRIGHT NOTICE AGREEMENT
This software and its documentation are copyright (2018) by Harvard University 
and The Broad Institute. All rights are reserved. This software is supplied 
without any warranty or guaranteed support whatsoever. Neither Harvard 
University nor The Broad Institute can be responsible for its use, misuse, or 
functionality. The software may be freely copied for non-commercial purposes, 
provided this copyright notice is retained.
