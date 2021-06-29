#!/usr/bin/env python

##### Code has not been tested on unsorted bam files, sort on barcode (CB):
##### samtools sort -t CB unsorted.bam > sorted_tags.bam
###
##### INPUT: .bam file to be sorted and output directory to place split BC
##### OUTPUT: .bam file for each unique barcode, best to make a new directory

### Python 3.6.8
import pysam
import os
pysam.set_verbosity(0) # warning of missing index is excluded

### Input varibles to set


if __name__ == "__main__":

    # input vbam file containing cells to split     
    unsplit_file = snakemake.input[0] # bam file    
    # where to place output files    
    out_dir = snakemake.output[0] # path to store split cells
    # variable to hold barcode index
    CB_hold = 'unset'
    itr = 0
    # read in upsplit file and loop reads by line
    samfile = pysam.AlignmentFile( unsplit_file, "rb")
    c = 0
    lines = 1
    for read in samfile.fetch( until_eof=True):
        # barcode itr for current read
        try:
            CB_itr = read.get_tag('CB')
            #if CB_itr in barcodes:
                # if change in barcode or first line; open new file
            if( CB_itr != CB_hold or itr == 0):                
                lines += 1
                if (lines % 500) == 0:
                    print(" {} cells spearated ...".format(lines))                    
                # close previous split file, only if not first read in file
                if( itr != 0):
                    split_file.close()
                CB_hold = CB_itr
                itr+=1                                                
                out_folder = "/".join(out_dir.split("/")[:-1]) + "/"
                split_file = pysam.AlignmentFile( out_folder +  "{}.bam".format( CB_itr), "wb", template=samfile)
            split_file.write( read)
        except KeyError:
            c += 1
            pass
    split_file.close()
    samfile.close()
        # write read with same barcode to file
    log_message = "\t-- Warning: {} reads without CB tag --".format(c)
    print(log_message)
    log_file = open(out_dir, "w")
    log_file.write(log_message)
    log_file.close()