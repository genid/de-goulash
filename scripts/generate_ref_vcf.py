#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  4 11:40:35 2021

@author: diego
"""

import os
import vcf
import subprocess
from pathlib import Path


def bcftools_filter(input_vcf, filter_vcf, output_vcf):
    """
    Parameters
    ----------
    input_vcf : FILE
        DESCRIPTION.
    filter_vcf : FILE
        DESCRIPTION.
    output_vcf : FILE
        return vcf positions from input_vcf where match in filter_vcf.

    Returns
    -------
    None.
    """        
    #output_vcf = input_vcf.replace(".vcf.gz",".filt.vcf")    
    #output_vcf = Path(input_vcf).stem + ".filt.vcf" # change extension
    f = open(output_vcf, "w") # change name of file per sample with highes match 
    out = subprocess.call(["./software/bcftools/bcftools", "filter", "-R", filter_vcf, input_vcf], 
                                                                  stdout = f)
    return output_vcf

if __name__ == "__main__":
        
    vcf_cell_filt = snakemake.input[0]
    vcf_ref_output = snakemake.output[0]
    dirpath_gen = snakemake.params[0]
    output_folder = os.path.dirname(vcf_ref_output) + "/"
            
    #vcf_cell_filt = "output/2/2.sample.vcf"    
    #vcf_ref_output = "output/2/2.ref.vcf"    
    #dirpath_gen = "/media/disk1/diego/git/Single-cell/input/1000G/"
    #output_folder = os.path.dirname(vcf_ref_output) + "/"

    vcf_cell = {}
    for cell_record in vcf.Reader(compressed=False, filename = vcf_cell_filt):    
        if cell_record.CHROM in vcf_cell:
            vcf_cell[cell_record.CHROM].append(cell_record.POS)        
        else:
            vcf_cell[cell_record.CHROM] = [cell_record.POS]
    
    chr_files = [file for file in os.listdir(dirpath_gen) if "genotypes" in file and \
                 file.endswith("gz")]
    for chr_, positions in vcf_cell.items():
        chr_look = "chr" + chr_ + "_"
        for chr_file in chr_files:
            if chr_look in chr_file:
                chr_file = dirpath_gen + chr_file
                chr_file_filt = output_folder + chr_look + ".vcf"      
                bcftools_filter(chr_file, vcf_cell_filt, chr_file_filt)                        
                            
                chr_file_filt_gz = output_folder + Path(chr_file_filt).stem + \
                                                    ".vcf.gz" # change extension            
                subprocess.check_output("./software/bgzip -c "+ chr_file_filt +" > "\
                                        +chr_file_filt_gz, shell=True)
                subprocess.check_output("./software/tabix -fp vcf " + chr_file_filt_gz, 
                                        shell=True)
    ## Concat filtered genotype by chromosome
    #out_1000G = output_folder + "REF_1000G.vcf"  
    out_1000G = vcf_ref_output
    cmd = "./software/bcftools/bcftools concat "+ output_folder +"*chr*.vcf.gz -o " + out_1000G    
    subprocess.call(cmd, shell=True)
    
    #vcf_ref_noheader = out_1000G.replace(".vcf", "-noheader.vcf")
    
    #cmd = "egrep -v '^##' " + out_1000G + " > " + vcf_ref_noheader
    #subprocess.call(cmd, shell=True)
    #vcf_cell_noheader = vcf_cell_filt.replace(".vcf","-noheader.vcf")
