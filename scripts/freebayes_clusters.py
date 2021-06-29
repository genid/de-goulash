#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  4 15:21:29 2021

@author: diego
"""
import subprocess
import os
import numpy as np

if __name__ == "__main__":                

    path_clusters   = snakemake.input[0]
    path_clusters   = "/".join(path_clusters.split("/")[:-1]) + "/"
    merge_vcf       = snakemake.output[0]
    ref_genome      = snakemake.params[0]
    regions         = snakemake.params[1]
    threads         = snakemake.params[2]    
    vcf_list        = []
    bam_files       = [path_clusters + bam for bam in os.listdir(path_clusters) if bam.endswith(".bam")]    
    if len(bam_files) > 0:
        for bam_file in bam_files:
            vcf_file = bam_file + ".vcf"
            vcf_list.append(vcf_file)
            cmd = "./scripts/freebayes-parallel.sh {} {} -f {} {} -C 2 -iXu > {}".format(
                                regions, threads, ref_genome, bam_file, vcf_file)            
            subprocess.call(cmd, shell = True)        
        if len(vcf_list) > 0:
            args_input = ""
            for vcf in vcf_list:
                args_input += "I="+vcf+" "        
            #cmd = "PicardCommandLine MergeVcfs {} O={}".format(args_input, merge_vcf)    
            cmd = "java -jar ./software/picard.jar MergeVcfs {} O={}".format(args_input, merge_vcf)    
            subprocess.call(cmd, shell = True)
