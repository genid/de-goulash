#!/usr/bin/env python

import os
import subprocess
import multiprocessing
from subprocess import Popen, PIPE
from functools import partial


def execute_mpileup(vcf_pos, out_folder, bam_file):
    mpileup = bam_file.split("/")[-1]
    mpileup = out_folder + mpileup.replace(".bam",".pu")        
    cmd = "samtools mpileup -d 1000 -l {} {} > {}".format(vcf_pos, bam_file, mpileup)
    #subprocess.run(cmd, shell=True) # default subprocess call
    timer_out = subprocess.Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE, universal_newlines=True)
    stout, sterror = timer_out.communicate()

if __name__ == "__main__":
    
    vcf_pos   = snakemake.input[0]
    log_file  = snakemake.output[0]
    dirpath   = snakemake.params[0]

    dirpath = os.walk(dirpath)
    samples = []
    for dirpath, dirnames, filenames in dirpath:
        for filename in [f for f in filenames if f.endswith(".bam")]:            
            path_file = os.path.join(dirpath, filename)                          
            samples.append(path_file)  
    out_folder = "/".join(log_file.split("/")[:-1]) + "/"    
    func = partial(execute_mpileup, vcf_pos, out_folder)    
    pool = multiprocessing.Pool()                                                                              
    pool.map(func, samples)
    pool.close()
    pool.join()
    tmp = open(log_file, "w")
    tmp.write("-- pileups generated --")
    tmp.close()
