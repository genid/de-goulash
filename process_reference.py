#!/usr/bin/env python3

import os
import subprocess


if __name__ == "__main__":

    
    N_size = 500000 # 10,000 bp for region
    reference_genome = os.sys.argv[1]

    path = os.path.dirname(reference_genome)
    print(path)
    
    cmd = "samtools faidx {}".format(reference_genome)
    subprocess.check_output(cmd, shell=True)
    dict_fai = reference_genome + ".dict"
    cmd = "samtools dict {} > {}".format(reference_genome, dict_fai)
    subprocess.check_output(cmd, shell=True)
    
    MT_reference_genome = path + "/MT_genome.fasta"
    cmd = "samtools faidx {} MT > {}".format(reference_genome, MT_reference_genome)
    subprocess.check_output(cmd, shell=True)
    
    cmd = "samtools faidx {}".format(MT_reference_genome)
    subprocess.check_output(cmd, shell=True)
    MT_dict_fai = MT_reference_genome + ".dict"
    cmd = "samtools dict {} > {}".format(MT_reference_genome, MT_dict_fai)
    subprocess.check_output(cmd, shell=True)
    
    regions = path + "/regions.txt"
    MT_regions = path + "/MT_regions.txt"

    cmd = "python fasta_generate_regions.py {} {} > {}".format(reference_genome,
                                                                    N_size, regions)
    subprocess.check_output(cmd, shell=True)
    cmd = "python fasta_generate_regions.py {} {} > {}".format(MT_reference_genome,
                                                                    N_size, MT_regions)
    subprocess.check_output(cmd, shell=True)


