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


def create_tmp_dirs(folder):

    flag = True
    if os.path.isdir(folder):    
        while(flag):
            print("WARNING! File "+folder+" already exists, \
                  \nWould you like to remove it?")
            choice = input("y/n: ")            
            if str(choice) == "y":                
                cmd = 'rm -r '+folder
                subprocess.call(cmd, shell=True)
                cmd = 'mkdir '+folder
                subprocess.call(cmd, shell=True)                
                flag = False
                return True                
            elif str(choice) == "n":                
                flag = False 
                return False                                  
            else:
                print("Please type y or n")                               
    else:
        cmd = 'mkdir '+folder
        subprocess.call(cmd, shell=True)        
        return True
    

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
    
    input_cells_vcf = snakemake.input[0]   
    input_exome_vcf = snakemake.input[1]    
    vcf_cell_noheader = snakemake.output[0]   
    output_folder = os.path.dirname(vcf_cell_noheader) + "/"

    #if not create_tmp_dirs(output_folder):
    #    os.sys.exit()

    #gz_input_cells = output_folder + Path(input_cells_vcf).stem + ".vcf.gz" # change extension
    gz_input_cells = output_folder + Path(input_cells_vcf).stem + ".vcf.gz" # change extension
    f = open(gz_input_cells, "w") # change name of file per sample with highes match 
    out = subprocess.call(["./software/bgzip", "-c", input_cells_vcf], stdout = f)
    out = subprocess.call(["./software/tabix", "-fp", "vcf", gz_input_cells])
    
    gz_input_cells_out = gz_input_cells.replace(".vcf.gz", ".filt.vcf")
    #input_exome_out = output_folder + input_exome_vcf.replace(".vcf.gz", ".filt.vcf")
    input_exome_out = input_exome_vcf.replace(".vcf.gz", ".filt.vcf")
    
    vcf_cells_file = bcftools_filter(gz_input_cells, input_exome_vcf, gz_input_cells_out)
    vcf_exome_file = bcftools_filter(input_exome_vcf, gz_input_cells, input_exome_out)
    
    ## Intersect of positions given the genotype and get the sample with highest match
    vcf_cells = vcf.Reader(compressed = False, filename = vcf_cells_file)
    #vcf_cells = vcf.Reader(compressed = False, filename = input_cells)
    vcf_exome = {}
    for exome_record in vcf.Reader(compressed = False, filename = vcf_exome_file):    
    #for exome_record in vcf.Reader(compressed = True, filename = input_exome):    
        vcf_exome[str(exome_record.CHROM) + "_" + str(exome_record.POS)] =\
                                                    exome_record.samples
    count = 0
    remove_pos = [] # stores positions that does not appear in the exome
    matches = {}
    for cell_record in vcf_cells:
        cell_id = str(cell_record.CHROM) + "_" + str(cell_record.POS)
        if cell_id in vcf_exome:
            count += 1
            cell_GT = cell_record.samples[0].data.GT                
            for sample_exome in vcf_exome[cell_id]: #loop throgh the samples
                #print(sample_exome)
                if cell_GT == sample_exome.data.GT:                
                    if sample_exome.sample in matches:                        
                        matches[sample_exome.sample] += 1
                    else:
                        matches[sample_exome.sample] = 1        
        else: # not found in dictionary so lets store and delete
            remove_pos.append((cell_record.CHROM, cell_record.POS))
            count = count
    # converts positions to delete as string
    remove_list = ""
    for i in range(len(remove_pos)):
        remove_list += str(remove_pos[i][0]) + "\t" + str(remove_pos[i][1]) + "\n"  
    remove_file = output_folder + "remove_list.txt"
    with open(remove_file, "w") as f:
        f.write(remove_list)
    matches_sort = sorted(matches.items(), key = lambda kv: kv[1], reverse = True)
    flag = False
    if len(matches_sort) >= 2: # two identical samples
        if matches_sort[0][1] == matches_sort[1][1]:
            flag = True
    
    if flag:
        print("{}".format("*"*20))
        output = "{}".format("Something unusual happened, there are two identical \
                             references, please manually checked")
        print(output)
        
    vcf_cell_filt = vcf_cell_noheader
        
    f = open(vcf_cell_filt, "w") # change name of file per sample with highes match 
    out = subprocess.call(["grep", "-vf", remove_file, vcf_cells_file], 
                              stdout = f)
        
    # generate out log file given match and sample name used in the vcf
    output = "Cell\tSample_REF\tSimilarity%\tMatches\tTotal\n"                        
    for match in matches_sort:        
        output += "{}\t{}\t{}\t{}\t{}\n".format(cell_record.samples[0].sample,
                                    match[0], round(match[1] / count,3),
                                    match[1], count)        
    with open (output_folder + "matches.tsv", "w") as f:
        f.write(output)
    
    '''
    if flag:
        print("{}".format("*"*20))
        output = "{}".format("Something unusual happened, there are two identical \
                             references, please manually checked")
        print(output)
        os.exit()
    else:
        #vcf_cell_filt = output_folder + matches_sort[0][0] + ".vcf"
        #vcf_cell_filt = output_folder + "tmp.vcf"
        vcf_cell_filt = vcf_cell_noheader
        
        f = open(vcf_cell_filt, "w") # change name of file per sample with highes match 
        out = subprocess.call(["grep", "-vf", remove_file, vcf_cells_file], 
                              stdout = f)
        
        # generate out log file given match and sample name used in the vcf
        output = "Cell\tSample_REF\tSimilarity%\tMatches\tTotal\n"                        
        for match in matches_sort:        
            output += "{}\t{}\t{}\t{}\t{}\n".format(cell_record.samples[0].sample,
                                       match[0], round(match[1] / count,3),
                                       match[1], count)        
        with open (output_folder + "matches.tsv", "w") as f:
            f.write(output)        
    
        #cmd = "egrep -v '^##' " + vcf_cell_filt + " > " + vcf_cell_noheader
        #subprocess.call(cmd, shell=True)
    '''
    
    
    
