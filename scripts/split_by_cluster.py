#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  4 15:21:29 2021
@author: diego
"""

import pandas as pd
import subprocess
import os


def trimm_cell_name(cell):
    cell = cell.split("/")[-1]
    cell = cell.split(".")[0]
    return cell


if __name__ == "__main__":

    
    cluster_file = snakemake.input[0]
    log_output = snakemake.output[0]
    cells_path = snakemake.params[0]
    threads = snakemake.params[1]
    output_path = "/".join(log_output.split("/")[:-1]) + "/"
        
    '''
    cluster_file = "../results_test/iter1/cells_clusters.txt"
    log_output = "../results_test/iter1/cells_clusters.log"    
    cells_path = "../results_test/cells/"
    threads = 10
    output_path = "/".join(log_output.split("/")[:-1]) + "/"
    
    #cluster_file = "iter2/cells_clusters.txt"
    #log_output = "iter2/cells_clusters.log"    
    #cells_path = "cells/"        
    #threads = 10
    #output_path = " /".join(log_output.split("/")[:-1]) + "/"
        
    print(cluster_file)
    print(log_output)
    print(cells_path)
    print(output_path)
    '''
    clusters = pd.read_csv(cluster_file, sep = "\t", index_col = 0, header = None)    
    clusters.index = list(map(trimm_cell_name, clusters.index))
    clusters.index = cells_path + clusters.index + ".bam"  
    out_message = ""
    #print(clusters)
    
    for i in set(clusters[4].values):        
        tmp_bams = []
        #chunk = 100 # merge every 100 cells
        #chunk = 50 # merge every 50 cells        
        
        start = 0
        fold_clusters = clusters.loc[clusters[4] == i]
        length_clusters = len(fold_clusters.index)        
        #print(fold_clusters)
        ratio = 0.1 # ratio
        chunk = round(length_clusters * ratio)
        #for j in range(chunk):

        for j in range(round(length_clusters/chunk)):
            #print(j)        
            if length_clusters <= chunk:
                bam_file = output_path + str(i) + "_" + str(j) + ".bam"                
                tmp_fold_clusters = fold_clusters.iloc[start:]
            else:
                bam_file = output_path + str(i) + "_" + str(j) + ".bam"                
                tmp_fold_clusters = fold_clusters.iloc[start:start + chunk]            
            length_clusters -= chunk
            start += chunk
            tmp_bams.append(bam_file)
            bam_files = " ".join(tmp_fold_clusters.index)                        
            cmd = "samtools merge -@{} -c {} {}".format(threads, bam_file, bam_files)            
            subprocess.call(cmd, shell=True)    
        out_message += "Merging bam files for cluster {} \n".format(str(i))
        bam_file = output_path + str(i) + ".bam"        
        bam_files = " ".join(tmp_bams)                
        cmd = "samtools merge -@{} -c {} {}".format(threads, bam_file, bam_files)
        subprocess.call(cmd, shell=True)
        cmd = "samtools index -@{} {}".format(threads, bam_file)
        subprocess.call(cmd, shell=True)
        if len(tmp_bams) > 0:
            cmd = "rm " + bam_files
            subprocess.call(cmd, shell=True)    
    log_file = open(log_output, "w")
    log_file.write(out_message)
    log_file.close()
    