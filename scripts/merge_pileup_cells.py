#!/usr/bin/env python
# coding: utf-8


import os
import subprocess
import numpy as np
import pandas as pd

def get_frequency_table(mpileup):    
    
    bases = ["A","T","G","C","+","-"]    
    frequency_table = {}    
    for i in mpileup.values:       
        fastadict = {"A":0,"T":0,"G":0,"C":0}                    
        sequence = i[4] #actual sequence                
        sequence = sequence.upper()         
        sequence = trimm_caret(sequence)                            
        sequence = sequence.replace("$", "")                   
        indel_pos = find_all_indels(sequence)
        ### Count number of indels
        indels = count_indels(sequence, indel_pos)        
        fastadict.update(indels)
        fastadict["-"] += sequence.count("*")
        ### Trimm Indels
        trimm_sequence = trimm_indels(sequence, indel_pos)        
        for seq in trimm_sequence:    
            if seq in fastadict:            
                fastadict[seq] +=1    
        #frequency_table.update({i[1]:list(fastadict.values())})    
        frequency_table.update({"chr" + str(i[0]) + "_" + str(i[1]) :list(fastadict.values())})    
    df_frequency_table = pd.DataFrame.from_dict(frequency_table, orient='index')
    df_frequency_table.columns = bases          
    
    return df_frequency_table


def find_all_indels(s):
    find_all = lambda c,s: [x for x in range(c.find(s), len(c)) if c[x] == s]
    list_pos = []
    for i in find_all(s,"-"):
        list_pos.append(i)
    for i in find_all(s,"+"):
        list_pos.append(i)    
    return sorted(list_pos)


def count_indels(s, pos):    
    dict_indel = {"+":0,"-":0}    
    if pos == []:
        return dict_indel            
    if len(pos) > 0:
        for i in range(0,len(pos)): 
            try: # in case it is not a number but a base pair e.g. A
                dict_indel[s[pos[i]]] += int(s[pos[i]+1])                                                                        
            except ValueError:                
                dict_indel[s[pos[i]]] += 1
                continue                
    return dict_indel


def trimm_indels(s, pos):    
    ## Receives a sequence and trimms indels            
    if pos == []:
        return s
    u_sequence = ""  
    start =  pos[0]
    count = (start+1)    
    try: # in case it is not a number but a base pair e.g. A
        end = count+int(s[count])+1
    except ValueError:                
        end = start+1            
    u_sequence = s[:start]    
    if len(pos) > 1:
        for i in range(1,len(pos)):                      
            start = end                
            u_sequence += s[start:pos[i]]
            start = pos[i]
            count = (start+1)            
            try: # in case it is not a number but a base pair e.g. A
                end = count+int(s[count])+1  
            except ValueError:
                end = start+1                
            if pos[-1] == pos[i]:    
                #print(s[end:])
                u_sequence += s[end:]
    else:        
        u_sequence += s[end:]        
    return u_sequence


def trimm_caret(s):           
    find_all = lambda c,s: [x for x in range(c.find(s), len(c)) if c[x] == s]
    list_pos = []
    for i in find_all(s,"^"):
        list_pos.append(i)    
    if list_pos == []:
        return s
    i = 0    
    start = 0
    end = 0
    sequence = ""
    while i < len(s):
        if s[i] == "^":        
            end = i
            sequence += (s[start:end])                    
            start = i+1
        elif i >= list_pos[-1]+1:
            sequence += (s[list_pos[-1]+1:])
            break
        i+=1        
    return sequence


def parse_vcf(vcf_file):    
    tmp_file = "tmp.txt"
    cmd = "grep -v '^#' "+vcf_file+ "> "+tmp_file
    subprocess.call(cmd, shell=True)
    if os.stat(tmp_file).st_size > 0:    
        df_vcf = pd.read_csv(tmp_file, sep="\t", header=None, low_memory=False)            
        chrom = df_vcf[0].values
        pos = df_vcf[1].values
        indexes = ["chr{}_{}".format(b_, a_) for a_, b_ in zip(pos, chrom)]    
        df_vcf["key"] = indexes    
        return df_vcf
    else:
        return pd.DataFrame()


def get_genotype(df_genotypes):                
    ref = df_genotypes[3].values
    dict_values = {}
    values = []
    for geno in df_genotypes.values:                
        if '0/0' in geno[9]:                  
            values.append(-1)
        elif '1/1' in geno[9]:                    
            values.append(1)
        elif '0/1' in geno[9]:                    
            #values.append(0)
            values.append(np.nan)        
        elif './.' in geno[9]: 
            values.append(np.nan)        
        dict_values[geno[-1]] = values[-1]
    return dict_values, ref


def get_informative_columns(df_all):
    # generate list with informative position at least 2 different genotypes
    columns_filter = []
    i = 0
    for column in df_all:    
        #list_values = [x for x in df_all[column].values if str(x) != 'nan' and str(x) != 'NaN']    
        list_values = [x for x in df_all[column].values if str(x) != 'nan']                
        if len(set(list_values)) > 1:
            #print(set(list_values))
            columns_filter.append(column)                
        #if i % 500000 == 0:
        #    print("{} {} ".format(i, len(columns_filter)))        
        i+=1
    print(len(columns_filter))
    return columns_filter


def filter_sample_name(s):
    return s.split("_")[0]


def trimm_names(name):
    name = name.split("/")[-1]
    name = name.split(".")[0]
    return name


def check_extension(filename,ext):
    flag = False
    for x in ext:
        if filename.endswith(x):
            flag = True
    return flag   


def check_if_folder(path,ext):    
    list_files = []
    if os.path.isdir(path):
        dirpath = os.walk(path)
        for dirpath, dirnames, filenames in dirpath:
            for filename in [f for f in filenames if check_extension(f, ext)]:
                files = os.path.join(dirpath, filename)
                list_files.append(files)
        return list_files
    else:
        return [path]


def base_calling(pu, threshold_coverage_pos, threshold_base_calling):
        
    dict_cell = {}
    tmp = {} 
    mpileup = pd.read_csv(pu, names=["chr","pos","ref","reads","seq","qual"], sep = "\t")    
    df_freq = get_frequency_table(mpileup)
    df_freq = df_freq.drop(['+','-'], axis=1)
    
    df_freq = df_freq[np.sum(df_freq.values, axis=1) >= threshold_coverage_pos]    
    
    list_col_indices = np.argmax(df_freq.values, axis=1)
    called_base = df_freq.columns[list_col_indices]
    
    total_count_bases = np.sum(df_freq.values, axis=1)
    max_count_bases = np.max(df_freq, axis=1)
    
    called_perc = round((max_count_bases/total_count_bases)*100,1)
    called_perc = np.nan_to_num(called_perc) #nan to 0
    
    df_freq["called_perc"] = np.array(called_perc, dtype=int)
    df_freq = df_freq[df_freq["called_perc"] >= threshold_base_calling]
    
    list_col_indices = np.argmax(df_freq[df_freq.columns[0:-1]].values, axis=1)    
    called_base = df_freq.columns[list_col_indices]                        
    for i in range(len(called_base)):
        tmp[df_freq.index[i]] = called_base[i]            
    dict_cell[pu] = tmp    
    dict_cell = pd.DataFrame(dict_cell).T            
    return dict_cell


if __name__ == "__main__":

    
    path_pu_cells = snakemake.input[0]    
    output_file = snakemake.output[0]
    dirpath_pu = "/".join(path_pu_cells.split("/")[:-1]) + "/"    
    threshold_coverage = int(snakemake.params[0]) # number of reads covered in total per sample
    threshold_coverage_pos = int(snakemake.params[1])  #number of reads covered in a given position
    threshold_base_calling = int(snakemake.params[0]) # Yleaf based calling parameter        
    extension = "pu"
    files_list = check_if_folder(dirpath_pu,extension)
    ## load pileup files in a list a make the basecalling for the positions
    j = 0
    i = 0
    list_cells = []
    for pu in files_list:
        try:        
            mpileup = pd.read_csv(pu, names=["chr","pos","ref","reads","seq","qual"], sep = "\t")                
            mpileup = mpileup[mpileup["reads"] >= 1]                
            if np.sum(mpileup["reads"]) >= (threshold_coverage): 
                df_freq = get_frequency_table(mpileup)                                             
                dict_cell = base_calling(pu, threshold_coverage_pos, threshold_base_calling)                        
                if dict_cell.empty:
                    continue            
                list_cells.append(dict_cell)
                i+=1
            else:
                continue
        except KeyError as e:
            #print("Not found in: "+pu)
            pass    
        if i%500 == 0:
            print("{}: {}".format(j,i))        
            print("#############################")    
        j += 1
    print("-- Start Concat --")    
    df = pd.concat(list_cells, axis=0, sort=False)
    print("-- Finish Concat --")

    df.index = list(map(trimm_names, df.index))
    df.to_csv(output_file)

