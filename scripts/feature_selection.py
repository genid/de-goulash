#!/usr/bin/env python
# coding: utf-8

import pandas as pd


if __name__ == "__main__":
               
    df = pd.read_csv(snakemake.input[0], index_col=0)
    output_file = snakemake.output[0]    
    thr_cell_i = int(snakemake.params[0])

    list_diffs = [list(set(df[df.columns[i]].values)) for i in range(len(df.columns))]
    # 1) Filter position by containin:  >= 3, must contain: allele A, allele B and nan (missing allele) 
    list_pos = [i for i in range(len(list_diffs)) if len(list_diffs[i]) > 2]
    df = df[df.columns[list_pos]]
    df.head()
    # 2) filter by SNPS
    snps = []
    for col in df.columns:
        # sum of missing values
        if (df[col].isna().sum() / len(df)) <= 0.99: #allows 99% of missing or less
            snps.append(col)
    df = df[snps]
    # filter cells
    thr_cell = len(snps) - thr_cell_i
    df[df.isna().sum(axis=1) <= thr_cell].shape
    df = df[df.isna().sum(axis=1) <= thr_cell]    
    print("-- Valid cells --")
    print("\t{} cells and  {} SNPs".format(df.shape[0], df.shape[1]))
    df.to_csv(output_file)
