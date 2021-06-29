## Snakemake start
import os
import subprocess
configfile: "config.yaml"

sample_bam=config["sample_bam"]
dirpath=os.sys.argv[1]
dirpath=dirpath.split("/")[0]
config["results"]=dirpath

rule reduced_cells:
    input:        
        barcodes=config["barcodes"],
        bam_file=expand("{sample_bam}", sample_bam=sample_bam)
    output:                
        temp(config["results"]+"/data/{sample_bam}_possorted.bam")
    params:
        log_level="debug",
        cores=config["cores"]
    log:
        config["results"]+"/logs/{sample_bam}_reduced_cells.log"
    shell:
        "(software/subset-bam_linux --bam {input.bam_file}"
        " --cell-barcodes {input.barcodes} --out-bam {output} --cores {params.cores} --log-level {params.log_level}) 2> {log}"


rule samtools_index:
    input:        
        config["results"]+"/data/{sample_bam}_possorted.bam"
    output:
        config["results"]+"/data/{sample_bam}_possorted.bam.bai"
    log:
        config["results"]+"/logs/{sample_bam}_samtools_index.log"
    shell:
        "(samtools index {input}) 2> {log}"


rule samtools_sort_tag:
    input:
        bam=config["results"]+"/data/{sample_bam}_possorted.bam",
        bai=config["results"]+"/data/{sample_bam}_possorted.bam.bai",        
    output:
        config["results"]+"/data/{sample_bam}_sorted_tags.bam"
    threads: 
        config["cores"]
    log:
        config["results"]+"/logs/{sample_bam}_samtools_sort_tag.log"
    shell:
        "(samtools sort -@{threads} -t CB {input.bam} > {output}) 2> {log}"
        

rule split_cells:
    input:
        bam=config["results"]+"/data/{sample_bam}_sorted_tags.bam"
    output:
        config["results"]+"/cells/{sample_bam}.log"
    log:
        config["results"]+"/logs/{sample_bam}_split_cells.log"
    script:
        "scripts/split_cells.py"


rule samtools_extract_bam_MT:
    input:
        bam_file=config["results"]+"/data/{sample_bam}_possorted.bam",
        bai=config["results"]+"/data/{sample_bam}_possorted.bam.bai",
        log_split=config["results"]+"/cells/{sample_bam}.log"
    output:
        config["results"]+"/data/{sample_bam}_MT.bam"
    log:
        config["results"]+"/logs/{sample_bam}_samtools_extract_bam_MT.log"
    shell:
        "(samtools view -b {input.bam_file} MT > {output}) 2> {log}"
        

rule samtools_index_bam_MT:
    input:
        bam=config["results"]+"/data/{sample_bam}_MT.bam"
    output:
        bai=config["results"]+"/data/{sample_bam}_MT.bam.bai" 
    log:
        config["results"]+"/logs/{sample_bam}_samtools_index_bam_MT.log"
    shell:
        "(samtools index {input.bam}) 2> {log}"
                

rule samtools_extract_ref_MT:
    input:                            
        ref=config["reference"],
        log=config["results"]+"/cells/{sample_bam}.log"
    output:
        config["results"]+"/data/{sample_bam}_MT.fasta"
    log:
        config["results"]+"/logs/{sample_bam}_samtools_extract_ref_MT.log"
    shell:
        "(samtools faidx {input.ref} MT > {output}) 2> {log}"


rule samtools_faidx_ref_MT:
    input:   
        config["results"]+"/data/{sample_bam}_MT.fasta"
    output:
        config["results"]+"/data/{sample_bam}_MT.fasta.fai"
    log:
        config["results"]+"/logs/{sample_bam}_samtools_faidx_ref_MT.log"
    shell:
        "(samtools faidx {input}) 2> {log}" 
        

rule freebayes_MT:
    input:
        bam_mt_file=config["results"]+"/data/{sample_bam}_MT.bam",
        bai=config["results"]+"/data/{sample_bam}_MT.bam.bai",
        fai=config["results"]+"/data/{sample_bam}_MT.fasta.fai"
    output:
        config["results"]+"/data/{sample_bam}_MT.vcf" 
    params:        
        mt_regions=config["regions_MT"],        
        mt_ref=config["reference_MT"],
        cores=config["cores"]
    log:
        config["results"]+"/logs/{sample_bam}_freebayes_MT.log"
    shell:
        "(./scripts/freebayes-parallel.sh "
        "{params.mt_regions} {params.cores} -f {params.mt_ref} "
        "-iXu -C 2 -q 1 {input.bam_mt_file} > {output}) 2> {log}"
        

rule qc_vcf_MT:
    input:
        config["results"]+"/data/{sample_bam}_MT.vcf"
    output:
        config["results"]+"/data/{sample_bam}_MT_qc.vcf"
    params:
        DP=config["dp"],
        QUAL=config["qual"]
    log:
        config["results"]+"/logs/{sample_bam}_qc_vcf_MT.log"
    shell:
        "(./software/vcffilter -f 'QUAL > {params.QUAL} & DP > {params.DP}' "
        "{input} > {output}) 2> {log}" 
        

rule generate_pileups_MT:
    input:        
        config["results"]+"/data/{sample_bam}_MT_qc.vcf",
        config["results"]+"/cells/{sample_bam}.log",
    output:
        config["results"]+"/iter1/pu/{sample_bam}.log"
    params:
        config["results"]+"/cells/"
    log:
        config["results"]+"/logs/{sample_bam}_generate_pileups_MT.log"
    script:
        "scripts/pileup.py"
    

rule merge_cells_MT:
    input:
        config["results"]+"/iter1/pu/{sample_bam}.log"
    output:
        config["results"]+"/iter1/{sample_bam}_merged_cells.csv"
    params:        
        threshold_coverage=config["threshold_coverage"],
        threshold_coverage_pos=config["threshold_coverage_pos"],
        threshold_base_calling=config["threshold_base_calling"]
    log:
        config["results"]+"/logs/{sample_bam}_merge_cells_MT.log"
    script:
        "scripts/merge_pileup_cells.py"
    

rule filter_cells_MT:
    input:
        config["results"]+"/iter1/{sample_bam}_merged_cells.csv"
    output:
        config["results"]+"/iter1/{sample_bam}_filt_cells.csv"
    params:
        thr_cell_1=config["thr_cell_1"]
    log:
        config["results"]+"/logs/{sample_bam}_filter_cells_MT.log"
    script:
        "scripts/feature_selection.py"
            

rule get_clusters_MT:
    input:
        config["results"]+"/iter1/{sample_bam}_filt_cells.csv"
    output:
        config["results"]+"/iter1/{sample_bam}_clusters.txt"
    params:
        n_neighbors=config["n_neighbors"],
        n_components=config["n_components"],
        output_plot=config["results"]+"/iter1/{sample_bam}_3dplot.eps",
	    clusters=config["clusters"]
    log:
        config["results"]+"/logs/{sample_bam}_get_clusters_MT.log"
    shell:
        "(Rscript --vanilla scripts/get_clusters.R -i {input} -o {output} "
        "-n {params.n_neighbors} -k {params.n_components} -p {params.output_plot} -c {params.clusters} 2> {log})"


rule merge_by_clusters_MT:
    input:
        config["results"]+"/iter1/{sample_bam}_clusters.txt"
    output:
        config["results"]+"/iter1/{sample_bam}_clusters.log"
    params:
        config["results"]+"/cells/",
        config["cores"]
    log:
        config["results"]+"/logs/{sample_bam}_merge_by_clusters_MT.log"    
    script:
        "scripts/split_by_cluster.py"
        

rule freebayes_clusters_MT:
    input:        
        config["results"]+"/iter1/{sample_bam}_clusters.log"
    output:
        config["results"]+"/iter1/{sample_bam}_merge_clusters.vcf"
    params:
        config["reference"],
        config["regions"],
        config["cores"]
    log:
        config["results"]+"/logs/{sample_bam}_freebayes_clusters_MT.log"
    script:
        "scripts/freebayes_clusters.py"
    

rule rm_dup_vcf_1:    
    input:
        config["results"]+"/iter1/{sample_bam}_merge_clusters.vcf"
    output:
        config["results"]+"/iter1/{sample_bam}_merge_clusters_uniq.vcf"
    log:
        config["results"]+"/logs/{sample_bam}_rm_dup_vcf_1.log"
    shell:
        "(./software/bcftools/bcftools norm -d all {input} > {output}) 2> {log}"


rule qc_vcf_1:
    input:
        config["results"]+"/iter1/{sample_bam}_merge_clusters_uniq.vcf"
    output:
        config["results"]+"/iter1/{sample_bam}_snps_iter1.vcf"
    params:
        DP=config["dp"],
        QUAL=config["qual"]
    log:
        config["results"]+"/logs/{sample_bam}_qc_vcf_1.log"
    shell:
        "(./software/vcffilter -f 'QUAL > {params.QUAL} & DP > {params.DP}' "
        "{input} > {output}) 2> {log}"


rule generate_pileups:
    input:        
        config["results"]+"/iter1/{sample_bam}_snps_iter1.vcf"        
    output:
        config["results"]+"/iter2/pu/{sample_bam}.log"
    params:
        config["results"]+"/cells/"
    log:
        config["results"]+"/logs/{sample_bam}_generate_pileups.log"
    script:
        "scripts/pileup.py"    


rule merge_cells:
    input:
        config["results"]+"/iter2/pu/{sample_bam}.log"
    output:
        config["results"]+"/iter2/{sample_bam}_merged_cells.csv"
    params:        
        threshold_coverage=config["threshold_coverage"],
        threshold_coverage_pos=config["threshold_coverage_pos"],
        threshold_base_calling=config["threshold_base_calling"]
    log:
        config["results"]+"/logs/{sample_bam}_merge_cells.log"
    script:
        "scripts/merge_pileup_cells.py"
    

rule filter_cells:
    input:
        config["results"]+"/iter2/{sample_bam}_merged_cells.csv"
    output:
        config["results"]+"/iter2/{sample_bam}_filt_cells.csv"    
    params:
        thr_cell_2=config["thr_cell_2"]
    log:
        config["results"]+"/logs/{sample_bam}_filter_cells.log"
    script:
        "scripts/feature_selection.py"
            

rule get_clusters:
    input:
        config["results"]+"/iter2/{sample_bam}_filt_cells.csv"
    output:
        config["results"]+"/iter2/{sample_bam}_clusters.txt"
    params:
        n_neighbors=config["n_neighbors"],
        n_components=config["n_components"],
        output_plot=config["results"]+"/iter2/{sample_bam}_3dplot.eps",
    	clusters=config["clusters"]	    
    log:
        config["results"]+"/logs/{sample_bam}_get_clusters.log"
    shell:
        "(Rscript --vanilla scripts/get_clusters.R -i {input} -o {output} "
        "-n {params.n_neighbors} -k {params.n_components} -p {params.output_plot} -c {params.clusters} 2> {log})"


rule merge_by_clusters:
    input:
        config["results"]+"/iter2/{sample_bam}_clusters.txt"
    output:
        config["results"]+"/iter2/{sample_bam}_clusters.log"
    params:
        config["results"]+"/cells/",
        config["cores"]
    log:
        config["results"]+"/logs/{sample_bam}_merge_by_clusters.log"
    script:
        "scripts/split_by_cluster.py"
        

rule freebayes_clusters:
    input:        
        config["results"]+"/iter2/{sample_bam}_clusters.log"
    output:
        config["results"]+"/iter2/{sample_bam}_merge_clusters.vcf"
    params:
        config["reference"],
        config["regions"],
        config["cores"]
    log:
        config["results"]+"/logs/{sample_bam}_freebayes_clusters.log"
    script:
        "scripts/freebayes_clusters.py"


