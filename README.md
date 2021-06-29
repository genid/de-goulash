# https://github.com/dkoppstein/ngsschool-snakemake-tutorial/blob/master/exercises/exercise3/Snakefile


# Installation 
# Note: Make sure you have installed 
    python>=3.6 & anaconda3  
    R language >=3.6.1 "Action of the toes"
    samtools
    freebayes

# 1) Snakemake and python libraries using conda environment
	#a) cretae new environment followed here
        https://snakemake.readthedocs.io/en/stable/tutorial/setup.html#requirements
	# b) Install vcfflib dependencies with conda                
        conda install -c bioconda vcflib
        pip install pandas
	pip install pyvcf
	 
# 2) R libraries
    # DINEOF
        install.packages("remote")
        remotes::install_github("marchtaylor/sinkr") # this is for dineof
    # NBclust
        install.packages("NbClust")
    # UMAP
        install.packages("umap")    
    # plot 3D
        install.packages("scatterplot3d")    
	install.packages("optparse")

# 3)Extract MT region from reference and indexed both genome.fasta and MT.fasta with
	"samtools faidx" command. For freebayes basedcalling this is needed
               
## Commands for snakemake
1) dry-run (no execute) print workflow
    snakemake -np iter1/{sample}_snps_iter1.vcf
    snakemake -np iter2/{sample}_snps_iter2.vcf
2) Create graph (no execute)
    snakemake --dag iter2/{sample}_snps_iter2.vcf | dot -Tsvg > dag_pipeline.svg
3) Execute pipeline

e.g. {sample} = "cells"
    snakemake --cores 20 iter2/{sample}_snps_iter2.vcf
    snakemake --cores 20 analysis_output/{sample}_forensic_params.txt
    snakemake --cores 10 haplogroups_{sample}/haplogrep/haplogroups.txt

3) Optional
    snakemake --cores 30 cells/cells.log --force --reason
    snakemake --cores 30 iter2/cells_snps_iter2.vcf --force --reason
    

# population indexes
population_names <- c("EUR_AF", "EAS_AF", "AMR_AF", "SAS_AF", "AFR_AF")


Commands
#snakemake --snakefile Snakefile --cores 20;

snakemake results_x/iter2/cells_merge_clusters.vcf --snakefile Snakefile --cores 40 --printshellcmds

snakemake results/iter2/cells_merge_clusters.vcf --snakefile Snakefile --cores 40 --printshellcmds;
snakemake --snakefile Snakefile_analysis --cores 10;

# optional
snakemake --snakefile Snakefile_analysis --cores 40 --printshellcmds -f


### Docker

sudo docker build -t single-cell-docker .


sudo docker run -it single-cell-docker sh # alpine
sudo docker run -it single-cell-docker bash # ubuntu



# Snakemake
sudo docker run -it -v /media/disk1/diego/genid/Single-cell-docker/:/single-cell single-cell-docker output/iter2/cells_merge_clusters.vcf --snakefile Snakefile --configfile config.yaml --cores 40 --printshellcmds

# Snakemake analysis
sudo docker run -it -v /media/disk1/diego/genid/Single-cell-docker/:/single-cell single-cell-docker --snakefile Snakefile_analysis --configfile config.yaml --cores 40 --printshellcmds




