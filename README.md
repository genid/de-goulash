# de-goulash

de-goulash is a bioinformatics pipeline build in [Snakemake](https://snakemake.readthedocs.io/en/stable/) which allows clustering mixed individuals using 10x single-cell RNA-seq.

The pipeline divides in two main steps. 

### 1) Deconvolution of mixed single-cell samples, including genotyping and clustering. 
The inputs needed includes the following: 
1. Possorted_genome_bam.bam
2. barcodes.tsv as output from 10x. this file contains the cells to use onwards.
3. genome.fasta (Human reference genome e.g. hg19 or hg38 in fasta format)
4. *MT.fasta (mitochondrial DNA sequence in fasta format, same build as genome)
5. *region.txt
6. *MT_regions.txt 

Input files with asterik * [4, 5, 6] can be generated with the python script. 
```
python process_reference.py [path/genome.fasta] 
```

### 2) Population assignation and statistical analysis. It needs the output variant calling for each assignated cluster from step 1 and it will calculate likelihood of Forensic Parameters, population assignation, execute haplogrep and finally Yleaf v.2.2. 
The inputs needed includes the following: 
1. Exone reference: exome_96_remmapedto38.vcf.gz
2. Reference population based on 1000G project: 100G_populations.txt
3. Path where the chromosomes for 1000G variant calling: /single-cell/input/1000G/ [https://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/]

### All input parameters need to be within the config.yaml file
### Snakemake step 1
#### Inputs
* sample_bam: /single-cell/input/1/possorted_genome_trimmed.bam
* barcodes: /single-cell/input/1/barcodes_reduced.txt
* reference: /single-cell/input/reference/genome.fasta
* regions: /single-cell/input/reference/regions.txt
* reference_MT: /single-cell/input/reference/MT.fasta
* regions_MT: /single-cell/input/reference/MT_regions.txt
#### Snakemake settings
* cores: 4
* dp: 50
* qual: 60
#### threshold for iteration 1
* thr_cell_1: 10
#### threshold for iteration 2
* thr_cell_2: 20
#### rule for merging cells python
* threshold_coverage: 10
* threshold_coverage_pos: 5
* threshold_base_calling: 90
#### rule for for clustering Rscript
* n_neighbors: 5
* n_components: 3
* clusters: 0 # if clusters > 1 then nBclust is executed to predict number of clusters to use
### Snakemake analysis step 2
#### Inputs
* ref_exome: /single-cell/input/exome_96_remmapedto38.vcf.gz
* ref_population: /single-cell/input/1000G/1000G_populations.txt
* dirpath_1000G: /single-cell/input/1000G/
* dirpath_analysis: output
* dp_2: 50
* qual_2: 60
#### Yleaf parameters
* read_depth: 1
* quality: 20
* base_calling: 90
* positions_file: /single-cell/software/Yleaf/Position_files/WGS_hg38.txt

## Usage in Docker (Linux) (recommended) 

We provided a docker image where you can run the pipeline without having to install any other dependency than docker. Although you need root permissions to proceed. 

Download docker image (2.03gb)
```
docker pull geniderasmusmc/de-goulash:1
```
Tested in Docker version 19.03.2, build 6a30dfc
```
docker --version
```

You can execute de-goulash Snakemake pipeline throught docker image-container. You have to manually mount the current directory where input files are located. 

## 1) de-goulash deconvolution

* Current directory where input files are located -> /current/directory/de-goulash/ 
* Default root location inside the container (do not change) -> :/single-cell 
* Container name -> geniderasmusmc/de-goulash:1
* Target file [only change output name e.g. output_test/iter2/cells_merge_clusters.vcf] -> output/iter2/cells_merge_clusters.vcf

```
docker run -it -v /current/directory/de-goulash/:/single-cell geniderasmusmc/de-goulash:1 output/iter2/cells_merge_clusters.vcf --snakefile Snakefile --configfile config.yaml --cores 1
```

## 2) de-goulash statistical analysis


```
docker run -it -v /media/disk1/diego/genid/Single-cell-docker/:/single-cell single-cell-docker --snakefile Snakefile_analysis --configfile config.yaml --cores 1
```

## Manual installation

Instead of using docker container you can install everything independently and run Snakemake directly

### Tested in:
* R 3.6.1 -- "Action of the toes"
* Python 3.7.3
* Linux Ubuntu 18.04
* Java Run time environment 8

Recommended use conda or Python3 venv 

### Install libraries
```
pip3 install requirements.txt
```
```
Rscript requirements.R
```
```
git clone https://github.com/genid/de-goulash.git
```


### To run through Snakemake pipeline

Step 1
```
snakemake output/iter2/cells_merge_clusters.vcf --snakefile Snakefile --configfile config.yaml --cores 1
```
Step 2
```
snakemake --snakefile Snakefile_analysis --configfile config.yaml  --cores 1
```
