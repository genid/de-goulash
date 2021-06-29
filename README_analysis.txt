## Convert vcfs to vcf.gz and indexed
# 1) Cells snps vs Exome snps
#1 filter exome vcf by cells_snps
bcftools filter -R cells_snps_iter2.vcf.gz exome_96_remmapedto38.vcf.gz > exome_snps_filt.vcf


#2 filter cells_snps by filtered exome vcf
bcftools filter -R exome_96_remmapedto38.vcf.gz cells_snps_iter2.vcf.gz > cells_snps_iter2_filt.vc

inputs:
    exome_snps_filt.vcf
    cells_snps_iter2_filt.vcf

## Compare common positions and see which one has the highest number of match genotypes
and store the value of match and the name of the sample that matches

output:
    cells_sample.vcf



# 2) Compare cells_sample.vcf to the 1000G reference and make a vcf containing all the positions
that matches with the genotypes
Get all valid chromosome from cells_sample.vcf and find the file from 1000G to filtered by position
bcftools filter -R cells_snps_iter2_filt.vcf 1000G/ALL.chr1_GRCh38.genotypes.20170504.vcf.gz > chr1.filt.vcf


# sort vcf\bcf
bcftools sort file.vcf.gz > file_sort.vcf.gz

# Merge 1000G vcfs
bcftools merge --threads 30 1000G/*sites.20170504.vcf.gz -Oz -o 1000G/merged.1000G.sites.vcf.gz
# Sort 1000G vcf
bcftools sort -m 5G -o 1000G/merged.1000G.sites.sort.vcf.gz -O z -T 1000G/ 1000G/merged.1000G.sites.vcf.gz
# Index
tabix -fp vcf 1000G/merged.1000G.sites.sort.vcf.gz

# Remove chr to the entire file
sed 's/chr//' exome_96_remmapedto38_sort.vcf > exome_96_remmapedto38.vcf

#convert vcf files to bcf and indexed (also sort in case of error)
bgzip -c file.vcf > file.vcf.gz
tabix -fp vcf file.vcf.gz

#This should be use in gatk SelectVariants to filtered out the biggest vcfs
# snps.list should be pass to make each vcf smaller
../software/gatk-4.1.9.0/gatk SelectVariants -R ../reference/genome.fasta -V exome_96_remmapedto38.vcf.gz -L snps.list -O exome_cells.vcf


a) Use cells_snps_iter2.vcf file generated from iter2 and compared against exome_96_remmapedto38.vcf
# Criterias:
    1) Compared exact genotype against each Sample (e.g S1, S2 etc) and select the Sample with highes percentage of identity.
    number_GT_matches/total_snps 
    2) Filtered and create S1.vcf based on the common ones from the vcf and exome vcf file. 

b) Compared S1.vcf with all vcfs from the 1000G vcf folders (per chromosome) and generate a new file
S3_1000G.vcf which matches by the chromosomical position.
    Use S1.vcf, S1_1000G.vcf in Structure part. (1000G_populations.txt is fixed input file)

c) In script from structure part store the last output in a text file


Inputs:
    cells_snps_iter2.vcf # file    
    Structure_part/1000G_populations.txt # file
    exome_96_remmapedto38.vcf # file
    1000G/ # folder including all vcf files per chromosome
    # S1.vcf # file_name base on sample name
    # S1_1000G.vcf # file_name base on sample and 100G reference posfix
Outputs:
    Familias_frequency_matrix417SNPs.txt
    STRUCTURE_417SNPs.txt
    forensic_parameters.txt


## command for extract information from 1000G via ftp

inputs:
1) list of all reference genomes
2) chromosome to look for
3) list of snps from specific chromosome

tabix -h ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/ALL.chrY.phase3_integrated_v1b.20130502.genotypes.vcf.gz | vcf-subset -c HG00096,HG00101,HG00103,HG00105,HG00107 ALL.chrY.phase3_integrated_v1b.20130502.genotypes.vcf.gz | bgzip -c | bcftools view | more
tabix -h ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/ALL.chrY.phase3_integrated_v1b.20130502.genotypes.vcf.gz | vcf-subset -c HG00096 ALL.chrY.phase3_integrated_v1b.20130502.genotypes.vcf.gz | bgzip -c | bcftools view | head
tabix -h ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/ALL.chrY.phase3_integrated_v1b.20130502.genotypes.vcf.gz  | vcf-subset -c HG00096,HG00101,HG00103,HG00105,HG00107 ALL.chrY.phase3_integrated_v1b.20130502.genotypes.vcf.gz | bgzip -c | bcftools view | grep "2659378\|2659347\|2662121"
