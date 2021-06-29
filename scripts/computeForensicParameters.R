
##################################
##
## Compute forensic parameters
## Version 1.3
## Author: Daniel Kling
##
#######################
## - Added STRUCTURE file preparations
## - Separated STRUCTURE analysis from rest
## - Seperated thinning procedure
##################################


###############################
##
## Function to thin data
##
###############################
thinData <- function(sample, mafs, maf, maf_difference, min_distance)
{
	# First based on MAF
	selectedmarkers <- array(FALSE, nrow(sample))

	######## PART 1 ###################
	## This is marker pruning
	## based on MAF. maf=0.2 sets 
	## minumum frequency to 0.2
	## maf_difference=0.3 sets the 
	## maximally allowed maf difference
	## between populations
	##
	## IFF maf=-1 then a different approach is used
	###################################

	for(i in 1:nrow(sample))
	{
		# Loop populations and make sure MAF>maf
		isincluded = TRUE
		for(j in 1:5)
		{
			if(mafs[i,j]>=maf & mafs[i,j]<=(1-maf))
			{
				# Do nothing
				
			}
			else
				isincluded = FALSE
			#else
			#	isincluded = FALSE
		}
		# Check difference
		if(mafs[i,6]>maf_difference)
			isincluded = FALSE

		if(maf==-1)
		{
			if(mafs[i,6]<maf_difference)
				isincluded = FALSE
			else
				isincluded = TRUE
		}
		selectedmarkers[i] <- isincluded
	}

	######## PART 2 ###################
	## This is marker pruning/thinning
	## based on bp distance.  
	## min_distance corresponds to inclusion
	##
	##
	###################################


	# Loop and thin if needed, should use cM in the future, for now stick to bp
	currchr = 0
	currpos = 10000000000000000000
	#selectedmarkers <- array(FALSE, nrow(sample))

	for (i in 1:nrow(sample))
	{
		if(selectedmarkers[i])
		{
			selectedmarkers[i] = FALSE
			currentrow = t(sample[i,])
			if(currchr!=currentrow[1])
			{
				# Select first marker on the chromosome
				selectedmarkers[i] = TRUE
				currchr = currentrow[1]
				currpos = currentrow[2]
		
			}
			else if( (as.numeric(currentrow[2])-as.numeric(currpos))>min_distance)
			{
				selectedmarkers[i] = TRUE
				currchr = currentrow[1]
				currpos = currentrow[2]
			}
				
		}
	}
	return(selectedmarkers)
}

library(optparse)
library(stringr)


option_list = list(
  make_option(c("-i", "--input"), type="character", help="input VCF file cells", metavar="character"),
  make_option(c("-r", "--reference"), type="character", help="input VCF file 100G reference", metavar="character"),
  make_option(c("-p", "--sample_population"), type="character", help="sanoke population", metavar="character"),
  make_option(c("-o", "--output"), type="character", help="output file name store forensics parameters", metavar="character"),
  make_option(c("-m", "--mainparams"), type="character", help="main parameters", metavar="character"),
  make_option(c("-e", "--extraparams"), type="character", help="extra parameters", metavar="character"),
  make_option(c("-f", "--folder"), type="character", help="extra parameters", metavar="character")
); 


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

sample_input <- opt$input
ref_input <- opt$reference
population_input <- opt$sample_population
output <- opt$output # name output file
main_params_file <- opt$mainparams
extra_params_file <- opt$extraparams 
folder <- opt$folder

#setwd("/media/disk1/diego/git/Single-cell/")
#sample_input <- "output/1/1.sample_noheader.vcf"
#output <- "output/1/1_forensic_params.txt"
#ref_input <- "output/1/1.ref_noheader.vcf"
#population_input <- "/media/disk1/diego/git/Single-cell/input/1000G/1000G_populations.txt"
#main_params_file <- "/media/disk1/diego/git/Single-cell/input/structure/mainparams"
#extra_params_file <- "/media/disk1/diego/git/Single-cell/input/structure/extraparams"

splitted_out <- strsplit(output, "/")[[1]]
dirpath <- paste0(splitted_out[1:length(splitted_out)-1], collapse = "/") # folder to store additional outputs

## Input files
ref <- as.matrix(read.table(ref_input, comment.char = "!"))
ind1 <- as.matrix(read.table(sample_input, sep="\t"))
sample_populations <- as.matrix(read.table(population_input, sep="\t"))

##################
### Settings #####
population <- "EUR"   # Used for forensic computations later - could also be based on maximum likelihood estimate
new_allele_freq <- 0.01  # Frequency given to alleles not observed in a population but detected in sample
min_distance <- 500000 # Let's use 500 kb for now, removes LD to a degree
maf <- 0.2 # Should not be too high (higher the better but stricter filtering)
maf_difference <- 0.3 # Maximum difference (absolute) between populations
pops <- c("EUR","EAS", "AMR", "SAS","AFR")
K <- 5 # STRUCTURE settings 

###########################
## Read reference data
###########################

ref_samples <- c("", "",t(ref[1,-c(1,2,3,4,5,6,7,8,9)]), "Sample")
ref <- ref[-1,]

#### Prepare STRUCTURE data frames ####
# ref_samples now contains the first column of the STRUCTURE file
# You need to create the second colum containing the known populations - to be used as locprior
known_populations <- array("0", length(ref_samples))
# sample_populations <- as.matrix(read.table("1000genomes_sample_infos.txt", sep="\t"))
#sample_populations <- as.matrix(read.table("1000G_populations.txt", sep="\t"))


# This information needs to pair HGxxxx with the population from some other file
populations <- unique(sample_populations[,2])
# population_ids <- 1:length(populations)
for(i in 3:(length(ref_samples)-1))
{
	matches <- match(ref_samples[i],sample_populations[,1])
	matches <- match(sample_populations[matches,2],populations)
	known_populations[i] <- matches
}

## Remove sample with unknown populations, i.e. does not have a row in the 1000G_populations file
ref_samples <- ref_samples[!is.na(known_populations)]
ref <- ref[,!is.na(known_populations)]
known_populations <- known_populations[!is.na(known_populations)] # Used for LOCPRIOR or POPINFO
popflags <- array(1,length(known_populations))
popflags[length(popflags)] <- 0

###########################
## Check input
###########################

# Sort data and remove markers not on the autosomes
ind1 <- ind1[!is.na(as.numeric(ind1[,1])),] 
ind1 <- ind1[order(as.numeric(ind1[,1]), as.numeric(ind1[,2])),]
ref <- ref[!is.na(as.numeric(ref[,1])),]
ref <- ref[order(as.numeric(ref[,1]), as.numeric(ref[,2])),]

if(nrow(ind1)!=nrow(ref))
{
	print("Error, we do not have a perfect overlap between reference and our profile!")
	matches <- match(paste(as.numeric(ind1[,1]),"-",as.numeric(ind1[,2]), sep=""), paste(as.numeric(ref[,1]),"-",as.numeric(ref[,2]), sep=""))
	ind1 <- ind1[!is.na(matches),]
	ref <- ref[matches[!is.na(matches)],]
}

# Assumes perfect overlapping input files and sorted


###############################
##
## Function to get MAFs
##
###############################
getMAFs <- function(ref, pops)
{
	selectedmarkers <- array(FALSE, nrow(ref))
	mafs <- matrix(0, nrow(ref), 6)
	colnames(mafs) <- c(pops, "Max diff")

	for(i in 1:nrow(ref))
	{
		info <- strsplit(as.character(ref[i,8]), ';')
		# Loop populations and make sure MAF>maf
		for(j in 1:5)
		{
			pos <- which(regexpr(pops[j], info[[1]])==1)
			if(!is.na(as.numeric(substr(info[[1]][pos],8,13))))
			{
				mafs[i,j] <- as.numeric(substr(info[[1]][pos],8,13))
				#if(mafs[i,j]>maf & mafs[i,j]<(1-maf))
				#{
					# Do nothing
				
				#}
				#else
				#	isincluded = FALSE
			}
			#else
			#	isincluded = FALSE
		}
		# Check difference
		max = 0
		min = 1
		for(j in 1:5)
		{
			if(as.numeric(mafs[i,j])>max)
				max = as.numeric(mafs[i,j])
			if(as.numeric(mafs[i,j])<min)
				min = as.numeric(mafs[i,j])
		}
		mafs[i,6] <- max-min
		#if((max-min)>maf_difference)
		#	isincluded = FALSE
		#selectedmarkers[i] <- isincluded
	}
	return(mafs)
}
mafs <- getMAFs(ref, pops)

length(which(mafs[,1]>0.3 & mafs[,1]<0.7))

selected_markers <- thinData(ind1, mafs, 0, 0.3, 0)
length(which(selected_markers))

##################################
## STRUCTURE
##################################

#selectedmarkers <- thinData(ind1, mafs, maf=-1, maf_difference=0.4, min_distance=500000)
selectedmarkers <- thinData(ind1, mafs, maf=0, maf_difference=1, min_distance=500000)

paste(length(which(selectedmarkers)), " markers kept after thinning", sep="")
ind1_pruned_thinned <- ind1[selectedmarkers,]
ref_pruned_thinned <- ref[selectedmarkers,]
mafs_pruned_thinned <- mafs[selectedmarkers,]

temp <- t(cbind(ref_pruned_thinned[,-c(1,4,5,6,7,8,9)],substr(ind1_pruned_thinned[,10],1,3)))
genotypes <- matrix("", nrow(temp), ncol(temp)*2)
for(i in 3:nrow(temp))
{
	for(j in 1:ncol(temp))
	{
		genotypes[i,1+2*(j-1)] <- substr(temp[i,j],1,1)
		genotypes[i,2+2*(j-1)] <- substr(temp[i,j],3,3)
	}
}

colnames(genotypes) <- rep(temp[2,],each=2) 

structure <- cbind(ref_samples, known_populations, popflags, genotypes)

## Remove two first rows unless we want this in the output
structure <- structure[-c(1,2),]
struct_file <- paste(dirpath, "/STRUCTURE_", nrow(ind1_pruned_thinned), "SNPs.txt", sep="")
write.table(file=struct_file, structure, quote=F, col.names=F, row.names=F, sep="\t")
## Run STRUCTURE with LOCPRIOR and PFROMPOPFLAGONLY active, see extraparams file
#command = paste("structure -K ", K ," -N ", nrow(structure) ," -m mainparams -e extraparams -L ",ncol(genotypes)/2 ," -i \"", struct_file, "\" -o \"", "structure_locprior_K",K, "\"", sep="")
command = paste("structure -K ", K ," -N ", nrow(structure) ," -m ", main_params_file, " -e ", main_params_file, " -L ",ncol(genotypes)/2 ," -i \"", struct_file, "\" -o \"", dirpath, "/structure_locprior_K",K, "\"", sep="")
system(command = command)

## Get results and plot
palette(colorRampPalette(c("#00AFBB", "#E7B800", "#FC4E07", "#0D0887FF", '#4daf4a','#984ea3', '#9C8847' ))(K) )
#tiff("/Figure_STRUCTURE.tiff", res=150, width=10, height=15, units='in')

setEPS()
postscript(file = paste(dirpath, "/Figure_STRUCTURE.eps", sep=""), width = 15, height = 30, paper = "a4")
#tiff(paste(dirpath, "/Figure_STRUCTURE.tiff", sep=""), res=150, width=10, height=15, units='in')
par(mfrow=c(1,1), mar=c(4,2,2,2))
#for(k in 5:12)
#{	
	results <- read.table(paste(dirpath, "/structure_locprior_K",K, "_f", sep=""), sep="\t", fill = T)
	# Find likelihood
	lik <- results[which(substr(as.character(results[,1]),1,12)=="Estimated Ln"),]
	pos <- regexpr("=", lik)		
	likelihood <- as.numeric(substr(lik, pos+1, 50))
	# Find sample start
	start <- which(substr(as.character(results[,1]),1,17)=="Inferred ancestry")+2
		
	plot_data <- matrix(0, 1, K+1)
	sample_data <- results[start:(start+(nrow(genotypes)-3)),1]
	for(s in 1:length(populations))
	{
		subdata <- sample_data[known_populations[-c(1,2)]==s]
		population_data <- matrix("", length(subdata), K+1)
		population_data[,1] <- s
		#	known_populations
		for(j in 1:length(subdata))
		{
			numbers <- unlist(strsplit(as.character(subdata[j]), split=" "))
			pos <- which(numbers==":")+2
			population_data[j,2:(K+1)] <- numbers[pos:(pos-1+K)]
		}
			
		# Sort individuals based on the clusters
		#if(s<7)
		#{
			max_col <- max.col(apply(population_data[,-1],2,as.numeric))[1]
			# next_max_col <- max.col(apply(population_data[,c(-1, -(max_col+1))], 2, as.numeric))
			population_data <- population_data[order(as.numeric(population_data[,max_col+1])),]
		#}
			
		plot_data <- rbind(plot_data, population_data)
	}

	## Add sample at the end
	subdata <- sample_data[length(sample_data)]
	population_data <- matrix("", 1, K+1)
	numbers <- unlist(strsplit(as.character(subdata), split=" "))
	pos <- which(numbers==":")+2
	population_data[1,2:(K+1)] <- numbers[pos:(pos-1+K)]
	plot_data <- rbind(plot_data, population_data)

	pop_sep <- array(0,1)
	plot_data <- plot_data[-1,]
	for(y in 2:nrow(plot_data))
	{
		if(plot_data[y,1]!=plot_data[y-1,1])
			pop_sep <- c(pop_sep, y)
	}
	bp <- barplot(t(plot_data[,-1]), col=1:K, border=NA, space=0, axes=F, cex.names=1.1, width=c(rep(1,each=nrow(genotypes)-1), 10))
	bp2 <- barplot(as.matrix(as.numeric(plot_data[nrow(plot_data),-1])), col=1:K, border=NA, axes=F, cex.names=1.1, width=50, beside=F, add=T, space=c(48,0), xlim=c(0,2400))
	
	# Create ablines and texts
	title(paste("K=", K,", logLik=", likelihood, sep=""))	
	abline(v=pop_sep, xpd=F, lwd=2)
	par(xpd=T)
	text(x=pop_sep+100, y=-0.05, labels = c(pops, "Sample"), srt=45, cex=1)
	par(xpd=F)
#}	
dev.off()

bp <- barplot(t(plot_data[,-1]), col=1:K, border=NA, space=0, axes=F, cex.names=1.1, width=c(rep(1,each=nrow(genotypes)-1), 10))
bp2 <- barplot(as.matrix(as.numeric(plot_data[nrow(plot_data),-1])), col=1:K, border=NA, axes=F, cex.names=1.1, width=50, beside=F, add=T, space=c(48,0), xlim=c(0,2400))
	
# Create ablines and texts
title(paste("K=", K,", logLik=", likelihood, sep=""))	
abline(v=pop_sep, xpd=F, lwd=2)
par(xpd=T)
text(x=pop_sep+100, y=-0.05, labels = c(pops, "Sample"), srt=45, cex=1)
	


###### END ##########



##################################
## Forensic parameters
##################################

## We thin the data for independent likelihood ratio computations

selectedmarkers <- thinData(ind1, mafs, maf, maf_difference, min_distance)

paste(length(which(selectedmarkers)), " markers kept after thinning", sep="")
ind1_pruned_thinned <- ind1[selectedmarkers,]
ref_pruned_thinned <- ref[selectedmarkers,]
mafs_pruned_thinned <- mafs[selectedmarkers,]


# Function to compute the LR comparing the hypotheses
# H1 : The person contributed to the crime scene sample
# H2 : A relative of the person contributed, related through the k0-k2 (defined for all non-inbred relationships)
computeLR <- function(a1, a2, k0, k1, k2, maf)
{
	if(a1==a2)
	{
		if(a1=="1")
			p = maf	
		else
			p = 1-maf
		return(p*p/(k0*p*p*p*p + k2*p*p + k1*p*p*p))
	}
	else
	{
		p = maf
		return(2*p*(1-p)/(k0*4*p*p*(1-p)*(1-p) + k2*2*p*(1-p) + k1*p*(1-p)))
	}
}


pm = 1 # Match probability, probability that two random individuals will have same genotypes
pe = 1 # Power of exclusion
pd = 1 # Power of discrimination

## Prepare frequency data for Familias - useful for subsequent simulations etc
familias_matrix <- matrix(0,3,nrow(ind1_pruned_thinned)+1)
familias_matrix[,1] <- c("Alleles", "1", "2") # We do not require correct allele names for simulations, simply frequencies, hence 1 and 2
familias_matrix[1,] <- c("Alleles", t(ind1_pruned_thinned[,3]))


rmp = 1
LR = 1
for(i in 1:nrow(ind1_pruned_thinned))
{
	g1 <- ind1_pruned_thinned[i,10]
	g1 <- substr(g1,1,3)
	a1 <- substr(g1, 1,1)
	a2 <- substr(g1, 3,3)
	# Should first check that the ref/alt allele is same for reference file
	if(ind1_pruned_thinned[i,4]==ref_pruned_thinned[i,4] && ind1_pruned_thinned[i,5]==ref_pruned_thinned[i,5])
	{
		#info <- strsplit(as.character(ref_pruned_thinned[i,8]), ';')
		# Find the position of the population we are interested in
		maf <- mafs_pruned_thinned[i,match(population, pops)] #as.numeric(substr(info[[1]][6],8,13))
		if(maf==0)
			maf = new_allele_freq
		familias_matrix[2,i+1] <- maf		
		familias_matrix[3,i+1] <- 1-maf
		pmi = (maf*maf*maf*maf+(1-maf)*(1-maf)*(1-maf)*(1-maf)+4*maf*maf*(1-maf)*(1-maf))	
		pm = pm * pmi
		pdi = 1 - pmi
		pd = pd * (1-pdi)
		h = maf*maf+(1-maf)*(1-maf)
		pei = (1-h)*(1-2*(1-h)*h*h)
		pe = pe * (1-pei)
		# Random match probability and LR
		if(a1==a2)
		{
			if(a1=="1")
				rmp = rmp * maf * maf
			else
				rmp = rmp * (1-maf) * (1-maf)
		}
		else
			rmp = rmp * 2 * maf * (1-maf)
		print(paste("RMP for marker", i, "is", rmp, sep=" "))
	
		# Compute LR for "it was my brother who did it" or any other relative, defined by the k0-k2 parameters
		LR = LR*computeLR(a1, a2, k0=0.25, k1=0.5, k2=0.25, maf)
		print(paste("LR for marker", i, "is", computeLR(a1, a2, k0=0.25, k1=0.5, k2=0.25, maf), sep=" "))
	}
	else
	{
		print(paste("Error in marker ", i,", possibly due to wrong ref/alt alleles", sep=""))
	}	
}

##### MISC ######
## This function writes a Familias compatible frequency database
write.table(file=paste(dirpath, "/Familias_frequency_matrix_", nrow(ind1_pruned_thinned), "SNPs.txt", sep=""), familias_matrix, quote=F, col.names=F, row.names=F, sep="\t")
#################


pd = 1 - pd
pe = 1 - pe

print(paste("Total number of markers used: ",nrow(ind1_pruned_thinned), sep=""))
print(paste("Combined power of exclusion (PE): ", pe, sep=""))
print(paste("Combined power of discrimination (PD): ", pd, sep=""))
print(paste("Combined match probability (PM): ", pm, sep=""))

print(paste("Total profile RMP: ", rmp, sep=""))
print(paste("Total profile 1/RMP: ", 1/rmp, sep=""))
print(paste("Total profile LR (it was my brother who did it): ", LR, sep=""))

print(output)
fileConn <- file(output)
writeLines(paste(
  "",
  paste("Total number of markers used: ",nrow(ind1_pruned_thinned),  "\n",  sep=""),
  paste("Combined power of exclusion (PE): ", pe, "\n", sep=""),
  paste("Combined power of discrimination (PD): ", pd, "\n", sep=""),
  paste("Combined match probability (PM): ", pm, "\n",  sep=""),
  
  paste("Total profile RMP: ", rmp, "\n", sep=""),
  paste("Total profile 1/RMP: ", 1/rmp, "\n", sep=""),
  paste("Total profile LR (it was my brother who did it): ", LR, "\n", sep="")), fileConn)
close(fileConn)

#files <- list.files(folder, pattern = "*.vcf")
#print(files)
#file.remove(paste(getwd(), folder, files, sep = "/"))
#warnings()


