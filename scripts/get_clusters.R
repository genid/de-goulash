#install.packages("manipulateWidget")
#install.packages("rgl", repos="http://R-Forge.R-project.org")
#install.packages("remotes")
#remotes::install_github("marchtaylor/sinkr") # this is for dineof
#library(data.table)
#library(irlba)
#library(ggplot2)
#library(dbscan)
#library(rgl)
#library(car)
#setwd("/media/disk1/diego/git/Single-cell/pipeline")

#process.beta.fread <- function(beta){
#  cpgs <- as.vector(beta$V1)
#  beta <- beta[,-1]
#  beta <- as.matrix(beta)
#  rownames(beta) <- cpgs
#  return(beta)
#}

library(NbClust)
library(sinkr)
library(umap)
library(scatterplot3d)
library(optparse)
library(data.table)


option_list = list(
    make_option(c("-i", "--input"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
    make_option(c("-o", "--output"), type="character", default="out.txt", 
              help="output file name [default= %default]", metavar="character"),
    make_option(c("-n", "--n_neighbors"), type="character", default = 300, 
              help="umap n_neighbors [default=300]", metavar="int"),    
    make_option(c("-k", "--n_components"), type="character", default=3, 
              help="umap n_components [default=3]", metavar="int"),
    make_option(c("-p", "--out_plot"), type="character", 
              help="output 3d plot filename", metavar="character"),
    make_option(c("-c", "--clusters"), type="integer", 
                help="0 runs NBclust if positive uses given cluster K", metavar="integer")
    
); 
 
#make_option(c("-d", "--min_dist"), type="character", default=10^-10,
#              help="umap min_dist  [default=10^-10]", metavar="int"),

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

namefile <- opt$input
output <- opt$output
n_neighbors <- as.integer(opt$n_neighbors)
min_dist <- as.integer(opt$min_dist)
n_components <- as.integer(opt$n_components)
out_plot <- opt$out_plot
clusters <- as.integer(opt$clusters)

## Umap parameters
#namefile <- "../results_test/iter2/cells_filt_cells.csv"

tmp <- fread(file=namefile, sep = ",", header = T)
tmp <- as.data.frame(tmp)
rownames(tmp) <- tmp$V1
tmp$V1 <- NULL
df.load <- tmp
df <- as.data.frame(tmp)

#tmp <- read.csv(file=namefile, sep = ",", row.names = 1)
#df.load <- read.csv(file=namefile, sep = ",", row.names = 1)
#df <- read.csv2(file = namefile, sep = ",", header=T, row.names = "X")

#n_neighbors <- length(df)
if(n_neighbors > length(df)){
  n_neighbors <- length(df)
}
min_dist <- 10^-10
#n_components <- 3
#clusters <- 0


################################
df.cells <- as.data.frame(lapply(df, as.factor))
df.cells <- as.data.frame(lapply(df.cells, as.numeric))
rownames(df.cells) <- rownames(tmp)
cols <- c()
for (col in colnames(df.cells)){
  val <- sum(df.cells[col] == 3)
  if(val > 0){
    cols <- c(cols, col)
    df.cells[col] <- df.cells[col] - 1
  }
}
df.cells[df.cells == 0] <- NA
#df.cells[df.cells == 1] <- -1
#df.cells[df.cells == 2] <- 1
for(i in 1:ncol(df.cells)){df.cells[,i] <- as.integer(df.cells[,i])}
df.cells <- as.matrix(df.cells)



################################

#n_neighbors <- length(df)
#min_dist <- 10^-10
#n_components <- 3
#n_neighbors <- 500

##########################################################################
#### Imputation Data INterpolating Empirical Orthogonal Functions (DINEOF)
##########################################################################

df.dineof <- dineof(Xo = df.cells, n.max=NULL, ref.pos=NULL, delta.rms=1e-5)

##########################################################################
#### UMAP
##########################################################################
df.tmp <- df.dineof$Xa[ , apply(df.dineof$Xa, 2, var) != 0]
custom.settings <- umap.defaults
custom.settings$n_neighbors <- n_neighbors
custom.settings$min_dist <- min_dist
custom.settings$n_components <- n_components
#print(custom.settings)
df.umap <- umap(df.tmp, config = custom.settings, method = c("naive"))
gr.umap <- df.umap$layout
rownames(gr.umap) <- rownames(tmp)

# NCclust
if(clusters == 0){
  nb.clust.data <- NbClust(data = gr.umap, diss = NULL, distance = "euclidean",
          min.nc = 2, max.nc = 15, method = "kmeans")
  
  out.clust <- as.data.frame(table(nb.clust.data$Best.nc[1,]))
  k <- as.integer(min(as.vector(out.clust[ out.clust$Freq == max(out.clust$Freq), ]$Var1)))
}else{
  k <- clusters
}
cl <- kmeans(gr.umap, centers = k, iter.max = 3000, nstart = 10)
df.plot <- data.frame(V1 = as.numeric(gr.umap[,1] ), V2 = as.numeric(gr.umap[,2]), V3 = as.numeric(gr.umap[,3]), cluster = as.factor(cl$cluster))
setEPS()
postscript(file = out_plot, width = 8, height = 6, title = "Umap genotype imputation", paper = "a4")
scatterplot3d(x = df.plot$V1, y = df.plot$V2, z = df.plot$V3, color = as.factor(df.plot$cluster),
              xlab = "umap 1", ylab = "umap 2", zlab = "umap 3")
dev.off()
write.table(cbind(rownames(df.plot), df.plot$V1, df.plot$V2, df.plot$V3, as.factor(cl$cluster)), file=output, row.names=FALSE, col.names=FALSE, sep="\t", quote = F)

