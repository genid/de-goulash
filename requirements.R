#pkgTest <- function(x)
#{
#  if (!require(x,character.only = TRUE))
#  {
#    install.packages(x,dep=TRUE)
#    if(!require(x,character.only = TRUE)) stop("Package not found")
#  }
#}
# loading libraries
install.packages("devtools")
devtools::install_github("marchtaylor/sinkr") # this is for dineof
install.packages("NbClust")
install.packages("umap")
install.packages("scatterplot3d")
install.packages("data.table")
install.packages("optparse")

