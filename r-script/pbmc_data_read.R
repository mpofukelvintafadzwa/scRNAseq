# Author: Nick Lytal
# Reading in PBMC data
library(Matrix)
data <- readMM("matrix.mtx") # Load the gene matrix
data.mat = as.data.frame(as.matrix(data))

# Now use barcodes.tsv and feature.tsv to fill in the row & column names
# Can convert to .csv in Excel first
cells <- read.table(file = 'barcodes.tsv', sep = '\t', header=F)
genes <- read.table(file = 'genes.tsv', sep = '\t', header = F)

colnames(data.mat) = cells[1:2700,]
rownames(data.mat) = genes[,2]
