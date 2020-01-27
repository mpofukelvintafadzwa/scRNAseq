# Load package
library(splatter)
library(scater)
path = "/home/ivory/VersionControl/Dropbox/SingleCellSequencing/Real Data Sets/Islam Data"
file = "GSE29087_original.csv"
filepath = paste(path, "/", "GSE29087_original.csv", sep = "")
dat = read.delim(filepath, sep=",", header=TRUE)

# Extract the row names which are gene names
rownames(dat) <- dat[,1]

# Remove the row names from dat
dat <- dat[,-1]

# save  the row names
genes <- rownames(dat)

# Get the coloumn names which is cell ID.
cellIDs <- colnames(dat)

# Build a scater object
library("SingleCellExperiment")
library("scater")