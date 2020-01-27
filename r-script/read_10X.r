# Author: Rahul Bhadani
# Email: rahulbhadani@email.arizona.edu

library("Matrix")
path = "/home/ivory/VersionControl/scRNASeq/data/tabula_muris/droplet/Kidney-10X_P4_5"
cellfile = "barcodes.tsv"
filepath = paste(path, "/", cellfile, sep = "")
cellbarcodes <- read.table(filepath)

genefile = "genes.tsv"
filepath = paste(path, "/", genefile, sep = "")
genenames <- read.table(filepath)

UMIfile = "matrix.mtx"
filepath = paste(path, "/", UMIfile, sep = "")
# readMM is used for reading sparse matrix file
molecules <- readMM(filepath)

head(cellbarcodes)

# Due to duplicate gene names, column 1 may be preferred over column 2
# since they will be uniquely named
rownames(molecules) <- genenames[,1]

colnames(molecules) <- paste("10X_P4_5", cellbarcodes[,1], sep="_")

# We need to get metadata and computational annotations
metadatapath = "/home/ivory/VersionControl/scRNASeq/data/tabula_muris"
metadatafile = "droplet_metadata.csv"
metafilepath = paste(metadatapath, "/", metadatafile, sep="")
metadata <- read.delim(metafilepath, sep=",", header = TRUE)

# Now we will read the droplet annotations:
anndatafile = "droplet_annotation.csv"
annfilepath = paste(metadatapath, "/", anndatafile, sep="")
anndata <- read.delim(annfilepath, sep=",", header = TRUE)

# If you look at the output of `head(cellbarcodes)`, 
# there is observable difference in formatting - 
# we want to make sure that cell type names match 
# from cell names and from annotations to easily 
# access metadata and annotations.
ann[,1] <- paste(anndata[,1], "-1", sep="")
ann_subset <- anndata[match(colnames(molecules), anndata[,1]),]
celltype <- ann_subset[,3]


# Now lets build the cell-metadata dataframe for specific 
# mouse ID "3_8_M" corresponding to 10X_P4_5, which you 
# can verify above from the output 
# of `metadata[metadata$channel == "10X_P4_5",]`
mouseID = "3_8_M"

cell_anns <- data.frame(mouse = rep(mouseID, times=ncol(molecules)), type=celltype)
rownames(cell_anns) <- colnames(molecules)

## At this point, we want to store sparse matrix `molecules` as a csv file
moleculesdata = as.data.frame(as.matrix(molecules))
rownames(moleculesdata) <- rownames(molecules)
colnames(moleculesdata) <- colnames(molecules)

# At this point, we want to store sparse matrix molecules as a csv file
write.csv(moleculesdata, paste(path, "/", "moleculatesdata.csv", sep=""))

#########################################################################
######  Synthetic Data Simulation from Real Data using Splatter #########
#########################################################################

## Do the count data estimation
library(splatter)

# First convert sparse matrix to regular dense matrix
m = as.matrix(molecules)
params <- splatEstimate(m)
# params is S4 class, hence, its attributes can be accessed using @
# see http://oshlacklab.com/splatter/reference/SplatParams.html for 
# list of all attibutes

nGenes <- params@nGenes # The number of genes to simulate.
nCells <- params@nCells # The number of cells to simulate.
seed <- params@seed #

nBatches <- params@nBatches #The number of batches to simulate.
batchCells <- params@batchCells #Vector giving the number of cells in each batch.
facLoc <- params@batch.facLoc #Location (meanlog) parameter for the batch effect factor log-normal distribution. Can be a vector.
facScale <- params@batch.facScale #Scale (sdlog) parameter for the batch effect factor log-normal distribution. Can be a vector.

meanshape <- params@mean.shape # Shape parameter for the mean gamma distribution.
meanrate <- params@mean.rate #Rate parameter for the mean gamma distribution.

libloc <- params@lib.loc #Location (meanlog) parameter for the library size log-normal distribution, or mean parameter if a normal distribution is used.
libscale <- params@lib.scale #Scale (sdlog) parameter for the library size log-normal distribution, or sd parameter if a normal distribution is used.
libnorm <- params@lib.norm #Logical. Whether to use a normal distribution for library sizes instead of a log-normal.

outprob <- params@out.prob #Probability that a gene is an expression outlier.
outfacLoc <- params@out.facLoc #Location (meanlog) parameter for the expression outlier factor log-normal distribution.
outfacScale <- params@out.facScale # Scale (sdlog) parameter for the expression outlier factor log-normal distribution.

nGroups <- params@nGroups # The number of groups or paths to simulate.
groupprob <- params@group.prob # Probability that a cell comes from a group.

deprob <- params@de.prob # Probability that a gene is differentially expressed in a group. Can be a vector.
dedownProb <- params@de.downProb # Probability that a differentially expressed gene is down-regulated. Can be a vector.
defacLoc <- params@de.facLoc # Location (meanlog) parameter for the differential expression factor log-normal distribution. Can be a vector.
defacScale <- params@de.facScale # Scale (sdlog) parameter for the differential expression factor log-normal distribution. Can be a vector.

bcvcommon <- params@bcv.common # Underlying common dispersion across all genes.
bcvdf <- params@bcv.df # Degrees of Freedom for the BCV inverse chi-squared distribution.

# The type of dropout to simulate. 
# "none" indicates no dropout, "experiment" is global dropout using the 
# same parameters for every cell, "batch" uses the same parameters 
# for every cell in each batch, "group" uses the same parameters for 
# every cell in each groups and "cell" uses a different 
# set of parameters for each cell.
dropouttype <- params@dropout.type 
dropoutmid <- params@dropout.mid # Midpoint parameter for the dropout logistic function.
dropoutshape <- params@dropout.shape # Shape parameter for the dropout logistic function.

# Vector giving the originating point of each path. 
# This allows path structure such as a cell type which 
# differentiates into an intermediate cell type that 
# then differentiates into two mature cell types.
# A path structure of this form would have a "from" 
# parameter of c(0, 1, 1) (where 0 is the origin). 
# If no vector is given all paths will start at the origin.
pathfrom <- params@path.from
# Vector giving the number of steps to simulate
# along each path. If a single value is given it
# will be applied to all paths. This parameter was
# previously called path.length.
pathnSteps <- params@path.nSteps
# Vector giving the skew of each path. Values closer 
# to 1 will give more cells towards the starting
# population, values closer to 0 will give more cells
# towards the final population. If a single value is
# given it will be applied to all paths.
pathskew <- params@path.skew

# Probability that a gene follows a non-linear path 
# along the differentiation path. This allows more 
# complex gene patterns such as a gene being equally 
# expressed at the beginning an end of a path but lowly 
# expressed in the middle.
pathnonlinearProb <- params@path.nonlinearProb
# Sigma factor for non-linear gene paths. 
# A higher value will result in more extreme non-linear
# variations along a path.
pathsigmaFac <- params@path.sigmaFac
params@dropout.type <- 'experiment'


