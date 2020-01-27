# Author: Nic Lytal
# Reading in PBMC data
library(Matrix)
data <- readMM("matrix.mtx") # Load the gene matrix
pbmc.data = as.data.frame(as.matrix(data))

# Now use barcodes.tsv and feature.tsv to fill in the row & column names
# Can convert to .csv in Excel first
cells <- read.table(file = 'barcodes.tsv', sep = '\t', header=F)
genes <- read.table(file = 'genes.tsv', sep = '\t', header = F)

colnames(pbmc.data) = cells[1:2700,]
rownames(pbmc.data) = genes[,1]
# Due to duplicate gene names, column 1 may be preferred over column 2
# since they will be uniquely named

#### READING THE 10X FILES DIRECTLY ####

# Rather than reading in the above all manually, you can use Seurat to read
# in the three provided files together as a single 10X object
# NOTE: Code from (https://satijalab.org/seurat/v3.1/pbmc3k_tutorial.html)
library(dplyr)
library(Seurat)

# Load the PBMC dataset
# Set working directory to the downloaded folder "filtered_gene_bc_matrices"
pbmc.data <- Read10X(data.dir = "./hg19/")
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)

# At this point, you can extract the count data as a matrix with the fosellowing:
data.matrix = as.matrix(pbmc[["RNA"]]@data)

# This resulting "data.matrix" object will have all the same rows/genes as
# the dataset in the Dropbox folder