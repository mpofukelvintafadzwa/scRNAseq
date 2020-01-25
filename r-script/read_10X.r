
library("Matrix")
path = "/home/ivory/VersionControl/scRNASeq/data/tabula_muris/droplet/Kidney-10X_P4_5"
cellfile = "barcodes.tsv"
filepath = paste(path, "/", cellfile, sep = "")
cellbarcodes <- read.table(filepath)

genefile = "genes.tsv"
filepath = paste(path, "/", genefile, sep = "")
genenames <- read.table(filepath)

genefile = "matrix.mtx"
filepath = paste(path, "/", genefile, sep = "")

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
