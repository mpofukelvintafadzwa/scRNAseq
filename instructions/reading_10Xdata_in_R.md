# Reading 10X Genomics data in R

## Introduction

10X Genomics data is obtained from Droplet-based single-cell RNA sequencing (scRNA-seq) technologies that allow researchers to obtain transcriptome-wide expression profiles for thousands of cells at once. Each cell is encapsulated in a droplet in a oil-water emulsion, along with a bead containing reverse transcription primers with a unique barcode sequence. After reverse transcription inside the droplet, each cell's cDNA is labelled with the barcode (called "cell barcode"). Bursting of the droplets yields a pool of cDNA for library preparation and sequencing. Reading the barcode data of the sequences can then be performed to obtain the expression profile for each cell.

10X data are stored as sparse matrix. Using bar-coding technology and Cell Ranger 3.0, the data format for count data is `.mtx` which stores this sparse matrix as a column of row coordinates, a column of column corodinates, and a column of expression values > 0. If we look at the .mtx file you will see two header lines followed by a line detailing the total number of rows, columns and counts for the full matrix. Since only the coordinates are stored in the `.mtx` file, the names of each row & column must be stored separately in the “genes.tsv” and “barcodes.tsv” files respectively.

## Downloading the dataset

For the purpose of following tutorials, we will use Tabula Muris Data. Use the [shell script]( https://github.com/rahulbhadani/scRNAseq/blob/master/data/tabula_muris/download_tabula_muris.sh) to download the dataset.


## Reading 10X data in R
1. To read the sparse matrix, we use Matrix package in R

```R
library("Matrix")
```
Now, we read the data. Note that since 10X data are stored as sparse matrix, row and coloumn information are stored separately.

```R
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
```

At this point, we want to store sparse matrix `molecules` as a csv file

```R
moleculesdata = as.data.frame(as.matrix(molecules))
rownames(moleculesdata) <- rownames(molecules)
colnames(moleculesdata) <- colnames(molecules)
```

We can save `moleculesdata` as csvn file in the following manner:

```R
write.csv(moleculesdata, paste(path, "/", "moleculatesdata.csv", sep=""))
```


Most of the time, we have only .mtx sparse matrix file that contains unique molecular identified (UMI) coun and .tsv file that has columns (containing cell names) and rows (containing gene names) as is the case with [PBMC3](https://s3-us-west-2.amazonaws.com/10x.files/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz) dataset. In that scenario, only above code is sufficient to read data. However, if we have datasets divided in .different pieces with annotations and metadata detailing the experiment under study, we might want to do something more.



We also need to get computational annotations and metadata for the data.

```R
metadatapath = "/home/ivory/VersionControl/scRNASeq/data/tabula_muris"
metadatafile = "droplet_metadata.csv"
metafilepath = paste(metadatapath, "/", metadatafile, sep="")
    
metadata <- read.delim(metafilepath, sep=",", header = TRUE)
```
Lets see how the metadata looks like:

```R
head(metadata)
```

Output:
```bash
 channel mouse.id  tissue   subtissue mouse.sex
1 10X_P4_0    3-M-8  Tongue        <NA>         M
2 10X_P4_1    3-M-9  Tongue        <NA>         M
3 10X_P4_2  3-M-8/9   Liver hepatocytes         M
4 10X_P4_3    3-M-8 Bladder        <NA>         M
5 10X_P4_4    3-M-9 Bladder        <NA>         M
6 10X_P4_5    3-M-8  Kidney        <NA>         M
```
Since we are using 10X_P4_5, we need to retreive metadata for 10X_P4_5

```R
metadata[metadata$channel == "10X_P4_5",]
```

```
channel mouse.id tissue subtissue mouse.sex
6 10X_P4_5    3-M-8 Kidney      <NA>         M
```


Now we will read the droplet annotations:

```R
anndatafile = "droplet_annotation.csv"
annfilepath = paste(metadatapath, "/", anndatafile, sep="")
annadata <- read.delim(annfilepath, sep=",", header = TRUE)
head(anndata)
```

Output:
```bash
                       cell  tissue cell_ontology_class                    cell_ontology_term_iri cell_ontology_id
1 10X_P4_3_AAAGTAGAGATGCCAG Bladder    mesenchymal cell http://purl.obolibrary.org/obo/CL_0008019       CL:0008019
2 10X_P4_3_AACCGCGTCCAACCAA Bladder    mesenchymal cell http://purl.obolibrary.org/obo/CL_0008019       CL:0008019
3 10X_P4_3_AACTCCCGTCGGGTCT Bladder    mesenchymal cell http://purl.obolibrary.org/obo/CL_0008019       CL:0008019
4 10X_P4_3_AACTCTTAGTTGCAGG Bladder        bladder cell http://purl.obolibrary.org/obo/CL_1001319       CL:1001319
5 10X_P4_3_AACTCTTTCATAACCG Bladder    mesenchymal cell http://purl.obolibrary.org/obo/CL_0008019       CL:0008019
6 10X_P4_3_AAGACCTAGATCCGAG Bladder    endothelial cell http://purl.obolibrary.org/obo/CL_0000115       CL:0000115
```
If you look at the output of `head(cellbarcodes)`:

```bash
                  V1
1 AAACCTGAGATGCCAG-1
2 AAACCTGAGTGTCCAT-1
3 AAACCTGCAAGGCTCC-1
4 AAACCTGTCCTTGCCA-1
5 AAACGGGAGCTGAACG-1
6 AAACGGGCAGGACCCT-1
```
there is observable difference in formatting - we want to make sure that cell type names match from cell names and from annotations to easily access metadata and annotations.

```R
ann[,1] <- paste(anndata[,1], "-1", sep="")
ann_subset <- anndata[match(colnames(molecules), anndata[,1]),]
celltype <- ann_subset[,3]
```

Now lets build the cell-metadata dataframe for specific mouse ID "3_8_M" corresponding to 10X_P4_5, which you can verify above from the output of `metadata[metadata$channel == "10X_P4_5",]`

```R
mouseID = "3_8_M"

cell_anns <- data.frame(mouse = rep(mouseID, times=ncol(molecules)), type=celltype)
rownames(cell_anns) <- colnames(molecules)
head(cell_ans)
```

Output:

```bash
                           mouse                        type
10X_P4_5_AAACCTGAGATGCCAG-1 3_8_M                   leukocyte
10X_P4_5_AAACCTGAGTGTCCAT-1 3_8_M            fenestrated cell
10X_P4_5_AAACCTGCAAGGCTCC-1 3_8_M          smooth muscle cell
10X_P4_5_AAACCTGTCCTTGCCA-1 3_8_M          kidney tubule cell
10X_P4_5_AAACGGGAGCTGAACG-1 3_8_M kidney collecting duct cell
10X_P4_5_AAACGGGCAGGACCCT-1 3_8_M          kidney tubule cell
```


--

<sup>*[web.archive.org](web.archive.org) version backup available.</sup>

## Author
- Rahul Bhadani (rahulbhadani@email.arizona.edu)
- With the help from Dr. Lingling An's Lab Members of Statistical and Bioinformatics Group
