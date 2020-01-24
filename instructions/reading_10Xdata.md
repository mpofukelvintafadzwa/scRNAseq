# Reading 10X Genomics data in R and Python

## Introduction

10X Genomics data is obtained from Droplet-based single-cell RNA sequencing (scRNA-seq) technologies that allow researchers to obtain transcriptome-wide expression profiles for thousands of cells at once. Each cell is encapsulated in a droplet in a oil-water emulsion, along with a bead containing reverse transcription primers with a unique barcode sequence. After reverse transcription inside the droplet, each cell's cDNA is labelled with the barcode (called "cell barcode"). Bursting of the droplets yields a pool of cDNA for library preparation and sequencing. Reading the barcode data of the sequences can then be performed to obtain the expression profile for each cell.

10X data are stored as sparse matrix. Using bar-coding technology and Cell Ranger 3.0, the data format for count data is .mtx which stores this sparse matrix as a column of row coordinates, a column of column corodinates, and a column of expression values > 0. If we look at the .mtx file you will see two header lines followed by a line detailing the total number of rows, columns and counts for the full matrix. Since only the coordinates are stored in the .mtx file, the names of each row & column must be stored separately in the “genes.tsv” and “barcodes.tsv” files respectively.

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
```

--
<sup>*[web.archive.org](web.archive.org) version backup available.</sup>

## Author
- Rahul Bhadani ( rahulbhadani@email.arizona.edu)

