# Reading 10X Genomics data in python

## Introduction

10X Genomics data is obtained from Droplet-based single-cell RNA sequencing (scRNA-seq) technologies that allow researchers to obtain transcriptome-wide expression profiles for thousands of cells at once. Each cell is encapsulated in a droplet in a oil-water emulsion, along with a bead containing reverse transcription primers with a unique barcode sequence. After reverse transcription inside the droplet, each cell's cDNA is labelled with the barcode (called "cell barcode"). Bursting of the droplets yields a pool of cDNA for library preparation and sequencing. Reading the barcode data of the sequences can then be performed to obtain the expression profile for each cell.

10X data are stored as sparse matrix. Using bar-coding technology and Cell Ranger 3.0, the data format for count data is `.mtx` which stores this sparse matrix as a column of row coordinates, a column of column corodinates, and a column of expression values > 0. If we look at the .mtx file you will see two header lines followed by a line detailing the total number of rows, columns and counts for the full matrix. Since only the coordinates are stored in the `.mtx` file, the names of each row & column must be stored separately in the “genes.tsv” and “barcodes.tsv” files respectively.

## Downloading the dataset

For the purpose of following tutorials, we will use Tabula Muris Data. Use the [shell script]( https://github.com/rahulbhadani/scRNAseq/blob/master/data/tabula_muris/download_tabula_muris.sh) to download the dataset.
