# Reading 10X Genomics data in python

## Introduction

10X Genomics data is obtained from Droplet-based single-cell RNA sequencing (scRNA-seq) technologies that allow researchers to obtain transcriptome-wide expression profiles for thousands of cells at once. Each cell is encapsulated in a droplet in a oil-water emulsion, along with a bead containing reverse transcription primers with a unique barcode sequence. After reverse transcription inside the droplet, each cell's cDNA is labelled with the barcode (called "cell barcode"). Bursting of the droplets yields a pool of cDNA for library preparation and sequencing. Reading the barcode data of the sequences can then be performed to obtain the expression profile for each cell.

10X data are stored as sparse matrix. Using bar-coding technology and Cell Ranger 3.0, the data format for count data is `.mtx` which stores this sparse matrix as a column of row coordinates, a column of column corodinates, and a column of expression values > 0. If we look at the .mtx file you will see two header lines followed by a line detailing the total number of rows, columns and counts for the full matrix. Since only the coordinates are stored in the `.mtx` file, the names of each row & column must be stored separately in the “genes.tsv” and “barcodes.tsv” files respectively.

## Downloading the dataset

For the purpose of following tutorials, we will use Tabula Muris Data. Use the [shell script]( https://github.com/rahulbhadani/scRNAseq/blob/master/data/tabula_muris/download_tabula_muris.sh) to download the dataset.

## Reading 10X data in R

We will use scanpy package to read 10X genomics data. You can install scanpy as

```python
pip install scanpy
```

```python
import scanpy
import pandas
import csv
path = "/home/ivory/VersionControl/scRNASeq/data/tabula_muris/droplet/Kidney-10X_P4_5" # Path where .mtx file and .csv files will be located
cellfile = "barcodes.tsv" # Cell files
genefile = "genes.tsv" # Gene file
UMIfile = "matrix.mtx" #Droplet based Unique Molecule Identifier count obtained from 10X Genomics Data

cell_file_path = path + "/" + cellfile
gene_file_path = path + "/" + genefile
UMI_file_path = path + "/" + UMIfile
```
Scanpy doesn't require pass cellfile, genefile and UMIfile separately as long as names of the files are barcodes.tsv, genes.tsv and matrix.mtx - pasing path of the folder where lie is sufficient. Hence, above information is rendudant for our case, but I have provided them here so that you can know what files we are reading.

```python
anndata = scanpy.read_10x_mtx(path, var_names='gene_symbols')

## get the UMI count data
count = anndata.X ## This is sparse matrix
```
Count data is stored in sparse matrix format that you can check by typing `type(count)`. We will convert sparse matrix to pandas data frame.  However, it should be noted that the dense matrix as a result of conversion is the transpose of the one you get by doing similar conversion in R.

```python
count =pandas.DataFrame(count.toarray())

## get the gene info
genes = anndata.var_names

## Get the cell info
cells = anndata.obs_names

count.columns =anndata.var_names
count.index = anndata.obs_names
```

 Now we transpose the dataframe so that gene names are same as row names and cell types are same  column names.
 
 ```python
count = count.T
count = count.astype(int)
```

Now we can save the dataframe as a csv file for further processing. Note we converted UMI entries to int as without conversion, pandas dataframe saves them as float which requires larger memory.

```python
count.to_csv(path + '/moleculercount.csv', quoting=csv.QUOTE_NONE)
```
