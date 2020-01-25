#!/usr/bin/env python

import scanpy
import pandas
import csv

# This python script lets read 10X Genomics using scanpy

path = "/home/ivory/VersionControl/scRNASeq/data/tabula_muris/droplet/Kidney-10X_P4_5" # Path where .mtx file and .csv files will be located
cellfile = "barcodes.tsv" # Cell files
genefile = "genes.tsv" # Gene file
UMIfile = "matrix.mtx" #Droplet based Unique Molecule Identifier count obtained from 10X Genomics Data

cell_file_path = path + "/" + cellfile
gene_file_path = path + "/" + genefile
UMI_file_path = path + "/" + UMIfile

# Scanpy doesn't require pass cellfile, genefile and UMIfile separately as long as names of the files are barcodes.tsv, genes.tsv and matrix.mtx - pasing path of the folder where lie is sufficient
anndata = scanpy.read_10x_mtx(path, var_names='gene_symbols')

## get the UMI count data
count = anndata.X ## This is sparse matrix

# Since Count data is sparse matrix, we need to convert to dataframe
# however, it should be noted that the dense matrix as a result of conversion
# is the transpose of the one you get by doing similar conversion in R.
count =pandas.DataFrame(count.toarray())

## get the gene info
genes = anndata.var_names

## Get the cell info
cells = anndata.obs_names

count.columns =anndata.var_names
count.index = anndata.obs_names

# Now transpose so that gene names are same as row names and cell types are same  column names
count = count.T

count = count.astype(int)


count.to_csv(path + '/count.csv', quoting=csv.QUOTE_NONE)

