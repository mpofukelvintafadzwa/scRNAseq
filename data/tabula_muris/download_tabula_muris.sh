#!/bin/bash

# Download FACS data
wget https://ndownloader.figshare.com/files/10038307

unzip 10038307

wget https://ndownloader.figshare.com/files/10038310

mv 10038310 FACS_metadata.csv

wget https://ndownloader.figshare.com/files/10039267

mv 10039267 FACS_annotations.csv

# Download 10X data
wget https://ndownloader.figshare.com/files/10038325
unzip 10038325
wget https://ndownloader.figshare.com/files/10038328
mv 10038328 droplet_metadata.csv
wget https://ndownloader.figshare.com/files/10039264
mv 10039264 droplet_annotation.csv

