#!/bin/sh
#This script will get Single Cell Explorer up and running

sudo apt-get update
sudo apt-get install -y build-essential
sudo apt-get install -y ssh libssl-dev libffi-dev libxml2-dev libxslt1-dev zlib1g-dev zip unzip libfftw3-dev libcurl3 openssl
sudo apt-get install -y python3 python3-pip python3-dev
cd ~
mkdir -p mongodb
cd mongodb
wget https://fastdl.mongodb.org/linux/mongodb-linux-x86_64-ubuntu1604-4.0.10.tgz 
tar -zxvf mongodb-linux-x86_64-ubuntu1604-4.0.10.tgz 
mkdir scdb
mkdir log
./mongodb-linux-x86_64-ubuntu1604-4.0.10/bin/mongod --dbpath "./scdb" --port 27017 --wiredTigerCacheSizeGB 1 --fork --logpath "./log/scdb.log"
cd ~/mongodb/
wget http://54.159.6.229/downloads/scDB.zip
unzip scDB.zip
mkdir dumpfiles 
mv scDB dumpfiles 
mongodb-linux-x86_64-ubuntu1604-4.0.10/bin/mongorestore dumpfiles

pip3 install django==2.2 gunicorn jupyter leidenalg numpy pandas pymongo sklearn scanpy torchvision  
cd ~ 
mkdir singleCell 
cd singleCell 
wget http://54.159.6.229/downloads/singleCellExplorer.zip 
unzip singleCellExplorer.zip 
cd singleCellExplorer 

wget http://54.159.6.229/downloads/scpipeline.py

#You can use this command to run things once the installation is complete
#sudo python3 manage.py runserver 0.0.0.0:80
