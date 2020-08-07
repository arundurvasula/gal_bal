# gal_bal
Pipeline for ABC simulations of multilocus balancing selection and introgression

1. Run `bash init.sh` to setup the directory stucture
2. Run `qsub scripts/pipeline_introgression.sh` to get the simulations for introgression

Note: make sure SLiM is available on your path. You can install it via the following:

```
# This script is for SLiM 3.2.1. Might not work for newer version of SLiM.
qrsh -l h_data=8G
cd ~
rm SLiM.zip

#download SLiM
wget http://benhaller.com/slim/SLiM.zip
unzip SLiM.zip

#make build directory
mkdir slim_build
cd slim_build

#load newer version of gcc
module load gcc/4.9.3

#tell cmake which one to use and compile
CC=gcc CXX=g++ cmake ../SLiM
make slim
```
