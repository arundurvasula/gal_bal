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

## Parameters for simulations ##

### Introgression simulation

4N until MRCA of all S .cerevisae. 

Ne = 10,000,000

migration rate - m uniform(0,1)

talpha = uniform(0,40,000,000)

### Fixed parameters

Selfing rate = 1/50,000 (Kruglyak NG)

Clonning rate = 1/50,000 (use instead of selfing rate)

Gene conversion rate = (this is probably what yeast uses for the most part)

Recombination rate =  3.133483e-06 (BYxRM linkage mapping)

mutation rate = 3.8e-10

number of generations = 3.2 billion.

#### Notes:

To set the generation time as a variable, see here: https://groups.google.com/forum/#!topic/slim-discuss/0jOCti_t380

Use `-t` and `-m` to benchmark time and memory from the command line.
