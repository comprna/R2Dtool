#!/bin/bash

# written by AJ Sethi
# to install genomicRanges on Gadi
# uninitialzie conda before using

module load R
module load intel-compiler/2021.6.0
module load boost/1.79.0

export R_LIBS=/g/data/lf10/as7425/apps/Rlibs

R
library(GenomicFeatures)

BiocManager::install("GenomicFeatures")
