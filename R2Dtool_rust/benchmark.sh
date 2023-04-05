#!/bin/bash

# variables
export input="/scratch/lf10/as7425/cheui.bed"
export in_gtf="/g/data/lf10/as7425/genomes/mouse_genome/GRCm39/Mus_musculus.GRCm39.104.chr.gtf"
export output="/scratch/lf10/as7425/R2Dtool_out.bed"

# set up path
module load R
export PATH="${PATH}:/home/150/as7425/R2Dtool/R2Dtool_rust/target/release"

##### liftover

# R
# time Rscript ~/R2Dtool/scripts/R2_lift.R ${input} ${in_gtf} ${output}

# R2Dtool
time R2Dtool_rust liftover -f gtf -H -g ${in_gtf} -i ${input} | tail -n +4 > ${output}
