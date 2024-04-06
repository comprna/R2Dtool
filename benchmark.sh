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
time Rscript ~/R2Dtool/scripts/R2_lift.R ${input} ${in_gtf} ${output}

# R2Dtool
time R2Dtool_rust liftover -f gtf -H -g ${in_gtf} -i ${input} | tail -n +4 > ${output}

########## test liftover

# Define the path to your script and input file
export SCRIPT_PATH="~/R2Dtool/scripts/R2_lift.R"
export INPUT_FILE="/scratch/lf10/as7425/cheui.bed"
export in_gtf="/g/data/lf10/as7425/genomes/mouse_genome/GRCm39/Mus_musculus.GRCm39.104.chr.gtf"
export output="/scratch/lf10/as7425/R2Dtool_out.bed"

module load R

# Define the line counts you want to test
LINE_COUNTS=("10000" "25000" "50000" "120000" "250000" "500000" "1000000" "1835555")

# Extract the header line from the input file
HEADER=$(head -n 1 $INPUT_FILE)

# Run the script on the full input file and record the runtime
FULL_RUNTIME=$(time $SCRIPT_PATH $INPUT_FILE 2>&1 >/dev/null | grep "real" | awk '{print $2}')

# Create a table header
echo "Line Count | Runtime (seconds)" > /scratch/lf10/as7425/liftover_R_times.txt

# Loop through each line count and run the script on that subset of the input file
for count in "${LINE_COUNTS[@]}"; do
  # Extract the header line and a random subset of lines
  echo "$HEADER" > ${scratch}/in.bed
  shuf $INPUT_FILE | head -n $(($count - 1))) >> ${scratch}/in.bed
  head ${scratch}/in.bed

  # Run the script on the subset file and record the runtime
  SUBSET_RUNTIME=$(time Rscript ~/R2Dtool/scripts/R2_lift.R "${scratch}/in.bed" "${in_gtf}" "${output}" | grep "real" | awk '{print $2}')

  # Output the line count and runtime for this subset
  echo "$count | $SUBSET_RUNTIME" >> /scratch/lf10/as7425/liftover_R_times.txt
done

# Output the runtime for the full input file
echo "Full | $FULL_RUNTIME" >> /scratch/lf10/as7425/liftover_R_times.txt


########## test liftover

# Define the path to your script and input file
export SCRIPT_PATH="~/R2Dtool/scripts/R2_lift.R"
export INPUT_FILE="/scratch/lf10/as7425/cheui.bed"
export in_gtf="/g/data/lf10/as7425/genomes/mouse_genome/GRCm39/Mus_musculus.GRCm39.104.chr.gtf"
export output="/scratch/lf10/as7425/R2Dtool_out.bed"
export PATH="${PATH}:/home/150/as7425/R2Dtool/R2Dtool_rust/target/release"

# Define the line counts you want to test
LINE_COUNTS=("10000" "25000" "50000" "120000" "250000" "500000" "1000000" "1835555")

# Extract the header line from the input file
HEADER=$(head -n 1 $INPUT_FILE)

# Run the script on the full input file and record the runtime
FULL_RUNTIME=$(time $SCRIPT_PATH $INPUT_FILE 2>&1 >/dev/null | grep "real" | awk '{print $2}')

# Create a table header
echo "Line Count | Runtime (seconds)" > /scratch/lf10/as7425/liftover_R_times.txt

# Loop through each line count and run the script on that subset of the input file
for count in "${LINE_COUNTS[@]}"; do
  # Extract the header line and a random subset of lines
  echo "$HEADER" > ${scratch}/in.bed
  shuf $INPUT_FILE | head -n $(($count - 1))) >> ${scratch}/in.bed
  head ${scratch}/in.bed

  # Run the script on the subset file and record the runtime
  SUBSET_RUNTIME=$(time time R2Dtool_rust liftover -f gtf -H -g ${in_gtf} -i "${scratch}/in.bed"  > "${output}" | grep "real" | awk '{print $2}')

  # Output the line count and runtime for this subset
  echo "$count | $SUBSET_RUNTIME" >> /scratch/lf10/as7425/liftover_R_times.txt
done

# Output the runtime for the full input file
echo "Full | $FULL_RUNTIME" >> /scratch/lf10/as7425/liftover_R_times.txt




#########

cd /home/150/as7425/R2Dtool

time ./target/release/r2d annotate -f gtf -H -g ../test/GRCm39_subset.gtf -i ../test/out_CHEUI_modelII.bed | head

time Rscript /scripts/R2_annotate.R "../test/out_CHEUI_modelII.bed" "../test/GRCm39_subset.gtf" "./test_out.txt"


cargo build --release 

#########

# testing 2024/04/05
cd /home/150/as7425/R2Dtool/

# run unit tests 
export RUST_BACKTRACE=1
export RUST_LOG=debug 
cargo test -- --nocapture 



cd /home/150/as7425/R2Dtool/
rm -rf target
cargo build --release --future-incompat-report


time ./target/release/r2d liftover -H -g ./test/GRCm39_subset.gtf -i ./test/out_CHEUI_modelII.bed 





time ./target/release/r2d annotate -H -g ./test/GRCm39_subset.gtf -i ./test/out_CHEUI_modelII.bed | head


# R2Dtool
time ./target/release/r2d liftover -H -g ./test/GRCm39_subset.gtf -i ./test/out_CHEUI_modelII.bed | more

# R2Dtool
export RUST_BACKTRACE=1
cargo build --release --future-incompat-report


time ./target/release/r2d annotate -H -g ./test/GRCm39_subset.gtf -i ./test/out_CHEUI_modelII.bed | more