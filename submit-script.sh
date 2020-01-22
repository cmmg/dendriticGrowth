#!/bin/sh
#This file is called submit-script.sh
#SBATCH --partition=univ # default "univ", if not specified
#SBATCH --time=2-12:00:00 # run time in days-hh:mm:ss
#SBATCH --nodes=2# require 2 nodes
#SBATCH --ntasks-per-node=16            # (by default, "ntasks"="cpus")
#SBATCH --mem-per-cpu=4000# RAM per CPU core, in MB (default 4 GB/core)
#SBATCH --job-name="cmake"
#SBATCH --job-name="takaki_explicit_grain"
#SBATCH --error=file%j.err
#SBATCH --output=file%j.out
#SBATCH --mail-user=kbhagat2@wisc.edu
#SBATCH --mail-type=ALL
#Make sure to change the above two lines to reflect your appropriate
# file locations for standard error and output

#Now list your executable command (or a string of them).
# Example for non-SLURM-compiled code:
source ~/.bashrc
cd /home/kbhagat2/workspace/takaki_explicit_grain/
rm -rf CMakeCache.txt cmake_install.cmake Makefile CMakeFiles 
cmake .
make release
mpirun -np 16 ./main
