#!/bin/bash

# Request 16 hours of runtime:
#SBATCH --time=16:00:00

# Use 60GB of Memory:
#SBATCH --mem=60GB

# Use 12 cores (can do up to 12 per matlab script; CCV license doesn't support more than that)
#SBATCH -n 12

# Request 1 and only 1 node (enables parfor to use all the cores on one machine)
#SBATCH --nodes=1-1
#SBATCH --contiguous

# Specify a job name:
#SBATCH -J Power_sim

# Specify an output file with an informative name that catches all the matlab output
#SBATCH -o Power_sim.out
#SBATCH -e Power_sim.out

# Email me if anything fun or interesting happens
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nickolai_riabov@brown.edu 

# Run a matlab script called 'sample_script' in the directory specified
cd /users/nriabov/McCloskey2014-Sims/Hausman
matlab-threaded -r "sample_script; exit"
date

