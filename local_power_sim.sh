#!/bin/bash

# Request 16 hours of runtime:
#SBATCH --time=16:00:00

# Use 46GB of Memory:
#SBATCH --mem=46GB

# Use 8 cores (can do up to 12 per matlab script; license doesn't support more than that)
#SBATCH -n 12

# Use contiguous cores
#SBATCH --contiguous

# Request 1 and only 1 node
#SBATCH --nodes=1-1

# Specify a job name:
#SBATCH -J Power_sim

# Specify an output file
#SBATCH -o Power_sim.out
#SBATCH -e Power_sim.out

# Email me if anything fun happens
#SBATCH --mail-type=ALL
#SBATCH --mail-user=nickolai_riabov@brown.edu 

# Run a matlab script called 'level_adjust_sim' in this same directory
cd /users/nriabov/McCloskey2014-Sims/Hausman
matlab-threaded -r "local_power_sim; exit"
date

