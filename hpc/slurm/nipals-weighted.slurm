#!/bin/bash

# Partition for the job:
#SBATCH -p mig 

# The name of the job:
#SBATCH --job-name="Weighted nipals"

# The project ID which this job should run under:
#SBATCH --account="punim0613"

# Maximum number of tasks/CPU cores used by the job:
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10

# The amount of memory in megabytes per process in the job:
#SBATCH --mem=50000

# Send yourself an email when the job:
# aborts abnormally (fails)
#SBATCH --mail-type=FAIL
# begins
#SBATCH --mail-type=BEGIN
# ends successfully
#SBATCH --mail-type=END

# Use this email address:
#SBATCH --mail-user=ajalalabadi@unimelb.edu.au

# The maximum running time of the job in days-hours:mins:sec
#SBATCH --time=0-20:00:00

# The modules to load:
module load R/3.5.0-GCC-6.2.0
# for rmd files
module load pandoc/2.5

# The job command(s):
R --vanilla < ./src/render/render-nipals-weighted.R