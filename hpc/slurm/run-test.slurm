#!/bin/bash
# Created by the Melbourne Bioinformatics job script generator for SLURM
# Fri Nov 16 2018 16:28:07 GMT+1100 (AEDT)

# Partition for the job:
#SBATCH -p mig 

# The name of the job:
#SBATCH --job-name="Run Test"

# The project ID which this job should run under:
#SBATCH --account="punim0613"

# Maximum number of tasks/CPU cores used by the job:
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4

# The amount of memory in megabytes per process in the job:
#SBATCH --mem=8000

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
#SBATCH --time=0-00:59:00

# check that the script is launched with sbatch
if [ "x$SLURM_JOB_ID" == "x" ]; then
   echo "You need to submit your job to the queuing system with sbatch"
   exit 1
fi

# Run the job from the directory where it was launched (default)

# The modules to load:
module load R/3.5.0-GCC-6.2.0
module load pandoc

# The job command(s):
Rscript "src/render_testRmd.R"