#!/bin/bash
#SBATCH --job-name="Singularity MCMC"
#SBATCH --output="mcmc.%j.%N.out"
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=9
#SBATCH -t 10:00:00

#Run the job
module load tacc-singularity

cd /work/04365/pablod/2_raj/BRAIN/SCFC-spectral-python/SCFC

singularity exec /work/04365/pablod/2_raj/ubuntu_mpi.img /bin/bash -c 'echo "NNN" | python3 run_MCMC_fromPipe.py'
