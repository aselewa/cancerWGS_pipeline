#!/bin/bash

# sbatch submission script to run main snakemake process. It then submits
# individual jobs from the compute node.

#SBATCH --job-name=snakemake
#SBATCH --output=snakelog.out
#SBATCH --time=24:00:00
#SBATCH --partition=xinhe
#SBATCH --account=pi-xinhe
#SBATCH --mem=2G
#SBATCH --tasks-per-node=1

source activate biotools

snakemake \
	-p \
	-j 12 \
	-r \
        --cluster-config cluster.json \
        --configfile configfile.yaml \
	-c "sbatch \
                --mem={cluster.mem} \
		--nodes={cluster.n} \
		--tasks-per-node={cluster.tasks} \
		--partition={cluster.partition} \
		--account={cluster.account} \
		--job-name={cluster.name} \
	        --output={cluster.logfile}"
