#!/bin/bash
#SBATCH --job-name=rDNA_bowtie2_index
#SBATCH --nodes=10
#SBATCH --output=output.out

module load Bowtie2/2.5.4-GCC-13.2.0

cd /home/rwang4/peteslab/Redd/proseq/hg38_rDNA/

bowtie2-build hg38_rDNA.fa hg38_rDNA_index

