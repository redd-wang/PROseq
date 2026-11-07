#!/bin/bash

#SBATCH --job-name=PROseq-R4
#SBATCH --nodes=1
#SBATCH --array=0-11
#SBATCH --cpus-per-task=8
#SBATCH --output=R4_%A_%a.txt
#SBATCh --mail-user=rwang4@fredhutch.org
#SBATCH --mail-type=END

###########################################################
# Step 2 of 3. Job array, one task per read pair (12 pairs -> 0-11).
# Run PROprep.sh first. Each task handles the Nth pair (sorted by name):
# trim + rRNA filter, spike-in align, genome align, UMI dedup, bigWigs.
# Run in the directory that has the "fastq" folder.
###########################################################

## Parameters
THREADS=${SLURM_CPUS_PER_TASK} # Threads follow what Slurm actually gave this task
UMI_LEN=6  # Length of UMI in basepairs

## UMI Flags (set to Y or N as appropriate)
FIVEP_UMI="Y" # Is there a UMI on the 5' end of the read?
THREEP_UMI="Y" # Is there a UMI on the 3' end of the read?

## Adaptor sequences to clip. Default = Tru-Seq small RNA
ADAPTOR_1="TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC"
ADAPTOR_2="GATCGTCGGACTGTAGAACTCTGAACGTGTAGATCTCGGTGGTCGCCGTATCATT"

## Genomes. Fill in paths.
GENOME_EXP="/shared/biodata/reference/iGenomes/Homo_sapiens/UCSC/hg38/Sequence/Bowtie2Index/genome"
GENOME_SPIKE="/home/rwang4/peteslab/Redd/reference/dm6hg38_spikegenome/spikedgenome" ## USE REPEAT MASKED VERSION!!
SPIKE_PREFIX="spike" ## This is the prefix you've used on your spike in chromosomes, ie >hg38chr1
RDNA="/home/rwang4/peteslab/Redd/reference/rDNA2/rDNA2"

## Mapq value for filtering multimappers
MAPQ=5

ml FastQC/0.12.1-Java-11
ml MAGeCK-VISPR/0.5.5-Python-3.7.4
ml deepTools/3.5.4.post1-gfbf-2022b
ml fhPython/3.9.6-foss-2021b
ml Bowtie2/2.5.1-GCC-12.2.0
ml soap-hla/1.0.0-GCCcore-12.2.0
ml QIIME2/2020.11
ml nullarbor/2.0.20191013
ml SAMtools/1.17-GCC-12.2.0
ml fastp

### Picking this task's pair
PAIR=$(ls fastq | sed 's/_R[1-2].*//' | uniq | sed -n "$((SLURM_ARRAY_TASK_ID+1))p")
if [ -z "$PAIR" ]
then
    echo "no pair for array index ${SLURM_ARRAY_TASK_ID}. exiting..."
    exit 1
fi
echo "task ${SLURM_ARRAY_TASK_ID} processing ${PAIR} with ${THREADS} threads"

mkdir -p logs/fastp logs/rRNA logs/spikeAlign logs/align logs/deDup logs/spikedeDup
mkdir -p trimmedFastq spikeBAM BAM BAMdeDuped spikeBAMdeDuped bw

### Trimming adapters and filtering rRNA reads
# fastp UMI location depends on which ends carry a UMI
if [[ $THREEP_UMI == "Y" && $FIVEP_UMI == "Y" ]]
then
    UMI_ARGS="--umi --umi_loc=per_read --umi_len=${UMI_LEN}"
elif [[ $THREEP_UMI == "Y" ]]
then
    UMI_ARGS="--umi --umi_loc=read1 --umi_len=${UMI_LEN}"
elif [[ $FIVEP_UMI == "Y" ]]
then
    UMI_ARGS="--umi --umi_loc=read2 --umi_len=${UMI_LEN}"
else
    UMI_ARGS=""
fi

if [ ! -s trimmedFastq/${PAIR}_R1.fastq ]
then
    echo "trimming adapters and filtering rRNA reads for ${PAIR}"
    (fastp \
    -i fastq/${PAIR}_R1.fastq \
    -I fastq/${PAIR}_R2.fastq \
    --adapter_sequence $ADAPTOR_1 \
    --adapter_sequence_r2 $ADAPTOR_2 \
    $UMI_ARGS \
    --stdout \
    --html logs/fastp/${PAIR}_fastp.html \
    -w $(echo ${THREADS}/3 | bc) \
    -c \
    --overlap_len_require 15 2> logs/fastp/${PAIR}_fastp.log) |
    (bowtie2 \
    --fast-local \
    --un-conc trimmedFastq/${PAIR}.fastq \
    --interleaved - \
    -x ${RDNA} \
    --threads $(echo ${THREADS}/3*2 | bc) 2> logs/rRNA/${PAIR}_rRNA_bowtie.log) > /dev/null

    # bowtie automatically names PE --un output ${PAIR}.1.fastq / ${PAIR}.2.fastq
    mv trimmedFastq/${PAIR}.1.fastq trimmedFastq/${PAIR}_R1.fastq
    mv trimmedFastq/${PAIR}.2.fastq trimmedFastq/${PAIR}_R2.fastq
fi


### Aligning to spike in genome to get normalization factors
module load Python/3.10.8-GCCcore-12.2.0
module load SciPy-bundle/2023.02-gfbf-2022b

export LD_LIBRARY_PATH=/app/software/Python/3.11.3-GCCcore-12.3.0/lib:$LD_LIBRARY_PATH

if [ ! -s "spikeBAM/${PAIR}.BAM" ]
then
    echo "aligning ${PAIR} to spike in genome"
    (bowtie2 \
    --local \
    --very-sensitive-local \
    --threads $(echo ${THREADS}/3*2 | bc) \
    --no-unal \
    --no-mixed \
    --no-discordant \
    -x "$GENOME_SPIKE" \
    -1 "trimmedFastq/${PAIR}_R1.fastq" \
    -2 "trimmedFastq/${PAIR}_R2.fastq" \
    2> logs/spikeAlign/${PAIR}_spikeAlign.log) |
    samtools view -hS -f 2 -q ${MAPQ} |
    perl -n -e 'print $_ if (/^\@/ || /'${SPIKE_PREFIX}'/ ) ' |
    samtools view -b |
    samtools sort -@ $(echo ${THREADS}/3 | bc) -o spikeBAM/${PAIR}.BAM
    samtools index spikeBAM/${PAIR}.BAM
fi


### Aligning to experimental genome
if [ ! -s "BAM/${PAIR}.BAM" ]
then
    echo "aligning ${PAIR} to experimental genome"
    (bowtie2 \
    --local \
    --sensitive-local \
    --threads $(echo ${THREADS}/3*2 | bc) \
    -x "$GENOME_EXP" \
    -1 "trimmedFastq/${PAIR}_R1.fastq" \
    -2 "trimmedFastq/${PAIR}_R2.fastq" \
    2> logs/align/${PAIR}_align.log) |
    samtools view -bS -f 2 -q ${MAPQ} |
    samtools sort -@ $(echo ${THREADS}/3 | bc) -o BAM/${PAIR}.BAM
    samtools index BAM/${PAIR}.BAM
fi


### deduplicating with UMIs
echo "deduplicating ${PAIR} with UMIs"
module unload SciPy-bundle
module load SciPy-bundle/2021.10-foss-2021b
module load SAMtools/1.19.2-GCC-13.2.0

export LD_LIBRARY_PATH=/app/software/Python/3.9.6-GCCcore-11.2.0/lib:$LD_LIBRARY_PATH

if [ ! -s "BAMdeDuped/${PAIR}_deDuped.BAM" ]
then
    (umi_tools dedup \
    -I "BAM/${PAIR}.BAM" \
    --umi-separator=":" \
    --paired \
    -S "BAMdeDuped/${PAIR}_deDuped.BAM" \
    )> "logs/deDup/${PAIR}_deDup.log" &&
    samtools index "BAMdeDuped/${PAIR}_deDuped.BAM"
fi

if [ ! -s "spikeBAMdeDuped/${PAIR}_deDuped.BAM" ]
then
    (
    umi_tools dedup \
    -I "spikeBAM/${PAIR}.BAM" \
    --paired \
    --umi-separator=":" \
    -S "spikeBAMdeDuped/${PAIR}_deDuped.BAM" \
    )> "logs/spikedeDup/${PAIR}_deDup.log" &&
    samtools index "spikeBAMdeDuped/${PAIR}_deDuped.BAM"
fi


# Making non-normalized bigWig files
FILE="BAMdeDuped/${PAIR}_deDuped.BAM"
if [ ! -s "bw/${PAIR}_deDuped_fwd.bw" ]
then
    bamCoverage \
    --bam $FILE \
    --skipNonCoveredRegions \
    --outFileName bw/${PAIR}_deDuped_fwd.bw \
    --binSize 1 \
    --numberOfProcessors ${THREADS} \
    --normalizeUsing None \
    --Offset 1 \
    --samFlagInclude 82
fi
if [ ! -s "bw/${PAIR}_deDuped_rev.bw" ]
then
    bamCoverage \
    --bam $FILE \
    --skipNonCoveredRegions \
    --outFileName bw/${PAIR}_deDuped_rev.bw \
    --binSize 1 \
    --numberOfProcessors ${THREADS} \
    --normalizeUsing None \
    --Offset 1 \
    --samFlagInclude 98
fi

echo "task ${SLURM_ARRAY_TASK_ID} (${PAIR}) done"
