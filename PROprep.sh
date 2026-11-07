#!/bin/bash

#SBATCH --job-name=PROseq-R4-prep
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --output=R4_prep.txt
#SBATCh --mail-user=rwang4@fredhutch.org
#SBATCH --mail-type=END

###########################################################
# Step 1 of 3. Run once, before the array job (PROalign_array.sh).
# Unzips and renames fastq files, then runs fastqc.
# Run in a directory that has a folder named "fastq".
###########################################################

ml FastQC/0.12.1-Java-11

# Unzipping if needed
echo "unzipping..."
for FILE in fastq/*
    do
        if [[ "$FILE" == *.gz ]]
        then
            gunzip $FILE &
        fi
    done
wait

# Removing extra info from filenames.
echo "renaming if needed..."
for FILE in $(ls fastq/)
    do
        NEW=fastq/"$(echo "$FILE" |
            sed 's/^[0-9]\+_[0-9]\+_[0-9]\+_[0-9A-Z]\+_//' |
            sed 's/_[ATCG]\{6,8\}_/_/')"
        if [ ! -s "$NEW" ]
        then
            mv fastq/"$FILE" "$NEW"
        fi
    done

mkdir -p logs
mkdir -p logs/fastqc

### Running fastqc on files
echo "running fastqc if needed..."
for FILE in fastq/*.fastq
    do
        if [ ! -s logs/fastqc/"$(basename ${FILE/.fastq/_fastqc.zip})" ]
        then
            fastqc "$FILE" -o logs/fastqc --quiet  &
        fi
    done
wait

echo "pairs found:"
ls fastq | sed 's/_R[1-2].*//' | uniq
echo "number of pairs: $(ls fastq | sed 's/_R[1-2].*//' | uniq | wc -l)"
