#!/bin/bash

#SBATCH --job-name=PROseq-R4-summary
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --output=R4_summary.txt
#SBATCH --mail-user=rwang4@fredhutch.org
#SBATCH --mail-type=END

###########################################################
# Step 3 of 3. Run after every array task has finished.
# Builds info/infoTable.tsv from the logs of all pairs.
###########################################################

### Generating infoTable
mkdir -p info

if [ ! -s info/infoTable.tsv ]
then
touch info/infoTable.tsv
echo -e Name'\t'RawReads'\t'NonDimerReads'\t'%dimer'\t'insertSize'\t'rRNAreads'\t'%rRNA'\t'passedFilters'\t'\
    bowtieConcordant'\t'bowtieMulti'\t'bowtieUnal'\t'bowtieOverallMap%'\t'bowtieConcordant%'\t'\
    bowtieMulti%'\t'bowtieUnal%'\t'uniqueMapped'\t'uniqueMappedNondup'\t'%PCRdups'\t'uniqueMappedSpikein'\t'\
    uniqueMappedSpikeinNondup'\t'spikeInPCRdups% >> info/infoTable.tsv

for SAMPLE in $(ls BAM/*.BAM | sed 's/.BAM//' | sed 's/BAM\///' )
do
    NAME=${SAMPLE}
    RAW_READS=$(cat logs/fastp/${SAMPLE}_fastp.log |
                grep "total reads:" | head -n 1 |
                awk '{print $3}')
    TRIMMED_READS=$(cat logs/fastp/${SAMPLE}_fastp.log |
                grep "total reads:" | tail -n 1 |
                awk '{print $3}')
    PER_DIMER=$(echo "(1-"${TRIMMED_READS}"/"${RAW_READS}")*100" | bc -l)%
    INSERT_SIZE=$(cat logs/fastp/${SAMPLE}_fastp.log |
                grep "Insert size peak" |
                awk '{print $8}')
    PASSED_FILTERS=$(cat logs/align/${SAMPLE}_align.log |
                grep "reads; of these:$" |
                awk '{print $1}')
    RRNA=$(echo ${TRIMMED_READS}"-"${PASSED_FILTERS} | bc )
    PER_RRNA=$(echo ${RRNA}"/"${RAW_READS}"*100" | bc -l)%
    B_CONC=$(cat logs/align/${SAMPLE}_align.log |
            grep "aligned concordantly exactly 1 time$" |
            awk '{print $1}')
    B_MULTI=$(cat logs/align/${SAMPLE}_align.log |
            grep "aligned concordantly >1 times$" |
            awk '{print $1}')
    B_UNAL=$(cat logs/align/${SAMPLE}_align.log |
            grep "aligned concordantly 0 times$" |
            awk '{print $1}')
    B_OAP=$(cat logs/align/${SAMPLE}_align.log |
            grep "overall alignment rate$" |
            awk '{print $1}')
    B_CONC_PER=$(echo ${B_CONC}"/"${PASSED_FILTERS}"*100" | bc -l)%
    B_MULTI_PER=$(echo ${B_MULTI}"/"${PASSED_FILTERS}"*100" | bc -l)%
    B_UNAL_PER=$(echo ${B_UNAL}"/"${PASSED_FILTERS}"*100" | bc -l)%
    UNIQ_MAPPED=$(cat logs/deDup/${SAMPLE}_deDup.log |
            grep "Input Reads:" | awk '{print $10}')
    UNIQ_MAPPED_DEDUP=$(cat logs/deDup/${SAMPLE}_deDup.log |
            grep "Number of reads out:" | awk '{print $8}')
    PER_DUPS=$(echo "(1-"${UNIQ_MAPPED_DEDUP}"/"${UNIQ_MAPPED}")*100" | bc -l)%
    UNIQ_MAPPED_SPIKE=$(cat logs/spikedeDup/${SAMPLE}_deDup.log |
            grep "Input Reads:" | awk '{print $10}')
    UNIQ_MAPPED_DEDUP_SPIKE=$(cat logs/spikedeDup/${SAMPLE}_deDup.log |
            grep "Number of reads out:" | awk '{print $8}')
    PER_DUPS_SPIKE=$(echo "(1-"${UNIQ_MAPPED_DEDUP_SPIKE}"/"${UNIQ_MAPPED_SPIKE}")*100" | bc -l)%

    echo -e $NAME'\t'\
    $RAW_READS'\t'\
    $TRIMMED_READS'\t'\
    $PER_DIMER'\t'\
    $INSERT_SIZE'\t'\
    $RRNA'\t'\
    $PER_RRNA'\t'\
    $PASSED_FILTERS'\t'\
    $B_CONC'\t'\
    $B_MULTI'\t'\
    $B_UNAL'\t'\
    $B_OAP'\t'\
    $B_CONC_PER'\t'\
    $B_MULTI_PER'\t'\
    $B_UNAL_PER'\t'\
    $UNIQ_MAPPED'\t'\
    $UNIQ_MAPPED_DEDUP'\t'\
    $PER_DUPS'\t'\
    $UNIQ_MAPPED_SPIKE'\t'\
    $UNIQ_MAPPED_DEDUP_SPIKE'\t'\
    $PER_DUPS_SPIKE  >> info/infoTable.tsv

done
fi
