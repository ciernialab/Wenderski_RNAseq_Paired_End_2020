#!/bin/bash
#
#SBATCH -c 4
#SBATCH --mem-per-cpu=2000
#SBATCH --job-name=gunzip
#SBATCH --output=trimmed_unzip.out
#SBATCH --time=2:00:00

#######################################################################################
mkdir -p trimmed_unzip

for sample in `cat SRR_Acc_List.txt`
do

echo ${sample} "starting unzip"

# -c keep original files unchanged
gunzip -c trimmed/${sample}_1.paired.fastq.gz > trimmed_unzip/${sample}_1.paired.fastq

echo ${sample} "finished unzip read 1"

gunzip -c trimmed/${sample}_2.paired.fastq.gz > trimmed_unzip/${sample}_2.paired.fastq

echo ${sample} "finished unzip read 2"

done

