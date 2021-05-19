#!/bin/bash
#
#SBATCH -c 2
#SBATCH --mem-per-cpu=2000
#SBATCH --job-name=SRAfetch
#SBATCH --output=SRAfetch.out
#SBATCH --time=12:00:00


#KCl_RNAseq analysis
#pull SRA files from GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE107436
#######################################################################################
# run from KCl_RNAseq directory
#puts data in shared data file for the experiment
#######################################################################################
mkdir -p SRA/
mkdir -p output/SRA_checksum

for sample in `cat SRR_Acc_List.txt`
do

cd ~/Wenderski_RNAseq_redo/SRA/

echo ${sample} "starting SRR pull"
prefetch ${sample}

#######################################################################################
#paired end sra > Fastq.gz
#######################################################################################

echo ${sample} "starting dump"

fastq-dump --origfmt --split-e --gzip ${sample}

echo ${sample} "done"

done
