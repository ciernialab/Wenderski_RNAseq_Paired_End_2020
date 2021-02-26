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

cd ~/KCl_RNAseq/SRA/

echo ${sample} "starting SRR pull"
prefetch ${sample}

#######################################################################################
#paired end sra > Fastq.gz
#######################################################################################
echo ${sample} "starting dump"

fastq-dump --origfmt --split-files --gzip ${sample}

echo ${sample} "done"
#######################################################################################
#validate each SRA file and save to output log
#https://reneshbedre.github.io/blog/fqutil.html
#######################################################################################

echo ${sample} "starting sra check"

vdb-validate /alder/home/xlum/KCl_RNAseq/SRA/${sample} 2> output/SRA_checksum/${sample}_SRAcheck.log

echo ${sample} "done"

done

#combine logs into one output file
cat output/SRA_checksum/*_SRAcheck.log > output/SRA_checksum/SRAcheck.log
