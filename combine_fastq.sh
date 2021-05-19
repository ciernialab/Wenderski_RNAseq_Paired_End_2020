#!/bin/bash
#
#SBATCH -c 1
#SBATCH --mem-per-cpu=1000
#SBATCH --job-name=CombineFastq
#SBATCH --output=CombineFastq.out
#SBATCH --time=1:00:00

#######################################################################################
#fastqc and then combine reports into one html and save to outputs folder
#######################################################################################
mkdir -p combined_fastq

xargs -n2 < ~/Wenderski_RNAseq/SRR_Acc_List.txt |
   while read a b ; do
      for c in 1 2 ; do

        echo "${a} read ${c} starting"

        cat SRA/${a}_${c}.fastq.gz SRA/${b}_${c}.fastq.gz > combined_fastq/${a}_${c}.fastq.gz

        echo "${a} read ${c} finished"
        
      done
   done

