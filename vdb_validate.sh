#!/bin/bash
#
#SBATCH -c 2
#SBATCH --mem-per-cpu=2000
#SBATCH --job-name=vdb_validate
#SBATCH --output=vdb_validate.out
#SBATCH --time=12:00:00



for sample in `cat SRR_Acc_List.txt`
do

echo ${sample} "starting vdb_validate"

#######################################################################################
#validate each SRA file and save to output log
#https://reneshbedre.github.io/blog/fqutil.html
#######################################################################################

vdb-validate /alder/home/xlum/Wenderski_RNAseq/SRA/${sample} 2> output/SRA_checksum/${sample}_SRAcheck.log

echo ${sample} "done"

done

#combine logs into one output file
cat output/SRA_checksum/*_SRAcheck.log > output/SRA_checksum/SRAcheck.log