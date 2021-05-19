#!/bin/bash
#
#SBATCH -c 1
#SBATCH --mem-per-cpu=1000
#SBATCH --job-name=CombineFastq2
#SBATCH --output=CombineFastq2.out
#SBATCH --time=00:10:00

#######################################################################################
#combines 2 fastq files of the same sequencing library together into a single file
#this must be run in the direction containing the fastq files
#######################################################################################

# make new directory
mkdir -p combinedFiles

num=1
# in a for loop, cat every 2 files together, and place into the new directory
for i in $(seq 892 2 929)
do

echo "SRR11313${i} starting"
sum=$(($i + $num))
cat "SRR11313${i}_1.fastq.gz" "SRR11313${sum}_1.fastq.gz" > "combinedFiles/SRR11313${i}_1.fastq.gz"
cat "SRR11313${i}_2.fastq.gz" "SRR11313${sum}_2.fastq.gz" > "combinedFiles/SRR11313${i}_2.fastq.gz"
echo "SRR11313${i} finished"
done
