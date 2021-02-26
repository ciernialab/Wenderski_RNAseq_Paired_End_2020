#!/bin/bash
#
#SBATCH -c 4
#SBATCH --mem-per-cpu=2000
#SBATCH --job-name=STARalignment
#SBATCH --output=STARalign.out
#SBATCH --time=6:00:00

#######################################################################################
#STAR alignment: https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf

#######################################################################################
# -p no error if existing, make parent directories as needed
mkdir -p star_out

module load gcc/6.2.0 star2.7.4a

for sample in `cat SRR_Acc_List.txt`
do

echo ${sample} "starting unzip fastq.gz"

gunzip -c {sample} > ${sample}.fastq

echo ${sample} "finished unzip fastq.gz"

echo ${sample} "starting STAR alignment"

# allow 10000 files to be open at once
ulimit -n 10000

STAR --runThreadN 12 \
--genomeDir $STAR_LIBS/mm10/star_indices/ \
--readFilesIn ${sample}_pass_1.fastq ${sample}_pass_2.fastq \
--outFilterType BySJout \
--outFilterMultimapNmax 20 \
--alignSJoverhangMin 8 \
--alignSJDBoverhangMin 1 \
--outFilterMismatchNmax 999 \
--outFilterMismatchNoverLmax 0.1 \
--alignIntronMin 20 --alignIntronMax 1000000 \
--alignMatesGapMax 1000000 \
--outSAMattributes NH HI nM AS MD \
--outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix star_out/${sample}

echo ${sample} "finished STARalignment"

done


#######################################################################################
#combine
#######################################################################################

multiqc output/STARlogs --filename output/STAR_alignment_report.html --ignore-samples Undetermined* --interactive