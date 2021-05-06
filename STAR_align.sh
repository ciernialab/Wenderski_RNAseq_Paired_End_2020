#!/bin/bash
#
#SBATCH -c 4
#SBATCH --mem-per-cpu=2000
#SBATCH --job-name=STARalignment
#SBATCH --output=STAR_align.out
#SBATCH --time=6:00:00

#######################################################################################
#STAR alignment: https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf

#######################################################################################
# -p no error if existing, make parent directories as needed
mkdir -p star_align
mkdir -p output/starlogs

for sample in `cat SRR_Acc_List.txt`
do

echo ${sample} "starting STAR align"

# allow 10000 files to be open at once
ulimit -n 10000

STAR --runThreadN 12 \
--genomeDir $STAR_LIBS/mm10/star_indices/ \
--readFilesIn trimmed_unzip/${sample}_1.paired.fastq trimmed_unzip/${sample}_2.paired.fastq \
--outFilterType BySJout \
--outFilterMultimapNmax 20 \
--alignSJoverhangMin 8 \
--alignSJDBoverhangMin 1 \
--outFilterMismatchNmax 999 \
--outFilterMismatchNoverLmax 0.1 \
--alignIntronMin 20 \
--alignIntronMax 1000000 \
--alignMatesGapMax 1000000 \
--outSAMattributes NH HI nM AS MD \
--outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix star_align/${sample}.sorted.bam \
2> output/starlogs/${sample}_starlog.log

echo ${sample} "finished STARalign"

done

