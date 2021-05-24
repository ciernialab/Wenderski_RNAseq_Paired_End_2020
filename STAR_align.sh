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

for sample in `cat SRR_Acc_List_2.txt`
do

echo ${sample} "starting STAR align"

# allow 10000 files to be open at once
ulimit -n 10000

# minimum intron size is 21 by default; genomic gap is considered splice junction if gap length is greater than 21bp.
# otherwise, it is considered a deletion.

STAR --runThreadN 12 \
--genomeDir $STAR_LIBS/mm10/star_indices/ \
--readFilesIn trimmed_unzip/${sample}_1.paired.fastq trimmed_unzip/${sample}_2.paired.fastq \
--outFilterType BySJout \ # reduce the number of spurious junctions
--outFilterMismatchNmax 999 \
--outFilterMismatchNoverLmax 0.04 \ # maximum number of mismatches per pair relative to read length (0.04 * 2 * read length)
--alignIntronMin 20 \
--alignIntronMax 1000000 \ # maximum intron length
--alignMatesGapMax 1000000 \ # maximum genomic distance between mates (read1 and read2) should be the same as alignIntronMax
--outSAMattributes NH HI nM AS MD \
--outSAMtype BAM SortedByCoordinate \
--quantMode GeneCounts \ # produces star matrix in out.tab file format
--outFileNamePrefix star_align/${sample}.sorted.bam \
2> output/starlogs/${sample}_starlog.log

echo ${sample} "finished STARalign"

done

