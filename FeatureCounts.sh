#!/bin/bash
#
#SBATCH --workdir /alder/home/xlum/Wenderski_RNAseq
#SBATCH -c 10                              # number of processors
#SBATCH -N 1                                # number of nodes
#SBATCH --output=featureCounts.out
#SBATCH --time=0-03:00:00

#######################################################################################
#RNA QuantSeq Analysis: https://www.lexogen.com/quantseq-data-analysis/

#######################################################################################
#load modules

#######################################################################################
mkdir -p featureCounts
mkdir -p output/countlogs

#it contains this:
for sample in `cat SRR_Acc_List.txt`
do


echo ${sample} "count started"

#forward strand, exon counts by gene_id, 
featureCounts -T 10 -s 1 -t exon -g gene_id -p \
-a /alder/data/cbh/ciernia-data/genomes/STAR_libs/mm10/annotation/Mus_musculus.GRCm38.92.gtf \
-o featureCounts/${sample}_featureCounts.txt \
star_align/${sample}_dedup.bam \
2> ~/Wenderski_RNAseq/output/countlogs/${sample}_featureCounts.log

#--quantMode GeneCounts 
echo ${sample} "count completed"
done



#######################################################################################
#collect all qc together with:
#module load multiqc/1.6


#######################################################################################
#######################################################################################