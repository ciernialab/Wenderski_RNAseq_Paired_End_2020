#!/bin/bash
#
#SBATCH -c 16
#SBATCH --mem-per-cpu=4000
#SBATCH --job-name=STARbuild
#SBATCH --output=STARbuildmm10.out
#SBATCH --time=10:00:00
#######################################################################################
#build STAR indexes for mm10 from ensemble
#https://leonjessen.wordpress.com/2014/12/01/how-do-i-create-star-indices-using-the-newest-grch38-version-of-the-human-genome/

#get fasta files mm10:
#tp://ftp.ensembl.org/pub/release-92/fasta/mus_musculus/dna/

mkdir -p STAR_libs/mm10/sequence
cd STAR_libs/mm10/sequence/
wget ftp://ftp.ensembl.org/pub/release-92/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.chromosome.{1..22}.fa.gz
wget ftp://ftp.ensembl.org/pub/release-92/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.chromosome.{MT,X,Y}.fa.gz
gunzip -c Mus_musculus.GRCm38.dna.chromosome.* > GRCm38_r92.all.fa
cd ../../../

#annotation: ftp://ftp.ensembl.org/pub/release-92/gtf/mus_musculus
mkdir -p STAR_libs/mm10/annotation
cd STAR_libs/mm10/annotation/
wget ftp://ftp.ensembl.org/pub/release-92/gtf/mus_musculus/Mus_musculus.GRCm38.92.gtf.gz
gunzip Mus_musculus.GRCm38.92.gtf
cd ../../../

#generate the STAR indices for 75bp reads
#hhttps://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf
mkdir -p STAR_libs/mm10/star_indices/

STAR --runThreadN 20 \
--runMode genomeGenerate \
--genomeDir STAR_libs/mm10/star_indices/ \
--genomeFastaFiles STAR_libs/mm10/sequence/GRCm38_r92.all.fa \
--sjdbGTFfile STAR_libs/mm10/annotation/Mus_musculus.GRCm38.92.gtf --sjdbOverhang 74
