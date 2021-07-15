#!/bin/bash
#
#SBATCH -c 12
#SBATCH --mem-per-cpu=4000
#SBATCH --job-name=UCSC
#SBATCH --output=HOMERUCSC.out
#SBATCH --time=12:00:00

#######################################################################################
#makeUCSCfile <tag directory-res 10 > -fragLength given -res 10 > UCSCbrowsertracks/
#10bp resolution
#######################################################################################
module load samtools/1.4
module load jre/1.8.0_121
module load R/3.6.1

mkdir -p UCSCbrowsertracks/
echo "making bedGraphs"


for sample in `cat SRR_Acc_List_2.txt`
do

#######################################################################################
# WT_2h
#######################################################################################

#make normalized bedgraphs:
#WT_2h
#TagDirectory/tag_SRR11313882
#TagDirectory/tag_SRR11313884
#TagDirectory/tag_SRR11313886
#TagDirectory/tag_SRR11313888
#TagDirectory/tag_SRR11313890

makeUCSCfile TagDirectory/tag_SRR11313882 -fragLength given -name WT_2h_Repl1 -res 10 > UCSCbrowsertracks/WT_2h_Repl1.bedGraph

makeUCSCfile TagDirectory/tag_SRR11313884 -fragLength given -name WT_2h_Repl2 -res 10 > UCSCbrowsertracks/WT_2h_Repl2.bedGraph

makeUCSCfile TagDirectory/tag_SRR11313886 -fragLength given -name WT_2h_Repl3 -res 10 > UCSCbrowsertracks/WT_2h_Repl3.bedGraph

makeUCSCfile TagDirectory/tag_SRR11313888 -fragLength given -name WT_2h_Repl4 -res 10 > UCSCbrowsertracks/WT_2h_Repl4.bedGraph

makeUCSCfile TagDirectory/tag_SRR11313890 -fragLength given -name WT_2h_Repl5 -res 10 > UCSCbrowsertracks/WT_2h_Repl5.bedGraph

#######################################################################################
# KO_2h
#######################################################################################
#make normalized bedgraphs:
#KO_2h
#TagDirectory/tag_SRR11313892
#TagDirectory/tag_SRR11313894
#TagDirectory/tag_SRR11313896
#TagDirectory/tag_SRR11313898
#TagDirectory/tag_SRR11313900
#TagDirectory/tag_SRR11313902
#TagDirectory/tag_SRR11313904

makeUCSCfile TagDirectory/tag_SRR11313892 -fragLength given -name KO_2h_Repl1 -res 10 > UCSCbrowsertracks/KO_2h_Repl1.bedGraph

makeUCSCfile TagDirectory/tag_SRR11313894 -fragLength given -name KO_2h_Repl2 -res 10 > UCSCbrowsertracks/KO_2h_Repl2.bedGraph

makeUCSCfile TagDirectory/tag_SRR11313896 -fragLength given -name KO_2h_Repl3 -res 10 > UCSCbrowsertracks/KO_2h_Repl3.bedGraph

makeUCSCfile TagDirectory/tag_SRR11313898 -fragLength given -name KO_2h_Repl4 -res 10 > UCSCbrowsertracks/KO_2h_Repl4.bedGraph

makeUCSCfile TagDirectory/tag_SRR11313900 -fragLength given -name KO_2h_Repl5 -res 10 > UCSCbrowsertracks/KO_2h_Repl5.bedGraph

makeUCSCfile TagDirectory/tag_SRR11313902 -fragLength given -name KO_2h_Repl6 -res 10 > UCSCbrowsertracks/KO_2h_Repl6.bedGraph

makeUCSCfile TagDirectory/tag_SRR11313904 -fragLength given -name KO_2h_Repl7 -res 10 > UCSCbrowsertracks/KO_2h_Repl7.bedGraph

#######################################################################################
# WT_1h_KCl
#######################################################################################
#make normalized bedgraphs:
#WT_1h_KCl
#TagDirectory/tag_SRR11313906
#TagDirectory/tag_SRR11313908
#TagDirectory/tag_SRR11313910
#TagDirectory/tag_SRR11313912
#TagDirectory/tag_SRR11313914

makeUCSCfile TagDirectory/tag_SRR11313906 -fragLength given -name WT_1h_KCl_Repl1 -res 10 > UCSCbrowsertracks/WT_1h_KCl_Repl1.bedGraph

makeUCSCfile TagDirectory/tag_SRR11313908 -fragLength given -name WT_1h_KCl_Repl2 -res 10 > UCSCbrowsertracks/WT_1h_KCl_Repl2.bedGraph

makeUCSCfile TagDirectory/tag_SRR11313910 -fragLength given -name WT_1h_KCl_Repl3 -res 10 > UCSCbrowsertracks/WT_1h_KCl_Repl3.bedGraph

makeUCSCfile TagDirectory/tag_SRR11313912 -fragLength given -name WT_1h_KCl_Repl4 -res 10 > UCSCbrowsertracks/WT_1h_KCl_Repl4.bedGraph

makeUCSCfile TagDirectory/tag_SRR11313914 -fragLength given -name WT_1h_KCl_Repl5 -res 10 > UCSCbrowsertracks/WT_1h_KCl_Repl5.bedGraph


#######################################################################################
# WT_7h
#######################################################################################
#make normalized bedgraphs:
#WT_7h
#TagDirectory/tag_SRR11313916
#TagDirectory/tag_SRR11313918
#TagDirectory/tag_SRR11313920
#TagDirectory/tag_SRR11313922

makeUCSCfile TagDirectory/tag_SRR11313916 -fragLength given -name WT_7h_Repl1 -res 10 > UCSCbrowsertracks/WT_7h_Repl1.bedGraph

makeUCSCfile TagDirectory/tag_SRR11313918 -fragLength given -name WT_7h_Repl2 -res 10 > UCSCbrowsertracks/WT_7h_Repl2.bedGraph

makeUCSCfile TagDirectory/tag_SRR11313920 -fragLength given -name WT_7h_Repl3 -res 10 > UCSCbrowsertracks/WT_7h_Repl3.bedGraph

makeUCSCfile TagDirectory/tag_SRR11313922 -fragLength given -name WT_7h_Repl4 -res 10 > UCSCbrowsertracks/WT_7h_Repl4.bedGraph

#######################################################################################
# WT_6h_KCl
#######################################################################################
#make normalized bedgraphs:
#WT_6h_KCl
#TagDirectory/tag_SRR11313924
#TagDirectory/tag_SRR11313926
#TagDirectory/tag_SRR11313928


makeUCSCfile TagDirectory/tag_SRR11313924 -fragLength given -name WT_6h_KCl_Repl1 -res 10 > UCSCbrowsertracks/WT_6h_KCl_Repl1.bedGraph

makeUCSCfile TagDirectory/tag_SRR11313926 -fragLength given -name WT_6h_KCl_Repl2 -res 10 > UCSCbrowsertracks/WT_6h_KCl_Repl2.bedGraph

makeUCSCfile TagDirectory/tag_SRR11313928 -fragLength given -name WT_6h_KCl_Repl3 -res 10 > UCSCbrowsertracks/WT_6h_KCl_Repl3.bedGraph

##########################################################################################
#######################################################################################
#make into ucsc format
#sed -i 's/old-text/new-text/g' input.txt
echo "converting to UCSC format"
#######################################################################################

#skip the first line, then add "chr" to each chromosome
sed -i "1n; s/^/chr/" UCSCbrowsertracks/WT_2h_Repl1.bedGraph
sed -i "1n; s/MT/M/g" UCSCbrowsertracks/WT_2h_Repl1.bedGraph

sed -i "1n; s/^/chr/" UCSCbrowsertracks/WT_2h_Repl2.bedGraph
sed -i "1n; s/MT/M/g" UCSCbrowsertracks/WT_2h_Repl2.bedGraph

sed -i "1n; s/^/chr/" UCSCbrowsertracks/WT_2h_Repl3.bedGraph
sed -i "1n; s/MT/M/g" UCSCbrowsertracks/WT_2h_Repl3.bedGraph

sed -i "1n; s/^/chr/" UCSCbrowsertracks/WT_2h_Repl4.bedGraph
sed -i "1n; s/MT/M/g" UCSCbrowsertracks/WT_2h_Repl4.bedGraph

sed -i "1n; s/^/chr/" UCSCbrowsertracks/WT_2h_Repl5.bedGraph
sed -i "1n; s/MT/M/g" UCSCbrowsertracks/WT_2h_Repl5.bedGraph

#skip the first line, then add "chr" to each chromosome

sed -i "1n; s/^/chr/" UCSCbrowsertracks/KO_2h_Repl1.bedGraph
sed -i "1n; s/MT/M/g" UCSCbrowsertracks/KO_2h_Repl1.bedGraph

sed -i "1n; s/^/chr/" UCSCbrowsertracks/KO_2h_Repl2.bedGraph
sed -i "1n; s/MT/M/g" UCSCbrowsertracks/KO_2h_Repl2.bedGraph

sed -i "1n; s/^/chr/" UCSCbrowsertracks/KO_2h_Repl3.bedGraph
sed -i "1n; s/MT/M/g" UCSCbrowsertracks/KO_2h_Repl3.bedGraph

sed -i "1n; s/^/chr/" UCSCbrowsertracks/KO_2h_Repl4.bedGraph
sed -i "1n; s/MT/M/g" UCSCbrowsertracks/KO_2h_Repl4.bedGraph

sed -i "1n; s/^/chr/" UCSCbrowsertracks/KO_2h_Repl5.bedGraph
sed -i "1n; s/MT/M/g" UCSCbrowsertracks/KO_2h_Repl5.bedGraph

sed -i "1n; s/^/chr/" UCSCbrowsertracks/KO_2h_Repl6.bedGraph
sed -i "1n; s/MT/M/g" UCSCbrowsertracks/KO_2h_Repl6.bedGraph

sed -i "1n; s/^/chr/" UCSCbrowsertracks/KO_2h_Repl7.bedGraph
sed -i "1n; s/MT/M/g" UCSCbrowsertracks/KO_2h_Repl7.bedGraph

#skip the first line, then add "chr" to each chromosome

sed -i "1n; s/^/chr/" UCSCbrowsertracks/WT_1h_KCl_Repl1.bedGraph
sed -i "1n; s/MT/M/g" UCSCbrowsertracks/WT_1h_KCl_Repl1.bedGraph

sed -i "1n; s/^/chr/" UCSCbrowsertracks/WT_1h_KCl_Repl2.bedGraph
sed -i "1n; s/MT/M/g" UCSCbrowsertracks/WT_1h_KCl_Repl2.bedGraph

sed -i "1n; s/^/chr/" UCSCbrowsertracks/WT_1h_KCl_Repl3.bedGraph
sed -i "1n; s/MT/M/g" UCSCbrowsertracks/WT_1h_KCl_Repl3.bedGraph

sed -i "1n; s/^/chr/" UCSCbrowsertracks/WT_1h_KCl_Repl4.bedGraph
sed -i "1n; s/MT/M/g" UCSCbrowsertracks/WT_1h_KCl_Repl4.bedGraph

sed -i "1n; s/^/chr/" UCSCbrowsertracks/WT_1h_KCl_Repl5.bedGraph
sed -i "1n; s/MT/M/g" UCSCbrowsertracks/WT_1h_KCl_Repl5.bedGraph

#skip the first line, then add "chr" to each chromosome

sed -i "1n; s/^/chr/" UCSCbrowsertracks/WT_7h_Repl1.bedGraph
sed -i "1n; s/MT/M/g" UCSCbrowsertracks/WT_7h_Repl1.bedGraph

sed -i "1n; s/^/chr/" UCSCbrowsertracks/WT_7h_Repl2.bedGraph
sed -i "1n; s/MT/M/g" UCSCbrowsertracks/WT_7h_Repl2.bedGraph

sed -i "1n; s/^/chr/" UCSCbrowsertracks/WT_7h_Repl3.bedGraph
sed -i "1n; s/MT/M/g" UCSCbrowsertracks/WT_7h_Repl3.bedGraph

sed -i "1n; s/^/chr/" UCSCbrowsertracks/WT_7h_Repl4.bedGraph
sed -i "1n; s/MT/M/g" UCSCbrowsertracks/WT_7h_Repl4.bedGraph

#skip the first line, then add "chr" to each chromosome

sed -i "1n; s/^/chr/" UCSCbrowsertracks/WT_6h_KCl_Repl1.bedGraph
sed -i "1n; s/MT/M/g" UCSCbrowsertracks/WT_6h_KCl_Repl1.bedGraph

sed -i "1n; s/^/chr/" UCSCbrowsertracks/WT_6h_KCl_Repl2.bedGraph
sed -i "1n; s/MT/M/g" UCSCbrowsertracks/WT_6h_KCl_Repl2.bedGraph

sed -i "1n; s/^/chr/" UCSCbrowsertracks/WT_6h_KCl_Repl3.bedGraph
sed -i "1n; s/MT/M/g" UCSCbrowsertracks/WT_6h_KCl_Repl3.bedGraph

#######################################################################################
#######################################################################################
#convert to bigwigs
echo "making bigwigs"
#######################################################################################
#convert to bigwig
#The input bedGraph file must be sort, use the unix sort command:
# WT_2h Repl1
sort -k1,1 -k2,2n UCSCbrowsertracks/WT_2h_Repl1.bedGraph > UCSCbrowsertracks/WT_2h_Repl1.sort.bedGraph

#remove track line (now at the end of sort files)
sed -i '$d' UCSCbrowsertracks/WT_2h_Repl1.sort.bedGraph

#fix extensions beyond chromosomes > removes entry
#bedClip [options] input.bed chrom.sizes output.bed
bedClip UCSCbrowsertracks/WT_2h_Repl1.sort.bedGraph $mm10chrsizes UCSCbrowsertracks/WT_2h_Repl1.sort2.bedGraph

#bedGraphToBigWig in.bedGraph chrom.sizes out.bw
bedGraphToBigWig UCSCbrowsertracks/WT_2h_Repl1.sort2.bedGraph $mm10chrsizes UCSCbrowsertracks/WT_2h_Repl1.bw

# WT_2h Repl2

sort -k1,1 -k2,2n UCSCbrowsertracks/WT_2h_Repl2.bedGraph > UCSCbrowsertracks/WT_2h_Repl2.sort.bedGraph

sed -i '$d' UCSCbrowsertracks/WT_2h_Repl2.sort.bedGraph

bedClip UCSCbrowsertracks/WT_2h_Repl2.sort.bedGraph $mm10chrsizes UCSCbrowsertracks/WT_2h_Repl2.sort2.bedGraph

bedGraphToBigWig UCSCbrowsertracks/WT_2h_Repl2.sort2.bedGraph $mm10chrsizes UCSCbrowsertracks/WT_2h_Repl2.bw

# WT_2h Repl3

sort -k1,1 -k2,2n UCSCbrowsertracks/WT_2h_Repl3.bedGraph > UCSCbrowsertracks/WT_2h_Repl3.sort.bedGraph

sed -i '$d' UCSCbrowsertracks/WT_2h_Repl3.sort.bedGraph

bedClip UCSCbrowsertracks/WT_2h_Repl3.sort.bedGraph $mm10chrsizes UCSCbrowsertracks/WT_2h_Repl3.sort2.bedGraph

bedGraphToBigWig UCSCbrowsertracks/WT_2h_Repl3.sort2.bedGraph $mm10chrsizes UCSCbrowsertracks/WT_2h_Repl3.bw

# WT_2h Repl4

sort -k1,1 -k2,2n UCSCbrowsertracks/WT_2h_Repl4.bedGraph > UCSCbrowsertracks/WT_2h_Repl4.sort.bedGraph

sed -i '$d' UCSCbrowsertracks/WT_2h_Repl4.sort.bedGraph

bedClip UCSCbrowsertracks/WT_2h_Repl4.sort.bedGraph $mm10chrsizes UCSCbrowsertracks/WT_2h_Repl4.sort2.bedGraph

bedGraphToBigWig UCSCbrowsertracks/WT_2h_Repl4.sort2.bedGraph $mm10chrsizes UCSCbrowsertracks/WT_2h_Repl4.bw

# WT_2h Repl5

sort -k1,1 -k2,2n UCSCbrowsertracks/WT_2h_Repl5.bedGraph > UCSCbrowsertracks/WT_2h_Repl5.sort.bedGraph

sed -i '$d' UCSCbrowsertracks/WT_2h_Repl5.sort.bedGraph

bedClip UCSCbrowsertracks/WT_2h_Repl5.sort.bedGraph $mm10chrsizes UCSCbrowsertracks/WT_2h_Repl5.sort2.bedGraph

bedGraphToBigWig UCSCbrowsertracks/WT_2h_Repl5.sort2.bedGraph $mm10chrsizes UCSCbrowsertracks/WT_2h_Repl5.bw


# KO_2h Repl1
sort -k1,1 -k2,2n UCSCbrowsertracks/KO_2h_Repl1.bedGraph > UCSCbrowsertracks/KO_2h_Repl1.sort.bedGraph

sed -i '$d' UCSCbrowsertracks/KO_2h_Repl1.sort.bedGraph

bedClip UCSCbrowsertracks/KO_2h_Repl1.sort.bedGraph $mm10chrsizes UCSCbrowsertracks/KO_2h_Repl1.sort2.bedGraph

bedGraphToBigWig UCSCbrowsertracks/KO_2h_Repl1.sort2.bedGraph $mm10chrsizes UCSCbrowsertracks/KO_2h_Repl1.bw

# KO_2h Repl2
sort -k1,1 -k2,2n UCSCbrowsertracks/KO_2h_Repl2.bedGraph > UCSCbrowsertracks/KO_2h_Repl2.sort.bedGraph

sed -i '$d' UCSCbrowsertracks/KO_2h_Repl2.sort.bedGraph

bedClip UCSCbrowsertracks/KO_2h_Repl2.sort.bedGraph $mm10chrsizes UCSCbrowsertracks/KO_2h_Repl2.sort2.bedGraph

bedGraphToBigWig UCSCbrowsertracks/KO_2h_Repl2.sort2.bedGraph $mm10chrsizes UCSCbrowsertracks/KO_2h_Repl2.bw

# KO_2h Repl3

sort -k1,1 -k2,2n UCSCbrowsertracks/KO_2h_Repl3.bedGraph > UCSCbrowsertracks/KO_2h_Repl3.sort.bedGraph

sed -i '$d' UCSCbrowsertracks/KO_2h_Repl3.sort.bedGraph

bedClip UCSCbrowsertracks/KO_2h_Repl3.sort.bedGraph $mm10chrsizes UCSCbrowsertracks/KO_2h_Repl3.sort2.bedGraph

bedGraphToBigWig UCSCbrowsertracks/KO_2h_Repl3.sort2.bedGraph $mm10chrsizes UCSCbrowsertracks/KO_2h_Repl3.bw

# KO_2h Repl4

sort -k1,1 -k2,2n UCSCbrowsertracks/KO_2h_Repl4.bedGraph > UCSCbrowsertracks/KO_2h_Repl4.sort.bedGraph

sed -i '$d' UCSCbrowsertracks/KO_2h_Repl4.sort.bedGraph

bedClip UCSCbrowsertracks/KO_2h_Repl4.sort.bedGraph $mm10chrsizes UCSCbrowsertracks/KO_2h_Repl4.sort2.bedGraph

bedGraphToBigWig UCSCbrowsertracks/KO_2h_Repl4.sort2.bedGraph $mm10chrsizes UCSCbrowsertracks/KO_2h_Repl4.bw

# WT_2h Repl5

sort -k1,1 -k2,2n UCSCbrowsertracks/KO_2h_Repl5.bedGraph > UCSCbrowsertracks/KO_2h_Repl5.sort.bedGraph

sed -i '$d' UCSCbrowsertracks/KO_2h_Repl5.sort.bedGraph

bedClip UCSCbrowsertracks/KO_2h_Repl5.sort.bedGraph $mm10chrsizes UCSCbrowsertracks/KO_2h_Repl5.sort2.bedGraph

bedGraphToBigWig UCSCbrowsertracks/KO_2h_Repl5.sort2.bedGraph $mm10chrsizes UCSCbrowsertracks/KO_2h_Repl5.bw

# KO_2h Repl6

sort -k1,1 -k2,2n UCSCbrowsertracks/KO_2h_Repl6.bedGraph > UCSCbrowsertracks/KO_2h_Repl6.sort.bedGraph

sed -i '$d' UCSCbrowsertracks/KO_2h_Repl6.sort.bedGraph

bedClip UCSCbrowsertracks/KO_2h_Repl6.sort.bedGraph $mm10chrsizes UCSCbrowsertracks/KO_2h_Repl6.sort2.bedGraph

bedGraphToBigWig UCSCbrowsertracks/KO_2h_Repl6.sort2.bedGraph $mm10chrsizes UCSCbrowsertracks/KO_2h_Repl6.bw

# WT_2h Repl7

sort -k1,1 -k2,2n UCSCbrowsertracks/KO_2h_Repl7.bedGraph > UCSCbrowsertracks/KO_2h_Repl7.sort.bedGraph

sed -i '$d' UCSCbrowsertracks/KO_2h_Repl7.sort.bedGraph

bedClip UCSCbrowsertracks/KO_2h_Repl7.sort.bedGraph $mm10chrsizes UCSCbrowsertracks/KO_2h_Repl7.sort2.bedGraph

bedGraphToBigWig UCSCbrowsertracks/KO_2h_Repl7.sort2.bedGraph $mm10chrsizes UCSCbrowsertracks/KO_2h_Repl7.bw


# WT_1h_KCl Repl1
sort -k1,1 -k2,2n UCSCbrowsertracks/WT_1h_KCl_Repl1.bedGraph > UCSCbrowsertracks/WT_1h_KCl_Repl1.sort.bedGraph

sed -i '$d' UCSCbrowsertracks/WT_1h_KCl_Repl1.sort.bedGraph

bedClip UCSCbrowsertracks/WT_1h_KCl_Repl1.sort.bedGraph $mm10chrsizes UCSCbrowsertracks/WT_1h_KCl_Repl1.sort2.bedGraph

bedGraphToBigWig UCSCbrowsertracks/WT_1h_KCl_Repl1.sort2.bedGraph $mm10chrsizes UCSCbrowsertracks/WT_1h_KCl_Repl1.bw

# WT_1h_KCl Repl2
sort -k1,1 -k2,2n UCSCbrowsertracks/WT_1h_KCl_Repl2.bedGraph > UCSCbrowsertracks/WT_1h_KCl_Repl2.sort.bedGraph

sed -i '$d' UCSCbrowsertracks/WT_1h_KCl_Repl2.sort.bedGraph

bedClip UCSCbrowsertracks/WT_1h_KCl_Repl2.sort.bedGraph $mm10chrsizes UCSCbrowsertracks/WT_1h_KCl_Repl2.sort2.bedGraph

bedGraphToBigWig UCSCbrowsertracks/WT_1h_KCl_Repl2.sort2.bedGraph $mm10chrsizes UCSCbrowsertracks/WT_1h_KCl_Repl2.bw

# WT_1h_KCl Repl3

sort -k1,1 -k2,2n UCSCbrowsertracks/WT_1h_KCl_Repl3.bedGraph > UCSCbrowsertracks/WT_1h_KCl_Repl3.sort.bedGraph

sed -i '$d' UCSCbrowsertracks/WT_6h_KCl_Repl3.sort.bedGraph

bedClip UCSCbrowsertracks/WT_1h_KCl_Repl3.sort.bedGraph $mm10chrsizes UCSCbrowsertracks/WT_1h_KCl_Repl3.sort2.bedGraph

bedGraphToBigWig UCSCbrowsertracks/WT_1h_KCl_Repl3.sort2.bedGraph $mm10chrsizes UCSCbrowsertracks/WT_1h_KCl_Repl3.bw

# WT_1h_KCl Repl4

sort -k1,1 -k2,2n UCSCbrowsertracks/WT_1h_KCl_Repl4.bedGraph > UCSCbrowsertracks/WT_1h_KCl_Repl4.sort.bedGraph

sed -i '$d' UCSCbrowsertracks/WT_1h_KCl_Repl4.sort.bedGraph

bedClip UCSCbrowsertracks/WT_1h_KCl_Repl4.sort.bedGraph $mm10chrsizes UCSCbrowsertracks/WT_1h_KCl_Repl4.sort2.bedGraph

bedGraphToBigWig UCSCbrowsertracks/WT_1h_KCl_Repl4.sort2.bedGraph $mm10chrsizes UCSCbrowsertracks/WT_1h_KCl_Repl4.bw

# WT_1h_KCl Repl5

sort -k1,1 -k2,2n UCSCbrowsertracks/WT_1h_KCl_Repl5.bedGraph > UCSCbrowsertracks/WT_1h_KCl_Repl5.sort.bedGraph

sed -i '$d' UCSCbrowsertracks/WT_1h_KCl_Repl5.sort.bedGraph

bedClip UCSCbrowsertracks/WT_1h_KCl_Repl5.sort.bedGraph $mm10chrsizes UCSCbrowsertracks/WT_1h_KCl_Repl5.sort2.bedGraph

bedGraphToBigWig UCSCbrowsertracks/WT_1h_KCl_Repl5.sort2.bedGraph $mm10chrsizes UCSCbrowsertracks/WT_1h_KCl_Repl5.bw


# WT_7h Repl1
sort -k1,1 -k2,2n UCSCbrowsertracks/WT_7h_Repl1.bedGraph > UCSCbrowsertracks/WT_7h_Repl1.sort.bedGraph

sed -i '$d' UCSCbrowsertracks/WT_7h_Repl1.sort.bedGraph

bedClip UCSCbrowsertracks/WT_7h_Repl1.sort.bedGraph $mm10chrsizes UCSCbrowsertracks/WT_7h_Repl1.sort2.bedGraph

bedGraphToBigWig UCSCbrowsertracks/WT_7h_Repl1.sort2.bedGraph $mm10chrsizes UCSCbrowsertracks/WT_7h_Repl1.bw

# WT_7h Repl2

sort -k1,1 -k2,2n UCSCbrowsertracks/WT_7h_Repl2.bedGraph > UCSCbrowsertracks/WT_7h_Repl2.sort.bedGraph

sed -i '$d' UCSCbrowsertracks/WT_7h_Repl2.sort.bedGraph

bedClip UCSCbrowsertracks/WT_7h_Repl2.sort.bedGraph $mm10chrsizes UCSCbrowsertracks/WT_7h_Repl2.sort2.bedGraph

bedGraphToBigWig UCSCbrowsertracks/WT_7h_Repl2.sort2.bedGraph $mm10chrsizes UCSCbrowsertracks/WT_7h_Repl2.bw

# WT_7h Repl3

sort -k1,1 -k2,2n UCSCbrowsertracks/WT_7h_Repl3.bedGraph > UCSCbrowsertracks/WT_7h_Repl3.sort.bedGraph

sed -i '$d' UCSCbrowsertracks/WT_7h_Repl3.sort.bedGraph

bedClip UCSCbrowsertracks/WT_7h_Repl3.sort.bedGraph $mm10chrsizes UCSCbrowsertracks/WT_7h_Repl3.sort2.bedGraph

bedGraphToBigWig UCSCbrowsertracks/WT_7h_Repl3.sort2.bedGraph $mm10chrsizes UCSCbrowsertracks/WT_7h_Repl3.bw

# WT_7h Repl4

sort -k1,1 -k2,2n UCSCbrowsertracks/WT_7h_Repl4.bedGraph > UCSCbrowsertracks/WT_7h_Repl4.sort.bedGraph

sed -i '$d' UCSCbrowsertracks/WT_7h_Repl4.sort.bedGraph

bedClip UCSCbrowsertracks/WT_7h_Repl4.sort.bedGraph $mm10chrsizes UCSCbrowsertracks/WT_7h_Repl4.sort2.bedGraph

bedGraphToBigWig UCSCbrowsertracks/WT_7h_Repl4.sort2.bedGraph $mm10chrsizes UCSCbrowsertracks/WT_7h_Repl4.bw



# WT_6h_KCl Repl1
sort -k1,1 -k2,2n UCSCbrowsertracks/WT_6h_KCl_Repl1.bedGraph > UCSCbrowsertracks/WT_6h_KCl_Repl1.sort.bedGraph

sed -i '$d' UCSCbrowsertracks/WT_6h_KCl_Repl1.sort.bedGraph

bedClip UCSCbrowsertracks/WT_6h_KCl_Repl1.sort.bedGraph $mm10chrsizes UCSCbrowsertracks/WT_6h_KCl_Repl1.sort2.bedGraph

bedGraphToBigWig UCSCbrowsertracks/WT_6h_KCl_Repl1.sort2.bedGraph $mm10chrsizes UCSCbrowsertracks/WT_6h_KCl_Repl1.bw

# WT_6h_KCl Repl2
sort -k1,1 -k2,2n UCSCbrowsertracks/WT_6h_KCl_Repl2.bedGraph > UCSCbrowsertracks/WT_6h_KCl_Repl2.sort.bedGraph

sed -i '$d' UCSCbrowsertracks/WT_6h_KCl_Repl2.sort.bedGraph

bedClip UCSCbrowsertracks/WT_6h_KCl_Repl2.sort.bedGraph $mm10chrsizes UCSCbrowsertracks/WT_6h_KCl_Repl2.sort2.bedGraph

bedGraphToBigWig UCSCbrowsertracks/WT_6h_KCl_Repl2.sort2.bedGraph $mm10chrsizes UCSCbrowsertracks/WT_6h_KCl_Repl2.bw

# WT_6h_KCl Repl3

sort -k1,1 -k2,2n UCSCbrowsertracks/WT_6h_KCl_Repl3.bedGraph > UCSCbrowsertracks/WT_6h_KCl_Repl3.sort.bedGraph

sed -i '$d' UCSCbrowsertracks/WT_6h_KCl_Repl3.sort.bedGraph

bedClip UCSCbrowsertracks/WT_6h_KCl_Repl3.sort.bedGraph $mm10chrsizes UCSCbrowsertracks/WT_6h_KCl_Repl3.sort2.bedGraph

bedGraphToBigWig UCSCbrowsertracks/WT_6h_KCl_Repl3.sort2.bedGraph $mm10chrsizes UCSCbrowsertracks/WT_6h_KCl_Repl3.bw

#zip bedgraphs
gzip UCSCbrowsertracks/*.bedGraph

echo "complete"

done