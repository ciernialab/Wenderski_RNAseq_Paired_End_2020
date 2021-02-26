# Wenderski RNAseq Analysis

Processing RNA-seq PE files sourced from the Crabtree Lab: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7211998/.

The RNAseq dataset is downloaded from published PE RNA-seq files from GEO and processing them for QC, aligning and analysis in Rstudio. The analysis is setup to run on the DMCBH Alder computing cluster using bash scripts.

## Step 1: Fetching Data from GEO
We will be analyzing RNAseq dataset from the following paper:https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7211998/.

The data is stored in NCBI GEO as SRR files which are highly compressed files that can be converted to fastq files using SRA Toolkit: https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP252982&o=acc_s%3Aa.


To download a text file of SRR files, use the Run Select function.
Navigate to the bottom of the page and select all the RNA-seq data. Careful to not select any ATAC-seq data. Download "Accession List".

### Make a directory for the experiment in your home directory
```
cd ~/ && mkdir KCl_RNAseq
```


Make a copy of the SRR_Acc.List.txt file in your home directory.
```
cd ~/KCl_RNAseq
nano SRR_Acc_List.txt


### Run the SRRpull.sh script

The script is set 


