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
```

### Run the SRRpull.sh script

The script is setup to run from the experiment directory inside home directory.

It makes a SRA directory in the experiment directory and uses a loop to pull each SRR file in SRR_Acc_List.txt using prefetch.

We added --split to fastq-dump command to get R1 and R2 fastq for each SRR entry because the files are PE.

Make a SRR_pull.sh script in the experiment directory and run the script as sbatch submission to Alder.
```
sbatch SRR.pull.sh
```

Check to see if the script is running.
```
squeue
```

Can also check to see if SRAfetch is running properly by checking the contents of the experiment directory to see if a SRAfetch.out has been generated.
```
head SRAfetch.out
```

This can also be used to monitor the progress of SRAfetch.

Once the script has finished running, make sure to check all the SRA files have been successfully copied over.

cd ~/KCl_RNAseq/SRA/SRA_checksum/SRAcheck.log

Make sure ALL files have 'OK' and "Database 'SRRNAME.sra' is consistent" listed. Need to rerun SRRpull.sh script if encountered any errors.

## Step 2: QC of Fastq Files
We need to check the quality of the fastq files before and after trimming. We are using FastQC from https://www.bioinformatics.babraham.ac.uk/projects/fastqc/. Refer to their tutorial for output file interpretations.
```
sbatch pretrim_fastqc.sh
```

Check PretrimFastQC_multiqc_report.html for details of sequence quality.

## Step 3: Trimming fastq files 


