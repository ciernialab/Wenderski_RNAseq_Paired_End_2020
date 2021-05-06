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
We need to remove adapters and poor quality reads before aligning.

Trimmomatic will look for seed matches of 16 bases with 2 mismatches allowed and will then extend and clip if a score of 30 for PE or 10 for SE is reached (~17 base match).

Minimum adapter length is 8 bases.

T = keeps both reads even if only one passes critieria.

Trims low quality bases at leading and trailing if quality score < 15.

Clipped the first 3 basepairs: HEADCROP 3 due to low quality and high G based on 1st round of trimming.

Sliding window: scans in a 4 base window, cuts when the average quality drops below 15.

Log outputs number of input reads, trimmed, and surviving reads in the trim_log_samplename.

It uses TruSeq3-PE.fa (comes wiht Trimmomatic download).

The path is set in bash_profile with $ADAPTERS
To check the content of the file:
```
less $ADAPTERS/TruSeq3-PE.fa
```

Run the script.
```
sbatch Trim.sh
```

## Repeat QC on post trim file
This step is the same as pretrim.
```
sbatch postrim_fastqc.sh
```

# Step 4: QC of Fastq Files: Contamination Screening
FastQ Screen is an application allowing us to search a large sequence dataset against a panel of different genomes to determine from where the data originate.

The program generates both text and graphical output to inform you what proportion of the library was able to map, either uniquely or to more than one location, against each of the specified reference genomes. 

Run script to check trimmed fastq files.
```
sbatch Fastqscreen.sh
```

The output is found in output/FastqScreen_multiqc_report.html

## Step 5: Align to mm10 genome using STAR
before aligning to STAR, we need to first unzip the fastq files in the trimmed directory.

Run the following script to unzip all paired trimmed fastq files
```
sbatch trimmed_unzip.sh
```

Run the following script to align trimmed fastq files to the mm10 genome using STAR. We specified the location of index files in bash_profile. 
```
sbatch STAR_align.sh
```

## Step 6: Filter aligned files
We will now convert sam files to bam and filter to remove PCR duplicates, remove unmapped reads and econdary alignments (multi-mappers), and remove unpaired reads.

Samtools is used to convert sam to bam.
Samtools fixmate is used to removed unmapped reads and 2ndary alignments.
Remove PCR duplicates using -F 0x400 and -f 0x2 to keep only propperly paired reads.
Index reads and collect additional QC metrics using picard tools and samtools flagstat.
QC metrics are then collected into a single report using multiqc.

Run the following script.
```
sbatch SamtoolsFiltering.sh
```

## Step 7: Counting Reads using featureCounts
Now that we have aligned reads to the mm10 genome, the next step is to count how many reads have been mapped to each gene.

The input files required are BAM files and an associated annotation file in GTF format. featureCounts (alternative htseq-count can be used instead) takes the alignment coordinates for each read and cross-references that to the coordinates for features described in the GTF file. featureCounts is best used for counting reads associated with gene but not splice isoforms and transcripts.

To set up featureCounts, download the latest subread package version from SourceForce website (https://sourceforge.net/projects/subread/files/subread-2.0.2/). Currently, the latest version is subread-2.0.2-source.tar.gz. Use FileZilla to transfer the download file to pipline-tools directory. Other useful links: http://bioinf.wehi.edu.au/subread-package/SubreadUsersGuide.pdf.

To unzip the subread package and delete the zip file,
```
tar zxvf subread-2.0.2-source.tar.gz && rm subread-2.0.2-source.tar.gz
```

Enter the src directory and then build it on Linux/Unix system,
```
make -f Makefile.Linux
```
Create path for the featureCounts and add to .bash_profile and source the changes made to .bash_profile before they take place in the current terminal session. The bash profile will load each time you start a new session.
```
#featureCounts
PATH=$PATH:/alder/home/xlum/pipline-tools/subread-2.0.2-source.//bin/
```

To use featureCounts, run the following script.
```
sbatch featureCounts.sh
```




