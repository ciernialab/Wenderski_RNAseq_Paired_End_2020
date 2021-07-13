# Wenderski RNAseq Analysis

Processing RNA-seq PE files sourced from the Crabtree Lab: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7211998/.

The RNAseq dataset is downloaded from published PE RNA-seq files from GEO Accession (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE147056) and processing them for QC, aligning and analysis in Rstudio. The analysis is setup to run on the DMCBH Alder computing cluster using bash scripts.

## Step 1: Fetching Data from GEO
We will be analyzing RNAseq dataset from the following paper:https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7211998/.

The data is stored in NCBI GEO as SRR files which are highly compressed files that can be converted to fastq files using SRA Toolkit: https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP252982&o=acc_s%3Aa.

To download a text file of SRR files, use the Run Select function.
Navigate to the bottom of the page and select all the RNA-seq data. Careful to not select any ATAC-seq data. Download "Accession List".

### Make a directory for the experiment in your home directory
```
cd ~/ && mkdir Wenderski_RNAseq
```

Make a copy of the SRR_Acc.List.txt file in your home directory.
```
cd ~/Wenderski_RNAseq
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

cd ~/Wenderski_RNAseq/SRA/SRA_checksum/SRAcheck.log

Make sure ALL files have 'OK' and "Database 'SRRNAME.sra' is consistent" listed. Need to rerun SRRpull.sh script if encountered any errors.

## Step 2: Concatenate fastq files
This paper sequenced the RNA library for a given sample on 2 separate sequencing runs. The paired-end sequencing data for each sample are stored in consecutive SRA files (i.e. SRR11313882 & SRR11313883 are from the same RNA library and so is SRR11313884 and SRR11313885...etc). These files needs to be contacaneted into a single file for both read 1 and read 2 of each sample. 

To do so, run the following script. The new fastq files are stored in a new directory - combined_fastqc under the first occurance of SRA file name for a given sample. For instance, SRR11313882 & SRR11313883 are combined together. The new fastq files will be stored under SRR11313882.
```
sbatch combine_fastq.sh
```

Alternatively can also run this script (which is much faster ~2 min) inside the directory containing the individual fastq files (i.e. SRA).
```
sbatch combine_fastq2.sh
```

## Step 3: QC of Fastq Files
We need to check the quality of the fastq files before and after trimming. We are using FastQC from https://www.bioinformatics.babraham.ac.uk/projects/fastqc/. Refer to their tutorial for output file interpretations.

We have combined all the fastq files from the same sample into one fastq file, so our SRR_Acc_List.txt also needs to be updated for downstream scripts that references the SRR_Acc_List.txt file. 

We will first ensure SRR_Acc_List.txt is sorted and make a copy of it.

```
sort -o SRR_Acc_List.txt SRR_Acc_List.txt
cp SRR_Acc_List.txt SRR_Acc_List_2.txt
```

Then we'll use the ex (text editor) command to match every line inside the file and deleting the next line.

```
ex SRR_Acc_List_2.txt <<\EX
:g/$/+d
:wq!
EX
```

Make sure to double check the txt file to ensure it contains every other SRA file name.

```
cat SRR_Acc_List_2.txt
```

Run the following script to perform quality check on the fastq files prior to trimming.
```
sbatch pretrim_fastqc.sh
```

Check PretrimFastQC_multiqc_report.html for details of sequence quality.

## Step 4: Trimming fastq files 
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
sbatch trim.sh
```

## Repeat QC on post trim file
This step is the same as pretrim.
```
sbatch postrim_fastqc.sh
```

# Step 5: QC of Fastq Files: Contamination Screening
FastQ Screen is an application allowing us to search a large sequence dataset against a panel of different genomes to determine from where the data originate.

The program generates both text and graphical output to inform you what proportion of the library was able to map, either uniquely or to more than one location, against each of the specified reference genomes. 

Run script to check trimmed fastq files.
```
sbatch fastqscreen.sh
```

The output is found in output/FastqScreen_multiqc_report.html

## Step 6: Align to mm10 genome using STAR

### Generating star index 
We first need to generate star indices for efficient mapping of RNA-seq fastq data to the indexed reference genome.

In desired directory, run the following script for setting up genome sequence, annotation, and star indices. The star index is optimized for 2 x 75bp sequencing data. Edit --sjdbOverhang for alternative read lengths.  Ideally, this length should be equal to the ReadLength-1, where ReadLength is the length of the reads. For more information: https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf.
```
sbatch mm10STARbuild.sh
```

### STAR Alignment
before aligning to STAR, we need to first unzip the fastq files in the trimmed directory.

Run the following script to unzip all paired trimmed fastq files
```
sbatch trimmed_unzip.sh
```

Run the following script to align trimmed fastq files to the mm10 genome using STAR. We specified the location of index files in bash_profile. 
```
sbatch STAR_align.sh
```

output star count matrix in the format: SRR11313882.ReadsPerGene.out.tab can be exported using FileZilla for expression analysis using Rstudio utilizing either edgeR or DESeq2.

## Step 7: Filter aligned files (Optional)
We will now convert sam files to bam and filter to remove PCR duplicates, remove unmapped reads and econdary alignments (multi-mappers), and remove unpaired reads. This step is optional because star count matrix provides the same gene count output as featureCounts which requires prior data filtering using SAMFilter.

Samtools is used to convert sam to bam.
Samtools fixmate is used to removed unmapped reads and 2ndary alignments.
-f 0x2 to keep only propperly paired reads.
PCR duplicates if NOT removed due to possibility of accidental removal of biological duplicates.
Index reads and collect additional QC metrics using picard tools and samtools flagstat.
QC metrics are then collected into a single report using multiqc.

Run the following script.
```
sbatch SamtoolsFiltering.sh
```

## Step 8: Counting Reads using featureCounts (Optional)
Now that we have aligned reads to the mm10 genome, the next step is to count how many reads have been mapped to each gene. This step is optional because star count matrix provides the same gene count output as featureCounts. 

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


## Step 9: Making TrackHubs using HOMOR

### Make Tag Directories
The first step to running HOMER is to make Tag Directories: http://homer.ucsd.edu/homer/ngs/tagDir.html We will make a folder Tag_Directories and then make tags for each sample WT and KO sample individually and all the input samples together. This is approach is specifically based for this experiment in which we have only 2-3 biological replicates per condition and input samples are not matched to individual samples. 

For WT 2h treatment with TTX/APV:
SRR11313882
SRR11313884
SRR11313886
SRR11313888
SRR11313890

For KO 2h treatment with TTX/APV:
SRR11313892
SRR11313894
SRR11313896
SRR11313898
SRR11313900
SRR11313902
SRR11313904

For WT 2h treatment with TTX/APV followed by 1h KCl stimulation:
SRR11313906
SRR11313908
SRR11313910
SRR11313912
SRR11313914

For WT 7h treatment with TTX/APV:
SRR11313916
SRR11313918
SRR11313920
SRR11313922

For WT 7h treatment with TTX/APV followed by 6h KCl stimulation:
SRR11313924
SRR11313926
SRR11313928

To make a tag directory for each sample run:
```
sbatch HOMER_MakeTags.sh
```

### UCSC Genome Browser Tracks
The basic strategy HOMER uses is to create a bedGraph formatted file that can then be uploaded as a custom track to the genome browser. This is accomplished using the makeUCSCfile program. To make a ucsc visualization file, type the following. To visualize the exact length of the reads, use "-fragLength given".

makeUCSCfile -o auto -fragLength given 

i.e. makeUCSCfile PU.1-RNA-Seq/ -o auto
(output file will be in the PU.1-RNA-Seq/ folder named PU.1-RNA-Seq.ucsc.bedGraph.gz)

The "-o auto" with make the program automatically generate an output file name (i.e. TagDirectory.ucsc.bedGraph.gz) and place it in the tag directory which helps with the organization of all these files. The output file can be named differently by specifying "-o outputfilename" or by simply omitting "-o", which will send the output of the program to stdout (i.e. add " > outputfile" to capture it in the file outputfile). It is recommended that you zip the file using gzip and directly upload the zipped file when loading custom tracks at UCSC.

To visualize the experiment in the UCSC Genome Browser, go to Genome Browser page and select the appropriate genome (i.e. the genome that the sequencing tags were mapped to). Then click on the "add custom tracks" button (this will read "manage custom tracks" once at least one custom track is loaded). Enter the file created earlier in the "Paste URLs or data" section and click "Submit".

There are two important parameters to consider during normalization of data. First, the total read depth of the experiment is important, which is obvious. The 2nd factor to consider is the length of the reads (this is new to v4.4). The problem is that if an experiment has longer fragment lengths, it will generate additional coverage than an experiment with shorter fragment lengths. In order to make sure there total area under the curve is the same for each experiment, experiments are normalized to a fixed number of reads as well as a 100 bp fragment length. If reads are longer than 100 bp, they are 'down-normalized' a fractional amount such that they produce the same relative coverage of a 100 bp fragment length. Experiments with shorter fragment lengths are 'up-normalized' a proportional amount (maximum of 4x or 25 bp). This allows experiments with different fragment lengths to be comparable along the genome browser.
Normalize the total number of reads to this number, default 1e7. This means that tags from an experiment with only 5 million mapped tags will count for 2 tags apiece. This script also fixes the chromosome names to include "chr" and change the MT chromosome to M. It also makes bigwig files needed for making a browser hub described below.

```
sbatch UCSCBrowserHOMER.sh
```

The bedGraph.gz files then then be loaded one at a time as custom tracks onto the UCSC genome browser. You can save the session and then look at them again later. However, this is very slow as you have to load each final individually. Instead, we can create a track hub, where all of our files can be loaded as a custom hub. 

We first have to convert our bedGraph and peak bed files into compressed formats: bigwig and bigBed. We can so this using the UCSC genome browser utilities. These tools are already installed in /alder/data/cbh/ciernia-data/pipeline-tools/UCSC/ using rsync -aP rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/linux.x86_64/ ./ 
The path was added to the bash_profile: PATH=$PATH:/alder/data/cbh/ciernia-data/pipeline-tools/UCSC and so the tools can be called by name. We also need the chromosome sizes for mm10. This can be retrieved from UCSC with: wget http://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes 
A copy of this file is in /alder/data/cbh/ciernia-data/pipeline-tools/UCSC/ and can be called using $mm10chrsizes. This was all done using the UCSCBrowserHOMER.sh script run.

Now that we have our compressed files we can setup our track hub:
Track hubs require a webserver to host the files. We can use github to host all files < 25MB. Github supports byte-range access to files when they are accessed via the raw.githubusercontent.com style URLs. To obtain a raw URL to a file already uploaded on Github, click on a file in your repository and click the Raw button. The bigDataUrl field (and any other statement pointing to a URL like bigDataIndex, refUrl, barChartMatrixUrl, etc.) of your trackDb.txt file should use the "raw.githubusercontent.com" style URL. 
https://genome.ucsc.edu/goldenPath/help/hgTrackHubHelp.html 

Track hubs can now be specified in a single text file: http://genome.ucsc.edu/goldenPath/help/hgTracksHelp.html#UseOneFile 
Once this text file is loaded onto github, you can get the RAW url for the text file and then check your hub works by pasting the url into the hub development took and clicking "Check Hub Settings". http://genome.ucsc.edu/cgi-bin/hgHubConnect?#hubDeveloper If your hub has no errors you can then click the "View on UCSC browser" to view your hub.

We have setup a hub using the TrackHubmm10.txt file. The link to this file is:
We have setup a hub using the TrackHubmm10.txt file. The link to this file is: http://microgliome.biochem.ubc.ca/BigWigs/Wenderski_RNAseq/TrackHubmm10_Wenderski_Master.txt
If you paste this link into the "My Hubs" entry area: http://genome.ucsc.edu/cgi-bin/hgTracks?db=mm10&hubUrl=http://microgliome.biochem.ubc.ca/BigWigs/Wenderski_RNAseq/TrackHubmm10_Wenderski_Master.txt you can view the hub. You can edit the hub, add tracks etc. Then just check the hub, reload it and go!
