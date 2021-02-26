# Wenderski RNAseq Analysis

Processing RNA-seq PE files sourced from the Crabtree Lab: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7211998/.

The RNAseq dataset is downloaded from published PE RNA-seq files from GEO and processing them for QC, aligning and analysis in Rstudio. The analysis is setup to run on the DMCBH Alder computing cluster using bash scripts.

## Step 1: Making a bash_profile

In your home directory make a .bash_profile in a text editor.

```
cd ~/ && nano .bash_profile
```
Remember to source the changes to your bash_profile.txt document before they take effet in the current terminal session.

Your bash profile will automatically load each time you start a new session.

```
source ~/.bash_profile
```

Check path variable.

```
echo $PATH
```
