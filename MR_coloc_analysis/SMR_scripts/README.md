# SMR analysis of gene expression in the plaque.

## Directory
- raw input files and preprocessed input files are stored in the top level folder SMR
- SMR results and log files are stored in each phenotype's eQTL_ folder
- alongside phenotype's SMR script smr_script_
Careful: This is a windows powershell script.
Usage: simply open powershell in folder of script and run .\smr_script.ps1

## History
First the analysis was run for depression only. Therefore, preprocessing steps were included in smr_script_DEP. There may be other traces of this.

## Preprocessing
- we use the first way described here https://yanglab.westlake.edu.cn/software/smr/#DataManagement
- create .flist files and the .esd files they refer to (preprocess_eQTL_for_SMR.R)
- from those files, generate .besd (and .epi and .esi) files
- those three file types are the input for SMR. 

## Files to have/Paths to adjust:
- LD reference files (tried and tested with 1000 genomes in plink file format)
- GWAS file in COJO format (COJO formatted by (preprocess_GWAS_for_SMR.R)
- probe file list named probe + counter + .flist (here contain only one probe each)
- SMR software for your OS (no such thing as 'obvious'. installation on Windows: download SMR https://yanglab.westlake.edu.cn/software/smr/#Download, unzip & arrange folder structure)

## Other prep required:
- decide on a pvalue threshold for each probe file list (here calculated by preprocess_eQTL_for_SMR.R)
