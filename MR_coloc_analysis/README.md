# Analysis scripts 
## Major depression and atherosclerotic disease: From shared genetics to shared pathways in blood, brain, heart, and atherosclerotic plaques

List of scripts:

setup.R
 - config script for this project (install and load packages, define custom functions)
 - called by most other scripts

add_rsid.R
- preprocessing script for GWAS files, adds rsid by matching chr and pos to EUR reference (further comments in file)
- Input: GWAS summary statistics as received
- Output: .tsv file of somewhat preprocessed summary statistics

HDL.R 
- script to run genetic correlations using High Definition Likelihood method

MR_DEP-CVD.Rmd
- MR with exposure DEP and outcomes CVD
- Input: .pgc file, CVD *_with_rsid.tsv files
- Output: list of outcomes .RData, list of results .RData
- .Rmd file also includes presentation of results (not saved)

MR_CVD-DEP.Rmd 
- reverse MR with exposures CVD and outcome DEP
- Input: .pgc file, CVD *_with_rsid.tsv files
- Output: list of outcomes .RData, list of results .RData
- .Rmd file also includes presentation of results (not saved)

forestplot_MR.R
- generate forestplot for MR results (forward and reverse)

coloc_abf, coloc_susie.R
- colocalisation analysis using coloc package. one version uses coloc.susie function, but as there were no meaningful results we switched back to coloc.adf

RACER_plots.Rmd
- generate RACER plots (mirror plots from RACER package) for regions with high PP.H4 in coloc


preprocess_eQTL_for_SMR.R
 1) create sum stats for each probe, one variant per row like GWAS
 2) create file list with one row per probe
 3) concatenate result files of all probes

preprocess_GWAS_for_SMR.R
- change GWAS format for SMR software
- format: COJO

postprocess_SMR_results.R
- check and visualize results downloaded from SMR portal
- results were from full GWAS, then we selected the coloc regions
- make heatmap figure across tissues
- load plaque SMR results and visualize z-scores (nominally-significant)
