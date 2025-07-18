rm(list = ls())

source("Scripts/setup.R")

# different ld reference panels
LD.path <- "ld-files/UKB_imputed_SVD_eigen99_extraction" # HapMap3


# phenotype
GWAS_list_full <- list("ALLSTROKE", "IS", "CES", "LAS", "SVD", "CAD", "DEP", "CAC", "CIMT")
combinations <- data.frame(combn(GWAS_list_full, 2))

#### Preprocess ####
NCAS_gigastroke <- c(ALLSTROKE = 73652, IS = 62100, CES = 10804, LAS = 6399, SVD = 6811)
NCON_gigastroke <- 1234808

# N_millionhearts <- 922788 # actually there is an N given in the sumstats. why did I put this one?


hdl_cols <- c("SNP", "A1", "A2", "b", "se", "N")

# Read and preprocess PGC depression sumstats
pgc_columns <- c("ID", "A1", "A2", "BETA", "SE", "NCAS", "NCON")
DEP <- as.data.frame(fread(here("Data", "DEP", "pgc-mdd2022-full-eur-v3.49.24.11.pgc"), select = pgc_columns))
DEP <- DEP |>
  mutate(N = NCAS + NCON, .keep = "unused")
colnames(DEP) <- hdl_cols

# Read CAC sumstats
CAC_columns <- c("variant_id", "effect_allele", "other_allele", "beta", "standard_error", "n")
CAC <- as.data.frame(fread(here("Data", "CAC", "CAC1000G_EA_FINAL_FULL.sumstats.gwascatalog.ssf.tsv"), select = CAC_columns))
colnames(CAC) <- hdl_cols

# Read CIMT sumstats
CIMT_columns <- c("refid", "Allele1", "Allele2", "Effect", "StdErr", "TotalSampleSize")
CIMT <- as.data.frame(fread(here("Data", "CIMT", "IMT.EA.META.MAF1.HetDF4_jun.csv"), select = CIMT_columns))
colnames(CIMT) <- hdl_cols

# Read and preprocess millionhearts CAD sumstats
cad_columns <- c("RSID", "effect_allele", "other_allele", "beta", "standard_error", "N")
CAD <- as.data.frame(fread(here("Data", "CAD", "CAD_with_RSID.tsv"), select = cad_columns))
colnames(CAD) <- hdl_cols

# Read gigastroke sumstats
dfs <- list()
for (i in c("ALLSTROKE", "IS", "CES", "LAS", "SVD")){
  filename <- paste0(i, "_with_RSID.tsv")
  path <- here("Data", i, filename)
  df <- as.data.frame(fread(path, select = c("RSID", "effect_allele", "other_allele","beta", "standard_error")))
  df["N"] <- rep(NCAS_gigastroke[i] + NCON_gigastroke, dim(df)[1])
  colnames(df) <- hdl_cols
  dfs[[i]] <- df
}

# Add DEP and CAD to the same list of dfs
dfs[["DEP"]] <- DEP
dfs[["CAD"]] <- CAD
dfs[["CAC"]] <- CAC
dfs[["CIMT"]] <- CIMT

save(dfs, file = here("Data", "all_preprocessed.RData"))
load(here("Data", "all_preprocessed.RData"))

#### HDL ####
GWAS_list_full <- list("ALLSTROKE", "IS", "CES", "LAS", "SVD", "CAD", "DEP", "CAC", "CIMT")
combinations <- data.frame(combn(GWAS_list_full, 2)) # 45 combinations

# loop the correlation using combinations[1,i] as df1 and combinations[2,i] as df2, LD.path stays the same
rgs <- list()
for (i in colnames(combinations)){
  name1 <- unlist(combinations[1, i])
  name2 <- unlist(combinations[2, i])
  gwas1 <- dfs[[name1]]
  gwas2 <- dfs[[name2]]
  
  label <- paste(name1, name2, sep = "_")
  rgs[[label]] <- HDL.rg(gwas1, gwas2, LD.path)
}

#save(rgs, file = here("Int-results", "rgs.RData"))
rgs

# Z-scores heritability
# DEP
0.018031040/0.0004806078
rgs$ALLSTROKE_DEP$estimates.df["Heritability_2", "Estimate"]/rgs$ALLSTROKE_DEP$estimates.df["Heritability_2", "se"]
# ALLSTROKE
0.004981529/0.0004448765
rgs$ALLSTROKE_DEP$estimates.df["Heritability_1", "Estimate"]/rgs$ALLSTROKE_DEP$estimates.df["Heritability_1", "se"]
# IS
0.004982404/0.000418881
rgs$IS_CES$estimates.df["Heritability_1", "Estimate"]/rgs$IS_CES$estimates.df["Heritability_1", "se"]
# CES
0.002681622/0.0004507071
rgs$IS_CES$estimates.df["Heritability_2", "Estimate"]/rgs$IS_CES$estimates.df["Heritability_2", "se"]
# LAS
0.001072099/0.0004265539
rgs$LAS_DEP$estimates.df["Heritability_1", "Estimate"]/rgs$LAS_DEP$estimates.df["Heritability_1", "se"]
# SVD
0.001171547/0.0002441437
rgs$SVD_DEP$estimates.df["Heritability_1", "Estimate"]/rgs$SVD_DEP$estimates.df["Heritability_1", "se"]
# CAD
0.003344898/0.0003793176
rgs$CAD_DEP$estimates.df["Heritability_1", "Estimate"]/rgs$CAD_DEP$estimates.df["Heritability_1", "se"]
# CAC
0.11045764/0.019356614
rgs$DEP_CAC$estimates.df["Heritability_2", "Estimate"]/rgs$DEP_CAC$estimates.df["Heritability_2", "se"]
# CIMT
0.04631130/0.007474673
rgs$DEP_CIMT$estimates.df["Heritability_2", "Estimate"]/rgs$DEP_CIMT$estimates.df["Heritability_2", "se"]


# Summarise
## Create matrix
var_names <- c("ALLSTROKE", "IS", "CES", "LAS", "SVD", "CAD", "CAC", "CIMT", "DEP")
rg_matrix <- matrix(rep(0), 100, ncol = 9, nrow = 9)
se_matrix <- matrix(rep(0), 100, ncol = 9, nrow = 9)
p_matrix <- matrix(rep(0), 100, ncol = 9, nrow = 9)
diag(rg_matrix) <- 1
diag(se_matrix) <- NA
diag(p_matrix) <- NA
colnames(rg_matrix) <- rownames(rg_matrix) <- colnames(se_matrix) <- rownames(se_matrix) <- colnames(p_matrix) <- rownames(p_matrix) <- var_names

## Populate matrix with values
for(j in colnames(rg_matrix)) {
  for(i in rownames(rg_matrix)){
    combined_name <- paste(j, i, sep = "_")
    reverse_combined_name <- paste(i, j, sep = "_")
    
    if(!is.null(rgs[[combined_name]])){
      rg_matrix[j, i] <- rgs[[combined_name]]$rg
      se_matrix[j, i] <- rgs[[combined_name]]$rg.se
      p_matrix[j, i] <- rgs[[combined_name]]$P
    }
    if (!is.null(rgs[[reverse_combined_name]])){
      rg_matrix[j, i] <- rgs[[reverse_combined_name]]$rg
      se_matrix[j, i] <- rgs[[reverse_combined_name]]$rg.se
      p_matrix[j, i] <- rgs[[reverse_combined_name]]$P
    }
  }
}
# Table of rgs and standard errors
rg_df <- as.data.frame(round(rg_matrix,3)) # correlations
se_df <- as.data.frame(round(se_matrix,3)) # standard error

rgse_df <- as.data.frame(do.call(cbind, lapply(1:ncol(rg_df), function(i) paste0(rg_df[ , i], " (", se_df[ , i], ")"  ) ))) # both in one
colnames(rgse_df) <- rownames(rgse_df) <- var_names # and nice names

#save(rgse_df, file = here("Results", "rgse_2decimals.RData")) # save
# load(here("Results", "rgse_2decimals.RData")) # load


# For manuscript: only correlations with depression
rg_matrix[,"DEP"]
se_matrix[,"DEP"]
p_matrix[,"DEP"]
rg_table <- cbind(round(rg_matrix[-9,"DEP"],3), 
      round(se_matrix[-9,"DEP"],3), 
      formatC(p_matrix[-9,"DEP"], format = "e", digits = 2)) # all in one table
colnames(rg_table) <- c("rg with DEP", "SE", "P-value")
rg_table
#write.csv(rg_table, file = here("Results", "rg_table.csv"), row.names = TRUE)
