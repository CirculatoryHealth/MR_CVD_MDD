
source("setup.R")

# different ld reference panels
LD.path <- "C:/Users/emma/Documents/git/ld_files/UKB_imputed_SVD_eigen99_extraction" # HapMap3


# phenotype
GWAS_list_full <- list("ALLSTROKE", "IS", "CES", "LAS", "SVD", "CAD", "DEP", "CAC", "CIMT_a", "CIMT_b")
combinations <- data.frame(combn(GWAS_list_full, 2))

#### Preprocess ####
NCAS_gigastroke <- c(ALLSTROKE = 73652, IS = 62100, CES = 10804, LAS = 6399, SVD = 6811)
NCON_gigastroke <- 1234808

# N_millionhearts <- 922788 # actually there is an N given in the sumstats. why did I put this one?


hdl_cols <- c("SNP", "A1", "A2", "b", "se", "N")

# Read and preprocess PGC depression sumstats
pgc_columns <- c("ID", "A1", "A2", "BETA", "SE", "NCAS", "NCON")
DEP <- as.data.frame(fread(here("sumstats", "DEP", "pgc-mdd2022-full-eur-v3.49.24.11.pgc"), select = pgc_columns))
DEP <- DEP |>
  mutate(N = NCAS + NCON, .keep = "unused")
colnames(DEP) <- hdl_cols

# Read CAC sumstats
CAC_columns <- c("variant_id", "effect_allele", "other_allele", "beta", "standard_error", "n")
CAC <- as.data.frame(fread(here("sumstats", "CAC", "CAC1000G_EA_FINAL_FULL.sumstats.gwascatalog.ssf.tsv"), select = CAC_columns))
colnames(CAC) <- hdl_cols

# Read CIMT sumstats
CIMT_columns <- c("refid", "Allele1", "Allele2", "Effect", "StdErr", "TotalSampleSize")
CIMT_a <- as.data.frame(fread(here("sumstats", "CIMT", "IMT.EA.META.MAF1.HetDF4_jun.csv"), select = CIMT_columns))
colnames(CIMT_a) <- hdl_cols

CIMT_b <- as.data.frame(fread(here("sumstats", "CIMT", "Plaque_meta_032218.csv"), select = CIMT_columns))
colnames(CIMT_b) <- hdl_cols

# Read and preprocess millionhearts CAD sumstats
cad_columns <- c("RSID", "effect_allele", "other_allele", "beta", "standard_error", "N")
CAD <- as.data.frame(fread(here("sumstats", "CAD", "CAD_with_RSID.tsv"), select = cad_columns))
colnames(CAD) <- hdl_cols

# Read gigastroke sumstats
dfs <- list()
for (i in c("ALLSTROKE", "IS", "CES", "LAS", "SVD")){
  filename <- paste0(i, "_with_RSID.tsv")
  path <- here("sumstats", i, filename)
  df <- as.data.frame(fread(path, select = c("RSID", "effect_allele", "other_allele","beta", "standard_error")))
  df["N"] <- rep(NCAS_gigastroke[i] + NCON_gigastroke, dim(df)[1])
  colnames(df) <- hdl_cols
  dfs[[i]] <- df
}

# Add DEP and CAD to the same list of dfs
dfs[["DEP"]] <- DEP
dfs[["CAD"]] <- CAD
dfs[["CAC"]] <- CAC
dfs[["CIMT_a"]] <- CIMT_a
dfs[["CIMT_b"]] <- CIMT_b

#save(dfs, file = here("sumstats", "all_preprocessed.RData"))
#load(here("sumstats", "all_preprocessed.RData"))

#### HDL ####

# For DEP, CAD, ALLSTROKE, IS, CES, LAS, SVD
# Run 15 genetic correlations
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

# Save that list and keep it safe
# save(rgs, file = here("rgs.RData"))
load(file = here("rgs.RData"))

# Run 6 genetic correlations for CAD
rgs_CAD <- list()
gwas1 <- dfs[["CAD"]]

for (i in c("DEP")){ # contains all phenotype names but CAD
  gwas2 <- dfs[[i]]
  label <- paste("CAD", i, sep = "_")
  rgs_CAD[[label]] <- HDL.rg(gwas1, gwas2, LD.path)
}

save(rgs_CAC, file = here("rgs_CAC.RData"))


# Run genetic correlations for CAC and CIMT
# CAC
rgs_CAC_2 <- list()
gwas1 <- CAC
for (i in c("CAD", "DEP", "CIMT_a", "CIMT_b")){ # contains all phenotype names but CAC
  gwas2 <- dfs[[i]]
  label <- paste("CAC", i, sep = "_")
  print(label)
  rgs_CAC_2[[label]] <- HDL.rg(gwas1, gwas2, LD.path)
}

save(rgs_CAC_2, file = here("rgs_CAC_2.RData"))

# CIMT_a
rgs_CIMT_a <- list()
gwas1 <- dfs[["CIMT_a"]]
for (i in c("ALLSTROKE", "IS", "CES", "LAS", "SVD", "CAD", "DEP", "CIMT_b")){ # contains all phenotype names but CAC & CIMT_a
  gwas2 <- dfs[[i]]
  label <- paste("CIMT_a", i, sep = "_")
  rgs_CIMT_a[[label]] <- HDL.rg(gwas1, gwas2, LD.path)
}

save(rgs_CIMT_a, file = here("rgs_CIMT_a.RData"))

# CIMT_b
rgs_CIMT_b <- list()
gwas1 <- dfs[["CIMT_b"]]
for (i in c("ALLSTROKE", "IS", "CES", "LAS", "SVD", "CAD", "DEP")){ # contains all phenotype names but CAC, CIMT_a, CIMT_b
  gwas2 <- dfs[[i]]
  label <- paste("CIMT_b", i, sep = "_")
  print(label)
  rgs_CIMT_b[[label]] <- HDL.rg(gwas1, gwas2, LD.path)
}

save(rgs_CIMT_b, file = here("rgs_CIMR_b.RData"))


# Make presentable output
load(here("Int-results", "rgs.RData"))
load(here("Int-results", "rgs_CAD.RData"))
load(here("Int-results", "rgs_CAC_1.RData"))
load(here("Int-results", "rgs_CAC_2.RData"))
load(here("Int-results", "rgs_CIMT.RData"))

## One structure for all correlations 
# this necessary due to messiness - there are two lists due to workflow and cors taking a long time to compute
all_rgs <- c(rgs, rgs_CAD, rgs_CAC, rgs_CAC_2, rgs_CIMT_a)

# Z-scores heritability
# DEP
0.018031040/0.0004806078
# ALLSTROKE
0.004981529/0.0004448765
# IS
0.004982404/0.000418881
# CES
0.002681622/0.0004507071
# LAS
0.001072099/0.0004265539
# SVD
0.001171547/0.0002441437
# CAD
0.003344898/0.0003793176
# CAC
0.11045764/0.019356614
# CIMT_A
0.04631130/0.007474673
# CIMT_B
0.01213623/0.006585895


# Prepare for plot
## Create matrix
var_names <- c("ALLSTROKE", "IS", "CES", "LAS", "SVD", "CAD", "CAC", "CIMT_a", "DEP")
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
    
    if(!is.null(all_rgs[[combined_name]])){
      rg_matrix[j, i] <- all_rgs[[combined_name]]$rg
      se_matrix[j, i] <- all_rgs[[combined_name]]$rg.se
      p_matrix[j, i] <- all_rgs[[combined_name]]$P
    }
    if (!is.null(all_rgs[[reverse_combined_name]])){
      rg_matrix[j, i] <- all_rgs[[reverse_combined_name]]$rg
      se_matrix[j, i] <- all_rgs[[reverse_combined_name]]$rg.se
      p_matrix[j, i] <- all_rgs[[reverse_combined_name]]$P
    }
  }
}

# Table of rgs and standard errors
rg_df <- as.data.frame(round(rg_matrix,2)) # correlations
se_df <- as.data.frame(round(se_matrix,2)) # standard error

rgse_df <- as.data.frame(do.call(cbind, lapply(1:ncol(rg_df), function(i) paste0(rg_df[ , i], " (", se_df[ , i], ")"  ) ))) # both in one
colnames(rgse_df) <- rownames(rgse_df) <- var_names # and nice names

save(rgse_df, file = here("Results", "rgse_2decimals.RData")) # save


# Plot - correlation matrix heatmap
# code from https://www.geeksforgeeks.org/how-to-create-correlation-heatmap-in-r/
# creating correlation matrix
corr_mat <- round(rg_matrix,2)
corr_mat <- corr_mat[nrow(corr_mat):1,] # flip it to appear with downward diagonal in the plot (trial and error don't know why this works)

# reduce the size of correlation matrix
melted_corr_mat <- melt(corr_mat, value.name = "rg")
head(melted_corr_mat)


# plotting the correlation heatmap
rg_heatmap <- ggplot(data = melted_corr_mat, aes(x=Var1, y=Var2, fill=rg)) + 
  geom_tile() +
  geom_text(aes(Var2, Var1, label = rg), color = "black", size = 4) +
  scale_fill_gradient(low = "#BBFFFF", high = "#00868B") + # check colours at https://r-charts.com/colors/
  labs(title = "Genetic correlations of CVD and depression",
       subtitle = "Calculated with R package HDL",
       x = NULL, y = NULL)
rg_heatmap

ggsave(filename = here("Results", "heatmap_plot.png"), plot = rg_heatmap)
