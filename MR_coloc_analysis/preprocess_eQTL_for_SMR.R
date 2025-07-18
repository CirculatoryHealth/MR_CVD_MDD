

# What does this script do?
# It prepares data for SMR analysis.

source("Scripts/setup.R")

# Function to round p-values as specified # mostly from ChatGPT
round_pvalue <- function(pvalue) {
  # Split into base and exponent
  parts <- strsplit(formatC(pvalue, format = "e", digits = 2), "e")[[1]]
  base <- as.numeric(parts[1])
  exponent <- as.integer(parts[2])
  
  # Round up the base to 2 decimal places
  rounded_base <- ceiling(base * 100 + 1) / 100  # Round up to the nearest hundredth + 1 to add a bit
  
  # Reassemble the number
  rounded_pvalue <- rounded_base * 10^exponent
  return(rounded_pvalue)
}

#### Regions for SMR analysis ####
# For the SMR analysis, we use eQTL and mQTL from certain regions only
# Namely, the regions around the lead SNPs identified by coloc for DEP and CVD comorbidity
# These regions are defined as 2,000,000 base pairs around the lead SNP
# Here, we make a list of these regions, to request the correct data
#
# I have a list of base pair positions:
leadSNP <- c("7:12263538", "19:4044579", "15:38843887", "7:114059156", "7:2760750", "7:1937261", "2:165045101", "2:164915279") 
# I need to provide a table of each of this positions minus 1,000,000 base pairs and plus 1,000,000 base pairs
# any digits before the colon should remain the same
lower_position <- as.numeric(sub(".*:", "", leadSNP)) - 1000000
upper_position <- as.numeric(sub(".*:", "", leadSNP)) + 1000000
# the digits before the colon represent chromosomes. print only those.
chromosome <- as.numeric(sub(":.*", "", leadSNP))
# now I need a table with columns leadSNP, chromosome, lower_position, upper_position
regions <- data.frame(leadSNP, chromosome, lower_position, upper_position)
# I need to save this table as a .csv file
write.csv(regions, here("Results", "regions.csv"), row.names = FALSE)


#### eQTL data ####
# Once cis eQTL data of these regions has been supplied, we want to explore it a bit
# Read in the file using fread
cis_eqtl <- fread(here("Data", "target_regions.nom_cis_eqtl_annotated.txt"))

# Check ZFAND2A in detail: Make scatterplot of effect sizes
ZFAND2A <- cis_eqtl |> 
  filter(symbol == "ZFAND2A") |> 
  arrange(pval_nominal)

GWAS <- fread(here("Data", "DEP", "DEP_COJO.txt"))

# Make dataframe for plot
GWAS_effects <- GWAS |> 
  select(SNP, b, se, p) |> 
  rename(b_GWAS = b, se_GWAS = se, p_GWAS = p)
eQTL_effects <- ZFAND2A |> 
  select(AltID, Beta, SE, pval_nominal) |>
  rename(SNP = AltID, b_eQTL = Beta, se_eQTL = SE, p_eQTL = pval_nominal)
effect_sizes <- inner_join(GWAS_effects, eQTL_effects, by = c("SNP")) |> 
  filter(p_eQTL < 0.000714)

plot(effect_sizes$b_eQTL, effect_sizes$b_GWAS,
     col = ifelse(effect_sizes$SNP == "rs73033363", "blue", "black"),
     pch = ifelse(effect_sizes$SNP == "rs73033363", 17, 1),
     cex = ifelse(effect_sizes$SNP == "rs73033363", 2, 1))

min(effect_sizes$p_eQTL)
effect_sizes |> 
  filter(p_eQTL < 2.172e-05)


# General investigation
unique(cis_eqtl$symbol) # which genes ended up in this database
dim(filter(cis_eqtl, pval_nominal < 0.05)) # how many SNPs are nominally significant
hist(filter(cis_eqtl, pval_nominal < 0.05)$pval_nominal) # distribution of p-values
dim(filter(cis_eqtl, pval_nominal < 5e-6)) # how many SNPs are genome-wide significant
# none. so now what?
table(cis_eqtl$chromosome) # all chromosomes are present

# Group cis_eqtl by EnsemblID and note number of rows as well as smallest p-value for each gene
cis_eqtl_summary <- cis_eqtl |> 
  group_by(symbol) |> 
  summarise(n = n(), min_pval = min(pval_nominal)) |> 
  arrange(min_pval) |> 
  mutate(pval_corrected = min_pval*n) # calculate Bonferroni corrected p-values


## Prepare files for SMR ##
# SMR flag will be used to create .befd file, see https://yanglab.westlake.edu.cn/software/smr/#MakeaBESDfile
# R is used to create input files in the correct format. 
# Steps:
# 1) create sum stats for each probe, one variant per row like GWAS
#   - columns: Chr, SNP, Bp, A1, A2, Freq, Beta, se, p
#   - remove multi-allelic SNPs. check for duplicates. if there are duplicates remove the rarer one and remove the same for each QTL
#   - FDR adjustment on p-values
# 
# 2) create file list with one row per probe
#   - columns: Chr, ProbeID, GeneticDistance, ProbeBp, Gene, Orientation, PathOfEsd
#   - get ProbeBp and Orientation from columns genomic_pos 
#
# Also, 
#   - check effect directions across QTLs
#
# 3) concatenate result files of all probes


# 1. Make myQTL.esd

# Filtering steps that are easier while it is still one dataframe
cis_eqtl <- cis_eqtl |> 
  filter(CAF_eQTL > 0.01 & CAF_eQTL < 0.99) |> # maf filter
  filter(nchar(CodedAlleleB) == 1 & nchar(OtherAlleleA) == 1) # remove multi-allelic SNPs

# how many distinct AltIDs are in cis_eqtl?
all_SNPs <- unique(cis_eqtl$AltID)

# Check whether the same effect allele was chosen across probes
for (i in all_SNPs){
  # if AltID = i, then CAF is the same for all rows
  those_rows <- filter(cis_eqtl, AltID == i)
  if( length(unique(those_rows$CAF_eQTL)) != 1){
    problems <- list(i, unique(those_rows$CAF_eQTL))
  }
}
# No problem!


# Make separate file for each probe
cis_eqtl_probefiles <- list()
# duplicated_snps <- list()
for (i in unique(cis_eqtl$symbol)) {
  cis_eqtl_probefiles[[i]] <- cis_eqtl |> 
    filter(symbol == i) |> 
    select(Chr = chromosome, SNP = AltID, Bp = position, A1 = CodedAlleleB, A2 = OtherAlleleA, Freq = CAF_eQTL, Beta = Beta, se = SE, p = pval_nominal) |> 
    # remove rows that are exact duplicates, i.e. the values in all columns are the same
    distinct() # remove duplicates
  
  # Check if there are (still) duplicate SNPs
  # duplicated_snps[[i]] <- sum(duplicated(cis_eqtl_probefiles[[i]]$SNP))
  # print(duplicated_snps)
  
  cis_eqtl_probefiles[[i]][,"FDR"] <- p.adjust(cis_eqtl_probefiles[[i]]$p, method = "BH")

}

# Preselect probes for SMR: There should be a cis eQTL with FDR < 0.1
# make a summary dataframe containing the minimum fdr found for each probe. 
fdr_summary <- data.frame()
for (i in unique(cis_eqtl$symbol)){
  FDR <- cis_eqtl_probefiles[[i]]$FDR
  if(!is.null(FDR)){
    fdr_summary[i,"min"] <- min(FDR, na.rm = TRUE)
  }
}
smr_genes <- rownames(fdr_summary |> filter(min <= 0.1)) # which probes have at least one SNP with FDR < 0.1?

rows_of_interest <- list()
temp_list <- list()
for(i in smr_genes){
   
  rows_of_interest[[i]] <- cis_eqtl_probefiles[[i]] |> 
    filter(FDR <= 0.1)
  
  cutoff <- rows_of_interest[[i]] |> 
    filter(FDR == max(FDR)) |>
    filter(p == max(p)) |> 
    select(p, FDR)
   
  print(cutoff)
  temp_list[[i]] <- cutoff[1,]
}

cutoff_df <- as.data.frame(do.call("rbind", temp_list))
rownames(cutoff_df) <- smr_genes

# Slightly rough process of rounding up the thresholds
rounded_pvalueThresholds <- sapply(cutoff_df$p, round_pvalue)


esd_prep <- list()
for(i in smr_genes){
  esd_prep[[i]] <- cis_eqtl_probefiles[[i]] |> 
    select(Chr, SNP, Bp, A1 , A2 , Freq , Beta ,se , p)
    
  # loop. for smr_genes only
  write.table(esd_prep[[i]], file = here("SMR", "esd_files", paste0(i, ".esd")), row.names = F, col.names = T, quote = F)
}


# 2. Make my.flist
file_list <- data.frame()
n <- 0
for (i in smr_genes) {
  n <- n + 1 # row counter of filelist
  info_row <- cis_eqtl |> 
    filter(symbol == i) |> 
    slice(1) # only first row
  file_list[1, "Chr"] <- info_row$chromosome
  file_list[1, "ProbeID"] <- paste0(info_row$symbol, "_", info_row$EnsemblID)
  file_list[1, "GeneticDistance"] <- 0
  file_list[1, "ProbeBp"] <- as.numeric(str_extract(info_row$genomic_pos, "(?<=\'start\': )[0-9]+"))
  file_list[1, "Gene"] <- info_row$symbol
  file_list[1, "Orientation"] <- as.numeric(str_extract(info_row$genomic_pos, "(?<=\'strand\': )[-+]?[0-9]+"))
  file_list[1, "PathOfEsd"] <- paste0("esd_files/", info_row$symbol, ".esd")
  
  # replace 1 with + and -1 with - just in case it is necessary
  file_list[1, "Orientation"] <- ifelse(file_list$Orientation == 1, "+", "-")
  
  # one file per gene
  write.table(file_list, file = here("SMR", paste0("probe" , n, ".flist")), row.names = F, col.names = T, quote = F)
}

# CONDUCT SMR (using command line script smr_script.ps1)

# 3. SMR Results
smr <- list()
for(i in 1:18){
  smr[[i]] <- fread(here("SMR", paste0("probe", i, ".smr.smr")))
}
smr_results <- as.data.frame(do.call("rbind", smr))

p.adjust(smr_results$p_SMR, method = "BH") # check FDR corrected pvalues

# add pvalue used as threshold for smr:
pvalue_thresholds <- as.data.frame(bind_cols("Gene"= smr_genes, "--peqtl-smr" = rounded_pvalueThresholds))
smr_results <- left_join(smr_results, pvalue_thresholds, by = c("Gene"))

# save
write.csv(smr_results, here("Results", "smr_results.csv"), row.names = F)
