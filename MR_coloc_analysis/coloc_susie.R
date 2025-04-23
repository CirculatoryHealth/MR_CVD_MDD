
source("C:/Users/emma/Documents/git/mr-coloc/setup.R")

# Script for colocalisation analysis
# Adapted from https://github.com/LN-Ceru/LN-Internship2/blob/main/MDD_coloc.r

# Read in top DEP SNPs for region selection

significant_SNPs <- fread(here("sumstats", "DEP", "DEP_with_rsid.tsv"))

# Read in exposure data once
exposure_data = read_exposure_data(filename = here("sumstats", "DEP", "pgc-mdd2022-full-eur-v3.49.24.11.pgc"),
                                   sep = "\t",
                                   snp_col = "ID",
                                   beta_col = "BETA",
                                   se_col = "SE",
                                   eaf_col = "FCON", # this is frequency in controls
                                   effect_allele_col = "A1",
                                   other_allele_col = "A2",
                                   pval_col = "PVAL",
                                   ncase_col = "NCAS",
                                   ncontrol_col = "NCON",
                                   #min_pval = 1e-200,
                                   log_pval = FALSE,
                                   chr_col= "#CHROM",
                                   pos_col= "POS")


exposure_data$exposure = 'DEP'
exposure_data$id.exposure = 'DEP'

print(paste('exposure read + dim:', dim(exposure_data)))

# exposure_data$beta.exposure <- log(exposure_data$beta.exposure) to transfer from OR if necessary
exposure_data$MAF = ifelse(exposure_data$eaf.exposure > 0.5, 1 - exposure_data$eaf.exposure, exposure_data$eaf.exposure)
exposure_data$varbeta = (exposure_data$se.exposure * exposure_data$se.exposure)
exposure_data$cc_ratio = max(exposure_data$ncase.exposure) / (max(exposure_data$ncase.exposure) + max(exposure_data$ncontrol.exposure))

head(exposure_data)


# Coloc in a loop
coloc_full_output =  list()

# second phenotype (here defined as outcome)
CVD_phenos <- c("ALLSTROKE") #, "IS", "CAD")

# sample sizes for gigastroke phenos and CAD
NCAS <- c(ALLSTROKE = 73652, IS = 62100, CES = 10804, LAS = 6399, SVD = 6811, CAD = 117010) # for gigastroke these numbers are from the GWAS catalog, for CAD from the readme
NCON <- c(ALLSTROKE = 1234808, IS = 1234808, CES = 1234808, LAS = 1234808, SVD = 1234808, CAD = 805778)


for (i in CVD_phenos) { # loop to read in outcome, preprocess, run coloc

#i <- "ALLSTROKE"
# read in outcome data
file <- paste0(i, "_with_RSID.tsv")

outcome_data <- read_outcome_data(
  filename = here("sumstats", i, file),
  sep = "\t",
  #snps = exposure_data$ID, # shouldn't this be exposure_data$SNP?
  snp_col = "RSID",
  beta_col = "beta",
  se_col = "standard_error",
  eaf_col = "effect_allele_frequency",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  pval_col = "p_value",
  #ncase_col = "NCAS",
  #ncontrol_col = "NCON",
  log_pval = FALSE,
  chr_col = "CHR",
  pos_col = "BP"
)

outcome_data$outcome <- outcome_data$id.outcome <- i

print("outcome read + dim:")
print(dim(outcome_data))

# outcome data wrangling
outcome_data$samplesize.outcome <- rep(NCAS[i] + NCON[i], dim(outcome_data)[1] ) # add column with total samplesize
outcome_data$MAF <- ifelse(outcome_data$eaf.outcome > 0.5, 1 - outcome_data$eaf.outcome, outcome_data$eaf.outcome)
outcome_data$varbeta <- (outcome_data$se.outcome * outcome_data$se.outcome)
outcome_data$cc_ratio <- rep((NCAS[i] / (NCAS[i] + NCON[i])), dim(outcome_data)[1]) ## create case control ratio from external info


#ensuring full overlap of SNPS
exposure_harmonised = exposure_data[exposure_data$SNP %in% outcome_data$SNP,]
outcome_harmonised = outcome_data[outcome_data$SNP %in% exposure_data$SNP,]

dim(exposure_harmonised)


# Colocalisation datawrangling
exposure_harmonised$pos.exposure <- as.numeric(exposure_harmonised$pos.exposure) # is this necessary?


# INITIALISE OUTPUT DATAFRAME
snp_output <- as.data.frame(matrix(nrow = 0, ncol = 10))
names(snp_output) = c("CV Phenotype", "lead snp", "CHR", "POS", "nsnps", "PP.H0", "PP.H1", "PP.H2", "PP.H3", "PP.H4")


# Make sure the SNPs we are going to loop through exist in both datasets
signif_snps <- significant_SNPs[significant_SNPs$RSID %in% outcome_harmonised$SNP, ]$RSID

# which(signif_snps == "rs10950392") 265

for (j in 265) { # loop to run colocalisation for each of the significant SNPs in exposure
  
  print(paste('snp number', j))
  
  #extract top SNP rowdata
  snp_list_chr <- c()
  lead_snp <- as.character(signif_snps[j])
  lead_pos <- exposure_harmonised[exposure_harmonised$SNP == lead_snp, ]$pos.exposure[1]
  lead_chr <- exposure_harmonised[exposure_harmonised$SNP == lead_snp, ]$chr.exposure[1]
  
  #subset by chromosome 
  chr_output <- outcome_harmonised[outcome_harmonised$chr.outcome == lead_chr, ]
  
  print(paste("chromosome", lead_chr))
  
  # Select window on this chromosome of 1e6 bp around the lead SNP
  window <- exposure_harmonised %>% filter(chr.exposure == lead_chr & pos.exposure >= lead_pos - 1e6 & pos.exposure <= lead_pos + 1e6)
  # snp_list_chr = unlist(append(snp_list_chr, window$pos.exposure))
  # snp_list_chr = unique(window$pos.exposure)
  
  chr_exposure_select <- window
  chr_output_select = subset(chr_output, chr_output$pos.outcome %in% chr_exposure_select$pos.exposure)
  
  #harmonising to be sure of direction of effect
  harmonised_data <- harmonise_data(chr_exposure_select, chr_output_select)
  chr_exposure_1 <- harmonised_data[, names(harmonised_data) %in% c(names(exposure_data), "remove", "varbeta.y", "MAF.y", "cc_ratio.y")]
  chr_output_1 <- harmonised_data[, names(harmonised_data) %in% c(names(outcome_data), "remove", "varbeta.x", "MAF.x", "cc_ratio.x")]
  
  
  # calculate LD matrix for the relevant SNPs
  # online version is limited to 500 variants so it's essential to use local version
  # note that ieugwasr uses plink 1.9
  LD <- ieugwasr::ld_matrix(
            variants = chr_exposure_1$SNP,
            with_alleles = FALSE,
            #pop = "EUR",
            #opengwas_jwt = get_opengwas_jwt(),
            bfile = "C:/Users/emma/Documents/git/ld-files/g1000_eur",
            plink_bin = "C:/Users/emma/AppData/Local/R/win-library/4.3/genetics.binaRies/bin/plink.exe"
  )
  
  # make sure N is a scalar not a vector
  dataset_susie_1 = list(LD = LD, snp = chr_exposure_1$SNP, position = as.numeric(chr_exposure_1$pos.exposure), N = as.numeric(chr_exposure_1$samplesize.exposure[1]), beta = as.numeric(chr_exposure_1$beta.exposure), MAF = as.numeric(chr_exposure_1$MAF.y), varbeta = as.numeric(chr_exposure_1$se.exposure)^2, pvalues = as.numeric(chr_exposure_1$pval.exposure), type = "cc", s = chr_exposure_1$cc_ratio[1])
  dataset_susie_2 = list(LD = LD, snp = chr_output_1$SNP, position = as.numeric(chr_output_1$pos.outcome), N = as.numeric(chr_output_1$samplesize.outcome[1]), pvalues = as.numeric(chr_output_1$pval.outcome), beta = as.numeric(chr_output_1$beta.outcome), varbeta = as.numeric(chr_output_1$se.outcome)^2, MAF = as.numeric(chr_output_1$MAF.x), type = "cc", s = chr_output_1$cc_ratio[1])
  
  D1 <- runsusie(dataset_susie_1)
  D2 <- runsusie(dataset_susie_2)
  
  susie.res <- coloc.susie(D1, D2)
  
  # #prep input datasets
  # dataset_1 = list(snp = chr_exposure_1$SNP, position = as.numeric(chr_exposure_1$pos.exposure), N = as.numeric(chr_exposure_1$samplesize.exposure), beta = as.numeric(chr_exposure_1$beta.exposure), MAF = as.numeric(chr_exposure_1$MAF.y), varbeta = as.numeric(chr_exposure_1$se.exposure)^2, pvalues = as.numeric(chr_exposure_1$pval.exposure), type = "cc", s = chr_exposure_1$cc_ratio[1])
  # dataset_2 = list(snp = chr_output_1$SNP, position = as.numeric(chr_output_1$pos.outcome), N = as.numeric(chr_output_1$samplesize.outcome), pvalues = as.numeric(chr_output_1$pval.outcome), beta = as.numeric(chr_output_1$beta.outcome), varbeta = as.numeric(chr_output_1$se.outcome)^2, MAF = as.numeric(chr_output_1$MAF.x), type = "cc", s = chr_output_1$cc_ratio[1])
  # print("number of SNPS datasets 1 and 2")
  # print(length(dataset_1$snp))
  # print(length(dataset_2$snp))
  # 
  # print("leadsnp")
  # print(lead_snp)
  # #print("check data 1")
  # #check_dataset(dataset_1)
  # #plot_dataset(dataset_1)
  # #print("check data 2")
  # #check_dataset(dataset_2)
  # #plot_dataset(dataset_2)
  # 
  # # Analysis and storage
  # coloc <- coloc.abf(dataset_1, dataset_2)
  # coloc$chrom <- lead_chr
  # 
  # snp_output[j, 1] <- i
  # snp_output[j, 2] <- lead_snp
  # snp_output[j, 3] <- lead_chr
  # snp_output[j, 4] <- lead_pos
  # snp_output[j, 5] <- coloc[[1]][[1]]
  # snp_output[j, 6] <- coloc[[1]][[2]]
  # snp_output[j, 7] <- coloc[[1]][[3]]
  # snp_output[j, 8] <- coloc[[1]][[4]]
  # snp_output[j, 9] <- coloc[[1]][[5]]
  # snp_output[j, 10] <- coloc[[1]][[6]]
  
}


## WRITE RESULTS FOR PROTEIN

# bind results
coloc_full_output[[i]] <- snp_output

}

coloc_full_output <- do.call(rbind, coloc_full_output)

coloc_full_output <- as.data.frame(coloc_full_output)
coloc_full_output <- coloc_full_output[rev(order(coloc_full_output$PP.H4)), ]

# coloc$results[rev(order(coloc$results$SNP.PP.H4)),]

# for ALLSTROKE, rs10950392 and rs7351050 have PP.H4 > 0.8

write.csv(coloc_full_output, here("Results", "coloc_DEP-CVD_output.csv"), row.names = F)

coloc_full_output <- read.csv(here("Results", "coloc_DEP-CVD_output.csv"))
sensitivity(coloc_full_output, "H4 > 0.8")




