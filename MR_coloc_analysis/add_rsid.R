#
# Non-automated script for sumstats formatting
# Issues addressed:
#   - add rsid (currently for EUR ancestry)
#   - standardise column names
#   - capitalise alleles 
#   - remove non-biallelic SNPs
#
# Overall process:
# identify chromosome and base pair columns
# link to the European ancestry reference file by chromosome and base pair thus matching with rsid column in the reference file
# return the original sumstats dataframe with additional rsid column
#
# Manual changes required: 
#   - pheno and filenames line 25-26
#   - i.a. path to input and output folder in line 28-29
#   - colnames of chr/bp/a1/a2 in line 42 
#   - if changing more colnames, adjust line 68
#
library(data.table)
library(here)
library(tidyverse)

# Name changing and checking section
phenotype <- "DEP"
original_filename <- "S3_table.csv"

readfrom <- here("Data", phenotype, original_filename)
writeto <- here("Data", phenotype, paste0(phenotype, "_with_RSID"))

# Read EUR file # Only run once
EUR <- fread(here(" ") ) # for rsid mapping add link to ld file (EUR). You can get them e.g. at https://vu.data.surfsara.nl/index.php/s/VZNByNwpD8qqINe

colnames(EUR) <- c("CHR", "RSID", "DUMMY", "BP", "A1", "A2")
coordinates <- c("CHR", "BP")

# Read sumstats
sumstats <- fread(readfrom)

head(sumstats) # pause to check

##### Necessary names for this script
sumstats <- sumstats %>% rename("CHR" = "CHR", "BP" = "BP", "effect_allele" = "A1", "other_allele" = "A2") # insert to the right the old names of sumstats columns


# Capitalize allele names if necessary
sumstats <- sumstats %>% mutate(
  effect_allele = toupper(effect_allele),
  other_allele = toupper(other_allele)
)

# More checks
table(sumstats$effect_allele)
table(sumstats$other_allele)

# Merge EUR and sumstats horizontally
df_temp <- EUR %>% inner_join(sumstats, by = coordinates, suffix = c("a", "b"), multiple = "all")
head(df_temp) # nervously check it

# Make sure the alleles are the same (taking both orders into account)
# This should also take care of not biallelic SNPs
df <- df_temp %>% filter((A1 == effect_allele & A2 == other_allele) | (A1 == other_allele & A2 == effect_allele))
table(df$effect_allele, df$other_allele) # reasonable allele combinations?

df_with_rsid <- as.data.frame(df[, - c("DUMMY", "A1", "A2")]) # save without extra EUR cols
head(df_with_rsid) 

# Optional colname adjustments. e.g. for MR scripts
df_with_rsid <- df_with_rsid %>% rename("effect_allele_frequency" = "Freq1", "beta" = "Effect", "standard_error" = "StdErr", "p_value" = "P") # insert to the right the old names of sumstats columns

# Store
# Maybe as R workspace
# save(df_with_rsid, file = paste0(writeto, ".RData"))
# Maybe as csv
write_tsv(df_with_rsid, file = paste0(writeto, ".tsv"))
