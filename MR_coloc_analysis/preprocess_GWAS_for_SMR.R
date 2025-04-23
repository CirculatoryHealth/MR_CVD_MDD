
# To COJO format for SMR
# https://yanglab.westlake.edu.cn/software/gcta/#COJO
# From the website: "Columns are SNP, the effect allele, the other allele, frequency of the effect allele, effect size, standard error, 
# p-value and sample size. The headers are not keywords and wi  ll be omitted by the program. 
# Important: "A1" needs to be the effect allele with "A2" being the other allele and "freq" should be the frequency of "A1".
# Note: 1) For a case-control study, the effect size should be log(odds ratio) with its corresponding standard error.

source("Scripts/setup.R")
DEP <- as.data.frame(fread(here("Data", "DEP", "pgc-mdd2022-full-eur-v3.49.24.11.pgc")))
DEP <- DEP |>
  mutate(SNP = ID,
         freq = (FCAS * NCAS + FCON * NCON) / (NCAS + NCON),
         b = BETA,
         se = SE,
         p = PVAL,
         N = NCAS + NCON, .keep = "unused") |>
  select(SNP, A1, A2, freq, b, se, p, N)
write.table(DEP, file = here("Data", "DEP", "DEP_COJO.txt"), sep = "\t", row.names = FALSE, quote = FALSE)



CAD <- as.data.frame(fread(here("Data", "CAD", "CAD_with_RSID.tsv")))
CAD <- CAD |>
  mutate(SNP = RSID,
         A1 = effect_allele,
         A2 = other_allele,
         freq = effect_allele_frequency,
         b = beta,
         se = standard_error,
         p = p_value,
         N = N, .keep = "unused") |>
  select(SNP, A1, A2, freq, b, se, p, N)
write.table(CAD, file = here("Data", "CAD", "CAD_COJO.txt"), sep = "\t", row.names = FALSE, quote = FALSE)


ALLSTROKE <- as.data.frame(fread(here("Data", "ALLSTROKE", "ALLSTROKE_with_RSID.tsv")))
ALLSTROKE <- ALLSTROKE |>
  mutate(SNP = RSID,
         A1 = effect_allele,
         A2 = other_allele,
         freq = effect_allele_frequency,
         b = beta,
         se = standard_error,
         p = p_value,
         N = rep(1308460, nrow(ALLSTROKE)), .keep = "unused") |>
  select(SNP, A1, A2, freq, b, se, p, N)
write.table(ALLSTROKE, file = here("Data", "ALLSTROKE", "ALLSTROKE_COJO.txt"), sep = "\t", row.names = FALSE, quote = FALSE)


CIMT <- as.data.frame(fread(here("Data", "CIMT", "CIMT_with_RSID.tsv")))
CIMT <- CIMT |>
  mutate(SNP = RSID,
         A1 = effect_allele,
         A2 = other_allele,
         freq = effect_allele_frequency,
         b = beta,
         se = standard_error,
         p = p_value,
         N = TotalSampleSize, .keep = "unused") |>
  select(SNP, A1, A2, freq, b, se, p, N)
write.table(CIMT, file = here("Data", "CIMT", "CIMT_COJO.txt"), sep = "\t", row.names = FALSE, quote = FALSE)





