#### PLot forunivariatre MR in the paper
### based on forestplot package https://cran.r-project.org/web/packages/forestplot/vignettes/forestplot.html

rm(list = ls())
source("Scripts/setup.R")

sci_10_format <- function(x, digits = 2) {
  ifelse(
    x == 0,
    "0",
    sprintf("%.2f × 10^%d", signif(x / 10^floor(log10(abs(x))), digits), floor(log10(abs(x))))
  )
}


load(here("Results", "MR", "MR_DEP-CVD_results.RData"))

mr_beta <- as.data.frame(do.call("rbind", map(main_results_mr, "mr_output"))) |>
  mutate(b = b * log(2))
special_or <- generate_odds_ratios(mr_beta)

ivw_forward <- special_or |>
  filter(method == "Inverse variance weighted") |>
  select(outcome, exposure, nsnp, b, se, or, or_lci95, or_uci95, pval) |>
  # FDR correction
  mutate(pval.fdr = p.adjust(as.numeric(pval), method = "BH")) |>
  # Rounding
  mutate(across(c(b, se), ~ round(.x, digits = 4))) |>
  mutate(across(c(or, or_lci95, or_uci95), ~ round(.x, digits = 2))) |>
  mutate(across(c(pval, pval.fdr), ~ formatC(.x, format = "e", digits = 2)))
ivw_forward

# sci_10_format(mr_beta$pval)



# Cochrane data from the 'rmeta'-package
cochrane_from_rmeta <-
  structure(list(
    mean  = c(NA, 1.1000, 1.1067, 1.0718, 1.0948, 1.2262, 1.1958, 1.0429, 1.0004),
    lower = c(NA, 1.0458, 1.0479, 0.9624, 0.9490, 1.0641, 1.1349, 0.9212, 0.9957),
    upper = c(NA, 1.1569, 1.1688, 1.1938, 1.2629, 1.4129, 1.2599, 1.1806, 1.0051)
  ))

ivw_forward$b
ivw_forward$b - ivw_forward$se
ivw_forward$b + ivw_forward$se
ivw_forward$b

cochrane_forward_beta <-
  structure(list(
    mean  = c(NA, 0.0953, 0.1014, 0.0694, 0.0905, 0.2039, 0.1788, 0.0420, 0.0004),
    lower = c(NA, 0.0695, 0.0735, 0.0144, 0.0176, 0.1316, 0.1521, -0.0213, -0.0020),
    upper = c(NA, 0.1211, 0.1293, 0.1244, 0.1634, 0.2762, 0.2055, 0.1053, 0.0028)
  ))

tabletext <- cbind(
  c("N SNPs", ivw_forward$nsnp),
  c("Exposure", ivw_forward$exposure),
  c("Outcome", ivw_forward$outcome),
  c("OR", "1.1", "1.11", "1.07", "1.09", "1.23", "1.2", "n.a.", "n.a."), # add note to CAC and CIMT!!
  c("p-value", "*9.47 x 10\u207B\u2078", "*1.52 x 10\u207B\u2077", "6.86 x 10\u207B\u00B2", "7.32 x 10\u207B\u00B2", "*4.76 x 10\u207B\u2075", "*3.76 x 10\u207B\u00B2\u00B2", "3.38 x 10\u207B\u00B9", "8.28 x 10\u207B\u00B9")
  # c("p-value (FDR)", ivw_forward$pval.fdr)
)

forestplot(tabletext,
  mean = cochrane_from_rmeta$mean,
  lower = cochrane_from_rmeta$lower,
  upper = cochrane_from_rmeta$upper,
  boxsize = 0.2,
  new_page = TRUE,
  xticks = c(0.9, 1.0, 1.1, 1.2, 1.3),
  xlog = TRUE,
  col = fpColors(box = "black", line = "black"),
  vertices = TRUE
)

# saved as pdf 12x8

# Same for reverse MR

load(here("Results", "MR", "CVD_ReverseMR_results.RData"))

reverse_beta <- as.data.frame(do.call("rbind", map(main_results_mr, "mr_output"))) |>
  mutate(b = b * log(2))
special_or <- generate_odds_ratios(reverse_beta)

ivw_reverse <- special_or |>
  filter(method == "Inverse variance weighted") |>
  select(outcome, exposure, nsnp, b, se, or, or_lci95, or_uci95, pval) |>
  # FDR correction
  mutate(pval.fdr = p.adjust(as.numeric(pval), method = "BH")) |>
  # Rounding
  mutate(across(c(b, se), ~ round(.x, digits = 4))) |>
  mutate(across(c(or, or_lci95, or_uci95), ~ round(.x, digits = 2))) |>
  mutate(across(c(pval, pval.fdr), ~ formatC(.x, format = "e", digits = 2)))

ivw_reverse

# Cochrane data from the 'rmeta'-package
cochrane_reverse <-
  structure(list(
    mean  = c(NA, 1.0172, 1.0122, 1.0058, 1.0091, NA, 1.0027, 1.0063, 0.8058),
    lower = c(NA, 0.9802, 0.9702, 0.9819, 0.9832, NA, 0.9870, 0.9958, 0.5536),
    upper = c(NA, 1.0555, 1.0560, 1.0303, 1.0356, NA, 1.0187, 1.0169, 1.1729)
  ))

ivw_reverse$b
ivw_reverse$b - ivw_reverse$se
ivw_reverse$b + ivw_reverse$se

cochrane_reverse_beta <-
  structure(list(
    mean = c(NA, 0.0170, 0.0121, 0.0058, 0.009, NA, 0.0027, 0.0063, -0.2159),
    lower = c(NA, -0.0019, -0.0095, -0.0065, -0.0043, NA, -0.0054, 0.0009, -0.2159),
    upper = c(NA, 0.0359, 0.0337, 0.0181, 0.0223, NA, 0.0108, 0.0117, -0.0244)
  ))

tabletext2 <- cbind(
  c("N SNPs", "23", "23", "7", "3", "0", "173", "6", "8"),
  c("Exposure", "ALLSTROKE", "IS", "CES", "LAS", "SVD", "CAD", "CAC", "CIMT"),
  c("Outcome", "DEP", "DEP", "DEP", "DEP", "DEP", "DEP", "DEP", "DEP"),
  c("OR", "1.02", "1.01", "1.01", "1.01", "n.a.", "1.0", "1.01", "0.81"),
  c("p-value", "1.93 x 10\u207B\u00B9", "4.18 x 10\u207B\u00B9", "4.96 x 10\u207B\u00B9", "3.26 x 10\u207B\u00B9", "n.a.", "6.31 x 10\u207B\u00B9", "9.18 x 10\u207B\u00B2", "1.04 x 10\u207B\u00B9")
  # c("p-value (FDR)", "4.51e-01", "5.78e-01", "5.78e-01", "5.70e-01", "NA", "6.31e-01", "3.63e-01", "3.63e-01")
)


forestplot(tabletext2,
  mean = cochrane_reverse$mean,
  lower = cochrane_reverse$lower,
  upper = cochrane_reverse$upper,
  boxsize = 0.2,
  new_page = TRUE,
  xticks = c(0.9, 1.0, 1.1, 1.2, 1.3),
  xlog = TRUE,
  col = fpColors(box = "black", line = "black"),
  vertices = TRUE
)



new_text <- rbind(tabletext, tabletext2)
new_text
forestplot(new_text,
  mean = c(cochrane_forward_beta$mean, cochrane_reverse_beta$mean),
  lower = c(cochrane_forward_beta$lower, cochrane_reverse_beta$lower),
  upper = c(cochrane_forward_beta$upper, cochrane_reverse_beta$upper),
  boxsize = 0.2,
  new_page = TRUE,
  xticks = c(-0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3),
  xlog = FALSE,
  col = fpColors(box = "black", line = "black"),
  vertices = TRUE
)

# save as pdf 12x8
