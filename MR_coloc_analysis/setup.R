
#### Install and load packages ####

## Special separate package installations
# 
# install.packages("devtools") # to install from github
# devtools::install_github("explodecomputer/genetics.binaRies") # to get a path to plink as described here https://mrcieu.github.io/ieugwasr/articles/local_ld.html
# devtools::install_github("oliviasabik/RACER") 
# 
# install.packages("remotes") # to install from github
# remotes::install_github("zhenin/HDL/HDL") # rg with high density likelihood method
# remotes::install_github("MRCIEU/TwoSampleMR")
# remotes::install_github("mrcieu/ieugwasr")
# remotes::install_github("chr1swallace/coloc",build_vignettes=TRUE)


# install.packages("BiocManager")
# BiocManager::install("IlluminaHumanMethylation450kmanifest")
# BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19")

## List of all packages
packages = c("coloc",
             "data.table",
             "flextable",
             "forestplot",
             "ggplot2",
             "here", 
             "HDL",
             "kableExtra",
             "ieugwasr",
             "IlluminaHumanMethylation450kmanifest",
             "IlluminaHumanMethylation450kanno.ilmn12.hg19",
             "data.table",
             "magrittr",
             "papaja",
             "patchwork",
             "png",
             "RACER",
             "readr",
             "reshape2",
             "tidyverse",
             "TwoSampleMR")

## Now load or install&load all
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)


# Own/Adapted Functions

# clump_data_locally()
# this is adapted from TwoSampleMR::clump_data
# purpose: enable use of locally stored ld reference files and plink version as opposed to accessing external server
# changes: 
#   - added input arguments bfile and plink_bin
#   - calling own function ld_clump_here rather than ieugwasr::ld_clump
clump_data_locally <- function (dat, clump_kb = 10000, clump_r2 = 0.001, clump_p1 = 1, 
                                clump_p2 = 1, pop = "EUR", bfile = NULL, plink_bin = plink_bin) 
{
  pval_column <- "pval.exposure"
  if (!is.data.frame(dat)) {
    stop("Expecting data frame returned from format_data")
  }
  if ("pval.exposure" %in% names(dat) & "pval.outcome" %in% 
      names(dat)) {
    message("pval.exposure and pval.outcome columns present. Using pval.exposure for clumping.")
  }
  else if (!"pval.exposure" %in% names(dat) & "pval.outcome" %in% 
           names(dat)) {
    message("pval.exposure column not present, using pval.outcome column for clumping.")
    pval_column <- "pval.outcome"
  }
  else if (!"pval.exposure" %in% names(dat)) {
    message("pval.exposure not present, setting clumping p-value to 0.99 for all variants")
    dat$pval.exposure <- 0.99
  }
  else {
    pval_column <- "pval.exposure"
  }
  if (!"id.exposure" %in% names(dat)) {
    dat$id.exposure <- random_string(1)
  }
  d <- data.frame(rsid = dat$SNP, pval = dat[[pval_column]], 
                  id = dat$id.exposure)
  
  out <- ld_clump_here(d, clump_kb = clump_kb, clump_r2 = clump_r2, 
                            clump_p = clump_p1, pop = pop, bfile = bfile, plink_bin = plink_bin)
  
  keep <- paste(dat$SNP, dat$id.exposure) %in% paste(out$rsid, out$id)
  return(dat[keep, ])
}


# ld_clump_here()
# based on ieugwasr::ld_clump
# call ld_clump_locally (my adapted function) rather than ld_clump_local (ieugwasr function)
ld_clump_here <- function (dat = NULL, clump_kb = 10000, clump_r2 = 0.001, clump_p = 0.99,
          pop = "EUR", access_token = NULL, bfile = NULL, plink_bin = NULL)
{
  stopifnot("rsid" %in% names(dat))
  stopifnot(is.data.frame(dat))
  if (is.null(bfile)) {
    message("Please look at vignettes for options on running this locally if you need to run many instances of this command.")
  }
  if (!"pval" %in% names(dat)) {
    if ("p" %in% names(dat)) {
      warning("No 'pval' column found in dat object. Using 'p' column.")
      dat[["pval"]] <- dat[["p"]]
    }
    else {
      warning("No 'pval' column found in dat object. Setting p-values for all SNPs to clump_p parameter.")
      dat[["pval"]] <- clump_p
    }
  }
  if (!"id" %in% names(dat)) {
    dat$id <- random_string(1)
  }
  if (is.null(bfile)) {
    access_token = check_access_token()
  }
  ids <- unique(dat[["id"]])
  res <- list()
  for (i in 1:length(ids)) {
    x <- subset(dat, dat[["id"]] == ids[i])
    if (nrow(x) == 1) {
      message("Only one SNP for ", ids[i])
      res[[i]] <- x
    }
    else {
      message("Clumping ", ids[i], ", ", nrow(x), " variants, using ",
              pop, " population reference")
      if (is.null(bfile)) {
        res[[i]] <- ld_clump_api(x, clump_kb = clump_kb,
                                 clump_r2 = clump_r2, clump_p = clump_p, pop = pop,
                                 access_token = access_token)
      }
      else {
        res[[i]] <- ld_clump_locally(x, clump_kb = clump_kb,
                                   clump_r2 = clump_r2, clump_p = clump_p, bfile = bfile,
                                   plink_bin = plink_bin)
      }
    }
  }
  res <- dplyr::bind_rows(res)
  return(res)
}


# ld_clump_locally
# adapted from original ieugwasr::ld_clump_local
# using plink2 the clumped file does not have a header. therefore read.table header = FALSE 
# and afterwards when subsetting from res use V3 as colname
ld_clump_locally <- function (dat, clump_kb, clump_r2, clump_p, bfile, plink_bin) 
{
  shell <- ifelse(Sys.info()["sysname"] == "Windows", "cmd", 
                  "sh")
  fn <- tempfile()
  write.table(data.frame(SNP = dat[["rsid"]], P = dat[["pval"]]), 
              file = fn, row.names = F, col.names = T, quote = F)
  fun2 <- paste0(shQuote(plink_bin, type = shell), " --bfile ", 
                 shQuote(bfile, type = shell), " --clump ", shQuote(fn, 
                                                                    type = shell), " --clump-p1 ", clump_p, " --clump-r2 ", 
                 clump_r2, " --clump-kb ", clump_kb, " --out ", shQuote(fn, 
                                                                        type = shell))
  system(fun2)
  res <- read.table(paste(fn, ".clumps", sep = ""), header = F)
  unlink(paste(fn, "*", sep = ""))
  y <- subset(dat, !dat[["rsid"]] %in% res[["V3"]])
  if (nrow(y) > 0) {
    message("Removing ", length(y[["rsid"]]), " of ", nrow(dat), 
            " variants due to LD with other variants or absence from LD reference panel")
  }
  return(subset(dat, dat[["rsid"]] %in% res[["V3"]]))
}
