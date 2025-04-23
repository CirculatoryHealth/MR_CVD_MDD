source("Scripts/setup.R")

# Read in the results that were downloaded from SMR portal
# SMR was run for a selection of 20 tissues and the MDD GWAS summary statistics
eSMR_MDD_output <- fread(here("Results", "SMR-portal", "trait_eSMR.merged.tsv"))
dim(eSMR_MDD_output)

# Load shared regions
regions <- fread(here("Results", "regions.csv"))
regions


eSMR_MDD_output |>
  filter(ProbeChr %in% regions$chromosome) |>
  dim()

# Filter for shared regions
# Select only those probes from eSMR_MDD_output that are in the shared regions, so those that have Probe_bp between lower_position and upper_position for the same chromosome
eSMR_shared_regions <- list()
for (i in 1:nrow(regions)) {
  region <- regions[i, ]
  leadSNP <- region$leadSNP
  eSMR_shared_regions[[leadSNP]] <- eSMR_MDD_output |>
    filter(ProbeChr == region$chromosome) |>
    filter(Probe_bp >= region$lower_position & Probe_bp <= region$upper_position)
  # fwrite(here("Results", "SMR-portal", paste0("trait_eSMR.merged_", region$chromosome, "_", region$lower_position, "_", region$upper_position, ".tsv")))
}
eSMR_shared_regions_df <- as.data.frame(do.call("rbind", eSMR_shared_regions))

eSMR_shared_regions_df <- eSMR_shared_regions_df |>
  filter(duplicated(eSMR_shared_regions_df) == FALSE) # There are quite some overlapping lines, I suspect because of the overlapping regions.


# Save the SMR results from shared regions
# write.csv(eSMR_shared_regions_df, file = here("Results", "SMR-portal", "eSMR_shared_regions.csv"))

# Explore more
eSMR_shared_regions_df |>
  filter(p_SMR <= 5e-8)

eSMR_shared_regions_df |>
  filter(p_SMR <= 0.05) |>
  # print only unique probe names
  distinct(index)

eSMR_shared_regions_df |>
  filter(p_SMR <= 5e-4) |>
  view()


# Plot
portal_SMR_results <- fread(file = here("Results", "SMR-portal", "eSMR_shared_regions.csv"))
# dim(portal_SMR_results)

display_portal_results <- portal_SMR_results |>
  select(index, qtl_name, ProbeChr, b_SMR, se_SMR, p_SMR, p_HEIDI, nsnp_HEIDI) |>
  arrange(index) |>
  # FDR correction
  mutate(p.fdr_SMR = p.adjust(as.numeric(p_SMR), method = "BH")) |>
  relocate(p.fdr_SMR, .after = p_SMR) |>
  filter(p.fdr_SMR < 0.1) |>
  # in the SMR paper they don't apply multiple testing correction for the HEIDI test
  # this is conservative for gene discovery: H0 of shared causal variant is rejected under stricter conditions
  # mutate(p.fdr_HEIDI = p.adjust(as.numeric(p_HEIDI), method = "BH")) |>
  # relocate(p.fdr_HEIDI, .after = p_HEIDI) |>
  # add variable of HEIDI result. In this column, if HEIDI > 0.05 then add "*" and if HEIDI <= 0.05 then add ""
  mutate(HEIDI.result = ifelse(p_HEIDI > 0.05, "*", "")) |>
  # change tissue labels by replacing eQTL_ with null
  mutate(qtl_name = gsub("eQTL_", "", qtl_name)) |>
  mutate(qtl_name = gsub("GTEx_", "", qtl_name)) |>
  mutate(nice_qtl_name = case_when(
    qtl_name == "Brain_Frontal_Cortex_BA9" ~ "Frontal Cortex",
    qtl_name == "Brain_Putamen_basal_ganglia" ~ "Putamen (Basal Ganglia)",
    qtl_name == "Brain_Caudate_basal_ganglia" ~ "Caudate (Basal Ganglia)",
    qtl_name == "Brain_Hippocampus" ~ "Hippocampus",
    qtl_name == "BrainMeta" ~ "Brain (All Regions)",
    qtl_name == "Brain_Nucleus_accumbens_basal_ganglia" ~ "Nucleus Accumbens (Basal Ganglia)",
    qtl_name == "Brain_Cortex" ~ "Cortex",
    qtl_name == "eQTLGen" ~ "Whole Blood",
    qtl_name == "Adipose_Visceral_Omentum" ~ "Adipose (Visceral Omentum)",
    qtl_name == "Heart_Atrial_Appendage" ~ "Atrial Appendage",
    qtl_name == "Heart_Left_Ventricle" ~ "Left Ventricle",
    qtl_name == "Brain_Cerebellum" ~ "Cerebellum",
    qtl_name == "Brain_Cerebellar_Hemisphere" ~ "Cerebellar Hemisphere",
    qtl_name == "Artery_Aorta" ~ "Aorta",
    qtl_name == "Brain_Spinal_cord_cervical_c-1" ~ "Spinal Cord",
    qtl_name == "Brain_Hypothalamus" ~ "Hypothalamus",
    qtl_name == "Brain_Anterior_cingulate_cortex_BA24" ~ "Anterior Cingulate Cortex"
  ))

# Order x axis
desired_order_tissues <- c(
  "Whole Blood", "Brain (All Regions)",
  "Anterior Cingulate Cortex", "Caudate (Basal Ganglia)", "Cerebellar Hemisphere",
  "Cerebellum", "Cortex", "Frontal Cortex", "Hippocampus", "Hypothalamus",
  "Nucleus Accumbens (Basal Ganglia)", "Putamen (Basal Ganglia)", "Spinal Cord",
  "Adipose (Visceral Omentum)",
  "Atrial Appendage", "Left Ventricle",
  "Aorta"
)

display_portal_results$nice_qtl_name <- factor(display_portal_results$nice_qtl_name, levels = desired_order_tissues)

desired_order_genes <- c(
  "AMZ1", "ARL4A", "ARRDC5", "C7orf50", "CHST12", "EEF2", "GNA12", "GPR146", "GPR85", "HDGFL2",
  "IQCE", "MAD1L1", "MAP2K2", "MDFIC", "MRM2", "NUDT1", "PIAS4", "PIP5K1C", "PLIN3", "RASGRP1",
  "SCN2A", "SEMA6B", "TMEM106B", "TTYH3", "UBXN6", "VWDE", "ZBTB7A", "ZFAND2A"
) # other way around :(

display_portal_results$index <- factor(display_portal_results$index, levels = desired_order_genes)


# widen dataframe. there should be 1 row per gene (index), and one column per tissue (qtl_name)
# portal_results_wide <- display_portal_results |>
#   # round b_SMR
#   mutate(b_SMR = round(b_SMR, digits = 4)) |>
#   pivot_wider(id_cols = index,
#     names_from = qtl_name,
#     values_from = c(b_SMR))

#### Visualize SMR results from the analysis of various tissues (blood, brain, heart) ####
plot <- ggplot(data = display_portal_results, aes(x = nice_qtl_name, y = index, fill = p.fdr_SMR)) +
  geom_tile() +
  geom_text(aes(label = HEIDI.result), color = "white", size = 4) +
  # adjust x-axis label rotation
  scale_x_discrete(position = "top", guide = guide_axis(angle = 90)) +
  labs(fill = "P-value (FDR)") +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12, face = "italic"),
    legend.title = element_text(size = 12),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid = element_line(colour = "grey92"),
    panel.grid.minor = element_line(linewidth = rel(0.5)),
    strip.background = element_rect(fill = "grey85", colour = "grey20"),
    axis.ticks = element_blank()
  )



plot

ggsave(here("Results", "SMR-portal", "SMR-HEIDI-results.pdf"), plot = plot, width = 8, height = 12)

# Brain: Frontal Cortex, Putamen (Basal Ganglia), Hippocampus, Nucleus Accumbens (Basal Ganglia), Cortex, Cerebellum,
# Cerebellar Hemisphere, Spinal Cord (Cervical), Hypothalamus, Anterior Cingulate Cortex
# Heart: Atrial Appendage, Left Ventricle
# Artery: Aorta
# Adipose: Visceral Omentum


#### Visualize SMR results from the atherosclerotic plaque analysis ####
smr_results <- read.csv(file = here("Results", "smr_results.csv"))
smr_results

smr_plot <- smr_results |>
  mutate(z_SMR = abs(b_SMR / se_SMR)) |>
  mutate(p.fdr_SMR = p.adjust(as.numeric(p_SMR), method = "BH"))

ggplot(data = smr_plot, aes(Gene)) +
  geom_bar(aes(weight = z_SMR)) +
  coord_flip() + # Flips the x and y axes
  geom_vline(xintercept = 1.96) +
  theme_minimal()

nominal_plaque_plot <- ggplot(data = smr_plot, aes(x = Gene, y = z_SMR)) +
  geom_col(fill = "#132B43") +
  geom_hline(aes(yintercept = 1.96, linetype = "Nominal significance (p < 0.05)"),
    color = "red", size = 0.8
  ) + # Vertical line at x = 1.96
  coord_flip() + # Flip axes for horizontal bars
  theme_minimal() +
  labs(x = "Gene", y = "Z-score", linetype = "Threshold") +
  scale_linetype_manual(values = c("Nominal significance (p < 0.05)" = "dashed")) + # Customize legend
  theme(
    legend.position = "bottom",
    title = element_text(face = "bold"), # Move legend to bottom
    # axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12, face = "italic"),
    legend.title = element_text(size = 12),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid = element_line(colour = "grey92"),
    panel.grid.minor = element_line(linewidth = rel(0.5)),
    strip.background = element_rect(fill = "grey85", colour = "grey20"),
    axis.ticks = element_blank()
  )
nominal_plaque_plot

ggsave(nominal_plaque_plot, file = here("Results", "nominal_plaque_plot_manuscript.pdf"), width = 7, height = 8)
