#==============================================================#
# Dose-response and empirical Points-of-Departure (PoD)        #
# RNA modality: Cytotoxicity Index (per replicate)             #
# Morphology modality: morph_dist (per condition)              #
#==============================================================#

# This script continues from main script: CellPainting.R

library(drc)

filter <- dplyr::filter; select <- dplyr::select   # guard against masking

#--- 1. Concentration map ----------------------------------------
conc_map <- ImageData %>%
  dplyr::select(compound, dose_level, conc_value) %>%
  dplyr::filter(!compound %in% c("control","water","dmso1","dmso2")) %>%
  dplyr::distinct() %>%
  dplyr::mutate(conc_value = as.numeric(conc_value))

#--- 2. RNA-side long table (per replicate) ----------------------
rna_long <- countsMetadata %>%
  tibble::rownames_to_column("sample_id") %>%
  dplyr::filter(dose %in% c("Low","Mid","High")) %>%
  dplyr::mutate(
    compound   = tolower(stringr::str_remove(Label, "_(Low|Mid|High)$")),
    dose_level = tolower(as.character(dose))
  ) %>%
  dplyr::left_join(conc_map, by = c("compound","dose_level")) %>%
  dplyr::transmute(sample_id, compound, conc_value,
                   response = Cytotoxicity_Index_scaled)

ctrl_rna <- countsMetadata %>%
  dplyr::filter(dose %in% c("control","water","dmso1","dmso2")) %>%
  dplyr::pull(Cytotoxicity_Index_scaled)

rna_long <- bind_rows(
  rna_long,
  expand.grid(compound = unique(rna_long$compound),
              k = seq_along(ctrl_rna)) %>%
    dplyr::mutate(sample_id = paste0("ctrl_", k),
                  conc_value = 0,
                  response = ctrl_rna[k]) %>%
    dplyr::select(sample_id, compound, conc_value, response)
)

#--- 3. Morphology-side long table (per condition) ---------------
morph_long <- QuestionA %>%
  dplyr::filter(!tolower(compound) %in% c("control","water","dmso1","dmso2")) %>%
  dplyr::transmute(compound   = tolower(compound),
                   dose_level = tolower(dose_level),
                   response   = morph_dist) %>%
  dplyr::left_join(conc_map, by = c("compound","dose_level"))

ctrl_morph <- QuestionA %>%
  dplyr::filter(tolower(compound) %in% c("control","water","dmso1","dmso2")) %>%
  dplyr::pull(morph_dist)

morph_long <- bind_rows(
  morph_long %>% dplyr::transmute(compound, conc_value, response),
  expand.grid(compound = unique(morph_long$compound),
              k = seq_along(ctrl_morph)) %>%
    dplyr::mutate(conc_value = 0, response = ctrl_morph[k]) %>%
    dplyr::select(compound, conc_value, response)
)

#--- 4. Benchmark response: control mean + 1 SD ------------------
bmr_rna   <- mean(ctrl_rna,   na.rm = TRUE) + sd(ctrl_rna,   na.rm = TRUE)
bmr_morph <- mean(ctrl_morph, na.rm = TRUE) + sd(ctrl_morph, na.rm = TRUE)

cat("BMR (RNA):       ", round(bmr_rna,   3), "\n")
cat("BMR (Morphology):", round(bmr_morph, 3), "\n")

#--- 5. Empirical PoD: lowest tested conc crossing BMR -----------
empirical_pod <- function(long_df, bmr) {
  long_df %>%
    dplyr::filter(conc_value > 0) %>%
    dplyr::group_by(compound, conc_value) %>%
    dplyr::summarise(mean_resp = mean(response, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(compound, conc_value) %>%
    dplyr::group_by(compound) %>%
    dplyr::summarise(
      ePoD = {
        crossed <- which(mean_resp > bmr)
        if (length(crossed) == 0) NA_real_ else conc_value[min(crossed)]
      },
      max_resp = max(mean_resp, na.rm = TRUE),
      .groups = "drop"
    )
}

rna_ePod   <- empirical_pod(rna_long,   bmr_rna)   %>% dplyr::rename(tPoD = ePoD, max_t = max_resp)
morph_ePod <- empirical_pod(morph_long, bmr_morph) %>% dplyr::rename(mPoD = ePoD, max_m = max_resp)

pod_summary <- dplyr::full_join(rna_ePod, morph_ePod, by = "compound") %>%
  dplyr::mutate(
    detected_t = !is.na(tPoD),
    detected_m = !is.na(mPoD),
    agree      = detected_t == detected_m
  )
print(pod_summary)

#--- 6. Detection cross-tab --------------------------------------
detect_tab <- table(
  RNA        = ifelse(pod_summary$detected_t, "detected", "no effect"),
  Morphology = ifelse(pod_summary$detected_m, "detected", "no effect")
)
print(detect_tab)

#--- 7. Per-compound dose-monotonicity Spearman ------------------
mono_test <- function(long_df) {
  long_df %>%
    dplyr::filter(conc_value > 0) %>%
    dplyr::group_by(compound) %>%
    dplyr::summarise(
      rho = suppressWarnings(cor(rank(conc_value), rank(response),
                                 method = "spearman")),
      .groups = "drop"
    )
}

mono_t <- mono_test(rna_long)   %>% dplyr::rename(rho_t = rho)
mono_m <- mono_test(morph_long) %>% dplyr::rename(rho_m = rho)
mono   <- dplyr::full_join(mono_t, mono_m, by = "compound")

mono_concord <- cor.test(mono$rho_t, mono$rho_m, method = "spearman", exact = FALSE)
print(mono_concord)

#--- 8. Concordance figure ---------------------------------------
p_mono <- ggplot(mono, aes(rho_t, rho_m, label = compound)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey70") +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "grey70") +
  geom_point(size = 3, colour = "#1F77B4") +
  ggrepel::geom_text_repel(size = 3, max.overlaps = Inf) +
  scale_x_continuous(limits = c(-1, 1)) +
  scale_y_continuous(limits = c(-1, 1)) +
  labs(x = "Transcriptional dose-monotonicity (Spearman ρ)",
       y = "Morphological dose-monotonicity (Spearman ρ)",
       caption = sprintf("Cross-modality Spearman ρ = %.2f, P = %.2g",
                         mono_concord$estimate, mono_concord$p.value)) +
  theme_bw(base_size = 12) +
  theme(plot.caption = element_text(size = 10))
print(p_mono)

#--- 9. Save outputs ---------------------------------------------
write.csv(pod_summary, file.path(output, "Empirical_PoD_summary.csv"),     row.names = FALSE)
write.csv(mono,        file.path(output, "Per_compound_monotonicity.csv"), row.names = FALSE)
ggsave(file.path(output, "Dose_monotonicity_concordance.pdf"),
       p_mono, width = 6, height = 5, dpi = 1000, device = "pdf")

# Dose responses per replicate, per-compound showing both modalities side-by-side at every dose
dose_panel <- bind_rows(
  rna_long   %>% dplyr::mutate(modality = "Transcriptional (CI)"),
  morph_long %>% dplyr::mutate(modality = "Morphological (distance)")
) %>% dplyr::filter(conc_value > 0)

p_dose_grid <- ggplot(dose_panel,
                      aes(x = conc_value, y = response, colour = modality)) +
  geom_point(alpha = 0.85, size = 1.6) +
  stat_summary(fun = mean, geom = "line", linewidth = 0.6) +
  scale_x_log10() +
  scale_colour_manual(values = c("Transcriptional (CI)"     = "#D62728",
                                 "Morphological (distance)" = "#1F77B4")) +
  facet_wrap(~ compound, scales = "free", ncol = 4) +
  labs(x = expression("Concentration ("*mu*"M, log"[10]*" scale)"),
       y = "Standardised response", colour = NULL) +
  theme_bw(base_size = 10) +
  theme(legend.position = "top",
        strip.text = element_text(size = 9))

print(p_dose_grid)

ggsave(file.path(output, "Per_compound_dose_response_RNA_vs_Morph.pdf"),
       p_dose_grid, width = 12, height = 10, dpi = 1000, device = "pdf")
