#==================================================================#
# Translational benchmarking: in vitro PoD vs Daston in vivo       #
# Continues from main script: CellPainting.R + PerCompound-Monotonicity.R
# Requires in workspace: rna_ePod, morph_ePod, output                 #
# Requires file: in_vivo_benchmark.csv                                #
#==================================================================#

library(pROC)
library(yardstick)
library(mltools)
library(caret)

filter <- dplyr::filter
select <- dplyr::select

#--- 1. Load in vivo benchmark and harmonise ----------------------
in_vivo <- read.csv("in_vivo_benchmark.csv", stringsAsFactors = FALSE) %>%
  dplyr::mutate(compound = tolower(stringr::str_replace_all(compound, " ", "_")))

#--- 2. Build the benchmark table ---------------------------------
benchmark <- rna_ePod %>%
  dplyr::full_join(morph_ePod, by = "compound") %>%
  dplyr::full_join(in_vivo,    by = "compound") %>%
  dplyr::mutate(
    truth = factor(
      ifelse(teratogen_class_binary == "positive", "teratogen", "non-teratogen"),
      levels = c("non-teratogen", "teratogen")
    )
  )
print(benchmark)

#--- 3. ROC, AUC and DeLong test ----------------------------------
roc_t <- pROC::roc(benchmark$truth, -log10(benchmark$tPoD), quiet = TRUE)
roc_m <- pROC::roc(benchmark$truth, -log10(benchmark$mPoD), quiet = TRUE)
auc_compare <- pROC::roc.test(roc_t, roc_m, method = "delong")

#--- 4. Threshold sweep for classification metrics ----------------
sweep_threshold <- function(score, truth) {
  ok    <- !is.na(score) & !is.na(truth)
  score <- score[ok]; truth <- truth[ok]
  thr   <- quantile(score, probs = seq(0.05, 0.95, by = 0.05), na.rm = TRUE)
  lapply(thr, function(t) {
    pred <- factor(ifelse(score <= t, "teratogen", "non-teratogen"),
                   levels = c("non-teratogen", "teratogen"))
    cm   <- caret::confusionMatrix(pred, truth, positive = "teratogen")
    data.frame(
      threshold = t,
      bacc = cm$byClass["Balanced Accuracy"],
      sens = cm$byClass["Sensitivity"],
      spec = cm$byClass["Specificity"],
      mcc  = mltools::mcc(actuals = truth, preds = pred)
    )
  }) %>% bind_rows()
}

best_t <- sweep_threshold(log10(benchmark$tPoD), benchmark$truth) %>%
  dplyr::slice_max(bacc, n = 1, with_ties = FALSE)
best_m <- sweep_threshold(log10(benchmark$mPoD), benchmark$truth) %>%
  dplyr::slice_max(bacc, n = 1, with_ties = FALSE)

cat("\nAUC tPoD =", round(pROC::auc(roc_t), 3),
    " | AUC mPoD =", round(pROC::auc(roc_m), 3), "\n")
print(auc_compare)
cat("\nBest tPoD metrics:\n"); print(best_t)
cat("\nBest mPoD metrics:\n"); print(best_m)

#--- 5. Daston Cmax concordance -----------------------------------
daston_df <- benchmark %>% dplyr::filter(!is.na(daston_pos_Cmax_uM))

cor_pearson_daston_t  <- cor.test(log10(daston_df$tPoD),
                                  log10(daston_df$daston_pos_Cmax_uM),
                                  method = "pearson")
cor_spearman_daston_t <- cor.test(log10(daston_df$tPoD),
                                  log10(daston_df$daston_pos_Cmax_uM),
                                  method = "spearman", exact = FALSE)
cor_pearson_daston_m  <- cor.test(log10(daston_df$mPoD),
                                  log10(daston_df$daston_pos_Cmax_uM),
                                  method = "pearson")
cor_spearman_daston_m <- cor.test(log10(daston_df$mPoD),
                                  log10(daston_df$daston_pos_Cmax_uM),
                                  method = "spearman", exact = FALSE)

cat("\n--- Daston Cmax concordance ---\n")
cat("tPoD vs Daston: Pearson r =", round(cor_pearson_daston_t$estimate, 3),
    " (P=", signif(cor_pearson_daston_t$p.value, 2), ")",
    " | Spearman rho =", round(cor_spearman_daston_t$estimate, 3),
    " (P=", signif(cor_spearman_daston_t$p.value, 2), ")\n")
cat("mPoD vs Daston: Pearson r =", round(cor_pearson_daston_m$estimate, 3),
    " (P=", signif(cor_pearson_daston_m$p.value, 2), ")",
    " | Spearman rho =", round(cor_spearman_daston_m$estimate, 3),
    " (P=", signif(cor_spearman_daston_m$p.value, 2), ")\n")

#--- 6. Figure A: ROC overlay -------------------------------------
roc_df <- bind_rows(
  data.frame(spec = roc_t$specificities, sens = roc_t$sensitivities,
             modality = sprintf("tPoD (AUC = %.2f)", pROC::auc(roc_t))),
  data.frame(spec = roc_m$specificities, sens = roc_m$sensitivities,
             modality = sprintf("mPoD (AUC = %.2f)", pROC::auc(roc_m)))
)

p_roc <- ggplot(roc_df, aes(1 - spec, sens, colour = modality)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
  geom_path(linewidth = 0.9) +
  scale_colour_manual(values = c("#1F77B4", "#D62728")) +
  labs(x = "False Positive Rate", y = "True Positive Rate", colour = NULL,
       caption = sprintf("DeLong P = %.3g", auc_compare$p.value)) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")
print(p_roc)

#--- 7. Figure B: in vitro PoD vs Daston Cmax (log-log) -----------
plot_df <- daston_df %>%
  tidyr::pivot_longer(cols = c(tPoD, mPoD),
                      names_to = "modality", values_to = "in_vitro_PoD") %>%
  dplyr::filter(!is.na(in_vitro_PoD))

p_daston <- ggplot(plot_df,
                   aes(daston_pos_Cmax_uM, in_vitro_PoD,
                       colour = modality, label = compound)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(size = 3) +
  ggrepel::geom_text_repel(size = 3, max.overlaps = Inf, show.legend = FALSE) +
  scale_x_log10() + scale_y_log10() +
  scale_colour_manual(values = c("tPoD" = "#D62728", "mPoD" = "#1F77B4")) +
  labs(x = expression("Daston positive Cmax ("*mu*"M, log"[10]*")"),
       y = expression("in vitro PoD ("*mu*"M, log"[10]*")"),
       colour = NULL,
       caption = sprintf("tPoD: r = %.2f (P=%.2g), rho = %.2f (P=%.2g)\nmPoD: r = %.2f (P=%.2g), rho = %.2f (P=%.2g)",
                         cor_pearson_daston_t$estimate, cor_pearson_daston_t$p.value,
                         cor_spearman_daston_t$estimate, cor_spearman_daston_t$p.value,
                         cor_pearson_daston_m$estimate, cor_pearson_daston_m$p.value,
                         cor_spearman_daston_m$estimate, cor_spearman_daston_m$p.value)) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")
print(p_daston)

#--- 8. Figure C: classification metrics --------------------------
metric_df <- bind_rows(
  best_t %>% dplyr::mutate(modality = "tPoD"),
  best_m %>% dplyr::mutate(modality = "mPoD")
) %>% tidyr::pivot_longer(c(bacc, sens, spec, mcc),
                          names_to = "metric", values_to = "value")

p_metrics <- ggplot(metric_df, aes(metric, value, fill = modality)) +
  geom_col(position = "dodge", width = 0.65) +
  geom_text(aes(label = sprintf("%.2f", value)),
            position = position_dodge(width = 0.65), vjust = -0.3, size = 3.2) +
  scale_fill_manual(values = c("tPoD" = "#D62728", "mPoD" = "#1F77B4")) +
  scale_y_continuous(limits = c(-0.1, 1.05), breaks = seq(0, 1, 0.2)) +
  labs(x = NULL, y = "Score", fill = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")
print(p_metrics)

#--- 9. Save outputs ----------------------------------------------
ggsave(file.path(output, "Translational_ROC.pdf"),
       p_roc, width = 5, height = 5, dpi = 1000, device = "pdf")
ggsave(file.path(output, "Translational_PoD_vs_Daston.pdf"),
       p_daston, width = 6.5, height = 5.5, dpi = 1000, device = "pdf")
ggsave(file.path(output, "Translational_Metrics.pdf"),
       p_metrics, width = 5, height = 4, dpi = 1000, device = "pdf")

write.csv(benchmark, file.path(output, "Translational_benchmark_table.csv"),
          row.names = FALSE)


# Re-theme to keep the panels visually consistent
panel_theme <- theme_bw(base_size = 11) +
  theme(plot.title       = element_text(face = "bold", size = 12),
        plot.caption     = element_text(size  = 9, hjust = 0),
        legend.position  = "bottom",
        legend.title     = element_blank(),
        panel.grid.minor = element_blank())

p_roc_panel     <- p_roc     + panel_theme
p_metrics_panel <- p_metrics + panel_theme
p_daston_panel  <- p_daston  + panel_theme

Figure4 <- (p_roc_panel | p_metrics_panel) / p_daston_panel +
  plot_annotation(tag_levels = "A") +
  plot_layout(heights = c(1, 1.15)) &
  theme(plot.tag          = element_text(size = 13, face = "bold"),
        plot.tag.position = c(0, 1))

print(Figure4)

ggsave(file.path(output, "Figure_4.pdf"),
       Figure4, width = 12, height = 11, dpi = 1000, device = "pdf")
ggsave(file.path(output, "Figure_4.png"),
       Figure4, width = 12, height = 11, dpi = 600, device = "png")