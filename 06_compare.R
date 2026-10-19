# Part 3: compare cell-level MAST with sample-level pseudobulk results.
library(dplyr)
library(ggplot2)
library(patchwork)
library(qs2)

result_dir <- "results/day2/de"
plot_dir <- "plots/day2"
de <- qs_read(file.path(result_dir, "05_de.qs2"))
mast_results <- de$mast
pbulk_results <- de$pseudobulk

######################################
# Resume Here
######################################
comparison <- dplyr::full_join(
  mast_results %>%
    dplyr::select(gene, mast_log2FC = avg_log2FC, mast_p = p_val_adj),
  pbulk_results %>%
    dplyr::select(gene, pbulk_log2FC = logFC, pbulk_FDR = FDR),
  by = "gene"
)
# utils::write.csv(comparison,
#                  file.path(result_dir, "teaching_only_method_comparison.csv"),
#                  row.names = FALSE)


volc1<-ggplot(comparison,aes(x = mast_log2FC,y = -log10(mast_p)))+geom_point()+theme_linedraw()
volc2<-ggplot(comparison,aes(x = pbulk_log2FC,y = -log10(pbulk_FDR)))+geom_point()+theme_linedraw()
volc1|volc2

plot(
  eulerr::euler(list(mast = mast_results$gene, pbulk = pbulk_results$gene)),
  quantities = TRUE
)

figure <- ggplot2::ggplot(comparison,
                          ggplot2::aes(x = mast_log2FC, y = pbulk_log2FC)) +
  ggplot2::geom_hline(yintercept = 0, color = "grey75") +
  ggplot2::geom_vline(xintercept = 0, color = "grey75") +
  ggplot2::geom_point(alpha = 0.35, size = 0.8) +
  ggplot2::labs(x = "MAST cell-level log2 fold change",
                y = "edgeR sample-level log2 fold change") +
  ggplot2::theme_classic()
# grDevices::png(file.path(plot_dir, "03_mast_vs_pseudobulk.png"),
#                width = 7, height = 6, units = "in", res = 160)
print(figure)

figure2 <- ggplot2::ggplot(comparison,
                          ggplot2::aes(x = -log10(mast_p), y = -log10(pbulk_FDR))) +
  ggplot2::geom_hline(yintercept = 0, color = "grey75") +
  ggplot2::geom_vline(xintercept = 0, color = "grey75") +
  ggplot2::geom_point(alpha = 0.35, size = 0.8) +
  ggplot2::labs(x = "MAST cell-level -log10(p)",
               y = "edgeR sample-level -log10(FDR)") +
  ggplot2::theme_classic()
# grDevices::png(file.path(plot_dir, "03_mast_vs_pseudobulk.png"),
#                width = 7, height = 6, units = "in", res = 160)
print(figure2)
# grDevices::dev.off()

message("Finished teaching-only comparison; do not interpret p-values as disease effects.")
