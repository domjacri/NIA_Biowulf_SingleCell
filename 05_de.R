# Part 3: compare MAST cell-level DE with edgeR sample-level pseudobulk.
# These five unmatched samples cannot identify a 5xFAD disease effect.
# Run: Rscript scripts/05_de.R

library(Seurat)
library(MAST)
library(edgeR)
library(dplyr)
library(ggplot2)
library(qs2)
library(patchwork)

project_dir <- "."
result_dir <- file.path(project_dir, "results/day2/de")
plot_dir <- file.path(project_dir, "plots/day2")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

obj<-qs_read("./results/day2/annotation/04_annotated.qs2")
DefaultAssay(obj) <- "RNA"


micro_cells <- rownames(obj@meta.data)[!is.na(obj$anno_subclass) & obj$anno_subclass == "334 Microglia NN"]

micro <- base::subset(obj, cells = micro_cells)

sample_info <- micro@meta.data |>
  dplyr::count(sample_id, study, condition, name = "n_microglia")

print(sample_info)

# MAST fits cell-level tests: cells from a sample are not independent replicates.
# The resulting p-values are not valid evidence of a disease association here.
Idents(micro) <- "condition"
mast_results <- Seurat::FindMarkers(
  object = micro, ident.1 = "5xFAD", ident.2 = "reference",
  assay = "RNA", test.use = "MAST", latent.vars = c("nCount_RNA", "percent.mt"),
  logfc.threshold = 0, min.pct = 0.1
)
mast_results <- tibble::rownames_to_column(as.data.frame(mast_results), var = "gene")
# utils::write.csv(mast_results,
#                  file.path(result_dir, "teaching_only_MAST_microglia.csv"),
#                  row.names = FALSE)

# Aggregate raw counts per biological sample; never aggregate per condition only.
# The s prefix avoids Seurat renaming a group that starts with a digit.
micro$sample_key <- paste0("s", match(micro$sample_id, sample_info$sample_id))
sample_info$sample_key <- paste0("s", seq_len(nrow(sample_info)))
pbulk_counts <- Seurat::AggregateExpression(
  micro, assays = "RNA", group.by = "sample_key", return.seurat = FALSE
)[["RNA"]]
pbulk_info <- sample_info[match(colnames(pbulk_counts), sample_info$sample_key), ]

group <- factor(pbulk_info$condition, levels = c("reference", "5xFAD"))
design <- stats::model.matrix(~ group)
pbulk_info <- as.data.frame(pbulk_info)
rownames(pbulk_info) <- pbulk_info$sample_key
y <- edgeR::DGEList(counts = pbulk_counts, samples = as.data.frame(pbulk_info))
keep <- edgeR::filterByExpr(y, design = design)
summary(keep)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- edgeR::calcNormFactors(y)
y <- edgeR::estimateDisp(y, design)
fit <- edgeR::glmQLFit(y, design)
qlf <- edgeR::glmQLFTest(fit, coef = "group5xFAD")
pbulk_results <- edgeR::topTags(qlf, n = Inf)$table |>
  tibble::rownames_to_column(var = "gene")

# utils::write.csv(pbulk_results,
#                  file.path(result_dir, "teaching_only_edgeR_pseudobulk_microglia.csv"),
#                  row.names = FALSE)


# Save both methods so the comparison can start in a fresh R session.
qs_save(list(mast = mast_results, pseudobulk = pbulk_results,
             counts = pbulk_counts, samples = pbulk_info),
        file.path(result_dir, "05_de.qs2"))
