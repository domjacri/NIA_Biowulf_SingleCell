# Part 2: import the five precomputed 10x filtered matrices and QC.
# Run: Rscript scripts/01_qc.R

library(Seurat)
library(ggplot2)
library(dplyr)
library(qs2)

project_dir <- "."
data_dir <- file.path(project_dir, "data")
output_dir <- file.path(project_dir, "results/day1")
plot_dir <- file.path(project_dir, "plots/day1")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)


####################
#the simplest load in possible
####################
test <- Seurat::Read10X(data.dir  = "./data/10x_mouse/filtered_feature_bc_matrix/")
test
obj_test <- CreateSeuratObject(test,assay = "RNA",project = "test_demo",min.cells = 3,min.features = 0)
obj_test
obj_test[["percent.mt"]] <- Seurat::PercentageFeatureSet(
  obj_test, pattern = "^mt-"
)
colnames(obj_test@meta.data)
VlnPlot(obj_test,features = c("nCount_RNA","nFeature_RNA","percent.mt"))

####################
#let's do a little more reproducibly
####################
# "Reference" means non-5xFAD source material, not a matched WT control.
samples <- tibble::tibble(
  sample_id = c("10x_mouse", "GSM8824189_S_1", "GSM8824190_S_2",
                "A1331_Idol-Camk2a-WT-5xFAD", "A2022_Idol-Camk2a-WT-5xFAD"),
  study = c("10x_mouse", "Kersey_2026", "Kersey_2026",
            "Karahan_2025", "Karahan_2025"),
  condition = c("reference", "reference", "reference", "5xFAD", "5xFAD"),
  archive = c(
    file.path(data_dir, "10x_mouse/5k_mouse_brain_CNIK_3pv3_filtered_feature_bc_matrix.tar.gz"),
    file.path(data_dir, "Kersey_2026/GSM8824189_S_1_filtered_feature_bc_matrix.tar.gz"),
    file.path(data_dir, "Kersey_2026/GSM8824190_S_2_filtered_feature_bc_matrix.tar.gz"),
    file.path(data_dir, "Karahan_2025/A1331_Idol-Camk2a-WT-5xFAD_filtered_feature_bc_matrix.tar.gz"),
    file.path(data_dir, "Karahan_2025/A2022_Idol-Camk2a-WT-5xFAD_filtered_feature_bc_matrix.tar.gz")
  )
)


seurat_samples <- vector("list", nrow(samples))
names(seurat_samples) <- samples$sample_id

for (i in seq_len(nrow(samples))) {
  sample_id <- samples$sample_id[i]
  unpack_dir <- file.path(data_dir, "unpacked", sample_id)
  dir.create(unpack_dir, recursive = TRUE, showWarnings = FALSE)

  utils::untar(samples$archive[i], exdir = unpack_dir)
  matrix_files <- list.files(unpack_dir, pattern = "^matrix[.]mtx([.]gz)?$",
                             recursive = TRUE, full.names = TRUE)
  counts <- Seurat::Read10X(data.dir = dirname(matrix_files))

  seurat_samples[[i]] <- Seurat::CreateSeuratObject(
    counts = counts, project = sample_id, min.cells = 3, min.features = 0
  )
  seurat_samples[[i]]$sample_id <- sample_id
  seurat_samples[[i]]$study <- samples$study[i]
  seurat_samples[[i]]$condition <- samples$condition[i]
  seurat_samples[[i]][["percent.mt"]] <- Seurat::PercentageFeatureSet(
    seurat_samples[[i]], pattern = "^mt-"
  )
  message(sample_id, ": ", ncol(seurat_samples[[i]]), " nuclei before QC")
}

raw_qc <- dplyr::bind_rows(lapply(seq_along(seurat_samples), function(i) {
  seurat_samples[[i]]@meta.data |>
    dplyr::mutate(sample_id = samples$sample_id[i]) |>
    dplyr::select(sample_id, nCount_RNA, nFeature_RNA, percent.mt)
}))
# utils::write.csv(raw_qc, file.path(output_dir, "01_nucleus_qc_before_filter.csv"), row.names = FALSE)

qc_plot <- ggplot2::ggplot(raw_qc, ggplot2::aes(x = sample_id, y = nFeature_RNA)) +
  ggplot2::geom_violin(fill = "#81b7be", scale = "width") +
  ggplot2::coord_flip() +
  ggplot2::labs(x = NULL, y = "Genes per nucleus", title = "Inspect QC by sample") +
  ggplot2::theme_classic()
# grDevices::png(file.path(plot_dir, "01_genes_per_nucleus.png"),
#                width = 9, height = 5, units = "in", res = 160)
print(qc_plot)
# grDevices::dev.off()

for (i in seq_along(seurat_samples)) {
  # Teaching thresholds; inspect sample-specific QC before keeping them.
  keep_cells <- colnames(seurat_samples[[i]])[
    seurat_samples[[i]]$nFeature_RNA >= 200 &
      seurat_samples[[i]]$percent.mt <= 5
  ]
  seurat_samples[[i]] <- base::subset(seurat_samples[[i]], cells = keep_cells)
}

cell_counts <- tibble::tibble(
  sample_id = samples$sample_id,
  study = samples$study,
  condition = samples$condition,
  n_cells_after_qc = vapply(seurat_samples, ncol, numeric(1))
)
# utils::write.csv(cell_counts, file.path(output_dir, "01_sample_summary.csv"), row.names = FALSE)
print(cell_counts)

obj <- base::merge(seurat_samples[[1]], y = seurat_samples[-1],
                   add.cell.ids = samples$sample_id, project = "BDMS_sc_2026")
DefaultAssay(obj) <- "RNA"
# Tip: qs2 is much quicker for saving/loading these objects: https://github.com/qsbase/qs2
qs_save(obj, file.path(output_dir, "01_merged.qs2"))

# Loading back in at the start of the next script:
# obj <- qs_read("results/day1/01_merged.qs2")
