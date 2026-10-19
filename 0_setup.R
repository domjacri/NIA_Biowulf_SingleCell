# Run once in the R console before the workshop, using the course RStudio session.
# Open BDMS_sc_2026_student.Rproj first; all paths start in this folder.

install.packages(c(
  "Seurat", "ggplot2", "dplyr", "tibble", "qs2", "harmony", "clustree",
  "patchwork", "eulerr", "remotes", "BiocManager", "rmarkdown", "knitr"
))
BiocManager::install(c("MAST", "edgeR", "ComplexHeatmap"), update = FALSE, ask = FALSE)
remotes::install_github("immunogenomics/presto", upgrade = "never")
remotes::install_github("neurogenetics/SAHA", upgrade = "never")
remotes::install_github("neurogenetics/SAHAdata", upgrade = "never")
