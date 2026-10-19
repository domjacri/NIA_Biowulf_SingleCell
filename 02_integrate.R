# Part 2: integrate study batches and inspect markers.

library(Seurat)
library(ggplot2)
library(dplyr)
library(qs2)
library(patchwork)

project_dir <- "."
result_dir <- file.path(project_dir, "results/day1")
plot_dir <- file.path(project_dir, "plots/day1")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

obj <- qs_read(file.path(result_dir, "01_merged.qs2"))
DefaultAssay(obj) <- "RNA"

# The import/merge may already have created one counts layer per sample.
obj[["RNA"]] <- SeuratObject::JoinLayers(obj[["RNA"]])#purposefully put together
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$study) #or split back out
obj <- Seurat::NormalizeData(obj)
obj <- Seurat::FindVariableFeatures(obj, nfeatures = 2000)
obj <- Seurat::ScaleData(obj)
obj <- Seurat::RunPCA(obj, npcs = 30)
obj <- Seurat::RunUMAP(obj, reduction = "pca", dims = 1:20,
                       reduction.name = "umap.unintegrated")
#let's take a short look
DimPlot(obj,reduction = "umap.unintegrated",group.by = "study")
#what if we blew through this stop sign...
obj <- Seurat::FindNeighbors(obj, reduction = "pca", dims = 1:20)
obj <- Seurat::FindClusters(obj,resolution = 2)
DimPlot(obj,reduction = "umap.unintegrated",group.by = "seurat_clusters",label = T)+NoLegend()

#so what would happen if we then explictly deal with the batches (studies)
obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj, nfeatures = 2000)
obj <- ScaleData(obj, vars.to.regress = c("percent.mt"))#regressing something like percent.mt will stop the percent.mt from driving clustering
obj <- RunPCA(obj)
library(harmony)
obj <- RunHarmony(object = obj, reduction = "pca", group.by.vars = c("orig.ident","study"),
                         reduction.save = 'harmony', plot_convergence = T, lambda = NULL)
print(ElbowPlot(obj, reduction = "harmony", ndims = 50))

n_dims <- 10
obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:n_dims)
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:n_dims, reduction.name = "umap_harmony")

DimPlot(obj,reduction = "umap_harmony",split.by = "study",group.by = "study")

obj<-FindClusters(obj,resolution = c(0.1,0.5,1.0))


p1<-DimPlot(obj,reduction = "umap_harmony",group.by = "RNA_snn_res.0.1",label = T)+NoLegend()
p2<-DimPlot(obj,reduction = "umap_harmony",group.by = "RNA_snn_res.0.5",label = T)+NoLegend()
p3<-DimPlot(obj,reduction = "umap_harmony",group.by = "RNA_snn_res.1",label = T)+NoLegend()

p1|p2|p3


obj$RNA_snn_res.2<-NULL

library(clustree)

clustree(obj,prefix = "RNA_snn_res.")

FeaturePlot(obj,features = c("Plp1","Aqp4","Map2"),reduction = "umap_harmony",ncol = 3)

Idents(obj)<-"RNA_snn_res.0.1"
DotPlot(obj,features = c("Plp1","Aqp4","Map2"))

Idents(obj) <- "RNA_snn_res.0.1"
markers <- Seurat::FindAllMarkers(#it will take whatever the Active idents are so need to run the line before to make sure it is the cluster you want!
  obj, assay = "RNA", only.pos = TRUE,
  min.pct = 0.2, logfc.threshold = 0.25)
#if you get a warning message about "immunogenomics/presto", please kill the command with the stop sign by the console below and run: remotes::install_github("immunogenomics/presto")
#then retry!
head(markers)
qs_save(obj,file = "./results/day1/02_integrated.qs2")
