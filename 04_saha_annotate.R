# Part 2.5: resume SAHA, finish resolution 0.1, then repeat at resolution 0.5.

library(Seurat)
library(ggplot2)
library(dplyr)
library(qs2)
library(SAHA)
library(SAHAdata)

obj <- qs_read("results/day1/02_integrated.qs2")
saha_inputs <- qs_read("results/day2/annotation/03_inputs_res01.qs2")
ann <- qs_read("results/day2/annotation/03_markers_res01.qs2")
data("ABC_meta", package = "SAHAdata")

data("ISOCTX_AvgExp", package = "SAHAdata")
ann <- SAHA::Create_SAHA_object(
  query = saha_inputs$avgexp, db = ISOCTX_AvgExp, data_type = "AvgExp",existing = ann
)
ann=Initialize_MarkerFree(ann = ann)
ann=Downsample(ann,custom_ds = saha_inputs$varfeat)
ann=NormalizeDS(ann,assay_query = "RNA")
ann=CorrelateDS(ann)
ann=Create_MarkerFree_Viz(ann,meta = ABC_meta)
call_SAHA_plots(ann, plot_type = "Marker-free",data_type = "AvgExp")

auto=AutoAnnotate(ann,data_type = "Both")
auto

ann <- Summarize_Params(ann)  # Run the function
ann@params$summary

Generate_SAHA_Report(ann, auto = auto, output_file = "./results/day2/annotation/saha_res01.html")

glimpse(obj@meta.data)
glimpse(auto)
glimpse(ABC_meta)

# Match auto clusters to the resolution 0.1 clustering.
obj$anno_subclass <- auto$marker_free[
  match(as.character(obj$RNA_snn_res.0.1), auto$cluster)
]

obj$anno_class <- ABC_meta$class[
  match(obj$anno_subclass, ABC_meta$subclass_spa)
]

DimPlot(obj,reduction = "umap_harmony",group.by = "anno_class")
DimPlot(obj,reduction = "umap_harmony",group.by = "anno_subclass")

obj$anno_subclass_res01 <- obj$anno_subclass
obj$anno_class_res01 <- obj$anno_class

##################do all over again....
Idents(obj) <- "RNA_snn_res.0.5"
saha_inputs<- Seurat2SAHA(obj)
qs_save(saha_inputs, file = "results/day2/annotation/03_inputs_res05.qs2")  

# Allen ISOCTX is an annotation aid. Review label validity for each tissue.
data("ISOCTX_Markers", package = "SAHAdata")
data("ABC_meta", package = "SAHAdata")
ann <- SAHA::Create_SAHA_object(
  query = saha_inputs$markers, db = ISOCTX_Markers, data_type = "Markers"
)
ann <- SAHA::Initialize_Markers(ann)
ann <- SAHA::Run_Marker_Based(ann)
ann <- SAHA::Create_MarkerBased_Viz(ann, meta = ABC_meta)
auto <- SAHA::AutoAnnotate(ann, data_type = "Markers")
print(auto)

print(ann@results$marker_based$dotplot_all)

data("ISOCTX_AvgExp", package = "SAHAdata")
ann <- SAHA::Create_SAHA_object(
  query = saha_inputs$avgexp, db = ISOCTX_AvgExp, data_type = "AvgExp",existing = ann
)
ann=Initialize_MarkerFree(ann = ann)
ann=Downsample(ann,custom_ds = saha_inputs$varfeat)
ann=NormalizeDS(ann,assay_query = "RNA")
ann=CorrelateDS(ann)
ann=Create_MarkerFree_Viz(ann,meta = ABC_meta)
call_SAHA_plots(ann, plot_type = "Marker-free",data_type = "AvgExp")

auto=AutoAnnotate(ann,data_type = "Both")
auto

ann <- Summarize_Params(ann)  # Run the function
ann@params$summary

Generate_SAHA_Report(ann, auto = auto, output_file = "./results/day2/annotation/saha_res05.html")

glimpse(obj@meta.data)
glimpse(auto)
glimpse(ABC_meta)

# Match auto clusters to the resolution 0.5 clustering.
obj$anno_subclass <- auto$marker_free[
  match(as.character(obj$RNA_snn_res.0.5), auto$cluster)
]

obj$anno_class <- ABC_meta$class[
  match(obj$anno_subclass, ABC_meta$subclass_spa)
]

DimPlot(obj,reduction = "umap_harmony",group.by = "anno_class")
DimPlot(obj,reduction = "umap_harmony",group.by = "anno_subclass")

qs_save(obj,"./results/day2/annotation/04_annotated.qs2")
