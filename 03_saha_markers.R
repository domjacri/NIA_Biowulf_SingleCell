# Part 2.5: begin SAHA annotation at resolution 0.1.
# Start here near the end of Day 1, or at the beginning of Day 2.

library(Seurat)
library(ggplot2)
library(dplyr)
library(qs2)
library(SAHA)
library(SAHAdata)

obj <- qs_read(file = "./results/day1/02_integrated.qs2")

Idents(obj) <- "RNA_snn_res.0.1"
saha_inputs<- Seurat2SAHA(obj)
qs_save(saha_inputs, file = "results/day2/annotation/03_inputs_res01.qs2")  

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

# Stop here if we have reached the end of Day 1.
qs_save(ann, file = "results/day2/annotation/03_markers_res01.qs2")
