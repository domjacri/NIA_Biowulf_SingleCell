# Starting Your Single-Cell Analysis on Biowulf

**NIA BDMS | October 19–20, 2026**  
**Instructor:** Dominic J. Acri, PhD, Laboratory of Neurogenetics, NIA/IRP

A two-day, hands-on workshop covering the steps from raw single-cell sequencing data to annotated clusters and differential expression. Materials include daily R Markdown documents and corresponding R scripts.

## Workshop topics

**Day 1**

- Run Cell Ranger on human PBMC FASTQs and explore its outputs.
- Import five mouse brain samples into Seurat and inspect quality metrics.
- Explore batch effects, Harmony integration, clustering resolutions, and markers.
- Begin cluster annotation using SAHA.

**Day 2**

- Complete SAHA annotation and compare clustering resolutions.
- Compare cell-level differential expression using MAST with sample-level pseudobulk analysis using edgeR.
- Discuss advanced single-cell analyses.

## Getting started

An NIH Biowulf account is required. Launch an RStudio session with **100 GB RAM**, as described in the slides.

In the RStudio Terminal, navigate to your preferred directory and run:

```bash
cp /data/CARD_datasets/acridj/BDMS_download.tar ./
tar -xf BDMS_download.tar
```

Open `BDMS_sc_2026_student/BDMS_sc_2026_student.Rproj`. Run `scripts/0_setup.R` once to install the required R packages.

## Using the materials

Follow `day1.Rmd` and `day2.Rmd`, or run the corresponding numbered scripts in `scripts/`. All paths are relative to the project folder, so the folder can be placed anywhere. Intermediate objects are saved as `.qs2` checkpoints for restarting between sections.

The course tarball includes input data and completed Cell Ranger output under `results/day1/cellranger_backup/`. Students’ Cell Ranger runs write to `results/day1/cellranger_run/`.

## Teaching datasets

The Cell Ranger demonstration uses human PBMC data. The Seurat exercises combine five mouse brain samples from 10x Genomics, Kersey 2026, and Karahan 2025.

These unmatched mouse samples support a comparison of analysis methods. Study and condition are confounded, so differential-expression results should not be interpreted as evidence of a 5xFAD disease effect.
