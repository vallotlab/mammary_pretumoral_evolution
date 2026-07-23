# mammary_pretumoral_evolution
This folder contains the code to reproduce the analysis from *Laisné, Schwartz 2026* paper. 

## 0. Data availability
Input and annotation data required to run the scripts are deposited to Zenodo (**DOI: XXX - private until publication**). 

## 1. Script description 
The **scripts** folder contains the following files, each correpsonding to a separate analysis block:

```
.
├── 0_Utils.R
├── 1.1_MiniBulk_DRAG_PreProcess.Rmd
├── 1.2_MiniBulk_DRAG_Figures.Rmd
├── 2.1_MiniBulk_RNA_analysis.Rmd
├── 2.1_MiniBulk_RNA_Figures.Rmd
├── 2.2_MiniBulkRNA_DRAG_Figures.Rmd
├── 2.3_MiniBulk_RNA_Tumor_vs_Juxta.Rmd
├── 3.1_sc_pseudotime.Rmd
├── 3.2_sc_Annotation_MutantAlv.Rmd
├── 3.3_sc_DRAG_PreProcess.Rmd
├── 3.4_sc_DRAG_Analyses.Rmd
├── 3.5_sc_DRAG_Analyses_Figures.Rmd
├── 3.6_BulkDNA_DRAG_Figure.Rmd
├── 4_sc_Comparison_MutantWT.Rmd
├── 5.1_scMultiome_PreProcess_Annotation.Rmd
├── 5.2_scMultiome_DiffusionPseudotime.Rmd
├── 5.3_scMultiome_NMF.rmd
└── 5.4_scMultiome_Figures.Rmd
```
Expected output of each script can be found in the **notebooks** folder.
\
\
To reproduce the full analysis, the scripts should be run in designated order to genereate the necessary intermediate files.
To reproduce the final figures, only the **_Figures** scripts can be run, starting from the final data object loaded to Zenodo (**DOI: XXX**). 

## 2. Dependencies
The code was run in Rstudio, using R (4.5.1) and Python (v3.10.19).
\
\
Required R packages are listed in the beginnig of each script.
Python dependencies and software versions are listed in *snap_py310_env.yml*.

Before running the scripts, the Python evironment *snap_py310_env.yml* should be created and all packages should be installed.
The installation time should not exceed 60 min.

## 3. Runtime
On a standard laptop, the runtime for all of the scripts does not exceed 30 min (MacBook Pro M4, 128Gb RAM).

