# Deposited code and data for 1/f E vs I glioma manuscript

This repository contains the analytic pipeline for the aperiodic 1/f E/I balance glioma manuscript, along with example input data. All scripts are functional and organized by analytical domain. All patient identifiers have been removed in accordance with peer review policy.

## Repository structure

### RestingState/
Scripts for preprocessing and analyzing resting-state ECoG data, computing the power spectral density (PSD) per channel and extracting the 1/f aperiodic slope. The two scripts form a sequential pipeline and are run in order.

- `RestingState_PSD_welch_method.py` — Loads a per-patient EEGLAB `.set` recording (via MNE), applies notch filters at 60/120/180/240 Hz to remove line noise and harmonics, and computes the PSD for every channel using Welch's method (3-second segments, 50% overlap). Per-channel frequency and power arrays are written to a `<ID>_output_psd/` directory as `.npy` files for downstream fitting.
- `RestingState_fooof_Computations.py` — Reads the per-channel frequency/power `.npy` arrays produced by the previous step, validates them (checks for NaN/inf and empty arrays), interpolates across the 120 Hz harmonic, and fits the aperiodic component with FOOOF/specparam (fixed aperiodic mode, 70–150 Hz). The aperiodic exponent (slope) and offset are extracted per electrode and saved to `<ID>_fooof_results.csv`.

**Example input data:** Resting-state ECoG is provided for 13 patients (see Releases: "Example Input Data_RS_ECoG"), allowing the pipeline to be run end to end. Files are provided per patient in EEGLAB format: `.set` (channel data and metadata) with the associated `.fdt` binary data arrays.

### snRNA/
Tools for integrating single-nucleus RNA sequencing data with electrophysiological metrics, covering Seurat-based preprocessing, cell-type annotation, differential expression, malignant-state scoring, and correlation with nearest-neighbor electrode slopes. Numbered notebooks are run in order; the `S*` series provides supplementary analyses and figures, and the `.R` files are helper scripts called during annotation.

- `0_Prepare_Main_Object.Rmd` — Builds the primary integrated Seurat object from raw snRNA-seq samples (QC, normalization, integration).
- `1-f_snRNA_analysis.Rmd` — Core analysis linking snRNA-seq-derived metrics to the 1/f aperiodic slope.
- `1_epilepsy_MapMyCells.Rmd` — Applies the Allen Institute MapMyCells classification to the epilepsy (normal-appearing) samples.
- `2_epilepsy_MapMyCells_Plots.Rmd` — Visualization of MapMyCells annotation results for the epilepsy samples.
- `3_DGE_Epilepsy.Rmd` — Differential gene expression analysis within the epilepsy cohort.
- `4_Neftel_states.Rmd` — Scores malignant cells along the Neftel et al. glioma cell-state axes (AC-, OPC-, NPC-, MES-like).
- `5_Comparative_DGE_Analysis.Rmd` — Comparative differential expression across conditions (glioma vs. normal-appearing).
- `ApplyMapMyCellsLabels.R` — Helper script that applies returned MapMyCells labels back onto the Seurat object.
- `ConvertSamplesToMapMyCellsInput.R` — Helper script that converts Seurat samples into the input format required by MapMyCells.
- `Nearest_Electrode_Avg.py` — Computes, for each nucleus/sample, the averaged 1/f slope of the nearest ECoG electrode(s) for snRNA–electrophysiology correlation.
- `S0_prepare_annotated_object.rmd` — Supplementary: prepares the annotated Seurat object.
- `S1_full_sample_UMAPs.rmd` — Supplementary: full-sample UMAP visualizations.
- `S3_sample_UMAPs_epilepsy.rmd` — Supplementary: per-sample UMAPs for the epilepsy samples.
- `S4_PCA.rmd` — Supplementary: PCA of the snRNA-seq data.

**Data access:** See `LINK_TO_DATA` in this folder for links to the example Seurat objects (glioma / normal-appearing and epilepsy), hosted on Box, which allow the scripts to be run and primary findings to be reproduced.

## Requirements
Code was written and tested using:
- Python (v3.13)
- R (v4.3.2)

## Notes
- All patient identifiers have been removed in accordance with peer review policy.
  
