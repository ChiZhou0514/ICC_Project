# Code Script for manuscript titled "Mature Tertiary Lymphoid Structures Fuel Stem-like CD8 T Cells to Boost Anti-Tumor Immunity in Intrahepatic Cholangiocarcinoma"

This repository is a collection of analysis scripts and notebooks used for genomic, single-cell, and spatial transcriptomics workflows.

**Important caveat:** most scripts are HPC-style pipelines with hard-coded absolute paths (commonly under `/home_2//` and `/NFS_home/...`). They are primarily intended to be run on the original environment where those paths, tools, and datasets exist.

## Contents

- [Directory Overview](#directory-overview)
- [Dependencies](#dependencies)
- [Detailed File Guide](#detailed-file-guide)
- [Portability Checklist](#portability-checklist)

## Directory Overview

- `Genomic analysis/`: WES-style mutation calling/post-processing, CNV calling, clonality, dN/dS, enrichment, and cohort comparison.
- `Single cell analysis/`: Seurat-based scRNA pipeline plus V(D)J (TCR/BCR) notebooks and metabolism scoring scripts.
- `Spatial Transcriptomics analysis/`: Spatial notebooks (public ST, DSP, 10X HD).

## Dependencies

### Python

- `numpy`, `pandas`, `tqdm`
- `scipy`, `seaborn`, `matplotlib`
- `skbio` (Shannon diversity in metrics script)
- `psutil` (Sentieon license server detection)

### R 

- Single-cell: `Seurat` [R1], `DoubletFinder` [R2], `monocle` [R3], `scMetabolism` [R4]
- Genomics: `maftools` [R5], `dndscv` [R6], `clonevol` [R7], `revolver` [R8], `ABSOLUTE` [R9]

### CLI tools

- `bcftools` [C1], `samtools` [C2]
- `picard` (LiftoverVcf, CollectWgsMetrics) [C3]
- `vcf2maf` [C4] + Ensembl VEP cache/data [C5]
- `CNVkit` (`cnvkit.py`) [C6]
- `msisensor-pro` [C7]
- `Sentieon` (license file + `licsrvr`) [C8]

## Detailed File Guide

Below is a per-file guide with **purpose**, **inputs**, **outputs**, **how to run**, and **notes**.

### `Genomic analysis/util.py`

- Purpose: shared helpers (cohort/outcome loading, TLS grouping join, patient/tissue enumeration, purity/ploidy lookup, multiprocessing wrapper).
- Inputs: cohort Excel/TSV references (hard-coded paths); TLS grouping workbook (`TLS_grouping.xlsx` or legacy name); local folder structure under `outcome/`.
- Outputs: returns pandas DataFrames and dictionaries; no fixed output file (unless called by other scripts).
- How to run: imported by almost all Python scripts in `Genomic analysis/`.
- Notes: includes column alias normalization (English + legacy Chinese column names) for robustness.

### `Genomic analysis/Call_Mutation/Mutect.py`

- Purpose: Sentieon-based WES pipeline for mapping + preprocessing + tumor-normal somatic calling; optional germline calling.
- Inputs:
  - FASTQs `A_R1.fastq.gz`, `A_R2.fastq.gz` (and optionally `B_*` lane).
  - Reference files and parameters from `RunMultiProcess` in `util.py`.
  - Matched normal BAM under `../NC/gatk/deduped.bam`.
- Outputs:
  - Preprocess: `gatk/deduped.bam`, `gatk/recal_data.table`, logs under `gatk/log.txt`.
  - Somatic: `gatk/output_filter.vcf.gz` (final filtered VCF).
  - Germline: `germline/germlinePASS.vcf`.
- How to run: execute as a script in the original HPC environment; entrypoints call `f_MappingAndProcess()`, `f_callMutation()`, `f_callMutation_G()`.
- Notes:
  - Starts Sentieon license server if missing.
  - Uses `shell=True`; review commands carefully when adapting.

### `Genomic analysis/Call_Mutation/VCFfilter.py`

- Purpose: VCF post-processing and VCF→MAF conversion; create pyclone-friendly SNV/indel VCF/MAF subsets.
- Inputs:
  - `output_filter.vcf.gz` produced by the caller.
  - LiftOver chain and reference fasta (for hg19 conversion path).
  - `vcf2maf.pl` + VEP cache/data.
- Outputs:
  - `lifted_overPASS.maf` (hg19 path) and/or `4pyclone_snpindels/output_pass2.maf` (hg38 path).
  - `4pyclone_snpindels/output_filterPro.vcf.gz` + index; `output_pass2.vcf`.
  - `4pyclone_snpindels/output_pass2absolute.maf` for ABSOLUTE (hg38 path).
- How to run: batch functions `f_vcfProcess()` (hg19+liftover) or `f_vcfProcess_hg38()` (hg38).
- Notes: filtering thresholds (AF/DP/AD) are project-specific; adjust for new cohorts.

### `Genomic analysis/Call_Mutation/mergeVCF.py`

- Purpose: merge region-level and patient-level MAFs; generate maftools clinical tables; add clonality labels from `clone_state9.tsv`.
- Inputs:
  - Per-region MAFs (e.g., `output_pass2.maf`, `output_pass3.maf`) and clone-state tables.
  - Cohort grouping from `util.getOutcome()` and TLS grouping via `util.getDataframe1()`.
- Outputs (examples):
  - Region-level: `analysisOutcome/maftools/regionLevel/clinical.tsv`, `clinical1.tsv`, `allSampleClonality.maf`, etc.
  - Patient-level: `analysisOutcome/maftools/patientLevel/short.maf`, `long.maf`, `allSample.maf`, plus `*Clonality.maf`.
- How to run: call the desired top-level function(s) in `__main__` according to the analysis you need.
- Notes: this script mixes multiple utilities; read function docstrings and run only the sections you need.

### `Genomic analysis/CNV_analysis/CNVkit.py`

- Purpose: orchestrate CNVkit CN calling and related downstream preparation (ABSOLUTE hook, GISTIC merges, CNV heatmaps).
- Inputs:
  - `gatk/deduped.bam` per sample.
  - CNVkit targets/antitargets + reference `.cnn` (hard-coded).
  - Purity/ploidy dictionaries from `ff_gettumorContent()`.
- Outputs:
  - CNVkit outputs under `cnv/` or `cnv1/` (coverage `.cnn`, `.cnr`, `.cns`, called `.cns`, exported `.seg`).
  - Optional: ABSOLUTE results (`absolute/`).
  - Optional: merged `.seg` files for GISTIC.
- How to run: enable the relevant calls in `__main__` (many lines are commented as “switches”).
- Notes: there are multiple CNV calling modes (`cnv` vs `cnv1`, segmetrics vs threshold); keep consistent with downstream scripts.

### `Genomic analysis/CNV_analysis/CNVkitGene.py`

- Purpose: gene-level CNV summaries, TLS/outcome group frequencies, merged CNV profiles, chi-square tests, maftools CN tables.
- Inputs:
  - CNVkit called tables such as `tumor.callsegmetrics.cns`.
  - TLS grouping workbook (English or legacy filename) with a sample ID column.
- Outputs:
  - `*_Status.tsv` per gene; merged CNV tables; `CNVfre.tsv`; maftools `cnTable.tsv` (patient-level and region-level).
- How to run: call `genCNVprofile()`, `mergeCNVprofile1()`, `f_getCN()`, `f_doTest()`, etc., depending on goal.
- Notes: `load_tls_grouping()` normalizes legacy column names; ensure your workbook contains the expected columns.

### `Genomic analysis/CNV_analysis/absolute.r`

- Purpose: run ABSOLUTE on CNV segments to infer purity/ploidy and absolute copy number.
- Inputs:
  - `cnv/tumor2absolute.seg`
  - `gatk/4pyclone_snpindels/output_pass2absolute.maf`
- Outputs: results under `absolute/` directory (ABSOLUTE output files).
- How to run: called by `CNVkit.py` (`Absolute(...)`) or run directly after setting the working directory.
- Notes: ABSOLUTE parameters are tuned for Illumina WES; adjust `min.ploidy/max.ploidy` if it fails.

### `Genomic analysis/microsatellite_analysis.py`

- Purpose: batch MSI calling using `msisensor-pro`.
- Inputs:
  - reference site list: `/home//software/msisensor/reference.site`
  - normal BAM `../../NC/gatk/deduped.bam` and tumor BAM `../gatk/deduped.bam`
- Outputs: `msi/msi.tsv` per sample.
- How to run: run `f_doMsi()` which uses multiprocessing to iterate patients.
- Notes: this uses `CNVkit.f_subprocess` for logging; output goes to `log.txt`.

### `Genomic analysis/Clone_analysis/pycloneMajorCopy.py`

- Purpose: prepare PyClone inputs using ref/alt counts and CNVkit minor/major CN; convert PyClone outputs to clonevol input; dispatch plotting.
- Inputs:
  - merged VCF: `analysisOutcome/SCpyclone_majorCopy/<patient>/merge.vcf`
  - CNVkit allele-specific CN calls: `cnv/tumor.callGermline.cns`
  - per-sample pyclone VCF/MAF in `gatk/4pyclone_snpindels/`
  - PyClone output tables under `pyclone_analysis/tables/`
- Outputs:
  - per-tissue PyClone tables (`PT.txt`, etc.)
  - `clonevol.txt` (clonevol-ready table)
  - clonevol plot artifacts produced by `clonevol.r`
- How to run: run the batch functions (`f_mergeVCF`, `f_prePyclone`, `f_preClonevol`, `f_plotPyclone`) as needed.
- Notes: several sections are project-specific (sample list overrides, parsing assumptions); validate on a small patient first.

### `Genomic analysis/Clone_analysis/clonevol.r`

- Purpose: visualize clonal evolution using clonevol from `clonevol.txt`.
- Inputs: `clonevol.txt` under the patient working directory.
- Outputs: clonevol plots (PDFs) written to the working directory (depending on clonevol defaults).
- How to run: set `setwd(...)` to a patient folder and run in R.
- Notes: uses bootstrap testing; can be slow for large variant sets.

### `Genomic analysis/Clone_analysis/REVOLVER.py`

- Purpose: convert per-patient clonevol tables into a REVOLVER cohort `input.tsv`.
- Inputs: `/home_2//analysis/pyclone_outcome_majorCopy/<patient>/clonevol.txt`.
- Outputs: `/home_2//analysis/revolver_outcome/input.tsv`.
- How to run: run `f_preFile()`; optionally run `analysisResult()` after REVOLVER produces `output.tsv`.
- Notes: cluster 1 is treated as clonal; CCFs are thresholded at 0.1 by default.

### `Genomic analysis/Clone_analysis/revolver.r`

- Purpose: run REVOLVER cohort model fitting, clustering, and plotting.
- Inputs: `input.tsv` produced by `REVOLVER.py`.
- Outputs: `output.tsv` plus multiple plot PDFs.
- How to run: set `setwd(...)` and run in R.
- Notes: requires the `revolver` R package and its dependencies.

### `Genomic analysis/dNds/dNdS.py`

- Purpose: export MAFs into `dndscv` TSV format for patient-level and region-level analyses.
- Inputs: `*Clonality.maf` files generated by maftools/merge steps plus clinical labels.
- Outputs: TSVs such as `short.tsv`, `long.tsv`, `NoTLS.tsv`, `ImmatureTLS.tsv`, etc.
- How to run: choose `MafOutcome1()` or `MafOutcome2()` in `__main__`.
- Notes: this is a preparatory step; the actual dN/dS runs are in the R scripts below.

### `Genomic analysis/dNds/dNdSOutCome.r`

- Purpose: run `dndscv` by outcome (long vs short), including clone/subclone splits.
- Inputs: patient-level TSVs produced by `dNdS.py`.
- Outputs: global dN/dS tables and gene-level significant tables under `GeneLevel/`.
- How to run: run in R after setting the correct working directory and `refdb`.
- Notes: q-value cutoffs are hard-coded; adjust if needed.

### `Genomic analysis/dNds/dNdSTLS.r`

- Purpose: run `dndscv` by TLS grouping (e.g., NoTLS vs Immature TLS), including clone/subclone splits.
- Inputs: region-level TSVs produced by `dNdS.py`.
- Outputs: global dN/dS summaries per group and clonality subset.
- How to run: run in R after setting the correct working directory and `refdb`.
- Notes: extend the script if you add more TLS categories.

### `Genomic analysis/maftools/maftools2.r`

- Purpose: maftools mutation landscape analysis for outcome grouping (short vs long).
- Inputs: `clinical1.tsv`, `cnTable.tsv`, `allSample.maf` (region level).
- Outputs: maftools plots (oncoplot, Ti/Tv, VAF plots, etc.) in the working directory.
- How to run: set `setwd(...)` and run in R.
- Notes: expects a project-specific `IMPACT1` field for filtering high-impact mutations.

### `Genomic analysis/maftools/maftools3.r`

- Purpose: maftools mutation landscape analysis for TLS grouping.
- Inputs: `clinical.tsv`, `cnTable.tsv`, `allSampleClonality.maf`.
- Outputs: maftools plots in the working directory.
- How to run: set `setwd(...)` and run in R.
- Notes: uses TLS labels from the clinical table; ensure your `clinical.tsv` matches the MAF sample IDs.

### `Genomic analysis/scripts_used_to_calculate_final_metrics.py`

- Purpose: aggregate multiple metrics (CNV, mutation, diversity, purity/ploidy, MSI) and produce comparison plots.
- Inputs: multiple upstream outputs (CNVkit `.cns`, MAFs, clonevol tables, MSI outputs, tumor content dicts).
- Outputs: PDFs and/or tables written to the configured output directories.
- How to run: this is typically executed as an analysis script after all upstream artifacts exist.
- Notes: hard-coded paths and cohort sizes are embedded in some functions; verify before reuse.

### `Genomic analysis/MutationEnrichment.py`

- Purpose: Hallmark gene-set mutation enrichment and statistical testing (LRT via MNLogit).
- Inputs: Hallmark gene set table, region-level analysis metadata (`Analysis1.tsv`), and MAF/clone-state tables.
- Outputs: `HallMarkCounts.tsv`, clone/subclone count tables, and enrichment result tables.
- How to run: run the `__main__` block or call specific functions based on your target label (TLS/outcome).
- Notes: normalizes legacy column names; p-value/FDR thresholds are project defaults.

### `Genomic analysis/cohort_compare/cohortCompare.py`

- Purpose: compare mutation frequencies for a curated gene list across multiple cohorts (our cohort vs CPTAC/TCGA/Zou).
- Inputs:
  - cohort-level MAFs (`ourCohort_patient.maf`, `TCGA.maf`, `NC.maf`, CPTAC file)
  - optional TCGA VCF folders if using the conversion helpers
- Outputs: `Gene_fre.tsv` and `Gene_fre.pdf`.
- How to run: run `calFre()` after preparing/placing the required cohort files.
- Notes: denominators are fixed cohort sizes; update them if your cohort sizes differ.

### `Single cell analysis/scRNAseq.r`

- Purpose: end-to-end scRNA pipeline in Seurat (load matrices → merge → doublet removal → integration/clustering → annotation → heatmaps → Monocle trajectory).
- Inputs: CellRanger `filtered_feature_bc_matrix` folders and multiple saved Seurat objects referenced by paths.
- Outputs: saved Seurat objects (`saveRDS`/`save`), plots (`ggsave`), heatmaps, Monocle objects.
- How to run: run top-to-bottom in an R session; confirm `setwd(...)` and all output directories.
- Notes: written as a sequential notebook-style pipeline; later sections depend heavily on earlier objects.

### `Single cell analysis/scMetabolism_metabolic_T.R`

- Purpose: compute metabolism scores (KEGG/REACTOME) for selected T cell subsets using `scMetabolism` (VISION).
- Inputs: saved integrated T cell Seurat object (`.Rdata`).
- Outputs: tab-delimited score matrices (e.g., `METABOLISM_score_KEGG_WithTLS.xls`).
- How to run: run in R; adjust subset logic and `ncores` for your environment.
- Notes: creates output directories if missing; heavy compute on large objects.

### `Single cell analysis/scMetabolism_metabolic_B.R`

- Purpose: compute metabolism scores (KEGG/REACTOME) for a B cell Seurat object.
- Inputs: saved integrated B cell Seurat object (`.Rdata`).
- Outputs: `METABOLISM_score_KEGG.xls`, `METABOLISM_score_REACTOME.xls`.
- How to run: run in R; adjust `ncores` as needed.
- Notes: assumes the loaded object contains `seurat.integrated.T_cell_sub` as in the original pipeline.

### `Single cell analysis/scTCR.ipynb`

- Purpose: scTCR clonotype analysis (diversity, expansion, overlap, distribution across samples/cell types).
- Inputs: Seurat metadata + TCR clonotype tables referenced by notebook paths.
- Outputs: plots saved by `ggsave(...)`, tables saved by `write.csv(...)`/`write.table(...)`.
- How to run: open in Jupyter and run sections in order (setup → summaries → overlap → subgroup analyses).
- Notes: an English overview cell is inserted at the top (`CODEX_NOTEBOOK_OVERVIEW_V1`) describing the intended run order.

### `Single cell analysis/scBCR.ipynb`

- Purpose: scBCR clonotype analysis (diversity, expansion, overlap; optional network/circos visualizations).
- Inputs: Seurat metadata + BCR clonotype tables referenced by notebook paths.
- Outputs: plots and exported tables.
- How to run: open in Jupyter; verify clonotype/barcode key columns before running.
- Notes: contains optional graph-based visualization code that may require additional packages.

### `Spatial Transcriptomics analysis/Public_ST_jupyter-notebook.ipynb`

- Purpose: spatial transcriptomics analysis on a public dataset (QC, clustering/visualization, ligand-receptor sections, ROC plots).
- Inputs: public ST objects and metadata referenced by notebook paths.
- Outputs: figures and exported tables depending on the notebook sections.
- How to run: open in Jupyter and run in order (import/QC → clustering → downstream).
- Notes: large objects may require substantial RAM.

### `Spatial Transcriptomics analysis/DSP_jupyter-notebook.ipynb`

- Purpose: DSP workflow with marker heatmaps and group comparisons across sample subsets.
- Inputs: DSP expression matrix + annotation tables referenced by notebook paths.
- Outputs: heatmaps and summary figures.
- How to run: open in Jupyter and run the relevant subset sections.
- Notes: heatmap scaling (`scale='row'` vs none) strongly affects interpretation; keep consistent.

### `Spatial Transcriptomics analysis/10XHD_jupyter-notebook.ipynb`

- Purpose: 10X HD spatial workflow: clustering/markers, region/subregion composition summaries, CellChat-based communication analysis blocks.
- Inputs: 10X HD objects, region/subregion annotation tables referenced by notebook paths.
- Outputs: composition CSVs and PDF figures, optional CellChat outputs.
- How to run: open in Jupyter; run heavy CellChat blocks last.
- Notes: output folders in `ggsave(...)`/`write.csv(...)` must exist or be created first.

## Portability Checklist

If you want to run this repo outside the original HPC environment:

1. Replace hard-coded `setwd(...)` and absolute file paths with your local paths.
2. Confirm the expected directory layout for `outcome/<patient>/WES/<tissue>/...`.
3. Ensure external tools are installed and on `PATH` (`bcftools`, `cnvkit.py`, `msisensor-pro`, etc.).
4. Verify cohort size denominators and hard-coded sample lists (they affect frequencies and plots).
5. Run one patient/sample end-to-end first to validate assumptions before batch processing.

## Citation

Wang P, Zhou C, Liu H, Wei Z, Gao Y, Wang Y, Meng F, Zhang Z, Zhou K, Rao K, Cao M, Guo W, Jin Q, Qiu S, Shi Y, Sun H, Gao Q, Zhou J, Hou Y, Peng DH, Fan J, Liu Q, Sun Y. Mature Tertiary Lymphoid Structures Fuel Stem-like CD8 T Cells to Boost Anti-Tumor Immunity in Intrahepatic Cholangiocarcinoma. Submit, 2026.

## Contacts

20310093@tongji.edu.cn or qiliu@tongji.edu.cn
