# Host Transcriptomic Analysis of Active vs. Latent Tuberculosis (GSE37250)

---

## 1. Project Overview

This repository contains a host transcriptomics analysis of the GEO microarray dataset **GSE37250**.

**Biological Question:**  
Which genes and biological pathways are differentially expressed between **HIV-negative participants with active tuberculosis (Active TB)** and **HIV-negative participants with latent TB infection (LTBI)**?

**Key Methodological Feature:**  
By including *Region* as a covariate in the linear model, we adjust for geographic differences between Malawi and South Africa when estimating the association between disease state and gene expression.

---

## 2. Dataset Summary

| Detail | Information |
| :--- | :--- |
| **Dataset** | GSE37250 (GEO) |
| **Platform** | Illumina HumanHT-12 V4.0 expression beadchip (**GPL10558**) |
| **Cohort** | HIV-negative adults only |
| **Comparison** | Active TB (n=97) vs. Latent TB infection (n=83) |
| **Geographic Distribution** | Malawi (n=86) - South Africa (n=94) |
| **Total Samples** | 180 |
| **Probes Tested** | 44,051 (mapped to 31,425 unique genes) |

---

## 3. Methods

This analysis used the **processed expression VALUE measurements** provided by GEO. According to the original study metadata, the data had already undergone background subtraction and quantile normalization in Genome Studio. We deliberately avoided blindly reapplying normalization or deleting negative values.

| Step | Method / Tool |
| :--- | :--- |
| **Data Retrieval** | GEOparse (Python) |
| **Data Processing** | Because the GEO-provided background-subtracted values included negative measurements, a global positive offset (`-min(expr) + 1`) was applied before log₂ transformation. This transformation was used consistently across all samples. |
| **Quality Control** | PCA, correlation matrix, and IQR-based filtering (removed only 2 invariant probes). |
| **Probe Count Reconciliation** | Starting from 47,323 probes, 3,270 probes without annotated gene symbols were removed. A subsequent IQR-based filter removed 2 invariant probes, resulting in **44,051 probes** taken forward for differential expression testing. |
| **Differential Expression** | `limma` (R/Bioconductor) with **empirical Bayes moderation**. Performed at the **probe level**; gene symbols were used for downstream visualization and biological interpretation. |
| **Statistical Model** | `~ 0 + disease + region` - The model estimates the Active TB vs. LTBI effect *while adjusting for geographical region*. |
| **Multiple Testing** | Benjamini-Hochberg False Discovery Rate (FDR < 0.05). |
| **Functional Enrichment** | Significant DE probes were mapped to gene symbols. Duplicate gene symbols were collapsed, and genes were filtered using absolute log2FC > 0.5. The resulting upregulated and downregulated gene lists contained **854** and **523** unique genes, respectively and were submitted separately to Enrichr using ChEA 2022 and GO databases. |

---

## 4. Results

### 4.1 Differential Expression Results

| Metric | Count |
| :--- | :--- |
| **Probes Tested** | 44,051 |
| **Significant DE Probes (FDR < 0.05)** | **12,049** |
| **Upregulated DE Probes** | 5,656 |
| **Downregulated DE Probes** | 6,393 |

### 4.2 Top Differentially Expressed Features

| Upregulated in Active TB (Top 5) | Downregulated in Active TB (Top 5) |
| :--- | :--- |
| **FCGR1A** | **FCGBP** |
| **ANKRD22** | **VPREB3** |
| **FCGR1B** | **CD79A** |
| **SERPING1** | **CXCR5** |
| **C1QB** | **ID3** |

> *Full results are available in `DE_results_limma_region.csv`.*

### 4.3 Biological Interpretation

**Active TB is associated with a hyper-inflammatory innate immune transcriptional profile alongside reduced expression of adaptive immunity genes.**

- **Enriched biological processes and immune programs associated with increased expression:** Interferon-γ-associated responses, neutrophil degranulation, Fc-gamma receptor-mediated phagocytosis and complement-related functions.
- **Enriched biological processes and immune programs associated with decreased expression:** B-cell receptor signaling, T-cell differentiation, and B-cell-associated transcriptional programs involving regulators such as **PAX5** and **TCF3/E2A**.

This transcriptional pattern is consistent with the immunopathogenesis of progressive TB.

---

## 5. Functional Enrichment (Enrichr)

To interpret the biological meaning of the DE probes, we performed functional enrichment using **Enrichr** on the collapsed, filtered gene lists. The upregulated list contained genes with FDR < 0.05 and log2FC > 0.5, while the downregulated list contained genes with FDR < 0.05 and log2FC < -0.5. Probe-level results were mapped to gene symbols and duplicate symbols were collapsed before enrichment. This resulted in **854 upregulated** and **523 downregulated** unique genes submitted to Enrichr.

### 5.1 Key Findings

**Upregulated in Active TB (Innate Immune Surge):**
- **Transcription Factors (ChEA 2022):** IRF8, IRF1, STAT3, STAT5A, NFE2L2
- **Top GO Terms:**
  - Defense Response to Bacterium (p=1.85e-17)
  - Response to Type II Interferon (p=5.73e-16)
  - Pattern Recognition Receptor Signaling (p=1.11e-15)
  - Cellular Response to Lipopolysaccharide (p=9.41e-16)

**Downregulated in Active TB (Reduced Adaptive Immune-Associated Expression):**
- **Transcription Factors (ChEA 2022):** KDM2B, STAT6, FOXO1, E2A, UTX
- **Enriched functional categories:** B-cell receptor signaling, T-cell differentiation, and primary immunodeficiency-related pathways.

### 5.2 Interactive Results & Reproducibility

The exact input gene lists are provided in the repository (`upregulated_genes.csv` and `downregulated_genes.csv`) for full reproducibility.

- **Upregulated genes:** [Enrichr Results](https://maayanlab.cloud/Enrichr/enrich?dataset=742b6508637ca6c91253abda833ef2e1)
- **Downregulated genes:** [Enrichr Results](https://maayanlab.cloud/Enrichr/enrich?dataset=adeffd18869c48290c6705e57b1e8b12)

---
## 6. Repository Structure

```text
TB-Host-Transcriptomics/
├── README.md                         
│
├── DATA & METADATA
│   ├── metadata.csv                  
│   ├── probe_annotation.csv           
│   └── GPL10558.soft.gz            
│
├── RESULTS
│   ├── DE_results_limma_region.csv   
│   ├── upregulated_genes.csv         
│   ├── downregulated_genes.csv        
│   ├── volcano_plot.png              
│   └── heatmap.png                    
│
└── CODE
    └── limma_analysis.R             
```                
## 7. How to Reproduce This Analysis
Prerequisites: R, RStudio, and the following packages:

limma
GEOquery
ggplot2
pheatmap
tidyr
dplyr

Steps:

Load expression_log2_filtered.csv and metadata.csv.
Build the design matrix: model.matrix(~ 0 + disease + region, data = meta).
Fit the linear model using lmFit(), apply the contrast ActiveTB - LatentTB, and run eBayes().
Extract results with topTable() and correct for multiple testing (Benjamini-Hochberg).
Generate volcano and heatmap plots using the code provided in limma_analysis.R.

## 8. Limitations
This analysis uses a publicly available cross-sectional cohort and therefore does not establish causality. The cohort consists of participants from Malawi and South Africa, and geographic region was included as a covariate to account for observed transcriptomic structure. The analysis does not independently validate the identified gene signatures in an external cohort. Differential expression was performed on microarray measurements and may be affected by probe-level annotation and multiple probes mapping to individual genes. Enrichment analyses provide biological associations rather than direct evidence of transcription factor activity or mechanism. Candidate molecular features require independent validation before clinical application.

## 9. Conclusion
Region-adjusted transcriptomic analysis of HIV-negative participants from GSE37250 identified extensive differences in whole-blood gene expression between active TB and LTBI. Active TB was characterized by increased expression of innate inflammatory, interferon-associated, Fc receptor, complement, and myeloid-associated genes, alongside reduced expression of genes associated with B-cell and adaptive immune functions.

These findings support the presence of a distinct systemic host transcriptional state associated with active TB compared with LTBI. The identified genes and pathways provide candidate molecular features for further biomarker development. Independent cohort validation will be required to determine their diagnostic and biological utility.

