# Spatial Transcriptomics Normalization & Scaling

## Introduction: From Quality Control to Expression Rescaling

During quality control, we established how spatial technical noise—such as tissue folds, air bubbles, staining artifacts, zero-count segments, and low-quality capture units—is identified and filtered using spatially-aware metrics (e.g., `SpotSweeper` local $Z$-scores and multi-factorial morphological gates).

The next crucial analytical step transforms raw transcript count matrices into statistically comparable expression values. By converting filtered raw counts (`counts`) into log-normalized expression matrices (`logcounts`) and centered/scaled features, we remove technical sampling variation while preserving genuine spatial biology prior to downstream dimensionality reduction (PCA/UMAP), spatially variable gene (SVG) detection, and spatial domain clustering.

## The Core Rationale: Why Normalize and Scale?

### The Fundamental Problem: Raw Counts Are Not Comparable

In both sequencing-based (10x Visium and Visium HD) and imaging-based (10x Xenium) spatial transcriptomics, raw UMI count totals vary substantially between physical capture units (spots, bins, or segmented cells). Two spatial capture units representing identical biological cell states can exhibit a ≥ 2× difference in total detected UMIs:

```text
=====================================================================
                RAW COUNTS ARE NOT DIRECTLY COMPARABLE               
=====================================================================

  Spot 1 (Deep Capture)                    Spot 2 (Shallow Capture) 
┌───────────────────────┐                ┌──────────────────────────┐
│ Gene A: 20 UMIs       │                │ Gene A: 10 UMIs          │
│ Gene B: 12 UMIs       │  Same Relative │ Gene B: 6 UMIs           │
│ Gene C: 6 UMIs        │ <── Biology ─> │ Gene C: 3 UMIs           │
│ Total : 38 UMIs       │                │ Total : 19 UMIs          │
└───────────────────────┘                └──────────────────────────┘
            ▲                                         ▲
            └─────────────────┬───────────────────────┘
        Technical Depth/Capture Efficiency Difference (≈2×)
=====================================================================
```

These raw count differences are driven primarily by three factors:

1. **Variable Capture Efficiency:** Local variations in spatial array permeabilization, probe hybridization/ligation, or optical detection efficiency across the tissue section.
2. **Sequencing Depth Variation:** Random sampling noise and variable read distribution during flowcell sequencing.
3. **Physical Cell Size & Density:** Natural differences in cellular volume, tangential sectioning plane, and nuclear packing density within spatial capture bins or segmented polygons.

### Normalization vs. Scaling: The Golden Rule

A common analytical pitfall is conflating **normalization** with **scaling**. They perform distinct mathematical operations across orthogonal matrix axes and serve different downstream purposes:

| Operation | Target / Axis | Primary Goal | Mathematical Action | Typical Downstream Use |
| :--- | :--- | :--- | :--- | :--- |
| **Normalization** | **Observations** (Cells / Spots / Bins) | Makes capture units comparable across the tissue slide. | Divides raw counts by unit-specific size factors ($s_i$) and applies a variance-stabilizing transformation (`log1p`). | Differential expression, marker visualization, spatial trend plotting. |
| **Scaling** | **Features** (Genes / Transcripts) | Makes genes comparable to one another across capture units. | Centers gene expression to mean $= 0$ and rescales variance to $\text{SD} = 1$ ($Z$-score transformation). | Highly Variable Gene (HVG) weighting, PCA, UMAP, $k$-means/Leiden clustering, heatmaps. |

> **The Golden Rule of Data Processing:** We always scale **after** normalization—never scale raw count matrices directly. Scaling unnormalized raw counts propagates technical library-depth artifacts directly into downstream principal components and spatial cluster assignments.

## Standard Log-Normalization Framework (`scuttle::logNormCounts`)

### Pipeline Workflow

The foundational method for single-cell and spatial transcriptomics normalization is **size-factor log-normalization**, implemented in Bioconductor via `scuttle::logNormCounts()` (and C+\+-accelerated in `scrapper::normalizeRnaCounts.se()`):

```text
==============================================================================
         STANDARD LOG-NORMALIZATION WORKFLOW (scuttle::logNormCounts)
==============================================================================
┌───────────────────┐     ┌───────────────────┐     ┌───────────────────┐     
│  Raw Counts       │     │ Calculate Size    │     │ Divide Counts     │     
│  Matrix           │ ──► │ Factors (s_i)     │ ──► │ by Size Factors   │ ──► 
│  counts(sfe)      │     │ (Mean = 1.0)      │     │ (X_ij / s_i)      │     
└───────────────────┘     └───────────────────┘     └───────────────────┘     
           ┌───────────────────┐     ┌───────────────────┐
           │ Apply log1p       │     │ Output logcounts  │
       ... │ Transform         │ ──► │ Assay             │
           │ log_e(X_ij/s_i +1)│     │ logcounts(sfe)    │
           └───────────────────┘     └───────────────────┘
==============================================================================
```

### Mathematical Formulation

Log-normalization computes the normalized expression value $Y_{ij}$ for gene $j$ in spatial observation $i$ using the transformation:

$$
Y_{ij} = \log_e \left( \frac{X_{ij}}{s_i} + 1 \right) = \log_e \left( \frac{X_{ij}}{\sigma_i \sum_{j} X_{ij}} + 1 \right)
$$

where:

* $X_{ij}$ is the raw UMI or transcript count for gene $j$ in observation $i$.
* $s_i$ is the observation-specific **size factor**, proportional to the total library depth $\sum_{j} X_{ij}$ and scaled by a constant $\sigma_i$ so that the mean size factor across all valid observations equals unity ($\text{mean}(s_i) = 1$).
* $+1$ is the **pseudocount**, ensuring that zero counts map cleanly to $\log_e(1) = 0$ and stabilizing variance across low-abundance transcripts.

### Statistical Properties & Diagnostic Inspection of Size Factors

Size factors capture technical library depth alongside physical capture volume. In well-curated spatial datasets, size factor distributions exhibit three characteristic properties:

* **Mean-Centered:** By construction, library size factors are scaled so that $\text{mean}(s_i) = 1.0$.
* **Right-Skewed Distribution:** Most capture units cluster below the mean near the median (e.g., median ≈ 0.77–0.78), whereas a long upper tail of large, transcriptionally active cells (such as proliferating crypt enterocytes or invasive adenocarcinoma glands) or dense multi-cell bins exhibits size factors reaching 8–10+.
* **Spatial Diagnostic Verification:** We routinely inspect size factors numerically via `summary(sizeFactors(sfe))` and project them into physical space using `Voyager::plotSpatialFeature()` (for `sfe` cell polygons) or `ggspavis::plotCoords()` (for `spe` binned grids) to verify that variations mirror underlying tissue morphology rather than technical edge defects.

```text
Size Factor Distribution (Visium HD & Xenium)
Mean = 1.0 | Right-skewed tail reflecting dense glands and large cells

   Frequency       |
      │       ┌──┐ | Mean (=1)
      │     ┌─┘  └─┐
      │   ┌─┘      |└──┐
      │ ┌─┘        |   └───┐
      │─┘          |       └──────────┐
      └─────────┬──|──────────────────┴────► Size Factor (s_i)
             Median                  Max
             (≈0.78)              (8.9–10.6)
```

### The Critical Prerequisite: Filter QC Failures Before Normalization

Log-normalization strictly depends on prior quality control and sequential 2D pruning. If zero-count observations ($s_i = 0$) or extreme technical dropouts are passed into `logNormCounts()`, the division $X_{ij} / s_i$ becomes mathematically undefined, generating non-finite values (`Inf` or `NaN`) that break downstream linear algebra routines. Consequently, empty capture units and QC-flagged cells must always be excluded before computing size factors.

## How Segmentation Quality Shapes Normalization Metrics

In segmented spatial architectures (such as 10x Xenium or Space Ranger v4+ segmented Visium HD), per-cell UMI totals—and therefore per-cell size factors—are directly governed by the fidelity of the upstream cell segmentation polygons (`cellSeg`):

```text
=============================================================
      HOW CELL SEGMENTATION ERRORS DISTORT SIZE FACTORS
=============================================================

  1. Accurate Segmentation    2. Under-Segmentation          
  ┌───────────────────────┐   ┌───────────────────────┐      
  │   Cell 1     Cell 2   │   │   Merged Polygon      │      
  │   (• • •)    (• • •)  │   │   (• • •  • • •)      │      
  └───────────────────────┘   └───────────────────────┘      
  Counts ≈ True Content       Counts = Inflated Sum          
  Size Factor (s_i) ≈ Normal  s_i = Artificially High

                3. Over-Segmentation          
                ┌────────────────────────┐    
                │ Split Cell 1a │ Cell 1b│    
                │   (• •)       │ (•)    │    
                └───────────────┴────────┘    
                Counts = Deflated Split       
                s_i = Artificially Low
=============================================================
```

* **Under-Segmentation (Merged Cells):** Combines multiple neighboring cells into a single polygon, artificially inflating per-cell UMI totals and generating oversized size factors ($s_i \gg 1$) that over-correct individual gene expression downward.
* **Over-Segmentation (Split Cells):** Fragments a single cell into multiple partial polygons, yielding deflated UMI totals and undersized size factors ($s_i \ll 1$) that artificially amplify background counts during normalization.

## The Spatial Normalization Dilemma & Spatially-Aware Methods

### The Confounding of Technical Depth and Tissue Architecture

In dissociated scRNA-seq, variation in total library size per droplet is predominantly driven by technical capture efficiency. In spatial transcriptomics, however, **library size is strongly confounded with true histological architecture**:

* High-density, metabolically active compartments (e.g., crowded colorectal adenocarcinoma glands or mucosal crypts of Lieberkühn) naturally yield higher UMI totals per unit area than sparse desmoplastic stroma, smooth muscle, or submucosal vasculature.
* Standard global library-size normalization assumes that all variation in total UMI counts is technical. Dividing every unit by its total library size risks flattening legitimate anatomical gradients in total RNA content.

```text
==========================================================================
                     THE SPATIAL NORMALIZATION DILEMMA
==========================================================================

Goal A: Remove Technical Depth     │ Goal B: Preserve Real Spatial Biology
        Variation                            Biology
────────────────────────────────── │ ──────────────────────────────────────
Eliminate local capture efficiency │ Retain natural cell density gradients 
gradients and flowcell tile depth  | and anatomical expression shifts.
differences.                       │

────────────► Conflict: Over-smoothing erases spatial signal! ◄────────────
===========================================================================
```

### SpaNorm: Joint Modeling of Spatial Anatomy and Library Depth

To resolve this confounding, spatially-aware normalization frameworks such as **`SpaNorm`** decompose spatial expression variation into biological and technical components:

* **Statistical Mechanism:** Fits a gene-wise Negative Binomial Generalized Linear Model (NB GLM) using spatial smoothing splines over physical coordinates $(x_c, y_c)$ to separate smooth biological spatial trends $f_g(x_c, y_c; \beta_g)$ from local library-size effects $\log(\text{LS}_c)$:

$$
\log(\mu_{c,g}) = \zeta_g + f_g(x_c, y_c; \beta_g) + \left\lbrace \alpha + h_g(x_c, y_c; \gamma_g) \right\rbrace \log(\text{LS}_c)
$$

* **Adjusted Output:** Computes Percentile-Adjusted Counts (PAC) or Pearson residuals by regressing out only the technical library-size term while preserving the spatial biological component $f_g(x_c, y_c; \beta_g)$.
* **Bioconductor Implementation:** Provided through `SpaNorm::SpaNorm()`, `SpaNormPCA()`, and `SpaNormSVG()`.

## Methodological Comparisons & Platform Pitfalls

### `sctransform`: Regularized Pearson Residuals vs. Spatial Signal Loss

`Seurat::SCTransform()` replaces `log1p` size-factor normalization with regularized Negative Binomial regression against log sequencing depth, returning Pearson residuals:

* **Advantages in scRNA-seq:** Eliminates the arbitrary $+1$ pseudocount, corrects the mean–variance dependency across expression dynamic ranges, and sharpens HVG identification.
* **Drawbacks in Spatial Omics:** Because sequencing depth correlates tightly with histological cell density, regressing out sequencing depth globally without spatial coordinates can inadvertently treat anatomical density gradients as technical noise, attenuating genuine regional expression signatures.

### Targeted Panel Bias in Imaging-Based Spatial Omics (10x Xenium)

Unlike whole-transcriptome Visium HD arrays, imaging-based platforms such as 10x Xenium profile targeted probe panels of 300–5,000 genes. This targeted design introduces a critical mathematical constraint for library-size normalization:

* **Panel Composition Skew:** Library-size normalization assumes that the sum of measured transcripts per cell ($\sum_j X_{ij}$) is proportional to total cellular mRNA content. If a targeted panel is heavily enriched for markers of Cell Type A (e.g., colorectal epithelial and carcinoma markers) but contains few markers for Cell Type B (e.g., quiescent fibroblasts or endothelial cells), Cell Type A will exhibit a much larger panel sum than Cell Type B, even if both cells contain identical total mRNA molecules.
* **Distortion of Relative Fold-Changes:** Dividing by the panel-restricted library size ($s_i$) compresses expression values in panel-enriched cells and inflates expression values in panel-sparse cells, distorting downstream log fold-changes, differential expression tests, and SVG rankings.
* **Analytical Rule:** Normalized values in targeted spatial panels must be interpreted strictly as **within-panel relative proportions**, and we always verify key markers across both raw `counts` and normalized `logcounts` (or use control-probe/cell-area scaling where appropriate).

```text
===========================================================================
              WHY TARGETED PANELS BIAS DEPTH NORMALIZATION
===========================================================================

 Cell Type A (Panel-Enriched)                Cell Type B (Panel-Sparse)
┌───────────────────────────┐               ┌───────────────────────────┐
│ ■■■ Panel Gene 1 (18 UMIs)│               │ ■   Panel Gene 1 (2 UMIs) │
│ ■■■ Panel Gene 2 (15 UMIs)│    Panel      │ ■   Panel Gene 2 (1 UMI)  │
│ ■   Off-Panel (2 UMIs)    │ ────────────► │ ■   Off-Panel (2 UMIs)    │
│ Total Panel UMIs = 33     │ Composition   │ Total Panel UMIs = 3      │
└───────────────────────────┘               └───────────────────────────┘
Higher total count driven by panel design, NOT true cellular RNA abundance!
===========================================================================
```

## Platform Comparison & Practical Implementation Workflow

### Cross-Platform Normalization Matrix

| Parameter / Feature | Visium HD (Binned `spe` / Segmented `sfe`) | 10x Xenium (In Situ Segmented `sfe`) |
| :--- | :--- | :--- |
| **Transcriptomic Scope** | Whole transcriptome (≈ 18,000 genes) | Targeted gene panel (300–5,000 genes; e.g., 422 CRC genes) |
| **Observation Units** | Binned squares (2 µm, 8 µm, 16 µm) or segmented cell polygons | Segmented single-cell polygons (`cellSeg`, `nucSeg`) |
| **Baseline Normalization** | `scuttle::logNormCounts()` or `SpaNorm::SpaNorm()` | `scuttle::logNormCounts()` (after QC & zero-count pruning) |
| **Zero-Count Handling** | Rare in 16 µm tissue bins; pruned during segmented cell QC | **Mandatory:** Exclude zero-count cells (`nCounts < 10`) prior to `logNormCounts()` |
| **Interpretation Scope** | Unbiased whole-transcriptome relative abundance | Within-panel relative abundance (subject to panel composition bias) |

### Practical Bioinformatic Implementation Steps

In our Bioconductor workflows across `SpatialExperiment` (`spe`) and `SpatialFeatureExperiment` (`sfe`) containers, normalization and scaling proceed through five structured steps:

1. **Enforce Sequential 2D QC Pruning:**
   Verify that compromised capture units (`!discard` or `qc_keep`) and zero-count observations have been removed, and unexpressed features (`rowSums(counts) > 0`) pruned, preventing $s_i = 0$ division errors.
2. **Compute Size Factors & Log-Normalized Counts:**
   Run `scuttle::logNormCounts()` (or C++\-accelerated `scrapper::normalizeRnaCounts.se()` for large out-of-memory matrices) to populate `sizeFactors()` and add the `logcounts` assay.
3. **Inspect Size Factor Distributions:**
   Evaluate the empirical distribution (confirming $\text{mean}(s_i) = 1.0$ and inspecting the right-skewed tail) and map `sizeFactor` onto tissue coordinates via `Voyager::plotSpatialFeature()` or `ggspavis::plotCoords()` to verify concordance with histological density.
4. **Compare Raw vs. Normalized Marker Topographies:**
   Contrast canonical histopathological markers (such as mucosal polymeric immunoglobulin receptor `PIGR` or goblet-cell `CLCA1`) across raw `counts`, `logcounts`, and spatially-aware `SpaNorm` assays to verify that technical depth gradients are dampened while biological crypt-to-carcinoma boundaries remain sharp.
5. **Center and Scale Features for Dimensionality Reduction:**
   Apply feature-level $Z$-score scaling (`mean = 0`, `SD = 1`) strictly **after** log-normalization when computing principal components (PCA) and downstream spatial embeddings.
