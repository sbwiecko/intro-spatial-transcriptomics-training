# Spatial Transcriptomics Feature Selection & Dimensionality Reduction

## From Normalized Matrices to Low-Dimensional Manifolds

In upstream modules, we performed quality control to eliminate compromised capture units and applied size-factor log-normalization (`logcounts`) to standardize sequencing depth.

While our normalized matrices accurately reflect relative gene expression, they present two analytical bottlenecks:

1. **High Dimensionality:** Whole-transcriptome datasets (such as 10x Visium HD) measure P ≈ 18,000 genes across tens of thousands of cells or bins, where the vast majority of features represent invariant housekeeping genes or stochastic technical noise.
2. **Spatial Autocorrelation:** Unlike dissociated single-cell RNA-seq (scRNA-seq), spatial capture units are physically anchored in histological space $(x, y)$. Neighboring cells share morphogen gradients, structural microenvironments, and cell–cell contacts.

**Intermediate processing** bridges raw matrix processing and biological interpretation through two tightly coupled steps: **Feature Selection** and **Dimensionality Reduction**.

```text
=================================================================================
                     SPATIAL TRANSCRIPTOMICS PIPELINE POSITION
=================================================================================

 [01–03 Data Struct.]  [04–06 QC/07–09 Norm]  [10 Intermediate]  [11+ Downstream]
 ┌──────────────────┐  ┌───────────────────┐  ┌───────────────┐  ┌──────────────┐
 │ Raw Files/Images │  │  Filter Outliers  │  │Feature Select │  │Spatial Domain│
 │ SCE / SPE / SFE  │─►│Size-Factor LogNorm│─►│ Dim Reduction │─►│  Clustering  │
 │ Matrix Objects   │  │ Feature Z-Scaling │  │ (PCA, BANKSY, │  │  Cell Type   │
 └──────────────────┘  └───────────────────┘  │     UMAP)     │  │  Annotation  │
                                              └───────────────┘  └──────────────┘
=================================================================================
```

## Feature Selection: Highly Variable Genes (HVGs) vs. Spatially Variable Genes (SVGs)

Feature selection reduces matrix dimensionality from ≈ 18,000+ genes down to a focused, highly informative subset (typically 1,000–2,000 genes), eliminating uninformative technical variation and accelerating downstream computations.

```text
=================================================================================
                    FEATURE SELECTION: MOLECULAR VS. SPATIAL
=================================================================================

     1. Highly Variable Genes (HVG)          2. Spatially Variable Genes (SVG)
┌─────────────────────────────────────┐  ┌───────────────────────────────────────┐
│ • Uses ONLY expression matrix       │  │ • Uses expression + (x,y) coordinates │
│ • Models mean-variance relationship │  │ • Identifies coherent spatial patterns│
│ • Ignores physical location (x,y)   │  │ • Ranks genes by spatial structure    │
└─────────────────────────────────────┘  └───────────────────────────────────────┘
                   │                                         │
                   ▼                                         ▼
        Captures major cell types              Captures spatial gradients &
      (e.g., Epithelium vs. Stroma)         tissue layers (e.g., Crypt L1–L6)
=================================================================================
```

### Highly Variable Genes (HVGs)

* **Mathematical Intuition:** Evaluates the mean-variance relationship across all cells via `scran::modelGeneVar()`. Decomposes total feature variance into a fitted technical trend (mean-dependent Poisson/negative binomial dispersion) and residual biological variance:
  $$\sigma^2_{\text{bio}, g} = \sigma^2_{\text{total}, g} - \sigma^2_{\text{tech}, g}$$
* **Limitation in Spatial Omics:** HVGs treat each cell as an isolated observation in unconstrained expression space. While highly effective at identifying dominant global lineages (e.g., adenocarcinoma epithelium vs. desmoplastic stroma), HVGs frequently miss localized spatial gradients, microvascular niches, or fine boundary markers whose total variance is modest across the entire slide.

### Spatially Variable Genes (SVGs)

* **Mathematical Intuition:** Incorporates the spatial coordinate matrix $(x, y)$ alongside expression to identify features whose spatial distribution deviates significantly from spatial randomness (Poisson point processes or spatial white noise).
* **Classification of SVG Algorithms** ([Yan, Hua & Li, *Nat Commun*, 2025](https://doi.org/10.1038/s41467-025-56080-w)):

```text
┌─────────────────────────────────────────────────────────────────────────────────┐
│                        3 MAIN CATEGORIES OF SVG METHODS                         │
├─────────────────────────────────────────────────────────────────────────────────┤
│ 1. Overall SVGs:                                                                │
│    Detects continuous spatial patterns across the entire tissue slide           │
│    • Gaussian Process Models: SpatialDE, nnSVG                                  │
│    • Spatial Autocorrelation Statistics: Moran's I, Geary's C                   │
│    • Non-Parametric Covariance Tests: SPARK, SPARK-X                            │
├─────────────────────────────────────────────────────────────────────────────────┤
│ 2. Spatial Domain-Marker SVGs:                                                  │
│    Identifies genes differentially expressed between discrete tissue domains    │
│    • Cluster-Based Testing: DESpace, SpaGCN                                     │
├─────────────────────────────────────────────────────────────────────────────────┤
│ 3. Cell-Type Specific SVGs:                                                     │
│    Detects spatial expression variation within a single cell type across niches │
│    • Niche-Specific Modeling: C-SIDE, CTSV                                      │
└─────────────────────────────────────────────────────────────────────────────────┘
```

### Benchmarking Key SVG Algorithms

| SVG Algorithm | Mathematical Basis | Computational Complexity | Gene-Specific Length Scale ($l_g$) | Best Application |
| :--- | :--- | :--- | :---: | :--- |
| **`nnSVG`** ([Weber 2023](https://doi.org/10.1038/s41467-023-39748-z)) | Nearest-Neighbor Gaussian Process (NNGP) | **Linear $O(n)$** | **Yes** | Unbiased discovery across fine and coarse spatial scales. |
| **`SPARK-X`** ([Zhu 2021](https://doi.org/10.1186/s13059-021-02404-0)) | Distance Covariance & Non-Parametric Kernels | **Linear $O(n)$** | **No** (Fixed grid) | Ultra-fast screening on massive datasets (>50k cells). |
| **`DESpace`** ([Cai 2024](https://doi.org/10.1093/bioinformatics/btae027)) | Spatial Cluster + NB GLM Testing | Dependent on clusters | N/A (Cluster-based) | Identifying marker genes for anatomical layers/domains. |
| **`Moran's I`** ([Anselin 1995](https://doi.org/10.1111/j.1538-4632.1995.tb00338.x)) | Spatial Autocorrelation Index ($W$ matrix) | $O(n \cdot k)$ on Graph | **No** (Fixed graph distance) | Rapid baseline ranking of spatial autocorrelation. |

> **Why Gene-Specific Length Scales Matter (`nnSVG`):** Biological structures operate across vastly different physical scales: invasive tumor sheets span hundreds of microns ($l_g$ > 200 µm), whereas microvascular capillaries or mucosal crypt stems operate across localized 15–30 µm niches ($l_g$ ≈ 20 µm). Fixed-scale methods (`SPARK-X`) smooth over fine-grained niches, whereas `nnSVG` estimates a dedicated length scale $l_g$ per gene.

## Dimensionality Reduction: Molecular vs. Spatially-Aware Embeddings

Dimensionality reduction compresses the selected features (1,000–2,000 HVGs or SVGs) into 20–50 orthogonal latent components (e.g., *Principal Components*) and 2D non-linear manifolds (UMAP), dampening stochastic gene dropouts and accelerating graph clustering.

### Standard Principal Component Analysis (PCA)

Standard PCA (`scater::runPCA`) rotates the normalized expression matrix along orthogonal axes of maximum molecular variance across cells:

* **Mechanism:** Computes linear combinations of features that maximize sample-to-sample variance without taking spatial coordinates $(x, y)$ into account.
* **Limitation in Spatial Omics:** Because standard PCA is unconstrained by tissue topology, top components can be dominated by cell-autonomous technical variance, cell-cycle phase, or stochastic expression spikes rather than histological architecture.

### Spatially-Aware Linear Embeddings (MULTISPATI PCA & SpatialPCA)

To recover genuine histological patterns, spatially-aware linear methods explicitly couple the gene expression matrix with spatial topology:

* **MULTISPATI PCA ([Dray 2008](https://doi.org/10.3170/2007-8-18312) / `Voyager`):** Introduces a spatial weight matrix $W$ (based on Moran's $I$) directly into the eigenvalue problem, maximizing spatial autocorrelation $\mathbf{v}^T X^T W X \mathbf{v}$ alongside sample variance. Positive eigenvalues isolate spatially continuous tissue compartments (e.g., crypt epithelium, desmoplastic stroma), while negative eigenvalues isolate high-frequency spatial noise or anti-correlated boundaries.
* **SpatialPCA ([Shang & Zhou 2022](https://doi.org/10.1038/s41467-022-30188-3)):** Employs a probabilistic PCA framework with a spatial kernel matrix to model spatial correlation across physical locations, enhancing cluster boundary preservation.

### The BANKSY Framework (Neighborhood Feature Expansion)

**BANKSY** (**B**uilding **A**ggregates with a **N**eighborhood **K**ernel and **S**patial **Y**ardstick; [Singhal et al., *Nat Genet*, 2024](https://doi.org/10.1038/s41588-024-01664-3)) embeds spatial context directly into feature space prior to dimensionality reduction by augmenting a cell's own transcriptomic vector with the mean of its physical neighbors.

#### Mathematical Formulation

1. **Spatial Weight Matrix ($W$):** From a spatial neighborhood graph ($k$-NN or Delaunay triangulation), we construct a row-standardized spatial weight matrix $W \in \mathbb{R}^{N \times N}$ where $\sum_j W_{ij} = 1$.
2. **Neighborhood Feature Aggregation:** For each cell $i$ with normalized expression vector $\mathbf{x}_i \in \mathbb{R}^P$, its local microenvironmental average $\bar{\mathbf{x}}_i$ is computed as:
   $$\bar{\mathbf{x}}_i = \sum_{j \in \mathcal{N}(i)} W_{ij} \mathbf{x}_j$$
3. **Feature Space Augmentation:** BANKSY concatenates the cell's own expression vector with its neighborhood mean, scaled by a spatial mixing parameter $\lambda \in [0, 1]$:
   $$\mathbf{b}_i = \begin{bmatrix} \sqrt{1 - \lambda} \cdot \mathbf{x}_i \\ \sqrt{\lambda} \cdot \bar{\mathbf{x}}_i \end{bmatrix}$$
4. **Spatio-Molecular PCA:** Standard PCA is computed directly on the augmented matrix $B \in \mathbb{R}^{N \times 2P}$.

```text
=============================================================================
    BANKSY: NEIGHBORHOOD FEATURE EXPANSION FOR SPATIO-MOLECULAR EMBEDDINGS
=============================================================================

 1. Physical Tissue Graph               2. Dual Feature Representation
 ────────────────────────               ──────────────────────────────
      (Neighbor j1)
           o                             Target Cell Vector (Cell-Autonomous)
            \ W_i,j1                     x_i = [Gene 1, Gene 2, ..., Gene P]
             \
   (j4) o ── Target Cell (i) ── o (j2)   Neighborhood Mean Vector
             /   x_i      \                 (Microenvironment)
            / W_i,j3       \ W_i,j2      x̄_i = ∑_{j ∈ N(i)} W_ij · x_j
           o                o
      (Neighbor j3)    (Neighbor j5)

 3. Matrix Augmentation with Spatial Mixing Parameter λ ∈ [0, 1]
 ──────────────────────────────────────────────────────────────
      ┌───────────────────────────┬──────────────────────────┐
  B = │      √(1 - λ) · x_i       │         √λ · x̄_i         │
      └───────────────────────────┴──────────────────────────┘
         Intrinsic Cell Identity    Microenvironmental Niche
              (P Features)               (P Features)

 4. Dimensionality Reduction & Analytical Trade-Off
 ─────────────────────────────────────────────────
             Standard PCA on Augmented Matrix B (2P Features)
                                  │
            ┌─────────────────────┴──────────────────────┐
            ▼                                            ▼
   λ ≈ 0.2: Fine Cell Typing                λ ≥ 0.7: Spatial Tissue Domains
 (Resists dropouts & noise,              (Delineates shared microenvironments
  preserves single-cell state)            e.g., Crypt Mucosa, Invasive Margin)
==============================================================================
```

#### The Analytical Role of $\lambda$

* **$\lambda = 0$ (Non-Spatial Baseline):** Identical to standard PCA. Groups cells strictly by intrinsic transcriptional identity regardless of their spatial location.
* **$\lambda \approx 0.2$ (Single-Cell Typing):** Recommended default. Smooths technical dropout noise through local neighbor aggregation while preserving sharp single-cell boundary transitions.
* **$\lambda \ge 0.7$ (Spatial Tissue Domains):** Microenvironment-dominated. Merges distinct interacting cell types that co-reside in the same anatomical niche (e.g., malignant epithelial cells, cancer-associated fibroblasts, and tumor-infiltrating lymphocytes within the invasive margin).

## Practical Implementation Workflow

In `10-feature_selection_dimensionality_reduction_VisiumHD_segmented.ipynb`, we execute this intermediate processing framework across our segmented colorectal carcinoma dataset (`sfe_visiumhd_segmented_normalized`):

1. **Model Gene Variance (HVGs):** Apply `scran::modelGeneVar()` on `logcounts` to select the top 1,000 highly variable features.
2. **Evaluate Spatial Autocorrelation (SVGs):** Compute Moran's $I$ and inspect spatial structure across canonical CRC markers (`PIGR`, `CEACAM5`, `MUC17`, `OLFM4`).
3. **Standard PCA:** Run `scater::runPCA()` on top HVGs, inspect variance explained via scree plots, and project top PC loadings onto spatial tissue coordinates with `plotSpatialFeature()`.
4. **BANKSY Spatio-Molecular Embedding:** Build the spatial neighbor graph, run `Banksy::computeBanksy()` and `Banksy::runBanksyPCA()`, and contrast non-spatial vs. spatial principal components.
5. **Non-Linear Manifold Visualization:** Compute UMAP on both standard PCA and BANKSY embeddings to compare cell-type separation and spatial domain continuity.
6. **Object Serialization:** Save the enriched `sfe` object containing dimensionality reduction slots (`PCA`, `BANKSY`, `UMAP`) via `alabaster.sfe::saveObject()` for downstream spatial clustering.
