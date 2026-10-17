# Cell Segmentation & Transcript Assignment in Spatial Transcriptomics

## The Assignment Problem

In spatial transcriptomics (SRT), measuring RNA abundance alongside physical coordinates helps us unravel tissue architecture and cell–cell communication. However, extracting biological meaning at single-cell resolution requires accurately assigning raw transcripts to individual cells.

Our approach to this "assignment problem" depends on the spatial capture architecture:

* **Imaging-Based ST (Single-Molecule, `SpatialFeatureExperiment`):** Generates subcellular $(x,y)$ transcript coordinates in physical microns that require **image and/or transcript-assisted segmentation** to draw polygon boundaries (`cellSeg`, `nucSeg`) around distinct single cells.
* **Sequencing-Based ST (Spot/Grid, `SpatialExperiment` vs. Segmented `SpatialFeatureExperiment`):** Captures transcripts on barcoded arrays (55 µm spots or 2/8/16 µm square bins) that either require mathematical **deconvolution** on multi-cell grids (`SPE`) or **H&E-guided nuclear segmentation and expansion** to group 2-µm capture reads into single-cell polygons (`SFE`).

## 1. True Image Segmentation (Imaging-Based Methods)

Because we have exact $(x,y)$ coordinates for every decoded RNA molecule in imaging-based ST (such as 10x Xenium or Vizgen MERFISH), our goal is to draw accurate cell boundary polygons on the registered morphology image and count the transcripts falling inside each polygon.

* **Evolution of Algorithms:** Image segmentation has evolved from classical morphological watersheds to deep learning and foundation models (e.g., U-Net, StarDist, Cellpose, and CellSAM).
* **Multi-Stain Approaches (e.g., 10x Xenium Onboard Analysis):** Modern pipelines combine multiple fluorescent channels to delineate cell borders accurately. For example, Xenium pairs DAPI (nuclear staining) with ATP1A1/CD45/E-Cadherin (membrane boundary staining) and 18S rRNA + αSMA/Vimentin (interior cytoplasmic staining).
* **Nuclear Expansion Fallback:** Where membrane stains are absent or locally weak, algorithms segment the DAPI-stained nucleus (`nucSeg`) first and isotropically expand the boundary outward by a fixed distance (e.g., 5–15 µm) or until colliding with a neighboring cell (`cellSeg`).

### Advanced Image & Transcript-Assisted Algorithms

While earlier methods relied solely on staining morphology, modern tools incorporate the spatial point cloud of transcripts itself to build or refine cell boundaries:

* **Baysor (Bayesian Segmentation):** A pioneering "segmentation-free" or prior-assisted method that clusters transcripts based on local transcriptional composition and spatial density, optionally integrating precalculated nuclear or membrane stains as Bayesian priors.
* **Segger (Graph Neural Network Segmentation):** Models nuclei/cell boundaries and individual transcripts as nodes in a heterogeneous spatial graph, using GPU-accelerated Graph Neural Networks (GNNs) to resolve ambiguous border transcripts significantly faster than Markov Chain Monte Carlo (MCMC) sampling in Baysor.
* **Proseg (Probabilistic Cell Segmentation):** Infers the specific 3D/2D cell shapes that best explain the observed physical transcript coordinates. Initialized from nuclear segmentation, Proseg iteratively refines cell boundaries while explicitly modeling "leaky" cells (extracellular transcript diffusion) via a spatial migration prior.
* **SMURF (Segmentation and Manifold Unrolling Framework):** Performs concentric expansion around segmented nuclei to assign high-resolution spatial capture spots to cells, iteratively maximizing intra-cluster cosine similarity and deconvolving shared border spots so that every transcript is assigned to at most one cell.

## 2. Spatial Binning, H&E Segmentation & Deconvolution (Sequencing-Based Methods)

Because sequencing-based technologies capture poly(A) or probe-hybridized transcripts on barcoded slide arrays rather than optically decoding molecules inside stained membranes, resolving single cells relies on three complementary strategies:

* **Spatial Binning (Visium HD Binned Output → `SpatialExperiment`):** Space Ranger aligns sequencing reads to the slide's barcoded array and aggregates native 2 × 2 µm capture squares into regular multi-scale grids (8 × 8 µm and 16 × 16 µm bins). While 8–16 µm bins provide robust UMI coverage for tissue-level profiling, individual squares frequently straddle adjacent cell borders.
* **Visium HD Cell Segmentation (Space Ranger v4.0+ & Bin2Cell → `SpatialFeatureExperiment`):** To bridge the gap between array grids and true single cells, Space Ranger v4+ applies a deep-learning **StarDist** model directly to the high-resolution H&E image to segment nuclei (`nucSeg`), expands those boundaries into whole-cell polygons (`cellSeg`), and aggregates the underlying 2 × 2 µm barcoded capture spots that fall within each cell footprint (similar in spirit to community tools like **Bin2Cell**).

    ![Visium HD data can be analyzed through cell segmentation, assigning 2-µm features to individual nuclei and cells based on the high-resolution microscope image.](https://cdn.10xgenomics.com/image/upload/v1765926999/analysis-guides/visium-hd-downstream-analysis/02_bins_seg.png)

  * *Limitations:* Because transcript capture is physically quantized into 2 × 2 µm array squares, a single 2 × 2 µm bin at a cell–cell interface can overlap two adjacent cells, creating subtle boundary spillovers. Furthermore, segmentation fidelity depends heavily on H&E staining clarity and tissue integrity.
* **Spot Deconvolution (e.g., RCTD — Robust Cell Type Decomposition):** When working with unsegmented multi-cell spots (55 µm Visium) or coarse bins (16 µm Visium HD), or when identifying doublets in segmented cells, we integrate an annotated single-cell RNA-seq reference atlas. Algorithms like **RCTD** learn cell-type-specific expression signatures from the reference and mathematically infer the cell-type proportions (or doublet status) within each spatial spot or segmented cell.

## 3. Correcting for Noise, Spillover, and Segmentation Errors

Even with state-of-the-art segmentation, spatial transcriptomics suffers from **"spillover" (molecular bleeding)**—caused by lateral RNA diffusion during permeabilization, optical crowding, or 3D cell overlap within a 2D histological section (e.g., transcripts from dense colonic epithelial crypts bleeding into adjacent lamina propria immune cells).

Specialized post-segmentation models are used to purify single-cell expression profiles:

### ResolVI

* **Mechanism:** Models spatial noise and lateral diffusion using a deep generative Variational Autoencoder (VAE).
* **Correction:** Decomposes each cell's observed count vector into a sum of three latent components:
    1. The *true endogenous* expression of the cell.
    2. *Spillover* expression diffused from physical neighbors.
    3. Unspecific *background* noise.

### SPLIT (Spatial Purification of Layered Intracellular Transcripts)

* **Mechanism:** An R framework that pairs RCTD reference decomposition with spatial neighborhood context to evaluate transcriptional purity per segmented cell.
* **Identification:** Flags heterotypic "doublet cells" whose transcriptomes reflect a mixture of a primary cell type and a contaminating neighbor.
* **Correction:** Mathematically purifies the count matrix by retaining the read fraction attributable to the primary cell identity and stripping away secondary spillover transcripts.

### cellAdmix

* **Mechanism:** Quantifies gene-level contamination via a Bayesian admixture model conditioned on spatial cell-type adjacency.
* **Factorization:** Builds a local $k$-NN transcript graph per cell and applies Non-negative Matrix Factorization (NMF) to separate endogenous programs from neighbor-derived contamination programs.
* **Correction:** Uses a Conditional Random Field (CRF) to assign individual transcript molecules to latent factors, removing molecules whose factor identity conflicts with the host cell type.

## 4. Production Workflow Example: The Snakemake Xenium Pipeline

In production settings, processing whole-slide imaging datasets (≈1 day of instrument runtime per Xenium slide) requires automated workflow orchestration (e.g., **Snakemake**) that chains resegmentation, format conversion, and spillover correction:

### Phase 1: Initial Output & Resegmentation

* **Resegmentation:** Ingests the primary Xenium Onboard Analysis bundle and executes containerized `xeniumranger` resegmentation or transcript-assisted refinement.
* **Image Registration:** Aligns post-run H&E pathology images to the DAPI/morphology coordinate frame so histological annotations map directly onto transcript coordinates.

### Phase 2: Advanced Segmentation & Container Formatting

From the aligned bundle, transcript coordinates and stains are routed into language-specific ecosystems:

* **Proseg Fork (R / Bioconductor):** Executes probabilistic transcript-assisted segmentation and exports directly to a **`SpatialFeatureExperiment` (`SFE`)** object.
* **Segger Fork (Python / scverse):** Executes GNN-based transcript segmentation and exports to the **`SpatialData`** container format.

### Phase 3: Downstream Processing Branches

Once encapsulated in `SFE` or `SpatialData`, downstream modules execute sequentially:

1. **ROI Cropping:** Trims slide-edge artifacts or isolates specific histological regions of interest.
2. **Multi-Sample Integration:** Embeds multiple patient sections into a shared latent space (e.g., via `scVI`).
3. **Annotation & Purification:** Runs **RCTD** against a reference scRNA-seq atlas to assign cell types, followed by **ResolVI** or **SPLIT** to strip lateral spillover and yield a purified single-cell spatial atlas.

## 5. Take-Home Messages

* **Foundational but Imperfect:** Cell segmentation sits at the very top of the single-cell spatial workflow and directly governs downstream clustering and marker specificity, yet every segmentation modality involves trade-offs between sensitivity (capturing cytoplasmic transcripts via expansion) and specificity (avoiding neighbor spillover).
* **Transcript-Assisted Refinement:** Combining optical nuclear/membrane stains with the spatial point cloud of transcripts (via Baysor, Proseg, or Segger) substantially improves boundary accuracy over morphology-only expansion.
* **Chaining Complementary Tools:** Robust single-cell spatial analysis rarely relies on a single algorithm; instead, we chain segmentation (StarDist / Proseg / Segger), reference-guided decomposition (RCTD), and spillover purification (ResolVI / SPLIT).
* **Bridging R and Python Ecosystems:** Because cutting-edge segmentation and deep-learning correction tools span both Python (`SpatialData`, `scverse`) and R/Bioconductor (`SpatialFeatureExperiment`, `Voyager`), reproducible pipelines rely on standardized on-disk formats (Parquet, GeoJSON, HDF5, `alabaster.sfe`) to pass geometries and matrices cleanly across environments.
