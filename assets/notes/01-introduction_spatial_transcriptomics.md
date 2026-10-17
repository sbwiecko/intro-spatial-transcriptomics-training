# Introduction to Spatial Transcriptomics (ST)

In traditional transcriptomics (bulk or dissociated single-cell RNA-seq), we analyze high-dimensional gene expression matrices across samples or individual cells. However, tissue dissociation destroys the physical architecture of the specimen, resulting in the complete loss of spatial context.

**Spatial Transcriptomics (ST)** techniques—named *Nature Methods'* "Method of the Year" in 2020—overcome this limitation by measuring gene expression while preserving the physical spatial coordinates $(x,y)$ of transcripts within an intact tissue section. Adding this spatial dimension allows us to investigate tissue microenvironments, cell–cell communication, and spatially variable gene expression across development and disease (such as the human colorectal adenocarcinoma sections we analyze throughout this course).

## 1. The Two Main ST Technologies

Current ST technologies broadly fall into two methodological categories, each generating distinct spatial representations and requiring tailored computational containers.

### A. Sequencing-Based (Spatial Barcoding)

These methods capture endogenous poly-adenylated mRNA (fresh-frozen) or hybridized whole-transcriptome ligation probes (FFPE) onto a spatially barcoded slide array before *ex vivo* next-generation sequencing (NGS).

* **Mechanism:** A tissue section is mounted on an array of barcoded capture features. Each feature carries oligonucleotides with a unique spatial barcode (determining the $(x,y)$ coordinate) and a Unique Molecular Identifier (UMI). Captured transcripts or ligated probes are barcoded, pooled, and sequenced.
* **Data Structure:** Analyzed conceptually as **Lattice Data** (stored in Bioconductor as a `SpatialExperiment` for regular grids or a `SpatialFeatureExperiment` when aggregated into segmented cell polygons).
* **Prominent Example — 10x Visium & Visium HD:**
  * *Visium V1/V2:* Uses 55 µm circular barcoded spots with a 100 µm center-to-center pitch. Because each 55 µm spot captures a mixture of 1–10+ cells, downstream pipelines rely on mathematical deconvolution to estimate cell-type proportions per spot.
  * *Visium HD:* Replaces sparse spots with a continuous lawn of **2 × 2 µm barcoded squares** (~18,000-gene human probe panel). Space Ranger aggregates these native 2-µm squares into multi-scale grids (**8 × 8 µm** and **16 × 16 µm** bins) and, in Space Ranger v4+, integrates H&E nuclear segmentation (StarDist) to group 2-µm squares into single-cell polygons.

### B. Imaging-Based (*In Situ*)

These methods detect RNA molecules directly within the intact tissue using high-resolution fluorescence microscopy, either via cyclic *in situ* hybridization (e.g., 10x Xenium, Vizgen MERFISH, NanoString CosMx) or *in situ* sequencing (e.g., Padlock/ISS, FISSEQ).

* **Mechanism:** Target transcripts are bound by specific padlock or primary probes and decoded across sequential rounds of fluorescent probe hybridization and high-resolution optical imaging.
* **Data Structure:** Analyzed conceptually as **Spatial Point Pattern Data** mapped to **Vector Polygons** (stored in Bioconductor as a `SpatialFeatureExperiment`). Every decoded transcript molecule receives an exact subcellular $(x,y)$ coordinate in physical microns (`µm`), requiring multi-stain cell segmentation (DAPI + membrane/cytoplasmic stains) to assign transcripts to individual cell boundaries.
* **Prominent Example — 10x Xenium:**
  * Supports targeted panels ranging from ≈300 to 5,000 genes alongside built-in negative control probes and optical codewords.
  * Achieves subcellular single-molecule resolution in native physical microns (`µm`) without an array grid.
  * Because capture is panel-targeted, biological interpretation is bounded by the pre-selected gene panel, and empty or severely depleted cell polygons must be filtered prior to normalization.

### Summary Comparison

| Feature | Sequencing-Based (Visium, Visium HD, Slide-seq) | Imaging-Based (Xenium, MERFISH, CosMx) |
| :--- | :--- | :--- |
| **Throughput** | Whole-transcriptome discovery (≈18,000+ genes via poly(A) or probe panels). | Targeted measurement requiring *a priori* panel design (≈100–5,000 genes). |
| **Resolution** | Array-bound resolution (55 µm spots down to 2/8/16 µm bins or H&E-segmented cells). | Subcellular, single-molecule point coordinates (`µm`) inside segmented polygons. |
| **Sensitivity** | Moderate per-molecule capture efficiency (array diffusion/hybridization limits). | High per-gene sensitivity and optical signal-to-noise ratio. |
| **Strengths** | Unbiased whole-transcriptome discovery; broad pathway and ligand–receptor screening. | Exact subcellular transcript localization; high sensitivity for low-abundance targets. |

## 2. Choosing the Right Technology

Selecting a spatial transcriptomics platform is driven by our biological question and experimental constraints:

* **Gene Throughput (Discovery vs. Targeted):** Does our hypothesis require unbiased whole-transcriptome screening across all signaling pathways (sequencing-based), or are we interrogating a defined panel of lineage and state markers at high sensitivity (imaging-based)?
* **Spatial Resolution:** Do we require exact subcellular transcript coordinates and multi-stain membrane boundaries (imaging-based), or is a continuous 8–16 µm array / H&E-segmented profile sufficient (Visium HD)?
* **Sensitivity & Low-Abundance Targets:** Are we profiling lowly expressed transcription factors or receptors that benefit from the high molecular capture efficiency of cyclic *in situ* imaging?
* **Sequence-Level Information:** Do we need direct cDNA sequence reads (e.g., poly(A)-captured fresh-frozen sequencing for somatic variants, TCR/BCR clonotypes, or viral transcripts)?
* **Practical Feasibility:** What is our tissue preservation format (FFPE vs. fresh-frozen), slide area, imaging turnaround time, and computational storage budget?

## 3. The Spatial Analytical Paradigm

In dissociated single-cell RNA-seq, the standard workflow proceeds linearly: **Quality Control → Normalization → Highly Variable Genes (HVGs) → PCA/UMAP → Clustering → Marker Annotation**.

In Spatial Transcriptomics, coupling the gene expression matrix to physical $(x,y)$ coordinates expands every stage of this paradigm:

* **Spatially Aware Quality Control:** Instead of applying rigid global thresholds that accidentally strip away low-UMI stromal or vascular compartments, we evaluate each cell or bin against its local physical neighborhood (e.g., via `SpotSweeper` local $Z$-scores).
* **Spatially Variable Genes (SVGs):** Beyond purely non-spatial variance (HVGs), we use spatial autocorrelation statistics (e.g., Moran's $I$, Geary's $C$) to identify genes whose expression forms coherent spatial gradients or anatomical structures.
* **Spatial Domains & Niches:** Rather than clustering cells solely by transcriptional similarity, spatially aware clustering algorithms combine gene expression with neighborhood graphs to delineate histological compartments (e.g., normal colonic mucosa, invasive adenocarcinoma glands, and desmoplastic stroma).
* **Cell–Cell Co-localization & Communication:** Using spatial neighborhood graphs, we statistically test whether specific cell types physically border or signal to one another more frequently than expected by chance.
* **Reference Integration & Deconvolution:** By pairing spatial datasets with annotated single-cell reference atlases (e.g., via `RCTD`), we map cell identities onto segmented polygons or deconvolve cell-type mixtures within multi-cell grid bins.

### Core Exploratory Operations

Throughout our hands-on notebooks, we systematically execute five core operations on our spatial objects:

1. **Select:** Cropping and zooming into specific histological regions of interest (ROIs).
2. **Score:** Quantifying quality metrics, gene set signatures, and spatial autocorrelation.
3. **Cluster:** Identifying transcriptomic cell states and coherent spatial domains.
4. **Characterize:** Annotating cell types and differential marker genes across tissue niches.
5. **Relate:** Testing for spatial neighborhood enrichment and cell–cell interactions.