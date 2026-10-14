# Introduction to Spatial Transcriptomics (ST)

As a researcher familiar with traditional transcriptomics (like bulk or single-cell RNA-seq), you already know how to analyze gene expression matrices. However, standard techniques require tissue dissociation, resulting in the complete loss of spatial context. 

**Spatial Transcriptomics (ST)** techniques, named Nature's "Method of the Year" in 2020, overcome this by measuring gene expression while preserving the exact spatial coordinates (x, y) of the transcripts within a tissue slice. This adds a physical dimension to the data, allowing us to study spatial tissue architecture, cell-cell communication, and spatially variable gene expression in development and disease.

## 1. The Two Main ST Technologies

Current ST technologies broadly fall into two methodological categories, each generating distinct types of data representations and requiring different analytical approaches.

### A. Sequencing-Based (Spatial Barcoding)

These methods capture poly-adenylated RNA on a spatially barcoded slide before reverse transcription and ex vivo next-generation sequencing (NGS).

* **Mechanism:** Tissues are placed on a lattice of spots. Each spot contains oligos with a unique spatial barcode. Transcripts are captured, barcoded, and sequenced.
* **Data Structure:** Analyzed conceptually as **Lattice Data**. The spatial locations are fixed grids or spots, and we study the expression of genes at these discrete locations.
* **Prominent Example — 10x Visium:**
  * *Visium V1/V2:* Uses 55µm barcoded spots. Because a 55µm spot often captures a mixture of multiple cells, analytical pipelines rely on "deconvolution" to estimate cell types within spots.
  * *Visium HD:* Dramatically increases resolution using a continuous lawn of oligos binned into a grid of 2x2µm squares (further binned to 8x8µm).

### B. Imaging-Based (_in situ_)

These methods detect transcripts directly within the intact tissue using high-resolution microscopy and fluorescent probes, either via *In Situ* Hybridization (e.g., MERFISH, smFISH) or *In Situ* Sequencing (e.g., FISSEQ).

* **Mechanism:** Uses sequential rounds of complementary fluorescent probes to decode specific transcripts directly inside the tissue.
* **Data Structure:** Analyzed conceptually as **Point Pattern Data**. We get exact (x, y) coordinates for individual transcript molecules, requiring image-based cell segmentation (using DAPI and boundary stains) to assign transcripts to individual cells.
* **Prominent Example — 10x Xenium:**
  * Supports targeted panels ranging from 500 to 5,000 genes.
  * Resolution is at the single-molecule level, mapping to physical microns without a spot grid.
  * Because it is targeted, data must be interpreted strictly within the panel's confines, and cells with zero counts are usually removed before downstream normalization.

### Summary Comparison

| Feature         | Sequencing-Based (Visium, Slide-seq)                                                  | Imaging-Based (Xenium, MERFISH)                                                 |
|:--------------- |:------------------------------------------------------------------------------------- |:------------------------------------------------------------------------------- |
| **Throughput**  | Unbiased, whole-transcriptome discovery (10,000+ genes).                              | Targeted measurement requiring *a priori* gene selection (100s to 5,000 genes). |
| **Resolution**  | Spot-level resolution (100µm down to 2µm) restricted to fixed grids.                  | Sub-cellular/single-molecule resolution with positional randomness.             |
| **Sensitivity** | Lower (capture-array limits).                                                         | High sensitivity and high signal-to-noise ratio.                                |
| **Strengths**   | Discovery of novel markers, unbiased exploration, sequence info (e.g., RNA velocity). | Precise transcript localization, mapping true cell-cell physical interactions.  |

## 2. Choosing the Right Technology

The selection of a spatial transcriptomics platform should be driven strictly by the core scientific question rather than methodological novelty. When designing an experiment, evaluate the following trade-offs:

* **Gene Throughput (Discovery vs. Targeted):** Does your hypothesis require unbiased whole-transcriptome discovery, or are you looking to map specific, known marker panels?
* **Resolution Needs:** Does the analysis require precise sub-cellular transcript localization and exact cell-cell boundaries (imaging), or is broad regional tissue architecture sufficient (sequencing)?
* **Sensitivity & Expression Levels:** Are you trying to capture low-abundance transcripts that require the high sensitivity of imaging methods?
* **Sequence Information:** Do you need exact sequence data for downstream analysis, such as identifying isoforms or splicing events (sequencing)?
* **Feasibility:** What is the tissue availability (FFPE vs. fresh frozen)? Do you have access to specialized microscopy, compute power, and sufficient budget?

## 3. The Analytical Paradigm

In single-cell RNA-seq, the standard workflow is: **Identify Highly Variable Genes (HVGs) → Cluster Cells → Identify Cluster Markers**.

In Spatial Transcriptomics, the addition of the $(x, y)$ coordinate matrix shifts this paradigm:

* **Spatially Variable Genes (SVGs):** Instead of just looking for high variance, we use spatial statistics (e.g., Moran's $I$, Geary's $C$) to find genes whose expression patterns strictly correlate with physical location.
* **Spatial Domains:** Instead of purely transcriptomic cell clusters, we cluster spots/cells by combining gene expression *and* spatial neighborhood information to define anatomical or functional "spatial domains" (e.g., cortex layers, tumor microenvironments).
* **Cell-Cell Communication & Co-localization:** Because we know where cells are, we can statistically test whether Cell Type A physically neighbors Cell Type B more often than by random chance (using point pattern analysis).
* **Integration / Deconvolution:** For sequencing-based methods with larger spots (like Visium V1), we computationally "deconvolve" the bulk-like ST spots using reference scRNA-seq datasets to estimate the proportions of different cell types residing within a single spot.

### Core Data Operations

During the upcoming exercises, we will routinely perform five classes of exploratory operations on the Spatial Gene Expression Matrix:

1. **Select:** Focusing on specific regions of interest.
2. **Score:** Evaluating spatial autocorrelation or gene set enrichment.
3. **Cluster:** Finding coherent spatial domains.
4. **Characterize:** Annotating the cell types making up the tissue.
5. **Relate:** Testing for spatial co-localization and communication between cells or domains.