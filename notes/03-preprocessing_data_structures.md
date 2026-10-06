# Spatial Transcriptomics Preprocessing & Data Structures

## Introduction

Here we provide a comprehensive synthesis of **Slide Deck 03: Preprocessing & Data Structures** from the SIB Training Group's Spatial Transcriptomics course. It serves as a pedagogical bridge between the physical/biochemical foundations of spatial omics technologies (covered in ["Introduction to spatial transcriptomics"](/notes/01-introduction_spatial_transcriptomics.md) & ["Cell segmentation for spatial transcriptomics"](/notes/02-segmentation_spatial_transcriptomics.md)) and the hands-on computational exercises in the accompanying Jupyter Notebook.

So far, we established the core technological split in spatial transcriptomics (SRT):

* **Sequencing-based technologies** (e.g., 10x Visium, Visium HD, Slide-seq) capture poly(A) or probe-hybridized transcripts on barcoded arrays. Data from these platforms feature structured coordinate grids (e.g., 55µm spots or 2µm binned squares) where spatial resolution is dictated by the array geometry rather than individual cell boundaries.
* **Imaging-based technologies** (e.g., 10x Xenium, Vizgen MERFISH, CosMx) detect RNA molecules _in situ_ using optical fluorescence and cyclic decoding. Data from these platforms yield subcellular transcript coordinates mapped directly to physical microns, requiring cell segmentation (via DAPI nuclear stains, boundary markers, or cell expansion) to define single cells.

Here we transitions from **hardware mechanics and slide chemistry** to **bioinformatic pipelines and data representations**. It details how raw sequencing reads (FASTQ) or decoded fluorescent spots (transcripts) are transformed into count matrices, spatial coordinate tables, polygon geometries, and S4 object structures in R/Bioconductor (specifically `SingleCellExperiment`, `SpatialExperiment`, and `SpatialFeatureExperiment`).

## Upstream preprocessing & file output structures

Primary data processing converts raw instrument signals into structured matrices and spatial metadata. The two major commercial pipelines—**Space Ranger** (for 10x Visium/Visium HD) and **Xenium Ranger / Onboard Analysis** (for 10x Xenium)—produce distinct file architectures reflecting their underlying spatial capture mechanics.

### Space Ranger pipeline (10x Visium & Visium HD)

Space Ranger processes demultiplexed FASTQ files and aligns sequencing reads to tissue coordinates.

#### FASTQ naming convention

Space Ranger inputs follow standard Illumina naming structures:

> `Library1_S1_L001_R1_001.fastq.gz`

| File Component | Meaning / Description |
| :--- | :--- |
| `Library1` | Library name derived from the tissue sample |
| `S1` | Sample number in the run |
| `L001` | Flowcell lane number |
| `R1` / `R2` | Read direction (`R1`: spatial barcode + UMI; `R2`: cDNA transcript read) |
| `I1` / `I2` | Dual-index reads (`i7` / `i5`) |
| `.001` | Chunk/split file number |
| `.fastq.gz` | Gzip-compressed FASTQ sequence file |

#### Space Ranger directory architecture (`binned_output`)

For Visium HD, Space Ranger aggregates 2µm barcoded capture spots into multi-resolution spatial bins:

```text
binned_output/
├── square_002um/    # 2µm×2µm native capture spots (sub-cellular, highly sparse)
├── square_008um/    # 8µm×8µm binned grid (4×4 pool of 16 spots; cell-scale)
├── square_016um/    # 16µm×16µm binned grid (8×8 pool of 64 spots; multi-cell)
│   ├── filtered_feature_bc_matrix/     # Count matrix in Matrix Market format
│   ├── filtered_feature_bc_matrix.h5   # Count matrix in HDF5 format
│   ├── spatial/                        # Tissue image, scalefactors, alignment JSONs
│   ├── analysis/                       # Clustering, PCA, UMAP, differential expression
│   └── cloupe.cloupe                   # Interactive Loupe Browser file
└── web_summary.html                    # QC report with sequencing & mapping metrics
```

#### Multi-resolution binning selection strategy

* **2µm (`square_002um`)**: Represents 1 native capture spot. Provides sub-cellular spatial resolution but yields highly sparse UMI counts. Recommended only when investigating sub-cellular transcript localization or fine anatomical structures.
* **8µm (`square_008um`)**: Aggregates a 4×4 grid of 16 spots. Approximates single-cell dimensions with balanced sequencing depth. Ideal for detailed cell-level profiling.
* **16µm (`square_016um`)**: Aggregates an 8×8 grid of 64 spots. Captures multi-cellular neighborhoods with deeper sequencing coverage. Recommended for regional tissue architecture and primary exploratory analyses.

#### Quality Control metrics (`web_summary.html`)

The HTML summary output provides critical quality gates across three categories:

1. **Sequencing Metrics**: Q30 bases in RNA reads, sequencing saturation, total read count.
2. **Mapping Metrics**: Percentage of reads mapped to genome (∼85%), mapped confidently to transcriptome (∼72-74%), intronic regions (∼1.6%), and intergenic regions (∼5.9%).
3. **Spot/Bin Metrics**: Total number of spots under tissue, mean/median reads per bin, and median UMI/gene counts per bin (e.g., mean 849.2 reads and 226.4 UMIs per 8µm bin).

#### Space Ranger cell segmentation outputs (Visium HD)

When cell segmentation algorithms (e.g., StarDist, Bin2Cell, SMURF) are integrated into Space Ranger, additional polygon and mapping files are generated:

* `cell_segmentations.geojson`: Vector polygon boundaries of detected cells.
* `nucleus_segmentations.geojson`: Vector polygon boundaries of detected nuclei.
* `graphclust_annotated_*_segmentations.geojson`: Cell/nucleus polygons annotated with graph clusters.
* `filtered_feature_cell_matrix*`: Cell-level (rather than bin-level) expression count matrix.
* `*_barcode_mappings.parquet`: Critical lookup table mapping 2µm binned spots to segmented cell IDs.

### Xenium Ranger / Onboard analysis (10x Xenium)

Unlike Visium HD, Xenium is an imaging-based platform that processes fluorescent imaging cycles onboard the instrument (XOA) or via `xeniumranger`.

* `transcripts.csv.gz`: Table containing molecule-level spatial coordinates ($x$, $y$ in physical microns), gene identity, quality score (Q-score), and assigned cell ID for every detected RNA molecule.
* `cells.csv.gz`: Metadata table for all segmented cells, recording centroid $x$/$y$ coordinates, cell area, nuclear area, and total UMI counts.
* `cell_feature_matrix.h5`: Sparse count matrix (cells × genes) in HDF5 format.
* `nucleus_boundaries.csv.gz` / `cell_boundaries.csv.gz`: Coordinate vertices defining cell and nuclear polygon boundaries.
* `morphology.ome.tif`: Pyramidal OME-TIFF containing full-resolution DAPI and morphology stains.

_Remind, in 10x Xenium, UMI counts are not derived by measuring raw fluorescence intensity; rather, they are generated by decoding and counting individual RNA molecules as discrete spatial points across automated cycles of imaging and cell segmentation._

### Key conceptual differences: Visium HD vs. Xenium

| Feature | Visium HD (Sequencing-Based) | Xenium (Imaging-Based) |
| :--- | :--- | :--- |
| **Spatial Coordinate System** | Binned square grid (2µm, 8µm, 16µm) or spot array | Physical microns ($x, y$ inµm) |
| **Raw Output Unit** | Sequencing reads mapped to barcoded spots | Individual decoded transcript locations |
| **Spatial Reference Folder** | Requires `spatial/` folder with `scalefactors_json.json` and alignment PNGs | No `spatial/` folder or scalefactors; native micron alignment |
| **Gene Coverage** | Whole-transcriptome unbiased capture | Targeted gene panel (∼280--5,000 genes) |
| **Cell Isolation Mode** | Spatial binning or post hoc image/expression segmentation | Direct onboard cell segmentation (DAPI + boundary stains) |

## Bioconductor spatial object hierarchy

In the R/Bioconductor ecosystem, spatial data representations have evolved by extending established single-cell object architecture.

Rather than creating rigid, static data containers, Bioconductor structures single-cell and spatial objects as modular, extendable containers. An object starts with raw measurements and accumulates metadata, normalized values, and low-dimensional embeddings as computational functions are executed.

```text
===============================================================================================
                      BIOCONDUCTOR OBJECT HIERARCHY & SLOT STRUCTURE
===============================================================================================
SingleCellExperiment (SCE) [Standard scRNA-seq Core]
│
├── Assays / Matrices     : assay(sce, "counts")     → Raw UMI count matrix
│                           assay(sce, "logcounts")  → Normalized expression matrix
├── Feature Annotations   : rowData(sce)             → Gene symbols, Ensembl IDs, HVG metrics
│                           rowRanges(sce)           → Genomic coordinates (GRanges)
├── Cell/Sample Metadata  : colData(sce)             → Observation (id, condition, cell type)
│                           colData(sce)             → QC metrics (sum/UMIs, detected, mito %)
├── Normalization Factors : sizeFactors(sce)         → Cell/spot library size factors
└── Dimension Reductions  : reducedDim(sce, "PCA")   → Principal Components
│                           reducedDim(sce, "UMAP")  → UMAP coordinates
│                           reducedDim(sce, "TSNE")  → t-SNE coordinates
│
│ (Inherits ALL SCE slots, assays, metadata, and reducedDims)
↓
SpatialExperiment (SPE) [SRT Grid & Centroid Extension (e.g., 10x Visium)]
│
├── Spatial Metadata      : colData(spe)             → barcode, `in_tissue`, `array_row`
├── Spatial Coordinates   : spatialCoords(spe)       → Matrix of (x, y) spot/cell centroids
└── Histology Image Data  : imgData(spe)             → Image paths, scale factors, resolution
│
│ (Inherits ALL SPE and SCE slots, assays, metadata, and images)
↓
SpatialFeatureExperiment (SFE) [Geospatial Vector & Subcellular Map (e.g., Xenium, MERFISH)]
│
├── Cell/Nucleus Polygons : colGeometries(sfe)       → Cell & nucleus boundary polygons (sf)
├── Transcript Locations  : rowGeometries(sfe)       → Exact (x, y) transcript spot locations
├── Spatial Neighborhood  : colGraphs(sfe)           → k-NN & Delaunay graphs (cells/bins)
├── Tissue Annotations    : annotGeometries(sfe)     → Anatomical regions & tissue borders
│
├── Substructure Graphs   : rowGraphs(sfe)           → Transcript interaction networks
│                           annotGraphs(sfe)         → Structural region network
└── Local Statistics      : localResults(sfe)        → Per-gene Local Moran's I, Getis-Ord Gi*
===============================================================================================
```

### `SingleCellExperiment` (SCE)

`SingleCellExperiment` is the foundational S4 class in Bioconductor for non-spatial single-cell RNA sequencing (scRNA-seq). Rather than remaining static, an SCE object **dynamically expands in place** as execution steps are run across a bioinformatics pipeline:

* **`assays()`**: List of feature × observation expression matrices.
  * **`assay(sce, "counts")`**: Holds raw, unnormalized UMI count data.
  * **`assay(sce, "logcounts")`**: Added during normalization (e.g., via `logNormCounts()`); stores log-transformed, size-factor-adjusted expression values.
* **`colData()`**: A `DataFrame` holding observation-level metadata.
  * **Observation Attributes**: Stores initial experimental descriptors such as sample `id`, experimental `condition`, `batch`, `author`, and annotated `cell_type`.
  * **Appended QC Metrics**: Post-processing tools (such as `calculateQCMetrics()`) append computed quality statistics as additional columns—including `sum` (total UMI count per cell), `detected` (number of expressed genes), and `subsets_mito_percent`.
* **`rowData()` / `rowRanges()`**: Metadata describing features.
  * **`rowData(sce)`**: Stores gene symbols, Ensembl IDs, and feature selection metrics (e.g., Highly Variable Gene variance).
  * **`rowRanges(sce)`**: Stores genomic coordinates as `GRanges` objects.
* **`sizeFactors()`**: A numeric vector storing cell-specific library scaling factors calculated during normalization (e.g., via `computeSumFactors()`).
* **`reducedDims()`**: A list storing low-dimensional coordinate matrices generated downstream—accessed as `reducedDim(sce, "PCA")`, `reducedDim(sce, "UMAP")`, or `reducedDim(sce, "TSNE")`.

> **Structural Limitation**: SCE contains no native slots for physical spatial coordinates, histological image alignments, spatial neighborhood graphs, or vector geometries.

### `SpatialExperiment` (SPE)

`SpatialExperiment` directly inherits all core slots, metadata, and functionality from `SingleCellExperiment` while extending the architecture for grid-based or centroid-based spatial transcriptomics (e.g., standard 10x Visium or unsegmented Visium HD).

* **Inherited SCE Core**: Preserves all `assays` (`counts`, `logcounts`), `colData`, `rowData`, `sizeFactors`, and `reducedDims`.
* **`colData(spe)` Extension**: In addition to cell/sample metadata and QC metrics, `colData(spe)` holds spot/grid spatial metadata such as `barcode`, `in_tissue` indicators, and array positions (`array_row`, `array_col`).
* **`spatialCoords()`**: A numeric matrix storing the physical $(x, y)$ centroid coordinates for every spot or cell centroid.
* **`imgData()`**: A `DataFrame` storing background histological images (H&E or fluorescence), image file paths, scale factors (mapping spatial coordinates to pixel dimensions), and resolution metadata.

> **Structural Limitation**: SPE treats spatial locations as point centroids or regular arrays; it lacks native support for multi-vertex vector shapes (cell/nucleus boundaries), individual subcellular transcript locations, or spatial neighborhood graphs.

### `SpatialFeatureExperiment` (SFE)

`SpatialFeatureExperiment` extends SPE and SCE by integrating **Geospatial Simple Features (`sf`)** and raster engines (`terra`/`EBImage`), turning the object into a multi-layered digital spatial map designed for subcellular imaging technologies (10x Xenium, MERFISH) and segmented Visium HD.

* **Inherited SPE/SCE Core**: Retains all expression matrices (`assays`), observation metadata (`colData`), feature annotations (`rowData`), centroid coordinates (`spatialCoords`), histology images (`imgData`), and low-dimensional embeddings (`reducedDims`).
* **Core Extended SFE Features**:
  - **`colGeometries()` (Cell/Nucleus Polygons)**: Stores multi-vertex vector polygons (`sf` data frames) representing precise cell (`cellSeg`) and nucleus (`nucSeg`) boundary geometries generated by cell segmentation algorithms (Cellpose, Baysor, Proseg).
  - **`rowGeometries()` (Transcript Locations)**: Stores the exact physical $(x, y)$ micron coordinates for every individual detected RNA transcript spot.
  - **`colGraphs()` (Spatial Neighbor Graphs)**: Stores spatial neighborhood topologies ($k$-NN, Delaunay triangulation, or polygon contiguity) as `spdep`/`Listw` objects mapping physical adjacencies between cells or spots.
  - **`annotGeometries()` (Tissue Annotations)**: Stores polygon boundaries for higher-level biological macro-structures independent of cell borders—such as histological layers, tumor cores, necrotic zones, or tissue slide boundaries.
* **Auxiliary SFE Extensions**:
  * **`rowGraphs()` & `annotGraphs()`**: Substructure neighborhood networks mapping transcript point clusters and anatomical region adjacencies.
  * **`localResults()`**: Stores per-gene local spatial statistics (such as Local Moran's $I$, Getis-Ord $G_i^*$, or LOSH) calculated directly on `logcounts` or `reducedDims`.

### Feature & capability comparison matrix

| Feature / Slot | `SingleCellExperiment` (SCE) | `SpatialExperiment` (SPE) | `SpatialFeatureExperiment` (SFE) |
| :--- | :---: | :---: | :---: |
| **Count Matrices (`assays`)** | ✓ | ✓ | ✓ |
| **Observation Metadata (`colData`)** | ✓ | ✓ | ✓ |
| **Feature Annotations (`rowData`)** | ✓ | ✓ | ✓ |
| **Reduced Dimensions (`reducedDims`)** | ✓ | ✓ | ✓ |
| **Centroid Coordinates (`spatialCoords`)** | — | ✓ | ✓ |
| **Tissue Histology Images (`imgData`)** | — | ✓ | ✓ |
| **Cell/Nucleus Boundaries (`colGeometries`)** | — | — | **✓** |
| **Transcript Locations (`rowGeometries`)** | — | — | **✓** |
| **Spatial Neighbor Graphs (`colGraphs`)** | — | — | **✓** |
| **Tissue Annotations (`annotGeometries`)** | — | — | **✓** |
| **Typical Data Application** | scRNA-seq | Visium HD (Unsegmented) | 10x Xenium & Segmented Visium HD |


## Conclusion & next steps

Before jumping into the code, here are the core conceptual takeaways from this module:

* **Pipeline Divergence Dictates Data Shape**: Space Ranger processes sequencing-based Visium and Visium HD data into multi-resolution spatial bins (e.g., 8µm or 16µm grids) that require alignment to histological reference images. Conversely, imaging-based platforms processed by Xenium Ranger output exact subcellular transcript coordinates and cell geometries mapped directly in physical microns.
* **Dynamic Object Hierarchy**: Bioconductor utilizes modular structures that expand based on spatial complexity. The `SpatialExperiment` (SPE) framework builds on standard single-cell architecture by adding centroid coordinate matrices and tissue image data, making it the standard choice for unsegmented Visium grids.
* **Vector Geometry Integration**: Analyzing high-resolution spatial data requires upgrading to the `SpatialFeatureExperiment` (SFE) framework. SFE natively incorporates geospatial simple features (`sf`), allowing it to store multi-vertex cell and nucleus polygons, individual transcript coordinate points, and spatial neighborhood topologies.
* **Spatial Context Drives Downstream Analysis**: Properly mapping hardware outputs into these structured R/Bioconductor slots is a prerequisite for robust analysis. Understanding where data resides—such as finding raw counts in `assays` or cell geometries in `colGeometries`—is essential for the analytical workflows you are about to execute.

Now that you understand the underlying data structures, the following workflow serves as a high-level roadmap for the practical Jupyter Notebook sessions spanning the rest of the course:

```text
┌─────────────────────────────────────────────────────────────────────────────────┐
│                        JUPYTER NOTEBOOK WORKFLOW ROADMAP                        │
├─────────────────────────────────────────────────────────────────────────────────┤
│ 1. Data Ingestion     → Read Space Ranger / Xenium outputs into SPE/SFE         │
│ 2. Slot Interrogation → Inspect counts, spatialCoords, colGeometries, imgData   │
│ 3. Spatial QC & Filt  → Filter low-count bins/cells & compute spatial metrics   │
│ 4. Spatial Norm       → Execute logNormCounts() and evaluate spatial trends     │
│ 5. Intermed. Process  → Perform feature selection & dimensionality reduction    │
│ 6. Spatial Clustering → Identify spatial domains and distinct tissue niches     │
│ 7. Spatial Statistics → Build spatial graphs, run Moran's I/LISA, find SVGs     │
│ 8. Multi-Scale/Sample → Integrate diverse datasets and compare resolutions      │
└─────────────────────────────────────────────────────────────────────────────────┘
```
