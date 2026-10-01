# Cell Segmentation & Transcript Assignment in ST

## The Assignment Problem

In spatial transcriptomics (SRT), measuring RNA abundance alongside physical coordinates helps unravel collaborative cell organizations and cell-cell communications. However, extracting biological meaning requires accurately assigning raw transcripts to individual cells. 

The approach to this "assignment problem" depends entirely on the technology used:

* **Imaging-Based ST (Single-Molecule):** Generates true transcript coordinates that require sophisticated **image segmentation** to draw boundaries around distinct single cells.
* **Sequencing-Based ST (Spot/Grid):** Captures transcripts inside larger areas (spots or bins) that often contain multiple cells, requiring mathematical **deconvolution** to estimate the underlying single-cell resolution.

## 1. True Image Segmentation (Imaging-Based Methods)

Because we have exact $(x,y)$ coordinates for every RNA molecule in imaging-based ST (like Xenium or MERFISH), the goal is to draw accurate cell boundaries (polygons) on the tissue image and count the transcripts inside them.

* **Evolution of Algorithms:** Image segmentation has evolved from classical morphological approaches (e.g., Watershed) to modern deep learning and foundation models (e.g., U-Net, StarDist, Cellpose, and CellSAM).
* **Multi-Stain Approaches (e.g., 10x Xenium Onboard Analysis):** Modern pipelines use multiple stains to define cell boundaries accurately. For example, Xenium utilizes DAPI for nuclear staining, ATP1A1/CD45/E-Cadherin for boundary staining, and 18S RNA + alphaSMA/vimentin for interior staining.
* **Cell Expansion:** When strict membrane markers are missing or noisy, algorithms often identify the nucleus first, then mathematically expand boundaries outwards by a set distance (e.g., 5µm or 15µm) to capture cytoplasmic transcripts.

### Advanced Image & Transcript-Assisted Algorithms

While earlier methods relied solely on morphology, modern tools incorporate the localized transcripts themselves to build or refine boundaries:

* **Baysor (Bayesian Segmentation):** Recognized as the first major method to utilize a "segmentation-free" approach based entirely on the transcriptional composition of cells. It can optionally integrate auxiliary stains by combining transcript colocalization with precalculated nuclear or membrane segmentation. *(Note: Newer unnamed GNN-based tools are emerging that resolve ambiguous border transcripts even faster than Baysor using GPU processing).*
* **Proseg (Probabilistic Cell Segmentation):** Operates by calculating the specific set of cell shapes that best explain the actual physical locations where observed transcripts landed. It initializes cell morphologies based on nuclear segmentation and iteratively expands them for the best fit. Crucially, it models "leaky" cells (a challenge particularly noted in Xenium data) by adding a mathematical prior on the distance transcripts can migrate.
* **SMURF (Segmentation and Manifold Unrolling Framework):** Increases Unique Molecular Identifier (UMI) counts by using concentric expansion to add nearby spatial spots to a nucleus's gene expression matrix. During each iteration, the framework checks whether the newly expanded cells increase the cosine similarity among cells in the same cluster, repeating until clustering stabilizes (measured by normalized mutual information). It deconvolves transcript counts from shared spots to ensure every transcript is assigned to at most one cell, requiring that the sum of cell-type proportions matches the normalized gene expression pattern for that specific spot.

## 2. Spatial Binning & Deconvolution (Sequencing-Based Methods)

Because sequencing-based technologies capture transcripts on barcoded arrays rather than imaging cell boundaries directly, "segmentation" takes the form of grid binning and mathematical deconvolution.

* **Spatial Binning (e.g., Space Ranger for Visium HD):** The standard Space Ranger pipeline aligns raw sequencing reads to physical tissue coordinates. For continuous-array technologies like Visium HD, it outputs expression matrices binned into arbitrary square grids (e.g., 2×2µm, 8×8µm, and 16×16µm squares) rather than perfect single-cell polygons.
* **Spot Deconvolution:** Because standard sequencing spots (like the 55µm Visium V1/V2 spots) contain mixtures of multiple cells, we must infer the exact cell-type composition of each spot.
* **Applying RCTD (Robust Cell Type Decomposition):** Deconvolution is achieved by integrating external single-cell reference databases. For example, using RCTD with a highly annotated single-cell atlas (e.g., Bone Marrow) allows the algorithm to learn the distinct transcriptional profiles of cell types, map those profiles back onto the spatial spots, and mathematically deduce the proportions of each cell type present in that mixed location.

## 3. Correcting for Noise, Spillover, and Errors

A major analytical hurdle in ST—even with perfect image segmentation—is "spillover" or "bleeding." This occurs due to the physical diffusion of RNA molecules to neighboring cells (e.g., transcriptionally active malignant cells contaminating the reads of adjacent healthy cells) or overlapping 3D cell boundaries in a 2D slice.

Advanced algorithms like ResolVI are used to model and correct these diffusion issues.

### ResolVI

* **Mechanism:** ResolVI addresses noise and diffusion bias using a variational autoencoder (VAE) latent variable model.
* **Correction:** It estimates the observed expression as a sum of three factors: 
    1. The *true* expression of the cell.
    2. The expression induced by *neighboring cells* (diffusion/spillover).
    3. Unspecific *background* expression.

### SPLIT (Spatial Purification of Layered Intracellular Transcripts)

* **Mechanism:** An R package that first utilizes RCTD with a reference scRNA-seq atlas for the initial decomposition of cells, evaluating the similarity of each cell's transcriptome to provided cell types.
* **Identification:** It specifically identifies "doublet-cells," which represent potential RNA spilling from neighbors.
* **Correction:** It corrects raw counts by "purifying" the cell. It keeps only the fraction of reads corresponding to the primary cell type (relative to the sum of the primary and secondary types), mathematically stripping away the contaminating secondary transcripts.

### cellAdmix

* **Mechanism:** Addresses segmentation errors by defining an "admixture score" using a Bayesian model. This estimates contamination for each gene while strictly taking cell type adjacency into account.
* **Factorization:** It computes neighborhood composition values by constructing a k-NN graph per cell and factorizes them using Non-negative Matrix Factorization (NMF) to recover pure expression profiles.
* **Correction:** Employs a Conditional Random Field (CRF) to assign individual molecules to the recovered factors. Through an automated "cell-bridging" approach, individual molecules whose factor assignments do not match their host cell's type are identified as spillover and removed.

## 4. The Snakemake Xenium Pipeline (Downstream Integration)

In practice, handling the massive outputs from machines like the 10x Xenium (which takes ~1 day of runtime for a whole slide) is orchestrated using workflow managers like **Snakemake**. Based on the standard Xenium analytical flowchart, the pipeline architecture follows these specific branching steps:

### Phase 1: Initial Output & Resegmentation

* **Resegmentation:** The pipeline begins with the initial Xenium Ranger output, which immediately undergoes a resegmentation step (often via a `xeniumranger` container) to produce a refined "resegmented bundle."
* **Image Alignment:** From this bundle, an independent branch handles the **Registration of H&E** images to ensure morphological stains perfectly overlay the transcript coordinates.

### Phase 2: Advanced Segmentation & Formatting

From the resegmented bundle, the data is pushed through advanced transcript-assisted algorithms and formatted for specific programming environments:

* **Proseg Fork:** Runs probabilistic segmentation, outputting directly to **SpatialFeatureExperiment (SFE)** objects for R-based downstream analysis.
* **Segger Fork:** Runs an alternative segmentation algorithm, outputting to the **spatialdata** format (a modern framework for Python-based spatial omics).

### Phase 3: Downstream Processing Branches

Once loaded into `spatialdata` or SFE objects, the workflow splits into specific functional modules:

1. **Sample Cropping:** Cropping samples to remove edge artifacts or isolate regions of interest (if the initial Xenium Explorer polygons were too loose).
2. **Integration:** Merging multiple tissue samples together using deep generative models like **scvi** to integrate the data into a unified latent space.
3. **Annotation & Correction:** A sequential chain where **RCTD** maps a matching reference atlas onto the spatial data, followed closely by **ResolVI** to correct for spillover and noise, resulting in the final annotated sample. *(Note: SPLIT correction is also sometimes evaluated at this stage).*

## 5. Take-Home Messages

* **Crucial but Imperfect:** Cell segmentation is vitally important but almost never perfect. It must be performed at the very beginning of the analytical pipeline (e.g., during Xenium onboard analysis) because it dictates the quality of all downstream steps.
* **Transcript-Assisted Segmentation:** Using the transcripts themselves can greatly enhance pure image-based segmentation, though one must account for technology-specific biases like RNA diffusion and spillover.
* **Tool Combinations are Key:** For the best results, you almost always need to combine multiple specialized tools rather than relying on a single algorithm (e.g., chaining SMURF + RCTD + ResolVI).
* **Ecosystem Bottlenecks:** Despite a comprehensive and growing list of segmentation tools, there is still a noticeable lack of smooth interoperability between the R and Python spatial ecosystems, requiring careful pipeline orchestration (like Snakemake).
