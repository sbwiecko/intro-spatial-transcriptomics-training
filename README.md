[![License: CC BY 4.0](https://img.shields.io/badge/License-CC_BY_4.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)
[![Forked from sib-swiss](https://img.shields.io/badge/Forked%20from-sib--swiss-blue)](https://github.com/sib-swiss/intro-spatial-transcriptomics-training)

# Intro to Spatial Transcriptomics: A Wet-Lab Scientist's Study Fork

Welcome! This repository is an annotated, restructured study fork of the [SIB Swiss Institute of Bioinformatics Introduction to Spatial Transcriptomics Data Analysis Course](https://github.com/sib-swiss/intro-spatial-transcriptomics-training).

I come from a **wet-lab background** and am investing significant effort into transitioning toward computational biology. To make the learning curve manageable and fully focus on the code and biology, **I simplified the project architecture**. I removed the heavy Docker containers in favor of standalone **Jupyter Notebooks (`.ipynb`)** so you can easily clone the repo and practice the code yourself.

---

## What You Will Find in This Fork

* **Streamlined Interactive Format:** Converted from the original `.qmd` files into pure Jupyter Notebooks (`.ipynb`) (via `quarto convert`, configured to use the R kernel via `jupyter: ir` in the YAML front matter) for direct, interactive execution.
* **Wet-Lab Context & Concept Deconstruction:** Line-by-line notes translating abstract computational steps into biological meaning.
* **Troubleshooting & Technical Deep Dives:** Detailed reflections and extra research on tricky parameters, spatial Bioconductor quirks, and statistical assumptions.
* **Refactored & Annotated Code:** Code cells with additional comments explaining the *why* behind data transformations.
* **Beginner-Friendly Focus:** Written specifically for researchers without formal bioinformatics degrees who want to understand spatial transcriptomics step by step.

### Key Code Optimizations in this Fork

Over the course of the notebooks, I adapted the original base R code to implement modern data science and tidyverse workflows. If you are learning R, you will find highly optimized, real-world solutions for:

* **Modern Biological Logic (`dplyr::case_when`):** Transitioned classic boolean matrix subsetting into elegant, readable `case_when()` pipelines. This makes logical filtering for spatial regions or cell types much more robust.
* **Spatial Object Wrangling:** Seamlessly manipulating `SpatialExperiment` and `SpatialFeatureExperiment` metadata directly, ensuring that spatial coordinates (e.g., Visium spots or Xenium segmented cells) and gene expression structures stay perfectly aligned.
* **Efficient I/O & Memory Management:** Refactored how large spatial matrices are stored and saved, significantly streamlining downstream loading and keeping the global environment clean when working with massive spatial datasets.
* **Parameterizing Shared Scripts & Local Wrappers:** Updated shared plotting functions to accept flexible aesthetic defaults (fonts, line widths), ensuring backward compatibility across all notebooks. Paired this with local function wrappers in notebooks to abstract away repetitive code when generating complex spatial plots.

---

## Analytical Workflow Covered

1. **Spatial Object Foundations:** Understanding and building `SpatialExperiment` and `SpatialFeatureExperiment` structures using Visium HD and Xenium spatial datasets.
2. **Quality Control & Preprocessing:** Executing spatial-specific QC using `SpotSweeper` and image-based QC metrics, followed by platform-specific normalization (Visium vs. Xenium).
3. **Feature Selection & Dimensionality Reduction:** Extracting informative features while accounting for spatial constraints, and performing dimensionality reduction.
4. **Clustering, Annotation & Spatial Statistics:** Applying spatial clustering, identifying marker genes, annotating clusters, and measuring spatial autocorrelation and differential spatial patterns.

---

## Learning Objectives

By working through these notebooks, you will be able to:

* **Understand SRT Principles:** Explain the applications and differences between sequencing-based and imaging-based spatially-resolved transcriptomics (SRT).
* **Navigate Limitations:** Identify potential pitfalls and limitations of SRT experiments and analysis workflows.
* **Process Spatial Data:** Assess and interpret raw spatial outputs and metadata files, understanding their structure and relevance for downstream analysis.
* **Master the Analytical Pipeline:** Apply crucial steps such as quality control, feature selection, dimensionality reduction, and differential gene expression/spatial pattern analysis to SRT data.
* **Apply Spatial Statistics:** Understand various spatial statistics and their application to biological questions.
* **Handle Complex Data:** Define applications for cell segmentation and use frequently-used methods to analyze multi-sample SRT experiments.

> **📝 Note on the Exercises:** The original course was designed for live delivery with lectures, polls, and group discussions. In this self-guided study fork, I have adapted the experience by integrating the core concepts directly into notes and notebooks. Each notebook begins with its own targeted learning objectives. I recommend reading through my explanatory notes and executing the code cells sequentially, thinking about the code before checking the provided outputs!

---

## Prerequisites & Scope

> **Important Scope Note:** This repository focuses specifically on **downstream spatial transcriptomics analysis and biological interpretation**. It does **not** cover raw sequencing/imaging preprocessing (e.g., SpaceRanger, FASTQ quality control, raw image cell segmentation algorithms) or feature engineering from scratch.

To get the most out of these notebooks, familiarity with the following concepts is recommended:

* **Basic R & Bioconductor:** Working with common data structures (`data.frame`, matrices) and core Bioconductor single-cell classes (`SingleCellExperiment`).
* **Single-Cell RNA-seq Foundations:** Understanding basic principles of cell/spot-based transcriptomics, count matrices, and dimensionality reduction (PCA/UMAP).
* **Exploratory Data Analysis (EDA):** Experience diagnosing batch effects, sample clustering, and data normalization before spatial integration.

### Looking to build foundational skills or explore more topics?

* **Biology-Informed Multi-Omics Integration:** Check out my companion remix of the SIB course on multi-omics integration at [sbwiecko/biology-informed-multiomics-training](https://github.com/sbwiecko/biology-informed-multiomics-training), which shares many of the underlying analytical methods and R/Bioconductor toolkits.
* **Statistical Foundations & Coding with Python:** For those wanting to master coding from the ground up using biostatistics as the vehicle, check out my book [*Intuitive Biostatistics with Python* (Oxford University Press)](https://global.oup.com/academic/product/intuitive-biostatistics-with-python-9780197845035?lang=en&cc=fr).
* **End-to-End RNA-seq & Preprocessing:** Check out my repository [sbwiecko/RNAseq_UPenn](https://github.com/sbwiecko/RNAseq_UPenn), rewritten alongside the excellent [UPenn DIY Transcriptomics course](https://diytranscriptomics.com/).
* **SIB Training Catalog:** Explore the full lineup of world-class courses designed by the exceptional SIB Training Team: [SIB Swiss Institute of Bioinformatics Courses & Training](https://www.sib.swiss/training/learning-paths)
* **EDA & Statistical Exploration:** [Modern Statistics for Modern Biology](https://www.huber.embl.de/msmb/) by Susan Holmes and Wolfgang Huber
* **Free Bioinformatics Resources:** Discover a curated collection of free bioinformatics training materials and resources at [GLITTR](https://glittr.org/)

---

## Setup & Running the Notebooks

**Prerequisite:** You must have **R** installed on your system.

No Docker or Quarto setup is required. The project relies on `renv` to manage dependencies. Since these are standard Jupyter Notebooks (`.ipynb`) with an R kernel, you can run them locally using **VS Code**.

### 1. Clone the repository

```bash
git clone https://github.com/sbwiecko/intro-spatial-transcriptomics-training.git
cd intro-spatial-transcriptomics-training
```

### 2. Restore the R environment

Open R in your terminal, console, or IDE (for example, by running `"C:\Program Files\R\R-4.6.1\bin\R.exe"` or `radian`) to restore dependencies via `renv`:

```R
install.packages("renv")
renv::restore()
```

> **💡 How was this environment built? (For the curious)**  
> I took the opportunity to upgrade the original course to a fresh version of R (4.6.x) and the latest Bioconductor packages. To do this without manually guessing dependencies, I read the original `renv.lock` JSON file, extracted the raw package list, and bulk-installed their newest versions via `BiocManager::install()`. Once everything was tested and running smoothly, I captured this modernized state into the new `renv.lock` file you see in the repository.  
> **Don't worry—you do not need to do this!** The hard work is done. Just execute `renv::restore()` as shown above to perfectly mirror this upgraded environment.

### 3. Running in VS Code

1. **Configure VS Code Working Directory (Crucial):** To guarantee that Jupyter always activates `.Rprofile` and locates your `renv` library, force VS Code to run kernels from the project root.
   - Open Settings (`Ctrl + ,` or `Cmd + ,`).
   - Search for: `Jupyter Notebook File Root`.
   - Set the value to: `${workspaceFolder}`.
2. **Register the Project Kernel:**
   ```R
   renv::install("IRkernel")
   IRkernel::installspec(name = "ir_spatial", displayname = "R (spatial-renv)")
   ```
3. Open any `.ipynb` notebook.
4. Click **Select Kernel** in the top-right -> **Jupyter Kernel...** -> **R (spatial-renv)**.

> **Optional check:** You can always run `.libPaths()` in a notebook's first cell to verify that packages are loading from the project's local `./renv/library`.

---

## Connect & Discuss

If you are also bridging the gap between wet-lab biology and bioinformatics, feel free to explore the notebooks, open an issue, or ask a question. Feedback, corrections, and discussions are always welcome!

---

## Credits & Upstream Authors

This repository is built upon the open course materials designed by the **SIB Swiss Institute of Bioinformatics**:

* **Original Course Site:** [https://sib-swiss.github.io/intro-spatial-transcriptomics-training/](https://sib-swiss.github.io/intro-spatial-transcriptomics-training/)
* **Original Repository:** [sib-swiss/intro-spatial-transcriptomics-training](https://github.com/sib-swiss/intro-spatial-transcriptomics-training)
* **Authors:** Deepak Tanwar, Joana Carlevaro-Fita, Martin Emons, Peiying Cai, Samuel Gunz, Mark Robinson, Ivan Berest, Julien Roux, Geert van Geest
