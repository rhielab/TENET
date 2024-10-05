# TENET
### R package for TENET (Tracing regulatory Element Networks using Epigenetic Traits) to identify key transcription factors

TENET identifies key transcription factors (TFs) and regulatory elements (REs) linked to a specific cell type by finding significantly correlated differences in gene expression and RE methylation between case and control input datasets, and identifying the top genes by number of significant RE DNA methylation site links. It also includes many additional tools to aid in visualization and analysis of the results, including plots displaying and comparing methylation and expression data and RE DNA methylation site link counts, survival analysis, TF motif searching in the vicinity of linked RE DNA methylation sites, custom TAD and peak overlap analysis, and UCSC Genome Browser track file generation. A utility function is also provided to download methylation, expression, and patient survival data from The Cancer Genome Atlas (TCGA) for use in TENET or other analyses.

# Acquiring and installing TENET

To use TENET, users will need to install the base package as well as its associated example experiment data package, [TENET.ExperimentHub](https://github.com/rhielab/TENET.ExperimentHub). **Note:** TENET.ExperimentHub will install automatically when TENET is installed from Bioconductor.

TENET also uses annotation datasets hosted in the Bioconductor AnnotationHub database. These datasets will be automatically loaded from AnnotationHub when necessary. They are also available separately via the [TENET.AnnotationHub](https://github.com/rhielab/TENET.AnnotationHub) package. It is not necessary to install the TENET.AnnotationHub package to use TENET's functions.

R 4.4 or a newer version is required.

On Ubuntu 22.04, installation was successful in a fresh R environment after [adding the official R Ubuntu repository](https://cran.r-project.org/bin/linux/ubuntu/) and running:

`sudo apt-get install r-base-core r-base-dev libcurl4-openssl-dev libfreetype6-dev libfribidi-dev libfontconfig1-dev libharfbuzz-dev libtiff5-dev libxml2-dev`

No dependencies other than R are required on macOS or Windows.

Two versions of this package are available.

To install the stable version from Bioconductor, start R and run:

```r
## Install BiocManager, which is required to install packages from Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

BiocManager::install(version = "devel")
BiocManager::install("TENET")
```

The development version containing the most recent updates is available from [our GitHub repository](https://github.com/rhielab/TENET).
To install the development version from GitHub, start R and run:

```r
## Install prerequisite packages to install the development version from GitHub
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
}

BiocManager::install(version = "devel")
BiocManager::install("rhielab/TENET.ExperimentHub")
BiocManager::install("rhielab/TENET")
```

# Loading TENET

To load the TENET package, start R and run:

```r
library(TENET)
```

To load the TENET.ExperimentHub package, start R and run:

```r
library(TENET.ExperimentHub)
```

To load the TENET.AnnotationHub package if it has been installed, start R and
run:

```r
library(TENET.AnnotationHub)
```

# Running TENET without internet access

Some TENET features and examples download datasets from the internet if they have not already been cached. You must run `TENETCacheAllData()` once while connected to the internet before using these TENET features or examples without internet access (for example, on HPC cluster nodes). See the documentation for `TENETCacheAllData` for more information.

## [Package documentation and vignette](https://bioconductor.org/packages/devel/bioc/html/TENET.html)
