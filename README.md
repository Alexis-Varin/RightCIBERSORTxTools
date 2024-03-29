
<!-- README.md is generated from README.Rmd. Please edit that file -->

# RightCIBERSORTxTools

<!-- badges: start -->
<!-- badges: end -->

This package provides a set of tools to generate CIBERSORTx input files
from RNA-seq data, namely a Reference Matrix (a Seurat object’s RNA
counts in a tab delimited txt file with cell identity as column name)
and a Mixture File (a tab delimited txt file with each column containing
CPM of a bulk RNA-seq or microarray experiment) all in a single
function.

## Installation

    devtools::install_github("Alexis-Varin/RightCIBERSORTxTools")

## Good Practice for CIBERSORTx use

- Your storage capacity on the CIBERSORTx web portal has a limit of 1 GB
  (please note that not just uploaded files take disk space but also job
  results). By experience, if your Seurat object has more than 25,000
  cells by 20,000 features (or any combination that would produce a
  matrix length over 500,000,000) your Reference Matrix will be bigger
  than 1 GB. This is the reason for the various subset and downsample
  arguments the function has to offer. You can use them to reduce the
  size of your Seurat object. I also included a failsafe option to stop
  the function if, before starting to extract the RNA counts, it detects
  that the Reference Matrix will be over the specified size limit. This
  prevents you from waiting a significant amount of time for nothing.

- While you can always downsample right before building the Reference
  Matrix with the provided parameters of the Reference Matrix Builder, I
  strongly recommend to do it on your Seurat object right after QC and
  before any downstream analysis (so before normalization, integration,
  UMAP, FindMarkers etc) as downsampling after will significantly alter
  your analysis, especially your differentially expressed genes. I
  therefore provide a parameter (check.size) to estimate the size of the
  Reference Matrix that would be built from your Seurat object and
  suggest a number of cells to downsample from it. You can then use the
  suggested number of cells to downsample your Seurat object before any
  downstream analysis or even let the function do it for you.

## Reference_Matrix_Builder

### Description

This function builds a Reference Matrix for CIBERSORTx from a Seurat
object.

### Dependencies

- Seurat version \>= 4.0.0
- SeuratObject version \>= 4.0.0
- data.table

### Usage

    Reference_Matrix_Builder = function(
        seurat_object,
        assay = "RNA",
        layer = "counts",
        ident.1 = NULL,
        ident.2 = NULL,
        double.ident = TRUE,
        sep.double.ident = "_",
        reverse.double.ident = FALSE,
        clusters.1 = NULL,
        clusters.1.invert = FALSE,
        clusters.2 = NULL,
        clusters.2.invert = FALSE,
        downsample.object.first = NULL,
        downsample.object.last = NULL,
        downsample.cluster = NULL,
        automatic.downsample = FALSE,
        check.size = FALSE,
        max.matrix.size = 900,
        cell.barcodes = FALSE,
        file.name = "Reference_Matrix",
        file.format = "txt",
        file.sep = "\t",
        path = NULL,
        write = TRUE,
        verbose = TRUE
    )

### Arguments

| Argument                    | Definition                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   |
|-----------------------------|:---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **seurat_object**           | A Seurat object, created with Seurat v4 or v5.                                                                                                                                                                                                                                                                                                                                                                                                                                               |
| **assay**                   | Character. If the Seurat object contains multiple RNA assays, you may specify which one to use (for example “RNA2” if you have created a second RNA assay you named “RNA2”. See Seurat v5 vignettes for more information). You may also use another assay such as SCT if you want to extract an expression matrix of this assay for projects other than CIBERSORTx. Must be RNA assay for CIBERSORTx.                                                                                        |
| **layer**                   | Character. Formerly known as slot. You can specify a layer such as “data” (normalized counts) or any other present in the Seurat object’s assay if you want to extract the expression matrix of this layer for projects other than CIBERSORTx. If you have split layers the function will always join them before extracting the expression matrix. Must be counts layer for CIBERSORTx.                                                                                                     |
| **ident.1**                 | Character. The first identity to use (for example “seurat_clusters” or “orig.ident”). Leave NULL if you want to use current identity (active.ident).                                                                                                                                                                                                                                                                                                                                         |
| **ident.2**                 | Character. The second identity to use. This will subset each ident.1 cluster with ident.2 identities (for example if ident.1 is “seurat_clusters” and ident.2 “orig.ident” cluster 0 will be divided into several clusters 0 each corresponding to an orig.ident, 0_sample1, 0_sample2 etc). Leave NULL if you don’t want to use a second identity.                                                                                                                                          |
| **double.ident**            | Logical. If TRUE, each cell barcode will be renamed with a double identity (ident.1_ident.2). If FALSE, each cell barcode will be renamed with ident.1. Ignored if ident.2 is NULL.                                                                                                                                                                                                                                                                                                          |
| **sep.double.ident**        | Character. The separator to use between ident.1 and ident.2 if double.ident = TRUE. Ignored if ident.2 is NULL.                                                                                                                                                                                                                                                                                                                                                                              |
| **reverse.double.ident**    | Logical. If TRUE, the order of ident.1 and ident.2 will be reversed when creating the double identity (ident.2_ident.1). Ignored if ident.2 is NULL.                                                                                                                                                                                                                                                                                                                                         |
| **clusters.1**              | Character. The clusters to subset from ident.1. Leave NULL if you want all clusters.                                                                                                                                                                                                                                                                                                                                                                                                         |
| **clusters.1.invert**       | Logical. If TRUE, invert the selection from clusters.1 (if clusters.1 = c(“0”,“1”,“2”) and clusters.1.invert = TRUE, all clusters except 0, 1 and 2 will be selected).                                                                                                                                                                                                                                                                                                                       |
| **clusters.2**              | Character. The clusters to subset from ident.2. Leave NULL if you want all clusters.                                                                                                                                                                                                                                                                                                                                                                                                         |
| **clusters.2.invert**       | Logical. If TRUE, invert the selection from clusters.2.                                                                                                                                                                                                                                                                                                                                                                                                                                      |
| **downsample.object.first** | Numeric. The number of cells to downsample from the entire Seurat object (not each cluster) before subsetting clusters.1 and clusters.2. Leave NULL if you don’t want to downsample.                                                                                                                                                                                                                                                                                                         |
| **downsample.object.last**  | Numeric. The number of cells to downsample from the entire Seurat object after subsetting clusters.1 and clusters.2. Leave NULL if you don’t want to downsample.                                                                                                                                                                                                                                                                                                                             |
| **downsample.cluster**      | Numeric. The number of cells to downsample from each cluster of ident.1. Will be performed after downsample.object.last. Leave NULL if you don’t want to downsample.                                                                                                                                                                                                                                                                                                                         |
| **automatic.downsample**    | Logical. If TRUE, automatically downsample the Seurat object so that the Reference Matrix file written to disk would be just under the max.matrix.size limit (empirical). Performed after all other subset and downsample parameters. Ignored if check.size = FALSE and write = FALSE. Please report an issue if you see a significant difference between the file size written to disk and max.matrix.size (for example max.matrix.size is set to 200 MB but the file is 400 MB or 100 MB). |
| **check.size**              | Logical. If TRUE, will return the Seurat object and print the estimated size of the Reference Matrix file written to disk and the number of cells to downsample. If automatic.downsample = TRUE will perform the downsample and return the Seurat object.                                                                                                                                                                                                                                    |
| **max.matrix.size**         | Numeric. The maximum size of the Reference Matrix file written to disk in MB. Will stop the function if the Reference Matrix file written to disk is estimated to be over this limit, or if automatic.downsample = TRUE will downsample the Seurat object instead so that the Reference Matrix output file is under the size limit. Ignored if check.size = FALSE and write = FALSE.                                                                                                         |
| **cell.barcodes**           | Logical. If TRUE, keep the cell barcodes and does not rename with cell identities, if you want to extract the expression matrix for projects other than CIBERSORTx. Must be FALSE for CIBERSORTx.                                                                                                                                                                                                                                                                                            |
| **file.name**               | Character. The name of the Reference Matrix file written to disk. Must not contain any space. Ignored if write = FALSE.                                                                                                                                                                                                                                                                                                                                                                      |
| **file.format**             | Character. The format of the Reference Matrix file written to disk. Must be txt or tsv for CIBERSORTx but you can also specify csv for example if you want to extract the expression matrix for projects other than CIBERSORTx. Accept any format the data.table::fwrite function would accept. Ignored if write = FALSE.                                                                                                                                                                    |
| **file.sep**                | Character. The separator to use in the Reference Matrix file written to disk. Must be Tabulation for CIBERSORTx but you can also specify a comma for example if you want to extract the expression matrix for projects other than CIBERSORTx. Accept any separator the data.table::fwrite function would accept. Ignored if write = FALSE.                                                                                                                                                   |
| **path**                    | Character. The path to write the Reference Matrix into. Leave NULL for current working directory. Ignored if write = FALSE.                                                                                                                                                                                                                                                                                                                                                                  |
| **write**                   | Logical. If TRUE, write to disk the Reference Matrix file.                                                                                                                                                                                                                                                                                                                                                                                                                                   |
| **verbose**                 | Logical. If FALSE, does not print progress messages and output, but warnings and errors will still be printed.                                                                                                                                                                                                                                                                                                                                                                               |

### Output

A data.table containing the Seurat object’s RNA counts or any other
specified assay layer, with cell identities or barcodes as column names
and feature names as first column. If write = TRUE, the data.table is
also written to disk. If check.size = TRUE, will instead return the
Seurat object.

### Example

    Reference_Matrix_Builder(
        seurat_object = pbmc1k,
        ident.1 = "seurat_clusters",
        ident.2 = "orig.ident",
        double.ident = FALSE,
        clusters.1 = c("Cluster.7","Cluster.14"),
        clusters.1.invert = TRUE,
        clusters.2 = "IMMUNE_CTRL",
        clusters.2.invert = TRUE,
        file.name = "STIM_Reference_Matrix",
        file.format = "tsv",
    )

    Reference_Matrix_Builder(
        seurat_object = pbmc1k,
        check.size = TRUE,
        max.matrix.size = 20,
    )

## Mixture_File_Builder

### Description

This function builds a Mixture File for CIBERSORTx from one or several
.txt or .csv files containing CPM of bulk RNA-seq or microarray
experiments.

TBD

## About the Author

Alexis Varin, PhD

Inserm Research Engineer in Bioinformatics Analyses

Inserm UMR Right

Bureau 228, 2ème étage

Bâtiment B3 Médecine

15 Bd Maréchal de Lattre de Tassigny 21000 Dijon, France

<alexis.varin@inserm.fr>

<https://linkedin.com/in/alexisvarin/>

<https://umr-right.com/>

## Citation

If you find this package useful and would like to cite it, please find
the citation information below :

TBD
