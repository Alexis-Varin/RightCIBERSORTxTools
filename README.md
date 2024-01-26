
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

## CIBERSORTx_Reference_Matrix_Builder

### Description

This function builds a Reference Matrix for CIBERSORTx from a Seurat
object.

### Dependencies

- Seurat
- SeuratObject

### Usage

    CIBERSORTx_Reference_Matrix_Builder(
      seurat_object,
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
      file.name = "Reference_Matrix.txt",
      path = NULL,
      write.table = FALSE
    )

### Arguments

| Argument                    | Definition                                                                                                                                                                                                                                                                                                                |
|-----------------------------|:--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **seurat_object**           | A Seurat object                                                                                                                                                                                                                                                                                                           |
| **ident.1**                 | Character. The first identity to use (for example “seurat_clusters” or “orig.ident”). Leave NULL if you want to use current identity (active.ident).                                                                                                                                                                      |
| **ident.2**                 | Character. The second identity to use. This will subset each ident.1 cluster with ident.2 identities (for example if ident.1 is “seurat_clusters” and ident.2 “orig.ident” cluster 0 will be divided into several clusters 0 each corresponding to an orig.ident). Leave NULL if you don’t want to use a second identity. |
| **double.ident**            | Logical. If TRUE, each cell barcode will be renamed with a double identity (ident.1_ident.2). If FALSE, each cell barcode will be renamed with ident.1. Ignored if ident.2 is NULL.                                                                                                                                       |
| **sep.double.ident**        | Character. The separator to use between ident.1 and ident.2 if double.ident = TRUE. Ignored if ident.2 is NULL.                                                                                                                                                                                                           |
| **reverse.double.ident**    | Logical. If TRUE, the order of ident.1 and ident.2 will be reversed when creating the double identity (ident.2_ident.1). Ignored if ident.2 is NULL.                                                                                                                                                                      |
| **clusters.1**              | Character. The clusters to subset from ident.1. Leave NULL if you want all clusters.                                                                                                                                                                                                                                      |
| **clusters.1.invert**       | Logical. If TRUE, invert the selection from clusters.1 (if clusters.1 = c(“0”,“1”,“2”) and clusters.1.invert = TRUE, all clusters except 0, 1 and 2 will be selected).                                                                                                                                                    |
| **clusters.2**              | Character. The clusters to subset from ident.2. Leave NULL if you want all clusters.                                                                                                                                                                                                                                      |
| **clusters.2.invert**       | Logical. If TRUE, invert the selection from clusters.2.                                                                                                                                                                                                                                                                   |
| **downsample.object.first** | Numeric. The number of cells to downsample from the Seurat object before subsetting clusters.1 and clusters.2. Leave NULL if you don’t want to downsample.                                                                                                                                                                |
| **downsample.object.last**  | Numeric. The number of cells to downsample from the Seurat object after subsetting clusters.1 and clusters.2. Leave NULL if you don’t want to downsample.                                                                                                                                                                 |
| **downsample.cluster**      | Numeric. The number of cells to downsample from each cluster. Will be performed after downsample.object.last. Leave NULL if you don’t want to downsample.                                                                                                                                                                 |
| **file.name**               | Character. The name of the Reference Matrix file to write. Must not contain any space. You can specify .txt or .tsv for the output file.                                                                                                                                                                                  |
| **path**                    | Character. The path to write the Reference Matrix into. Leave NULL for current working directory.                                                                                                                                                                                                                         |
| **write.table**             | Logical. If TRUE, write to disk the Reference Matrix as a .txt or .tsv file. If FALSE, only returns the Reference Matrix as a data.frame.                                                                                                                                                                                 |

### Output

A data.frame containing the Seurat object’s RNA counts, with cell
identities as column names and gene names as row names. If write.table =
TRUE, a .txt or .tsv file of the data.frame is also written to disk.

### Example

    CIBERSORTx_Reference_Matrix_Builder(
    seurat_object = pbmc1k,
    ident.1 = "seurat_clusters",
    ident.2 = "orig.ident",
    double.ident = FALSE,
    clusters.1 = c("Cluster.7","Cluster.14"),
    clusters.1.invert = TRUE,
    clusters.2 = "IMMUNE_CTRL",
    clusters.2.invert = TRUE,
    file.name = "STIM_Reference_Matrix.tsv",
    write.table = TRUE
    )

## CIBERSORTx_Mixture_File_Builder

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
