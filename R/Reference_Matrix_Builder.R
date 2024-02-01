#' @title Reference_Matrix_Builder
#'
#' @description This function builds a Reference Matrix for CIBERSORTx from a Seurat object.
#'
#' @param seurat_object A Seurat object.
#' @param ident.1 Character. The first identity to use (for example "seurat_clusters" or "orig.ident"). Leave NULL if you want to use current identity (active.ident).
#' @param ident.2 Character. The second identity to use. This will subset each ident.1 cluster with ident.2 identities (for example if ident.1 is "seurat_clusters" and ident.2 "orig.ident" cluster 0 will be divided into several clusters 0 each corresponding to an orig.ident). Leave NULL if you don't want to use a second identity.
#' @param double.ident Logical. If TRUE, each cell barcode will be renamed with a double identity (ident.1_ident.2). If FALSE, each cell barcode will be renamed with ident.1. Ignored if ident.2 is NULL.
#' @param sep.double.ident Character. The separator to use between ident.1 and ident.2 if double.ident = TRUE. Ignored if ident.2 is NULL.
#' @param reverse.double.ident Logical. If TRUE, the order of ident.1 and ident.2 will be reversed when creating the double identity (ident.2_ident.1). Ignored if ident.2 is NULL.
#' @param clusters.1 Character. The clusters to subset from ident.1. Leave NULL if you want all clusters.
#' @param clusters.1.invert Logical. If TRUE, invert the selection from clusters.1 (if clusters.1 = c("0","1","2") and clusters.1.invert = TRUE, all clusters except 0, 1 and 2 will be selected).
#' @param clusters.2 Character. The clusters to subset from ident.2. Leave NULL if you want all clusters.
#' @param clusters.2.invert Logical. If TRUE, invert the selection from clusters.2.
#' @param downsample.object.first Numeric. The number of cells to downsample from the Seurat object before subsetting clusters.1 and clusters.2. Leave NULL if you don't want to downsample.
#' @param downsample.object.last Numeric. The number of cells to downsample from the Seurat object after subsetting clusters.1 and clusters.2. Leave NULL if you don't want to downsample.
#' @param downsample.cluster Numeric. The number of cells to downsample from each cluster of ident.1. Will be performed after downsample.object.last. Leave NULL if you don't want to downsample.
#' @param automatic.downsample Logical. If TRUE, automatically downsample the Seurat object so that the Reference Matrix file written to disk would be just under the max.matrix.size limit (empirical). Performed after all other subset and downsample parameters. Ignored if write.table = FALSE. Please report an issue if you see a significant difference between the .txt or .tsv file written to disk and max.matrix.size (for example max.matrix.size is set to 200 MB but the .txt file is 400 MB or 100 MB).
#' @param max.matrix.size Numeric. The maximum size of the Reference Matrix file written to disk in MB. Will stop the function if the Reference Matrix file written to disk is estimated to be over this limit, or if automatic.downsample = TRUE will downsample the Seurat object instead so that the Reference Matrix output file is under the size limit. Ignored if write.table = FALSE.
#' @param file.name Character. The name of the Reference Matrix file written to disk. Must not contain any space. You can specify .txt or .tsv for the output file.
#' @param path Character. The path to write the Reference Matrix into. Leave NULL for current working directory.
#' @param write.table Logical. If TRUE, write to disk the Reference Matrix as a .txt or .tsv file.
#'
#' @return A data.frame containing the Seurat object's RNA counts, with cell identities as column names and feature names as row names. If write.table = TRUE, a .txt or .tsv file of the data.frame is also written to disk.
#'
#' @examples
#' Reference_Matrix_Builder(
#' seurat_object = pbmc1k,
#' ident.1 = "seurat_clusters",
#' clusters.1 = c("Cluster.0","Cluster.4","Cluster.5"),
#' file.name = "Reference_Matrix_1.txt",
#' write.table = TRUE
#' )
#'
#' Reference_Matrix_Builder(
#' seurat_object = pbmc1k,
#' ident.2 = "seurat_clusters",
#' downsample.object.first = 300,
#' write.table = TRUE
#' )
#'
#' Reference_Matrix_Builder(
#' seurat_object = pbmc1k,
#' ident.1 = "seurat_clusters",
#' ident.2 = "orig.ident",
#' double.ident = FALSE,
#' clusters.1 = c("Cluster.7","Cluster.14"),
#' clusters.1.invert = TRUE,
#' clusters.2 = "IMMUNE_CTRL",
#' clusters.2.invert = TRUE,
#' file.name = "STIM_Reference_Matrix.tsv",
#' write.table = TRUE
#' )
#'
#' @import Seurat
#' @import SeuratObject
#' @export

Reference_Matrix_Builder = function(
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
    automatic.downsample = FALSE,
    max.matrix.size = 900,
    file.name = "Reference_Matrix.txt",
    path = NULL,
    write.table = FALSE) {

  if (is.null(ident.2) & is.character(clusters.2))
    stop("You must provide an ident.2 to subset clusters.2 from")

  if (is.character(ident.1))
    Idents(seurat_object) = ident.1

  seurat_object@meta.data$ident.1 = Idents(seurat_object)

  if (is.character(ident.2)) {
    checkident.1 = Idents(seurat_object)
    Idents(seurat_object) = ident.2
    seurat_object@meta.data$ident.2 = Idents(seurat_object)
    checkident.2 = Idents(seurat_object)
    if (checkident.1[1] == checkident.2[1])
      stop("ident.1 and ident.2 or seurat_object@active.ident and ident.2 are the same, please choose different identities")
    if (isTRUE(reverse.double.ident))
      seurat_object@meta.data$double.ident = paste(checkident.2,checkident.1,sep = sep.double.ident)
    else seurat_object@meta.data$double.ident = paste(checkident.1,checkident.2,sep = sep.double.ident)
    Idents(seurat_object) = "double.ident"
  }

  if(is.numeric(downsample.object.first)) {
    if (downsample.object.first > length(colnames(seurat_object)))
      warning("downsample.object.first is greater than the number of cells in seurat_object, no downsampling was done")
    current.ident = Idents(seurat_object)
    seurat_object@meta.data$seurat_object = "seurat_object"
    Idents(seurat_object) = "seurat_object"
    cell.list = WhichCells(seurat_object, downsample = downsample.object.first)
    Idents(seurat_object) = current.ident
    seurat_object = seurat_object[,cell.list]
  }

  if(is.character(clusters.1)) {
    Idents(seurat_object) = "ident.1"
    seurat_object = subset(seurat_object, idents = clusters.1, invert = clusters.1.invert)
    if (is.character(ident.2))
      Idents(seurat_object) = "double.ident"
  }

  if (is.character(clusters.2)) {
    Idents(seurat_object) = "ident.2"
    seurat_object = subset(seurat_object, idents = clusters.2, invert = clusters.2.invert)
    Idents(seurat_object) = "double.ident"
  }

  if(is.numeric(downsample.object.last)) {
    if (downsample.object.last > length(colnames(seurat_object)))
      warning("downsample.object.last is greater than the number of cells in seurat_object, no downsampling was done")
    current.ident = Idents(seurat_object)
    seurat_object@meta.data$seurat_object2 = "seurat_object2"
    Idents(seurat_object) = "seurat_object2"
    cell.list = WhichCells(seurat_object, downsample = downsample.object.last)
    Idents(seurat_object) = current.ident
    seurat_object = seurat_object[,cell.list]
  }

  if(is.numeric(downsample.cluster)) {
    Idents(seurat_object) = "ident.1"
    cell.list = WhichCells(seurat_object, downsample = downsample.cluster)
    seurat_object2 = seurat_object[,cell.list]
    if (length(colnames(seurat_object2)) == length(colnames(seurat_object)))
      warning("downsample.cluster is greater than the number of cells in each cluster, no downsampling was done")
    seurat_object = seurat_object2
    if (is.character(ident.2))
      Idents(seurat_object) = "double.ident"
  }

  if (isFALSE(double.ident))
    Idents(seurat_object) = "ident.1"

  if (length(seurat_object[["RNA"]]$counts) > 490000*max.matrix.size/1.024 & isTRUE(write.table))
    if (isFALSE(automatic.downsample))
      stop(sprintf("The Reference Matrix output file is projected to be over your size limit of %s MB on CIBERSORTx web portal (%s cells by %s features), please subset clusters, downsample the number of cells or set automatic.downsample = TRUE",max.matrix.size,length(colnames(seurat_object[["RNA"]]$counts)),length(rownames(seurat_object[["RNA"]]$counts))))
    else {
      projected.cell.number = trunc(490000/1.024*max.matrix.size/length(rownames(seurat_object[["RNA"]]$counts)))
      warning(sprintf("The Reference Matrix output file is projected to be over your size limit of %s MB on CIBERSORTx web portal (%s cells by %s features), automatically downsampling to %s cells to be under the size limit",max.matrix.size,length(colnames(seurat_object[["RNA"]]$counts)),length(rownames(seurat_object[["RNA"]]$counts)),projected.cell.number))
      current.ident = Idents(seurat_object)
      seurat_object@meta.data$seurat_object3 = "seurat_object3"
      Idents(seurat_object) = "seurat_object3"
      cell.list = WhichCells(seurat_object, downsample = projected.cell.number)
      Idents(seurat_object) = current.ident
      seurat_object = seurat_object[,cell.list]
    }

  refmat = as.data.frame(seurat_object[["RNA"]]$counts)
  colnames(refmat) = seurat_object@active.ident

  if (isTRUE(write.table)) {
    if(is.null(path))
      path = getwd()
    if (isTRUE(grepl(" ",file.name))) {
      warning("Your file name contains one or several spaces, renaming with underscores as CIBERSORTx will otherwise report an error with your Reference Matrix")
      file.name = gsub(" ","_",file.name)
    }
    refmat = cbind(Gene = rownames(refmat),refmat)
    write.table(refmat, file = paste0(path,"/",file.name), sep = "\t", quote = F, row.names = F)
    refmat[[1]] = NULL
  }

  return(refmat)
}
