#' @title Reference_Matrix_Builder
#'
#' @description This function builds a Reference Matrix for CIBERSORTx from a Seurat object.
#'
#' @param seurat_object A Seurat object, created with Seurat v4 or v5.
#' @param assay Character. If the Seurat object contains multiple RNA assays, you may specify which one to use (for example "RNA2" if you have created a second RNA assay you named "RNA2". See Seurat v5 vignettes for more information). You may also use another assay such as SCT if you want to extract an expression matrix of this assay for projects other than CIBERSORTx. Must be RNA assay for CIBERSORTx.
#' @param layer Character. Formerly known as slot. You can specify a layer such as "data" (normalized counts) or any other present in the Seurat object's assay if you want to extract the expression matrix of this layer for projects other than CIBERSORTx. If you have split layers the function will always join them before extracting the expression matrix. Must be counts layer for CIBERSORTx.
#' @param ident.1 Character. The first identity to use (for example "seurat_clusters" or "orig.ident"). Leave NULL if you want to use current identity (active.ident).
#' @param ident.2 Character. The second identity to use. This will subset each ident.1 cluster with ident.2 identities (for example if ident.1 is "seurat_clusters" and ident.2 "orig.ident" cluster 0 will be divided into several clusters 0 each corresponding to an orig.ident, 0_sample1, 0_sample2 etc). Leave NULL if you don't want to use a second identity.
#' @param double.ident Logical. If TRUE, each cell barcode will be renamed with a double identity (ident.1_ident.2). If FALSE, each cell barcode will be renamed with ident.1. Ignored if ident.2 is NULL.
#' @param sep.double.ident Character. The separator to use between ident.1 and ident.2 if double.ident = TRUE. Ignored if ident.2 is NULL.
#' @param reverse.double.ident Logical. If TRUE, the order of ident.1 and ident.2 will be reversed when creating the double identity (ident.2_ident.1). Ignored if ident.2 is NULL.
#' @param clusters.1 Character. The clusters to subset from ident.1. Leave NULL if you want all clusters.
#' @param clusters.1.invert Logical. If TRUE, invert the selection from clusters.1 (if clusters.1 = c("0","1","2") and clusters.1.invert = TRUE, all clusters except 0, 1 and 2 will be selected).
#' @param clusters.2 Character. The clusters to subset from ident.2. Leave NULL if you want all clusters.
#' @param clusters.2.invert Logical. If TRUE, invert the selection from clusters.2.
#' @param downsample.object.first Numeric. The number of cells to downsample from the entire Seurat object (not each cluster) before subsetting clusters.1 and clusters.2. Leave NULL if you don't want to downsample.
#' @param downsample.object.last Numeric. The number of cells to downsample from the entire Seurat object after subsetting clusters.1 and clusters.2. Leave NULL if you don't want to downsample.
#' @param downsample.cluster Numeric. The number of cells to downsample from each cluster of ident.1. Will be performed after downsample.object.last. Leave NULL if you don't want to downsample.
#' @param automatic.downsample Logical. If TRUE, automatically downsample the Seurat object so that the Reference Matrix file written to disk would be just under the max.matrix.size limit (empirical). Performed after all other subset and downsample parameters. Ignored if check.size = FALSE and write = FALSE. Please report an issue if you see a significant difference between the file size written to disk and max.matrix.size (for example max.matrix.size is set to 200 MB but the file is 400 MB or 100 MB).
#' @param check.size Logical. If TRUE, will return the Seurat object and print the estimated size of the Reference Matrix file written to disk and the number of cells to downsample. If automatic.downsample = TRUE will perform the downsample and return the Seurat object.
#' @param max.matrix.size Numeric. The maximum size of the Reference Matrix file written to disk in MB. Will stop the function if the Reference Matrix file written to disk is estimated to be over this limit, or if automatic.downsample = TRUE will downsample the Seurat object instead so that the Reference Matrix output file is under the size limit. Ignored if check.size = FALSE and write = FALSE.
#' @param cell.barcodes Logical. If TRUE, keep the cell barcodes and does not rename with cell identities, if you want to extract the expression matrix for projects other than CIBERSORTx. Must be FALSE for CIBERSORTx.
#' @param file.name Character. The name of the Reference Matrix file written to disk. Must not contain any space. Ignored if write = FALSE.
#' @param file.format Character. The format of the Reference Matrix file written to disk. Must be txt or tsv for CIBERSORTx but you can also specify csv for example if you want to extract the expression matrix for projects other than CIBERSORTx. Accept any format the data.table::fwrite function would accept. Ignored if write = FALSE.
#' @param file.sep Character. The separator to use in the Reference Matrix file written to disk. Must be Tabulation for CIBERSORTx but you can also specify a comma for example if you want to extract the expression matrix for projects other than CIBERSORTx. Accept any separator the data.table::fwrite function would accept. Ignored if write = FALSE.
#' @param path Character. The path to write the Reference Matrix into. Leave NULL for current working directory. Ignored if write = FALSE.
#' @param write Logical. If TRUE, write to disk the Reference Matrix file.
#' @param verbose Logical. If FALSE, does not print progress messages and output, but warnings and errors will still be printed.
#'
#' @return A data.table containing the Seurat object's RNA counts or any other specified assay layer, with cell identities or barcodes as column names and feature names as first column. If write = TRUE, the data.table is also written to disk. If check.size = TRUE, will instead return the Seurat object.
#'
#' @examples
#' Reference_Matrix_Builder(
#' seurat_object = pbmc1k,
#' ident.1 = "seurat_clusters",
#' clusters.1 = c("Cluster.0","Cluster.4","Cluster.5"),
#' file.name = "Reference_Matrix_1",
#' )
#'
#' Reference_Matrix_Builder(
#' seurat_object = pbmc1k,
#' ident.2 = "seurat_clusters",
#' downsample.object.first = 300,
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
#' file.name = "STIM_Reference_Matrix",
#' file.format = "tsv",
#' )
#'
#' Reference_Matrix_Builder(
#' seurat_object = pbmc1k,
#' check.size = TRUE,
#' max.matrix.size = 20,
#' )
#'
#' @import Seurat
#' @import SeuratObject
#' @import data.table
#' @export

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
    verbose = TRUE) {

  if (isTRUE(verbose))
    cat("Starting...","\n")

  if (is.null(ident.2) & is.character(clusters.2))
    stop("You must provide an ident.2 to subset clusters.2 from")

  to.put.back = Idents(seurat_object)

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

  gc(verbose = F)

  if(is.numeric(downsample.object.first)) {
    if (isTRUE(verbose))
      cat("Downsampling...","\n")
    if (downsample.object.first > length(colnames(seurat_object)))
      warning("downsample.object.first is greater than the number of cells in seurat_object, no downsampling was done",immediate. = T)
    current.ident = Idents(seurat_object)
    seurat_object@meta.data$seurat_object = "seurat_object"
    Idents(seurat_object) = "seurat_object"
    cell.list = WhichCells(seurat_object, downsample = downsample.object.first)
    Idents(seurat_object) = current.ident
    seurat_object = seurat_object[,cell.list]
  }

  if(is.character(clusters.1)) {
    if (isTRUE(verbose))
      cat("Subsetting...","\n")
    Idents(seurat_object) = "ident.1"
    seurat_object = subset(seurat_object, idents = clusters.1, invert = clusters.1.invert)
    if (is.character(ident.2))
      Idents(seurat_object) = "double.ident"
  }

  if (is.character(clusters.2)) {
    if (isTRUE(verbose))
      cat("Subsetting...","\n")
    Idents(seurat_object) = "ident.2"
    seurat_object = subset(seurat_object, idents = clusters.2, invert = clusters.2.invert)
    Idents(seurat_object) = "double.ident"
  }

  if(is.numeric(downsample.object.last)) {
    if (isTRUE(verbose))
      cat("Downsampling...","\n")
    if (downsample.object.last > length(colnames(seurat_object)))
      warning("downsample.object.last is greater than the number of cells in seurat_object, no downsampling was done",immediate. = T)
    current.ident = Idents(seurat_object)
    seurat_object@meta.data$seurat_object2 = "seurat_object2"
    Idents(seurat_object) = "seurat_object2"
    cell.list = WhichCells(seurat_object, downsample = downsample.object.last)
    Idents(seurat_object) = current.ident
    seurat_object = seurat_object[,cell.list]
  }

  if(is.numeric(downsample.cluster)) {
    if (isTRUE(verbose))
      cat("Downsampling...","\n")
    Idents(seurat_object) = "ident.1"
    cell.list = WhichCells(seurat_object, downsample = downsample.cluster)
    seurat_object2 = seurat_object[,cell.list]
    if (length(colnames(seurat_object2)) == length(colnames(seurat_object)))
      warning("downsample.cluster is greater than the number of cells in each cluster, no downsampling was done",immediate. = T)
    seurat_object = seurat_object2
    if (is.character(ident.2))
      Idents(seurat_object) = "double.ident"
  }

  if (isFALSE(double.ident))
    Idents(seurat_object) = "ident.1"

  if (inherits(seurat_object[[assay]],"Assay5")) {
    if (isTRUE(verbose))
      cat("Joining layers...","\n")
    seurat_object[[assay]] = JoinLayers(seurat_object[[assay]])
    if (isTRUE(verbose))
      cat("Extracting the expression matrix...","\n")
    refmat = LayerData(seurat_object, assay = assay, layer = layer)
  }
  else {
    if (isTRUE(verbose))
      cat("Extracting the expression matrix...","\n")
    refmat = GetAssayData(seurat_object, assay = assay, slot = layer)
  }

  projected.cell.number = trunc(500000/1.02*max.matrix.size/length(rownames(refmat)))

  if (isTRUE(check.size)) {
    refmat.size = length(refmat)/500000/1.02

    if (projected.cell.number >= length(colnames(refmat)))
      cat("Current estimated Reference Matrix size on CIBERSORTx web portal is between ",
          trunc(refmat.size/1.01),
          " and ",
          trunc(refmat.size*1.02),
          " MB :",
          "\n",
          "Matrix of ",
          length(colnames(refmat)),
          " cells by ",
          length(rownames(refmat)),
          " features",
          "\n",
          "You do not need to downsample your Seurat object","\n",sep="")

    else {
      if (isTRUE(automatic.downsample)) {
        cat("Current estimated Reference Matrix size on CIBERSORTx web portal is between ",
            trunc(refmat.size/1.01),
            " and ",
            trunc(refmat.size*1.02),
            " MB :",
            "\n",
            "Matrix of ",
            length(colnames(refmat)),
            " cells by ",
            length(rownames(refmat)),
            " features",
            "\n",
            "Downsampling the Seurat object to ",
            projected.cell.number,
            " cells for a Reference Matrix under ",
            max.matrix.size,
            " MB...","\n",sep="")
        current.ident = Idents(seurat_object)
        seurat_object@meta.data$seurat_object3 = "seurat_object3"
        Idents(seurat_object) = "seurat_object3"
        cell.list = WhichCells(seurat_object, downsample = projected.cell.number)
        Idents(seurat_object) = current.ident
        seurat_object = seurat_object[,cell.list]
      }
      else
        cat("Current estimated Reference Matrix size on CIBERSORTx web portal is between ",
            trunc(refmat.size/1.01),
            " and ",
            trunc(refmat.size*1.02),
            " MB :",
            "\n",
            "Matrix of ",
            length(colnames(refmat)),
            " cells by ",
            length(rownames(refmat)),
            " features",
            "\n",
            "You may want to downsample your Seurat object to ",
            projected.cell.number,
            " cells for a Reference Matrix under ",
            max.matrix.size,
            " MB","\n",sep="")
    }

    if (isTRUE(verbose))
      cat("Cleaning...","\n")
    seurat_object$ident.1 = NULL
    if (is.character(ident.2)) {
      seurat_object$ident.2 = NULL
      seurat_object$double.ident = NULL
    }
    if(is.numeric(downsample.object.first))
      seurat_object$seurat_object = NULL
    if(is.numeric(downsample.object.last))
      seurat_object$seurat_object2 = NULL
    if(isTRUE(automatic.downsample) & projected.cell.number < length(colnames(refmat)))
      seurat_object$seurat_object3 = NULL
    Idents(seurat_object) = to.put.back
    gc(verbose = F)
    if (isTRUE(verbose))
      cat("Done.","\n")
    return(seurat_object)
  }

  if (length(refmat) > 500000*max.matrix.size/1.02 & isTRUE(write))
    if (isFALSE(automatic.downsample))
      stop(paste0("The Reference Matrix file is projected to be over the size limit of ",
               max.matrix.size,
               " MB on CIBERSORTx web portal :",
               "\n",
               " Matrix of ",
               length(colnames(refmat)),
               " cells by ",
               length(rownames(refmat)),
               " features",
               "\n",
               " Please subset clusters, downsample the number of cells or set automatic.downsample = TRUE"))
    else {
      warning(paste0("The Reference Matrix file is projected to be over the size limit of ",
      max.matrix.size,
      " MB on CIBERSORTx web portal :",
      "\n",
      " Matrix of ",
      length(colnames(refmat)),
      " cells by ",
      length(rownames(refmat)),
      " features",
      "\n",
      " Automatically downsampling to ",
      projected.cell.number,
      " cells to be under the size limit..."),immediate. = T)
      current.ident = Idents(seurat_object)
      seurat_object@meta.data$seurat_object3 = "seurat_object3"
      Idents(seurat_object) = "seurat_object3"
      cell.list = WhichCells(seurat_object, downsample = projected.cell.number)
      Idents(seurat_object) = current.ident
      seurat_object = seurat_object[,cell.list]
      if (inherits(seurat_object[[assay]],"Assay5")) {
        seurat_object[[assay]] = JoinLayers(seurat_object[[assay]])
        refmat = LayerData(seurat_object, assay = assay, layer = layer)
      }
      else
        refmat = GetAssayData(seurat_object, assay = assay, slot = layer)
    }

  if (isTRUE(verbose))
    cat("Building the data.table...","\n")
  refmat = as.data.frame(refmat)
  refmat = cbind(Gene = rownames(refmat),refmat)
  refmat = setDT(refmat)
  if (isFALSE(cell.barcodes))
    colnames(refmat) = c("Gene", as.character(seurat_object@active.ident))

  if (isTRUE(write)) {
    if(is.null(path))
      path = getwd()
    if (isTRUE(grepl(" ",file.name))) {
      warning("The file name contains one or several spaces, renaming with underscores as CIBERSORTx will otherwise report an error with the Reference Matrix...",immediate. = T)
      file.name = gsub(" ","_",file.name)
    }
    if (isTRUE(verbose))
      cat("Writing to disk...","\n")
    fwrite(refmat, file = paste0(path,"/",file.name,".",file.format), sep = file.sep, quote = F, row.names = F)
  }

  if (isTRUE(verbose))
    cat("Cleaning...","\n")
  gc(verbose = F)
  if (isTRUE(verbose))
    cat("Done.","\n")
  return(refmat)
}
