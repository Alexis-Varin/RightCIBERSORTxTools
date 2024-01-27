#' @title Reference_Matrix_Sizer
#'
#' @description This function calculates the size of a Reference Matrix which would be built from a Seurat object and suggests a number of cells to downsample from the object or performs that downsample to be under a specified size limit.
#'
#' @param seurat_object A Seurat object.
#' @param max.matrix.size Numeric. The maximum size of the Reference Matrix file written to disk in MB.
#' @param downsample Logical. If TRUE, downsample the Seurat object so that the Reference Matrix file written to disk would be just under the max.matrix.size limit (empirical). Please report an issue if you see a significant difference between the .txt or .tsv file written to disk and max.matrix.size (for example max.matrix.size is set to 200 MB but the .txt file is 400 MB or 100 MB).
#'
#' @return The Seurat object, an estimation of the Reference Matrix size in MB and a suggestion of the number of cells to downsample. If downsample = TRUE, the Seurat object will instead be downsampled to the suggested number of cells.
#'
#' @examples
#' Reference_Matrix_Sizer(
#' seurat_object = pbmc1k,
#' max.matrix.size = 20,
#' downsample = FALSE
#' )
#'
#' Reference_Matrix_Sizer(
#' seurat_object = pbmc1k,
#' max.matrix.size = 50,
#' downsample = FALSE
#' )
#'
#' @import Seurat
#' @import SeuratObject
#' @export

Reference_Matrix_Sizer = function(
    seurat_object,
    max.matrix.size = 900,
    downsample = FALSE) {

  projected.cell.number = trunc(500000/1.024*max.matrix.size/length(rownames(seurat_object$RNA@counts)))
  refmat.size = round(length(seurat_object$RNA@counts)/500000/1.024,2)

  if (projected.cell.number > length(colnames(seurat_object)))
    warning(sprintf("Current estimated Reference Matrix size on CIBERSORTx web portal is %s MB (%s cells by %s features). You do not need to downsample your Seurat object.",refmat.size, length(colnames(seurat_object)), length(rownames(seurat_object))))

  else {
    if (isTRUE(downsample)) {
      current.ident = Idents(seurat_object)
      seurat_object@meta.data$seurat_object = "seurat_object"
      Idents(seurat_object) = "seurat_object"
      cell.list = WhichCells(seurat_object, downsample = projected.cell.number)
      Idents(seurat_object) = current.ident
      seurat_object = seurat_object[,cell.list]
      warning(sprintf("Current estimated Reference Matrix size on CIBERSORTx web portal is %s MB (%s cells by %s features). The Seurat object has been downsampled to %s cells for a Reference Matrix under %s MB.",refmat.size, length(colnames(seurat_object)), length(rownames(seurat_object)), projected.cell.number, max.matrix.size))
    }
    else
      warning(sprintf("Current estimated Reference Matrix size on CIBERSORTx web portal is %s MB (%s cells by %s features). You may want to downsample your Seurat object to %s cells for a Reference Matrix under %s MB.",refmat.size, length(colnames(seurat_object)), length(rownames(seurat_object)), projected.cell.number, max.matrix.size))
  }

  return(seurat_object)
}
