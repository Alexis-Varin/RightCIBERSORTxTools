#' @title Mixture_File_Builder
#'
#' @description This function builds a Mixture File for CIBERSORTx from bulk RNA-seq and/or microarray files and/or R objects.
#'
#' @param objects A data.frame, data.table or matrix object or a mixed list of data.frame, data.table or matrix objects corresponding to bulk RNA-seq or microarray experiments, with genes as rows and gene names as row names or first column. A single object may contain multiple samples and/or experiments as columns as long as it contains a single gene names column as first column or gene names as row names (e.g. | gene | sample1 | sample2 | sample3 | etc). You may mix objects with different gene order and number (e.g. an object containing 13567 genes and another with 18458 genes).
#' @param files.path Character. The path to read txt, csv and/or tsv files corresponding to bulk RNA-seq or microarray experiments, with genes as rows and gene names as first column. A single file may contain multiple samples and/or experiments as columns as long as it contains a single gene names column as first column. You may mix files with different gene order and number.
#' @param file.name Character. The name of the Mixture File written to disk. Must not contain any space. Ignored if write = FALSE.
#' @param file.format Character. The format of the Mixture File written to disk. Must be txt or tsv for CIBERSORTx but you can also specify csv for example for projects other than CIBERSORTx. Accept any format the data.table::fwrite function would accept. Ignored if write = FALSE.
#' @param file.sep Character. The separator to use in the Mixture File written to disk. Must be Tabulation for CIBERSORTx but you can also specify a comma for example for projects other than CIBERSORTx. Accept any separator the data.table::fwrite function would accept. Ignored if write = FALSE.
#' @param write.path Character. The path to write the Mixture File into. Leave NULL for current working directory. Ignored if write = FALSE.
#' @param write Logical. If TRUE, write to disk the Mixture File.
#' @param verbose Logical. If FALSE, does not print progress messages and output, but warnings and errors will still be printed.
#'
#' @return A Mixture File for CIBERSORTx
#'
#' @examples
#'
#' # Example 1: Using objects
#'
#'
#' @import data.table
#' @export

Mixture_File_Builder = function(
    objects = NULL,
    files.path = NULL,
    file.name = "Mixture_File",
    file.format = "txt",
    file.sep = "\t",
    write.path = NULL,
    write = TRUE,
    verbose = TRUE) {

  if (is.null(objects) & is.null(files.path))
    stop("Please provide at least one object or a path to your files")

  if (is.character(files.path)) {
    if (isTRUE(verbose))
      cat("Reading files from ",files.path,"...","\n",sep="")
    dt.list = lapply(list.files(path=files.path, pattern="*.txt|*.csv|*.tsv", full.names = T),
                     function(x) fread(file = x, header = T))
  }

  if (!is.null(objects)) {
    if (isTRUE(verbose))
      cat("Adding objects...","\n")
    if (!is.list(objects)) objects = list(objects)
    objects = lapply(objects, function(x) {
      if (is.matrix(x)) x = as.data.frame(x)
      x = setDT(x)
      })
    if (is.character(files.path))
      dt.list = c(dt.list,objects)
    else
      dt.list = objects
  }

  if (isTRUE(verbose))
    cat("Building the data.table...","\n")
  dt.list = lapply(dt.list, function(x) {names(x)[1] = "Gene"; x})
  merged.dt = Reduce(function(x, y) merge(x, y, all = TRUE), dt.list)
  merged.dt[is.na(merged.dt)] = 0

  if (isTRUE(write)) {
    if(is.null(write.path))
      write.path = getwd()
    if (isTRUE(grepl(" ",file.name))) {
      warning("The file name contains one or several spaces, renaming with underscores as CIBERSORTx will otherwise report an error with the Mixture File...",immediate. = T)
      file.name = gsub(" ","_",file.name)
    }
    if (isTRUE(verbose))
      cat("Writing to disk...","\n")
    fwrite(x=merged.dt, file = paste0(write.path,"/",file.name,".",file.format), sep = file.sep, quote = F)
  }

  if (isTRUE(verbose))
    cat("Cleaning...","\n")
  gc(verbose = F)
  if (isTRUE(verbose))
    cat("Done.","\n")
  return(merged.dt)
}
