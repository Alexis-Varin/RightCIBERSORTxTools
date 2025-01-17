% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Mixture_File_Builder.R
\name{Mixture_File_Builder}
\alias{Mixture_File_Builder}
\title{Mixture_File_Builder}
\usage{
Mixture_File_Builder(
  objects = NULL,
  files.path = NULL,
  file.name = "Mixture_File",
  file.format = "txt",
  file.sep = "\\t",
  write.path = NULL,
  write = TRUE,
  verbose = TRUE
)
}
\arguments{
\item{objects}{A data.frame, data.table or matrix object or a mixed list of data.frame, data.table and/or matrix objects corresponding to bulk RNA-seq and/or microarray experiments, with genes as rows and gene names as row names or first column. A single object may contain multiple samples and/or experiments as columns as long as it contains a single gene names column as first column or gene names as row names (e.g. | gene | sample1 | sample2 | sample3 | etc). You may mix objects with different gene order and number (e.g. an object containing 13567 genes and another with 18458 genes).}

\item{files.path}{Character. The path to read txt, csv and/or tsv files corresponding to bulk RNA-seq or microarray experiments, with genes as rows and gene names as first column. A single file may contain multiple samples and/or experiments as columns as long as it contains a single gene names column as first column. You may mix files with different gene order and number.}

\item{file.name}{Character. The name of the Mixture File written to disk. Must not contain any space. Ignored if write = FALSE.}

\item{file.format}{Character. The format of the Mixture File written to disk. Must be txt or tsv for CIBERSORTx but you can also specify csv for example for projects other than CIBERSORTx. Accept any format the data.table::fwrite function would accept. Ignored if write = FALSE.}

\item{file.sep}{Character. The separator to use in the Mixture File written to disk. Must be Tabulation for CIBERSORTx but you can also specify a comma for example for projects other than CIBERSORTx. Accept any separator the data.table::fwrite function would accept. Ignored if write = FALSE.}

\item{write.path}{Character. The path to write the Mixture File into. Leave NULL for current working directory. Ignored if write = FALSE.}

\item{write}{Logical. If TRUE, write to disk the Mixture File.}

\item{verbose}{Logical. If FALSE, does not print progress messages and output, but warnings and errors will still be printed.}
}
\value{
A data.table containing the counts of the objects and/or files provided, with each column being a sample. If write = TRUE, the data.table is also written to disk.
}
\description{
This function builds a Mixture File for CIBERSORTx from bulk RNA-seq and/or microarray files and/or R objects.
}
