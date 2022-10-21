#' NetActivity: compute gene set scores from a deep learning framework
#'
#' NetActivity enables to compute gene set scores from previously trained
#' sparsely-connected autoencoders. The package contains a function to prepare
#' the data (`prepareSummarizedExperiment`) and a function to compute the gene
#' set scores (`computeGeneSetScores`). The package `NetActivityData` contains
#' different pre-trained models to be directly applied to the data.
#' Alternatively, the users might use the package to compute gene set scores
#' using custom models.
#'
#' @docType package
#' @name NetActivity
#'
#' @importFrom methods is
#' @importFrom utils data
#' @importFrom DelayedMatrixStats rowSds
#' @importFrom DelayedArray DelayedArray
#' @importFrom SummarizedExperiment assay SummarizedExperiment rowData
#' @import airway
#' @importFrom DESeq2 DESeqDataSet makeExampleDESeqDataSet varianceStabilizingTransformation
#' @import NetActivityData
NULL
