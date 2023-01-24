#' Prepare a SummarizedExperiment for computing gene set scores computation
#'
#' This function will prepare the data for the computation of gene set scores.
#' The function will perform two steps. First, the function will check whether
#' the genes present in the trained model are present in the input
#' `SummarizedExperiment`. Missing genes will be set to 0 for all samples.
#' Second, the function will standardize the gene expression values, so gene
#' values have a mean of 0 and a standard deviation of 1.
#'
#' Notice that the function will not accept raw count data. We recommend to
#' convert count data to continuous values using the Variant Stabilization
#' Transformation from [DESeq2::varianceStabilizingTransformation].
#'
#' This function can also prepare the data for a model not present in
#' `NetActivityData`. In this case, `param` should be a matrix where the columns
#' are the genes and the rows the gene sets.
#'
#' @param SE A `SummarizedExperiment`
#' @param model A string matching a model in `NetActivityData` or a custom
#' matrix.
#' @return A `SummarizedExperiment` with the data prepared for gene set score
#' computation with `computeGeneSetScores`
#'
#' @export
#' @examples
#' library(airway)
#' data(airway)
#' ddsSE <- DESeq2::DESeqDataSet(airway, design = ~ cell + dex)
#' vst <- DESeq2::varianceStabilizingTransformation(ddsSE)
#' out <- prepareSummarizedExperiment(vst, "gtex_gokegg")
prepareSummarizedExperiment <- function(SE, model){
    ## Check SE
    if (!is(SE, "SummarizedExperiment")){
        stop("SE should be a SummarizedExperiment")
    }
    if (DelayedArray::type(SummarizedExperiment::assay(SE)) == "integer"){
        stop("Count data is not accepted by the model. Data should be normalized before passing it to the model.")
    }

    ## Load model
    if (is(model, "character")){
        av.mods <- data(package = "NetActivityData")
        av.mods <- av.mods$results[, 3]
        if (!model %in% av.mods){
            stop("model name not present in NetActivityData.")
        } else {
            data(list = model, package = "NetActivityData", envir = environment())
            model <- get(model)
        }
    } else if (!is(model, "matrix")){
        stop("model should be a string with a model name present in NetActivityData or a matrix with the gene set weights.")
    }

    ## Check genes
    genes <- colnames(model)
    com.genes <- intersect(genes, rownames(SE))
    n.com.genes <- sum(!genes %in% com.genes)
    if (n.com.genes == length(genes)){
        stop("None of the genes of the model were present in the input data. Check the genes annotation.")
    }
    if (n.com.genes > 0){
        warning(sprintf("%i genes present in the model not found in input data. The expression of all samples will be set to 0.", n.com.genes))
        out_probes <- setdiff(genes, rownames(SE))
        out <- matrix(0, ncol(SE),
            nrow = length(out_probes), ncol = ncol(SE),
            dimnames = list(out_probes, colnames(SE)))
        if (is(SummarizedExperiment::assay(SE), "DelayedArray")){
            out <- DelayedArray::DelayedArray(out)
        }
        new_assay <- rbind(SummarizedExperiment::assay(SE), out)
        SE <- SummarizedExperiment::SummarizedExperiment(new_assay, colData = SummarizedExperiment::colData(SE))
    }

    ## Standardize gene expression
    vals <- SummarizedExperiment::assay(SE)
    mat <- (vals - DelayedArray::rowMeans(vals, na.rm = TRUE))/DelayedMatrixStats::rowSds(vals, na.rm = TRUE)
    mat[is.nan(mat)] <- 0
    SummarizedExperiment::assay(SE) <- mat

    SE
}
