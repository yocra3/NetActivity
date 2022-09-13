#' Compute gene set scores
#'
#' This function will compute the gene set scores using gene weights previously
#' computed. The package `NetActivityData` contains different pre-trained models
#' that can be used to compute the gene set scores.
#'
#' This function can also compute the gene set scores for a model not present in
#' `NetActivityData`. In this case, `param` should be a matrix where the columns
#' are the genes and the rows the gene sets.
#'
#'
#' Notice that the function will not accept raw count data. We recommend to
#' convert count data to continuous values using the Variant Stabilization
#' Transformation from [DESeq2::varianceStabilizingTransformation].
#'

#' @param SE A `SummarizedExperiment`
#' @param model A string matching a model in `NetActivityData` or a custom
#' matrix.
#' @return A `SummarizedExperiment` with the gene set scores.
#'
#' @export
#' @examples
#' library(airway)
#' data(airway)
#' ddsSE <- DESeq2::DESeqDataSet(airway, design = ~ cell + dex)
#' vst <- DESeq2::varianceStabilizingTransformation(ddsSE)
#' out <- prepareSummarizedExperiment(vst, "gtex_gokegg")
#' scores <- computeGeneSetScores(out, "gtex_gokegg")
computeGeneSetScores <- function(SE, model){

    ## Check SE
    if (!is(SE, "SummarizedExperiment")){
        stop("SE should be a SummarizedExperiment")
    }

    if (DelayedArray::type(SummarizedExperiment::assay(SE)) == "integer"){
        stop("Count data is not accepted by the model. Data should be normalized before passing it to the model.")
    }

    if (sum(duplicated(rownames(SE))) > 0){
        stop("SE has features with the same gene name. Remove duplicated genes before proceeding.")
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

    genes <- colnames(model)
    n.com.genes <- sum(!genes %in% rownames(SE))

    ## Check genes
    if (n.com.genes > 0){
        stop("All genes present in the model should be present in the SE. The function prepareSummarizedExperiment can solve this issue.")
    }

    mat <- SummarizedExperiment::assay(SE[genes, ])

    scores <- t(mat) %*% t(model)

    res <- SummarizedExperiment::SummarizedExperiment(t(scores), colData =  SummarizedExperiment::colData(SE))

}
