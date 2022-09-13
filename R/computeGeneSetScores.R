#' Compute gene set scores
#'
#' This function will compute the gene set scores using gene weights previously
#' computed. The package `NetActivityData` contains different pre-trained models
#' that can be used to compute the gene set scores. Models included in
#' `NetActivityData` also includes gene set annotation.
#'
#' This function can also compute the gene set scores for a model not present in
#' `NetActivityData`. In this case, `model` should be a matrix where the columns
#' are the genes and the rows the gene sets. When using a custom model, we can
#' add the gene set annotation using the `annot` paramenter. `annot` parameter
#' should contain a column named `GeneSet` matching the gene set ids from the
#' weights matrix (rownames of weights matrix).
#'
#' Notice that the function will not accept raw count data. We recommend to
#' convert count data to continuous values using the Variant Stabilization
#' Transformation from [DESeq2::varianceStabilizingTransformation].
#'

#' @param SE A `SummarizedExperiment`
#' @param model A string matching a model in `NetActivityData` or a custom
#' matrix.
#' @param annot A `data.frame` with the gene set annotation, only when using a
#' custom model.
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
computeGeneSetScores <- function(SE, model, annot = NULL){

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
        annot <- paste0(model, "_annot")

        av.mods <- data(package = "NetActivityData")
        av.mods <- av.mods$results[, 3]
        if (!model %in% av.mods){
            stop("model name not present in NetActivityData.")
        } else {
            data(list = model, package = "NetActivityData", envir = environment())
            model <- get(model)

            data(list = annot, package = "NetActivityData", envir = environment())
            annot <- get(annot)
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

    if (!is.null(annot)){
        if (!is(annot, "data.frame")){
            warning("annot should be a data.frame. Output will not contain gene set annotation.")
        }
        else if (!"GeneSet" %in% colnames(annot)){
            warning("annot should contain the column GeneSet with the gene set ids. Output will not contain gene set annotation.")
        }
        else if (all(!rownames(res) %in% annot$GeneSet)){
            warning("None of the gene sets present in annot GeneSet column are present in the model. Output will not contain gene set annotation.")
        }
        else {
            rownames(annot) <- annot$GeneSet
            SummarizedExperiment::rowData(res) <- annot[rownames(res), ]
            if (!all(rownames(res) %in% annot$GeneSet)){
                warning("Some gene sets might not contain annotation.")
            }
        }
    }
    return(res)
}
