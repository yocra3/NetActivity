data(list = "gtex_gokegg", package = "NetActivityData", envir = environment())
model <- get("gtex_gokegg")

se <- DESeq2::makeExampleDESeqDataSet(n = 50, m = 5)
rownames(se) <- colnames(model)[1:50]
se$sex <- factor(c("M", "M", "F", "M", "F"))
ddsSE <- DESeq2::DESeqDataSet(se, design = ~ sex)
vst <- DESeq2::varianceStabilizingTransformation(ddsSE)
input <- prepareSummarizedExperiment(vst, "gtex_gokegg")


test_that("Computation of gene set scores", {

  out <- computeGeneSetScores(input, "gtex_gokegg" )
  expect_s4_class(out, "SummarizedExperiment")
  expect_equal(SummarizedExperiment::colData(out), SummarizedExperiment::colData(input))

})

test_that("Custom matrix", {
  mat <- matrix(1:6, nrow = 2, dimnames = list(c("GS1", "GS2"), c("ENSG00000000003", "ENSG00000000419", "COT")))

  expect_error(computeGeneSetScores(input, mat), "All genes present in the model should be present in the SE. The function prepareSummarizedExperiment can solve this issue.")


  mat2 <- matrix(c(1, 0, 0, 1, -1, 0), nrow = 2, dimnames = list(c("GS1", "GS2"), c("ENSG00000000003", "ENSG00000000938", "ENSG00000000419")))
  out2 <- computeGeneSetScores(input, mat2)
  expect_equal(SummarizedExperiment::colData(out2), SummarizedExperiment::colData(input))

  gene_mat <- SummarizedExperiment::assay(input[colnames(mat2), ])
  gs1 <- gene_mat["ENSG00000000003", ]*mat2[1, 1] + gene_mat["ENSG00000000419", ]*mat2[1, 3]
  expect_equal(SummarizedExperiment::assay(out2["GS1", ]), gs1, ignore_attr = TRUE)
  expect_equal(SummarizedExperiment::assay(out2["GS2", ]), gene_mat["ENSG00000000938", ], ignore_attr = TRUE)


})


test_that("Check inputs", {

  expect_error(computeGeneSetScores(SummarizedExperiment::assay(input), "gtex_gokegg"), "SE should be a SummarizedExperiment")
  expect_error(computeGeneSetScores(input, "pac"), "model name not present in NetActivityData.")
  expect_error(computeGeneSetScores(input, c(1, 4, 5)), "model should be a string with a model name present in NetActivityData or a matrix with the gene set weights.")
  expect_error(computeGeneSetScores(se, "gtex_gokegg"), "Count data is not accepted by the model. Data should be normalized before passing it to the model.")

  rownames(input) <- "A"
  expect_error(computeGeneSetScores(input, "gtex_gokegg"), "SE has features with the same gene name. Remove duplicated genes before proceeding.")

})
