data(list = "gtex_gokegg", package = "NetActivityData", envir = environment())
model <- get("gtex_gokegg")

se <- DESeq2::makeExampleDESeqDataSet(n = 50, m = 5)
rownames(se) <- colnames(model)[1:50]
se$sex <- factor(c("M", "M", "F", "M", "F"))
ddsSE <- DESeq2::DESeqDataSet(se, design = ~ sex)
vst <- DESeq2::varianceStabilizingTransformation(ddsSE)

test_that("Function standardizes and adds missing genes", {

  out <- prepareSummarizedExperiment(vst, "gtex_gokegg")
  expect_s4_class(out, "SummarizedExperiment")
  expect_equal(SummarizedExperiment::colData(out), SummarizedExperiment::colData(vst))

  ## Check standardize
  val <- SummarizedExperiment::assay(vst["ENSG00000000003", ])
  expect_equal( SummarizedExperiment::assay(out["ENSG00000000003", ]), (val - mean(val))/sd(val),
               ignore_attr = TRUE
  )
  ## Check adding genes
  expect_true("ENSG00000273540" %in% rownames(out))
  expect_equal( SummarizedExperiment::assay(out["ENSG00000273540", ]), rep(0, ncol(out)), ignore_attr = TRUE)

})

test_that("Custom matrix", {
  mat <- matrix(1:6, nrow = 2, dimnames = list(c("GS1", "GS2"), c("ENSG00000000003", "ENSG00000000419", "COT")))

  expect_warning(out2 <- prepareSummarizedExperiment(vst, mat), "1 genes present in the model not found in input data. The expression of all samples will be set to 0.")
  expect_s4_class(out2, "SummarizedExperiment")
  expect_equal(SummarizedExperiment::colData(out2), SummarizedExperiment::colData(vst))

  ## Check standardize
  val <- SummarizedExperiment::assay(vst["ENSG00000000003", ])
  expect_equal( SummarizedExperiment::assay(out2["ENSG00000000003", ]), (val - mean(val))/sd(val),
               ignore_attr = TRUE
  )

  ## Check adding genes
  expect_true("COT" %in% rownames(out2))
  expect_equal(SummarizedExperiment::assay(out2["COT", ]), rep(0, ncol(vst)), ignore_attr = TRUE)

})


test_that("Check inputs", {

  expect_error(prepareSummarizedExperiment(SummarizedExperiment::assay(vst), "gtex_gokegg"), "SE should be a SummarizedExperiment")
  expect_error(prepareSummarizedExperiment(vst, "pac"), "model name not present in NetActivityData.")
  expect_error(prepareSummarizedExperiment(vst, c(1, 4, 5)), "model should be a string with a model name present in NetActivityData or a matrix with the gene set weights.")
  expect_error(prepareSummarizedExperiment(se, "gtex_gokegg"), "Count data is not accepted by the model. Data should be normalized before passing it to the model.")

})



test_that("Custom matrix without genes", {
  mat2 <- matrix(1:6, nrow = 2, dimnames = list(c("GS1", "GS2"), c("PAC", "MAH", "COT")))

  expect_error(prepareSummarizedExperiment(vst, mat2), "None of the genes of the model were present in the input data. Check the genes annotation.")
})
