# code to generate tests

#context("bayesianStateSpace")

options <- analysisOptions("bayesianStateSpace")
options$dependent <- "contNormal"
options$checkboxLocalLinearTrend <- TRUE
options$checkboxPlotAggregatedStates <- TRUE
options$checkBoxResidual <- TRUE
set.seed(1)
results <- runAnalysis("bayesianStateSpace", "bstsTest.csv", options)


test_that("Residual Plot matches", {
  plotName <- results[["results"]][["bstsMainContainer"]][["collection"]][["bstsMainContainer_bstsErrorPlots"]][["collection"]][["bstsMainContainer_bstsErrorPlots_bstsResidualPlot"]][["data"]]
  testPlot <- results[["state"]][["figures"]][[plotName]][["obj"]]
  jaspTools::expect_equal_plots(testPlot, "residual-plot")
})

test_that("Model Summary table results match", {
  table <- results[["results"]][["bstsMainContainer"]][["collection"]][["bstsMainContainer_bstsModelSummaryTable"]][["data"]]
  jaspTools::expect_equal_tables(table,
                                 list(-0.0323841668445357, 1.19950377915337, 0.424655892444242, 1.07541445856843
                                 ))
})

test_that("Aggregated State plot matches", {
  plotName <- results[["results"]][["bstsMainContainer"]][["collection"]][["bstsMainContainer_bstsStatePlots"]][["collection"]][["bstsMainContainer_bstsStatePlots_bstsAggregatedStatePlot"]][["data"]]
  testPlot <- results[["state"]][["figures"]][[plotName]][["obj"]]
  jaspTools::expect_equal_plots(testPlot, "aggregated-state")
})



