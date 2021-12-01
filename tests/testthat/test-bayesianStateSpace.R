# code to generate tests

#context("bayesianStateSpace")

dateTimeVarNames <- c("dateYear", "dateWeek", "dateDay", "timeHour", "timeMinute", "timeSecond")


test_that("Model Summary table results match (AR - manual)", {
  options <- jaspTools::analysisOptions("bayesianStateSpace")
  options$dependent <- "contNormal"
  options$checkboxAr <- TRUE
  options$mcmcDraws <- 10

  set.seed(1)
  for (i in 1:6){
    options$dates <- dateTimeVarNames[i]
    results <- jaspTools::runAnalysis("bayesianStateSpace", "bstsTest.csv", options)
    table <- results[["results"]][["bstsMainContainer"]][["collection"]][["bstsMainContainer_bstsModelSummaryTable"]][["data"]]
    jaspTools::expect_equal_tables(table,
                                   list(0.411857605442053, 1.05018027614687, 0.542967391044393, 0.811701682421661
                                   ))
  }

})

test_that("Model Summary table results match (AR - automatic)", {
  options <- jaspTools::analysisOptions("bayesianStateSpace")
  options$dependent <- "contNormal"
  options$checkboxAr <- TRUE
  options$lagSelectionMethod <- "autoAR"
  options$maxNoLags <- 4
  options$mcmcDraws <- 10

  set.seed(1)
  for (i in 1:6){
    options$dates <- dateTimeVarNames[i]
    results <- jaspTools::runAnalysis("bayesianStateSpace", "bstsTest.csv", options)
    table <- results[["results"]][["bstsMainContainer"]][["collection"]][["bstsMainContainer_bstsModelSummaryTable"]][["data"]]
    jaspTools::expect_equal_tables(table,
                                   list(0.535488123460269, 1.05773803700731, 0.539404673068745, 0.72136258779103
                                   ))
  }

})


test_that("Model Summary table results match (local level)", {
  options <- jaspTools::analysisOptions("bayesianStateSpace")
  options$dependent <- "contNormal"
  options$checkboxLocalLevel <- TRUE
  options$mcmcDraws <- 10

  set.seed(1)
  for (i in 1:6){
    options$dates <- dateTimeVarNames[i]
    results <- jaspTools::runAnalysis("bayesianStateSpace", "bstsTest.csv", options)
    table <- results[["results"]][["bstsMainContainer"]][["collection"]][["bstsMainContainer_bstsModelSummaryTable"]][["data"]]
    jaspTools::expect_equal_tables(table,
                                   list(-0.0432055353493266, 1.09380783381569, 0.520147149047757, 1.08103597020621
                                   ))
  }

})


test_that("Model Summary table results match (local linear trend)", {
  options <- jaspTools::analysisOptions("bayesianStateSpace")
  options$dependent <- "contNormal"
  options$checkboxLocalLinearTrend <- TRUE
  options$mcmcDraws <- 10

  set.seed(1)
  for (i in 1:6){
    options$dates <- dateTimeVarNames[i]
    results <- jaspTools::runAnalysis("bayesianStateSpace", "bstsTest.csv", options)
    table <- results[["results"]][["bstsMainContainer"]][["collection"]][["bstsMainContainer_bstsModelSummaryTable"]][["data"]]
    jaspTools::expect_equal_tables(table,
                                   list(-0.0777667265276332, 1.20104889746718, 0.423311739669594, 1.09879731436253
                                   ))
  }

})

test_that("Model Summary table results match (semi-local linear trend)", {
  options <- jaspTools::analysisOptions("bayesianStateSpace")
  options$dependent <- "contNormal"
  options$checkboxLocalLevel <- TRUE
  options$mcmcDraws <- 10

  set.seed(1)
  for (i in 1:6){
    options$dates <- dateTimeVarNames[i]
    results <- jaspTools::runAnalysis("bayesianStateSpace", "bstsTest.csv", options)
    table <- results[["results"]][["bstsMainContainer"]][["collection"]][["bstsMainContainer_bstsModelSummaryTable"]][["data"]]
    jaspTools::expect_equal_tables(table,
                                   list(-0.0432055353493266, 1.09380783381569, 0.520147149047757, 1.08103597020621
                                   ))
  }

})

test_that("Model Summary table results match (seasonal)", {
  options <- jaspTools::analysisOptions("bayesianStateSpace")
  options$dependent <- "contNormal"
  options$checkboxLocalLevel <- TRUE
  options$seasonalities <- list(list(mu = 0,
                                     nSeasons = 2,
                                     name = "season",
                                     sample.size = 0.01,
                                     seasonDuration = 1,
                                     sigma="",
                                     sigma.guess=""
  ))
  options$mcmcDraws <- 10

  set.seed(1)
  for (i in 1:6){
    options$dates <- dateTimeVarNames[i]
    results <- jaspTools::runAnalysis("bayesianStateSpace", "bstsTest.csv", options)
    table <- results[["results"]][["bstsMainContainer"]][["collection"]][["bstsMainContainer_bstsModelSummaryTable"]][["data"]]
    jaspTools::expect_equal_tables(table,
                                   list(-0.0394752039178554, 1.12753310072425, 0.490267461631722, 1.07910143584074
                                   ))
  }

})

options <- jaspTools::analysisOptions("bayesianStateSpace")
options$dependent <- "contNormal"
options$covariates <- c("contcor1", "contcor2")
options$dates <- "dateYear"
options$postSummaryTable <- TRUE
options$checkboxLocalLevel <- TRUE
options$mcmcDraws <- 10
options$expectedModelSize <- 3
options$modelTerms <- list(list(components = "contcor1", isNuisance = FALSE), list(
  components = "contcor2", isNuisance = FALSE))
set.seed(1)


test_that("Model Summary table results match (local level - covariates)", {

  for (i in 1:6){
    options$dates <- dateTimeVarNames[i]
    results <- jaspTools::runAnalysis("bayesianStateSpace", "bstsTest.csv", options)

    table <- results[["results"]][["bstsMainContainer"]][["collection"]][["bstsMainContainer_bstsModelSummaryTable"]][["data"]]
    jaspTools::expect_equal_tables(table,
                                   list(0.0446411563379834, 1.07342874539745, 0.538811093349419, 1.03451899091434
                                   ))
  }

})

test_that("Posterior Summary of Coefficients table results match (local level - covariates)", {
  for (i in 1:6){
    options$dates <- dateTimeVarNames[i]
    results <- jaspTools::runAnalysis("bayesianStateSpace", "bstsTest.csv", options)
    table <- results[["results"]][["bstsMainContainer"]][["collection"]][["bstsMainContainer_bstsCoefficientSummaryTable"]][["data"]]
    jaspTools::expect_equal_tables(table,
                                   list("<unicode><unicode><unicode>", "contcor1", 0.0710805371191776,
                                        0.314637544172676, 1, 1, "", 0.375148071900377, "<unicode><unicode><unicode>",
                                        "contcor2", -0.287511723102937, -0.11088682615349, 1, 1, "",
                                        0.0421351582659255, 0, "(Intercept)", 0, 0, 0, 1, "", 0))

  }

})


options <- jaspTools::analysisOptions("bayesianStateSpace")
options$dependent <- "contNormal"
options$factors <- c("facGender", "facExperim")
options$dates <- "dateYear"
options$postSummaryTable <- TRUE
options$checkboxLocalLevel <- TRUE
options$mcmcDraws <- 10
options$expectedModelSize <- 3
options$modelTerms <- list(list(components = "facGender", isNuisance = FALSE), list(
  components = "facExperim", isNuisance = FALSE))
set.seed(1)


test_that("Model Summary table results match (local level - factors)", {

  for (i in 1:6){
    options$dates <- dateTimeVarNames[i]
    results <- jaspTools::runAnalysis("bayesianStateSpace", "bstsTest.csv", options)

    table <- results[["results"]][["bstsMainContainer"]][["collection"]][["bstsMainContainer_bstsModelSummaryTable"]][["data"]]
    jaspTools::expect_equal_tables(table,
                                   list(0.0579395417976514, 1.06599469439761, 0.546337972019031, 1.02729362023781
                                   ))
  }

})

test_that("Posterior Summary of Coefficients table results match (local level - factors)", {
  for (i in 1:6){
    options$dates <- dateTimeVarNames[i]
    results <- jaspTools::runAnalysis("bayesianStateSpace", "bstsTest.csv", options)
    table <- results[["results"]][["bstsMainContainer"]][["collection"]][["bstsMainContainer_bstsCoefficientSummaryTable"]][["data"]]
    jaspTools::expect_equal_tables(table,
                                   list("<unicode><unicode><unicode>", "facGenderm", 0.253231187396135,
                                        0.609442717859108, 1, 1, "", 0.638874711991056, "<unicode><unicode><unicode>",
                                        "facExperimexperimental", -0.222080709900353, 0.0399657905691665,
                                        1, 1, "", 0.271697756233028, 0, "(Intercept)", 0, 0, 0, 1, "",
                                        0))

  }

})

options$checkboxPlotAggregatedStates <- TRUE
options$checkboxPlotComponentStates <- TRUE
options$checkBoxForecastError<- TRUE
options$checkBoxResidual <- TRUE
options$checkControlChart <- TRUE
options$checkControlProbPlot <- TRUE
for (i in 1:6){
  options$dates <- dateTimeVarNames[i]

  results <- jaspTools::runAnalysis("bayesianStateSpace", "bstsTest.csv", options)

  test_that("State plots matches", {
    plotName <- results[["results"]][["bstsMainContainer"]][["collection"]][["bstsMainContainer_bstsStatePlots"]][["collection"]][["bstsMainContainer_bstsStatePlots_bstsComponentStatePlot"]][["data"]]
    testPlot <- results[["state"]][["figures"]][[plotName]][["obj"]]
    jaspTools::expect_equal_plots(testPlot, paste0("component-states-",i))

    plotName <- results[["results"]][["bstsMainContainer"]][["collection"]][["bstsMainContainer_bstsStatePlots"]][["collection"]][["bstsMainContainer_bstsStatePlots_bstsAggregatedStatePlot"]][["data"]]
    testPlot <- results[["state"]][["figures"]][[plotName]][["obj"]]
    jaspTools::expect_equal_plots(testPlot, paste0("aggregated-state-",i))
  })
  test_that("Residual Plots match", {
    plotName <- results[["results"]][["bstsMainContainer"]][["collection"]][["bstsMainContainer_bstsErrorPlots"]][["collection"]][["bstsMainContainer_bstsErrorPlots_bstsResidualPlot"]][["data"]]
    testPlot <- results[["state"]][["figures"]][[plotName]][["obj"]]
    jaspTools::expect_equal_plots(testPlot, paste0("residual-plot-",i))

    plotName <- results[["results"]][["bstsMainContainer"]][["collection"]][["bstsMainContainer_bstsErrorPlots"]][["collection"]][["bstsMainContainer_bstsErrorPlots_bstsForecastErrorPlot"]][["data"]]
    testPlot <- results[["state"]][["figures"]][[plotName]][["obj"]]
    jaspTools::expect_equal_plots(testPlot, paste0("forecast-error-plot-",i))
  })
  test_that("Control chart plots match", {
    plotName <- results[["results"]][["bstsMainContainer"]][["collection"]][["bstsMainContainer_bstsControlPlots"]][["collection"]][["bstsMainContainer_bstsControlPlots_bstsControlPlotThreshold"]][["data"]]
    testPlot <- results[["state"]][["figures"]][[plotName]][["obj"]]
    jaspTools::expect_equal_plots(testPlot, paste0("threshold-plot-",i))

    plotName <- results[["results"]][["bstsMainContainer"]][["collection"]][["bstsMainContainer_bstsControlPlots"]][["collection"]][["bstsMainContainer_bstsControlPlots_bstsControlPlotProbability"]][["data"]]
    testPlot <- results[["state"]][["figures"]][[plotName]][["obj"]]
    jaspTools::expect_equal_plots(testPlot, paste0("probability-plot-",i))
  })
}

test_that("Prediction plot matches", {
  options <- jaspTools::analysisOptions("bayesianStateSpace")
  options$dependent <- "contNormal"
  options$mcmcDraws <- 10
  options$checkboxLocalLevel <- TRUE
  options$dates <- "dateYear"
  options$predictionHorizon <- 12


  results <- jaspTools::runAnalysis("bayesianStateSpace", "bstsTest.csv", options)
  plotName <- results[["results"]][["bstsMainContainer"]][["collection"]][["bstsMainContainer_bstsPredictionPlot"]][["data"]]
  testPlot <- results[["state"]][["figures"]][[plotName]][["obj"]]
  jaspTools::expect_equal_plots(testPlot, "prediction-plot")
})






