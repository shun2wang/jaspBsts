library(jaspTools)
library(testthat)

jaspTools::runTestsTravis(module = getwd())
testthat::snapshot_review()
