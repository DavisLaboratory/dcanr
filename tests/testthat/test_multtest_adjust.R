library(dcanr)
library(Biobase)
library(edgeR)
library(SummarizedExperiment)

context('Statistical tests ')

#select methods that can be run
getInfMethods <- function() {
  infmethods = c('zscore', 'diffcoex', 'dicer', 'ebcoexpress')
  if (!require('EBcoexpress')) {
    infmethods = setdiff(infmethods, 'ebcoexpress')
  }
  if (!require('COSINE')) {
    infmethods = setdiff(infmethods, 'ecf')
  }
  if (!require('GeneNet')) {
    infmethods = setdiff(infmethods, 'ggm-based')
  }

  return(infmethods)
}

getCondition <- function() {
  return(rep(1:2, each = 30))
}

getData <- function() {
  #data to test methods
  set.seed(360)
  x = matrix(rnorm(240), 4, 60)
  colnames(x) = 1:ncol(x)
  rownames(x) = 1:nrow(x)
  return(x)
}

#compute dcScores for each method
getScoreList <- function(){
  x = getData()
  cond = getCondition()

  scorelist = lapply(getInfMethods(), function (m) dcScore(x, cond, m))
  names(scorelist) = getInfMethods()

  return(scorelist)
}

#compute matrices holding the results of tests (p-values, probs or orig scores)
getTestMatrices <- function() {
  #generate test matrices
  testmats <- lapply(getScoreList(), function(s) dcTest(s, getData(), getCondition()))

  return(testmats)
}

test_that('Testing calls work', {
  for (m in setdiff(getInfMethods(), c('ebcoexpress', 'diffcoex'))) {
    expect_is(dcAdjust(getTestMatrices()[[!!m]]), 'matrix')
  }

  if ('ebcoexpress' %in% getInfMethods())
    expect_is(dcAdjust(getTestMatrices()[['ebcoexpress']]), 'matrix')
  expect_is(dcAdjust(getTestMatrices()[['diffcoex']]), 'matrix')
})

test_that('Testing attribute changes', {
  for (m in setdiff(getInfMethods(), c('ebcoexpress', 'diffcoex'))) {
    expect_output(str(attr(dcAdjust(getTestMatrices()[[!!m]]), 'dc.test')), regexp = 'adj')
  }
})
