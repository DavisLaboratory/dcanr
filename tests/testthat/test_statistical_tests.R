library(dcanr)

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

#compute dcScores for each method
getScoreList <- function(){
  #data to test methods
  set.seed(360)
  x = matrix(rnorm(240), 4, 60)
  colnames(x) = 1:ncol(x)
  rownames(x) = 1:nrow(x)
  cond = rep(1:2, each = 30)

  scorelist = lapply(getInfMethods(), function (m) dcScore(x, cond, m))
  names(scorelist) = getInfMethods()

  return(scorelist)
}

#compute matrices holding the results of tests (p-values, probs or orig scores)
getTestMatrices <- function() {
  #generate test matrices
  testmats <- lapply(scorelist, function(s) dcTest(s, x, cond))

  return(testmats)
}

test_that('Testing calls work', {
  for (m in getInfMethods()) {
    expect_is(dcTest(getScoreList[[!!m]], x, cond), 'matrix')
  }
})

test_that('Correct dimensions of results', {
  for (m in getInfMethods()) {
    expect_equal(dim(getTestMatrices[[!!m]]), c(4, 4))
  }
})

test_that('Attributes attached to results', {
  for (m in getInfMethods()) {
    expect_output(str(attr(getTestMatrices[[!!m]], 'dc.test')), regexp = '.')
  }
})

test_that('Diagonals are NAs', {
  for (m in getInfMethods()) {
    expect_equal(sum(is.na(diag(getTestMatrices[[!!m]]))), nrow(x))
  }
})

test_that('Values are in range for statistical tests', {
  for (m in setdiff(getInfMethods(), 'diffcoex')) {
    expect_gte(min(getTestMatrices[[!!m]], na.rm = TRUE), 0)
    expect_lte(max(getTestMatrices[[!!m]], na.rm = TRUE), 1)
  }
})
