library(dcanr)
library(Biobase)
library(edgeR)
library(SummarizedExperiment)

context('Statistical tests ')

#data to test methods
set.seed(360)
x <- matrix(rnorm(40), 4, 10)
colnames(x) = 1:ncol(x)
rownames(x) = 1:nrow(x)
cond <- rep(1:2, each = ncol(x) / 2)

#pick a representative set - exclude any method that uses perm test during check
getInfMethods <- function() {
  infmethods = c('zscore', 'diffcoex', 'ebcoexpress')
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

scorelist = lapply(getInfMethods(), function (m) dcScore(x, cond, m))
names(scorelist) = getInfMethods()

#generate test matrices
testmats <- lapply(scorelist, function(s) dcTest(s, x, cond))

test_that('Testing calls work', {
  for (m in getInfMethods()) {
    expect_is(dcAdjust(testmats[[!!m]]), 'matrix')
  }
})

test_that('Testing attribute changes', {
  for (m in setdiff(getInfMethods(), c('diffcoex', 'ebcoexpress'))) {
    expect_output(str(attr(dcAdjust(testmats[[!!m]]), 'dc.test')), regexp = 'adj')
  }
  expect_output(str(attr(dcAdjust(testmats[['diffcoex']]), 'dc.test')), regexp = 'none')
  if ('ebcoexpress' %in% getInfMethods())
    expect_output(str(attr(dcAdjust(testmats[['ebcoexpress']]), 'dc.test')), regexp = 'none')
})
