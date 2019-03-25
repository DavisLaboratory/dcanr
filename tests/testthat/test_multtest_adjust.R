library(dcanr)
library(Biobase)
library(edgeR)
library(SummarizedExperiment)

context('Statistical tests ')

#data to test methods
set.seed(360)
x = matrix(rnorm(240), 4, 60)
colnames(x) = 1:ncol(x)
rownames(x) = 1:nrow(x)
cond = rep(1:2, each = 30)

#run all methods and store results
infmethods = dcMethods()
if (!require('EBcoexpress')) {
  infmethods = setdiff(infmethods, 'ebcoexpress')
}
if (!require('COSINE')) {
  infmethods = setdiff(infmethods, 'ecf')
}
if (!require('GeneNet')) {
  infmethods = setdiff(infmethods, 'ggm-based')
}

scorelist = lapply(infmethods, function (m) dcScore(x, cond, m))
names(scorelist) = infmethods

#generate test matrices
testmats <- lapply(scorelist, function(s) dcTest(s, x, cond))

test_that('Testing calls work', {
  for (m in setdiff(infmethods, c('ebcoexpress', 'diffcoex'))) {
    expect_is(dcAdjust(testmats[[!!m]]), 'matrix')
  }

  if ('ebcoexpress' %in% infmethods)
    expect_is(dcAdjust(testmats[['ebcoexpress']]), 'matrix')
  expect_is(dcAdjust(testmats[['diffcoex']]), 'matrix')
})

test_that('Testing attribute changes', {
  for (m in setdiff(infmethods, c('ebcoexpress', 'diffcoex'))) {
    expect_output(str(attr(dcAdjust(testmats[[!!m]]), 'dc.test')), regexp = 'adj')
  }
})
