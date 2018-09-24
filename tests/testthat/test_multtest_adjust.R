library(dcevalr)
library(Biobase)
library(edgeR)
library(SummarizedExperiment)
library(EBcoexpress)

context('Statistical tests ')

#data to test methods
set.seed(360)
x <- matrix(rnorm(240), 4, 60)
colnames(x) = 1:ncol(x)
rownames(x) = 1:nrow(x)
cond <- rep(1:2, each = 30)

#run all methods and store results
scorelist <- lapply(dcMethods(), function (m) dcScore(x, cond, m))
names(scorelist) <- dcMethods()

#generate test matrices
testmats <- lapply(scorelist, function(s) dcTest(s, x, cond))

test_that('Testing calls work', {
  for (m in setdiff(dcMethods(), c('ebcoexpress', 'diffcoex'))) {
    expect_is(dcAdjust(testmats[[!!m]]), 'matrix')
  }

  expect_is(dcAdjust(testmats[['ebcoexpress']]), 'matrix')
  expect_is(dcAdjust(testmats[['diffcoex']]), 'matrix')
})

test_that('Testing attribute changes', {
  for (m in setdiff(dcMethods(), c('ebcoexpress', 'diffcoex'))) {
    expect_output(str(attr(dcAdjust(testmats[[!!m]]), 'dc.test')), regexp = 'adj')
  }
})
