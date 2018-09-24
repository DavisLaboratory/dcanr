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
scorelist <- lapply(dcMethods(), function (m) diffScore(x, cond, m))
names(scorelist) <- dcMethods()

test_that('Testing calls work', {
  for (m in dcMethods()) {
    expect_is(diffTest(scorelist[[!!m]], x, cond), 'matrix')
  }
})

#generate test matrices
testmats <- lapply(scorelist, function(s) diffTest(s, x, cond))

test_that('Correct dimensions of results', {
  for (m in dcMethods()) {
    expect_equal(dim(testmats[[!!m]]), c(4, 4))
  }
})

test_that('Attributes attached to results', {
  for (m in dcMethods()) {
    expect_output(str(attr(testmats[[!!m]], 'dc.test')), regexp = '.')
  }
})

test_that('Diagonals are NAs', {
  for (m in dcMethods()) {
    expect_equal(sum(is.na(diag(testmats[[!!m]]))), nrow(x))
  }
})

test_that('Values are in range for statistical tests', {
  for (m in setdiff(dcMethods(), 'diffcoex')) {
    expect_gte(min(testmats[[!!m]], na.rm = TRUE), 0)
    expect_lte(max(testmats[[!!m]], na.rm = TRUE), 1)
  }
})
