library(dcanr)
library(Biobase)
library(edgeR)
library(SummarizedExperiment)

context('Statistical tests ')

#data to test methods
set.seed(360)
x <- matrix(rnorm(240), 4, 60)
colnames(x) = 1:ncol(x)
rownames(x) = 1:nrow(x)
cond <- rep(1:2, each = 30)

#run all methods and store results
infmethods = dcMethods()
if (!require('EBcoexpress')) {
  infmethods = setdiff(infmethods, 'ebcoexpress')
}
if (!require('COSINE')) {
  infmethods = setdiff(infmethods, 'ecf')
}

scorelist <- lapply(infmethods, function (m) dcScore(x, cond, m))
names(scorelist) <- infmethods

test_that('Testing calls work', {
  for (m in infmethods) {
    expect_is(dcTest(scorelist[[!!m]], x, cond), 'matrix')
  }
})

#generate test matrices
testmats <- lapply(scorelist, function(s) dcTest(s, x, cond))

test_that('Correct dimensions of results', {
  for (m in infmethods) {
    expect_equal(dim(testmats[[!!m]]), c(4, 4))
  }
})

test_that('Attributes attached to results', {
  for (m in infmethods) {
    expect_output(str(attr(testmats[[!!m]], 'dc.test')), regexp = '.')
  }
})

test_that('Diagonals are NAs', {
  for (m in infmethods) {
    expect_equal(sum(is.na(diag(testmats[[!!m]]))), nrow(x))
  }
})

test_that('Values are in range for statistical tests', {
  for (m in setdiff(infmethods, 'diffcoex')) {
    expect_gte(min(testmats[[!!m]], na.rm = TRUE), 0)
    expect_lte(max(testmats[[!!m]], na.rm = TRUE), 1)
  }
})
