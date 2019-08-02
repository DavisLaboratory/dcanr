library(dcanr)
library(Biobase)
library(edgeR)
library(SummarizedExperiment)
library(Matrix)

context('Inference methods ')

#data to test methods
set.seed(360)
x <- matrix(rnorm(40), 4, 10)
colnames(x) = 1:ncol(x)
rownames(x) = 1:nrow(x)
cond <- rep(1:2, each = ncol(x) / 2)

getInfMethods <- function() {
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

  return(infmethods)
}

test_that('Inference method calls work', {
  expect_is(dcScore(x, cond), 'matrix')
  expect_is(dcScore(x, cond, dc.method = 'zscore'), 'matrix')
  for (m in getInfMethods()) {
    expect_is(dcScore(x, cond, !!m, lambda = 0.5), 'matrix')
  }

  expect_error(dcScore(x, cond, 'ldgm'), 'Need to specify either lambda')

  if (!'ebcoexpress' %in% getInfMethods()) {
    expect_error(dcScore(x, cond, 'ebcoexpress'), 'needed for this function to work')
  }
  if (!'ecf' %in% getInfMethods()) {
    expect_error(dcScore(x, cond, 'ecf'), 'needed for this function to work')
  }
  if (!'ggm-based' %in% getInfMethods()) {
    expect_error(dcScore(x, cond, 'ggm-based'), 'needed for this function to work')
  }

  expect_error(dcScore(x, cond, 'fakemethod'), 'not TRUE')
  expect_error(dcScore(x, cond, cor.method = 'kendall'), 'should be one of ')
  expect_error(dcScore(x, c(cond, 3)), 'not TRUE')
  expect_error(dcScore(x, c(cond, 2)), 'not TRUE')
})

test_that('Correct dimensions of results', {
  for (m in getInfMethods()) {
    expect_equal(dim(dcScore(x, cond, !!m, lambda = 0.5)), c(4, 4))
  }
})

test_that('Attributes attached to results', {
  for (m in getInfMethods()) {
    expect_equal(attr(dcScore(x, cond, !!m, lambda = 0.5), 'dc.method'), m)
  }
})

test_that('Diagonals are NAs', {
  for (m in getInfMethods()) {
    expect_equal(sum(is.na(diag(dcScore(x, cond, !!m, lambda = 0.5)))), nrow(x))
  }
})

test_that('Row and column names are the same', {
  for (m in getInfMethods()) {
    expect_equal(rownames(dcScore(x, cond, !!m, lambda = 0.5)), colnames(dcScore(x, cond, !!m, lambda = 0.5)))
  }
})

test_that('Different condition types', {
  expect_is(dcScore(x, cond), 'matrix')
  condchar = LETTERS[1:2][cond]
  expect_is(dcScore(x, condchar), 'matrix')
  expect_is(dcScore(x, as.factor(condchar)), 'matrix')
  expect_error(dcScore(x, as.matrix(cond)), 'not TRUE')
})

test_that('Different matrix types', {
  dge = DGEList(counts = round(x * 10) + 100)
  eset = ExpressionSet(assayData = x)
  se = SummarizedExperiment(x)
  Mat = Matrix(x)

  expect_is(dcScore(x, cond), 'matrix')
  expect_is(dcScore(as.data.frame(x), cond), 'matrix')
  expect_is(dcScore(dge, cond), 'matrix')
  expect_is(dcScore(eset, cond), 'matrix')
  expect_is(dcScore(se, cond), 'matrix')
})
