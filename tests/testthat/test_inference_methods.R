library(dcanr)
library(Biobase)
library(edgeR)
library(SummarizedExperiment)
library(Matrix)

context('Inference methods ')

#data to test methods
set.seed(360)
x <- matrix(rnorm(240), 4, 60)
colnames(x) = 1:ncol(x)
rownames(x) = 1:nrow(x)
cond <- rep(1:2, each = 30)

#define methods
infmethods = dcMethods()
if (!require('EBcoexpress')) {
  infmethods = setdiff(infmethods, 'ebcoexpress')
}
if (!require('COSINE')) {
  infmethods = setdiff(infmethods, 'ecf')
}

test_that('Inference method calls work', {
  expect_is(dcScore(x, cond), 'matrix')
  expect_is(dcScore(x, cond, dc.method = 'zscore'), 'matrix')
  for (m in infmethods) {
    expect_is(dcScore(x, cond, !!m), 'matrix')
  }

  try(detach('package:EBcoexpress', unload = T), silent = TRUE)
  try(detach('package:mclust', unload = T), silent = TRUE)
  try(detach('package:minqa', unload = T), silent = TRUE)
  if(require(EBcoexpress))
    expect_is(dcScore(x, cond, 'ebcoexpress'), 'matrix')
  else
    expect_error(dcScore(x, cond, 'ebcoexpress'), '\'EBcoexpress\' needed for this function to work')

  expect_error(dcScore(x, cond, 'fakemethod'), 'not TRUE')
  expect_error(dcScore(x, cond, cor.method = 'kendall'), 'not TRUE')
  expect_error(dcScore(x, c(cond, 3)), 'not TRUE')
  expect_error(dcScore(x, c(cond, 2)), 'not TRUE')
})

test_that('Correct dimensions of results', {
  for (m in infmethods) {
    expect_equal(dim(dcScore(x, cond, !!m)), c(4, 4))
  }
})

test_that('Attributes attached to results', {
  for (m in infmethods) {
    expect_equal(attr(dcScore(x, cond, !!m), 'dc.method'), m)
  }
})

test_that('Diagonals are NAs', {
  for (m in infmethods) {
    expect_equal(sum(is.na(diag(dcScore(x, cond, !!m)))), nrow(x))
  }
})

test_that('Row and column names are the same', {
  for (m in infmethods) {
    expect_equal(rownames(dcScore(x, cond, !!m)), colnames(dcScore(x, cond, !!m)))
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

