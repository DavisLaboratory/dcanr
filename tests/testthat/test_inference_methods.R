library(dcevalr)
library(Biobase)
library(edgeR)
library(SummarizedExperiment)

context('Inference methods ')

#data to test methods
set.seed(360)
x <- matrix(rnorm(240), 4, 60)
colnames(x) = 1:ncol(x)
rownames(x) = 1:nrow(x)
cond <- rep(1:2, each = 30)

test_that('Inference method calls work', {
  expect_is(diffScore(x, cond), 'matrix')
  expect_is(diffScore(x, cond, dc.method = 'zscore'), 'matrix')
  for (m in setdiff(dcMethods(), 'ebcoexpress')) {
    expect_is(diffScore(x, cond, !!m), 'matrix')
  }

  try(detach('package:EBcoexpress', unload = T), silent = TRUE)
  try(detach('package:mclust', unload = T), silent = TRUE)
  try(detach('package:minqa', unload = T), silent = TRUE)
  expect_error(diffScore(x, cond, 'ebcoexpress'), 'loading the EBcoexpress library')
  library(EBcoexpress)
  expect_is(diffScore(x, cond, 'ebcoexpress'), 'matrix')

  expect_error(diffScore(x, cond, 'fakemethod'), 'not TRUE')
  expect_error(diffScore(x, cond, cor.method = 'kendall'), 'not TRUE')
  expect_error(diffScore(x, c(cond, 3)), 'not TRUE')
  expect_error(diffScore(x, c(cond, 2)), 'not TRUE')
})

test_that('Correct dimensions of results', {
  for (m in dcMethods()) {
    expect_equal(dim(diffScore(x, cond, !!m)), c(4, 4))
  }
})

test_that('Attributes attached to results', {
  for (m in dcMethods()) {
    expect_equal(attr(diffScore(x, cond, !!m), 'dc.method'), m)
  }
})

test_that('Diagonals are NAs', {
  for (m in dcMethods()) {
    expect_equal(sum(is.na(diag(diffScore(x, cond, !!m)))), nrow(x))
  }
})

test_that('Row and column names are the same', {
  for (m in dcMethods()) {
    expect_equal(rownames(diffScore(x, cond, !!m)), colnames(diffScore(x, cond, !!m)))
  }
})

test_that('Different condition types', {
  expect_is(diffScore(x, cond), 'matrix')
  condchar = LETTERS[1:2][cond]
  expect_is(diffScore(x, condchar), 'matrix')
  expect_is(diffScore(x, as.factor(condchar)), 'matrix')
  expect_error(diffScore(x, as.matrix(cond)), 'not TRUE')
})

test_that('Different matrix types', {
  dge = DGEList(counts = round(x * 10) + 100)
  eset = ExpressionSet(assayData = x)
  se = SummarizedExperiment(x)

  expect_is(diffScore(x, cond), 'matrix')
  expect_is(diffScore(as.data.frame(x), cond), 'matrix')
  expect_is(diffScore(dge, cond), 'matrix')
  expect_is(diffScore(eset, cond), 'matrix')
  expect_is(diffScore(se, cond), 'matrix')
})

