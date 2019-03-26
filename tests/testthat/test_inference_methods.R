library(dcanr)

context('Inference methods ')

#select methods that can be run
getInfMethods <- function() {
  infmethods = c('zscore', 'diffcoex', 'ebcoexpress', 'ecf', 'ggm-based')
  # infmethods = dcMethods()
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

#compute dcScores for each method
getData <- function(){
  #data to test methods
  set.seed(360)
  x = matrix(rnorm(240), 4, 60)
  colnames(x) = 1:ncol(x)
  rownames(x) = 1:nrow(x)

  return(x)
}

test_that('Inference method calls work', {
  expect_is(dcScore(getData(), getCondition()), 'matrix')
  expect_is(dcScore(getData(), getCondition(), dc.method = 'zscore'), 'matrix')
  for (m in getInfMethods()) {
    expect_is(dcScore(getData(), getCondition(), !!m), 'matrix')
  }

  try(detach('package:EBcoexpress', unload = TRUE), silent = TRUE)
  try(detach('package:mclust', unload = TRUE), silent = TRUE)
  try(detach('package:minqa', unload = TRUE), silent = TRUE)
  if(require('EBcoexpress'))
    expect_is(dcScore(getData(), getCondition(), 'ebcoexpress'), 'matrix')
  else
    expect_error(dcScore(getData(), getCondition(), 'ebcoexpress'), '\'EBcoexpress\' needed for this function to work')

  expect_error(dcScore(getData(), getCondition(), 'fakemethod'), 'not TRUE')
  expect_error(dcScore(getData(), getCondition(), cor.method = 'kendall'), 'not TRUE')
  expect_error(dcScore(getData(), c(getCondition(), 3)), 'not TRUE')
  expect_error(dcScore(getData(), c(getCondition(), 2)), 'not TRUE')
})

test_that('Correct dimensions of results', {
  for (m in getInfMethods()) {
    expect_equal(dim(dcScore(getData(), getCondition(), !!m)), c(4, 4))
  }
})

test_that('Attributes attached to results', {
  for (m in getInfMethods()) {
    expect_equal(attr(dcScore(getData(), getCondition(), !!m), 'dc.method'), m)
  }
})

test_that('Diagonals are NAs', {
  for (m in getInfMethods()) {
    expect_equal(sum(is.na(diag(dcScore(getData(), getCondition(), !!m)))), nrow(getData()))
  }
})

test_that('Row and column names are the same', {
  for (m in getInfMethods()) {
    expect_equal(rownames(dcScore(getData(), getCondition(), !!m)), colnames(dcScore(getData(), getCondition(), !!m)))
  }
})

test_that('Different condition types', {
  expect_is(dcScore(getData(), getCondition()), 'matrix')
  cond_char = LETTERS[1:2][getCondition()]
  expect_is(dcScore(getData(), cond_char), 'matrix')
  expect_is(dcScore(getData(), as.factor(cond_char)), 'matrix')
  expect_error(dcScore(getData(), as.matrix(getCondition())), 'not TRUE')
})

test_that('Different matrix types', {
  dge = edgeR::DGEList(counts = round(getData() * 10) + 100)
  eset = Biobase::ExpressionSet(assayData = getData())
  se = SummarizedExperiment::SummarizedExperiment(getData())
  Mat = Matrix::Matrix(getData())

  expect_is(dcScore(getData(), getCondition()), 'matrix')
  expect_is(dcScore(as.data.frame(getData()), getCondition()), 'matrix')
  expect_is(dcScore(dge, getCondition()), 'matrix')
  expect_is(dcScore(eset, getCondition()), 'matrix')
  expect_is(dcScore(se, getCondition()), 'matrix')
})

