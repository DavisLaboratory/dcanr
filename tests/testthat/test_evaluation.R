library(dcanr)
library(igraph)
library(stringr)

context('Simulation accessors ')

data("sim102")
#pick a representative set
infmethods = c('zscore', 'diffcoex', 'dicer', 'ebcoexpress')
if (!require('EBcoexpress')) {
  infmethods = setdiff(infmethods, 'ebcoexpress')
}
if (!require('COSINE')) {
  infmethods = setdiff(infmethods, 'ecf')
}

test_that('Result of pipeline ', {
  for(m in infmethods) {
    res = dcPipeline(sim102, dc.func = m)
    expect_equal(length(res), 2, info = m)
    expect_equal(as.numeric(lapply(res, function(x) length(V(x)))), rep(nrow(sim102$data), 2), info = m)
    expect_equal(as.character(lapply(res, class)), rep('igraph', 2), info = m)
  }
})

test_that('Evaluation of pipeline results ', {
  for(m in infmethods) {
    nets = dcPipeline(sim102, dc.func = m)
    res = dcEvaluate(sim102, dclist = nets, perf.method = 'f.measure', combine = TRUE)
    expect_equal(length(res), 1, info = m)
    res = dcEvaluate(sim102, dclist = nets, perf.method = 'f.measure', combine = FALSE)
    expect_equal(length(res), length(nets), info = m)
    expect_equal(names(res), names(nets), info = m)
  }
})

test_that('Retrieve precomputed results ', {
  for(m in dcMethods()) {
    res = dcPipeline(sim102, dc.func = m, precomputed = TRUE)
    expect_equal(length(res), 2, info = m)
    expect_equal(as.numeric(lapply(res, function(x) length(V(x)))), rep(nrow(sim102$data), 2), info = m)
    expect_equal(as.character(lapply(res, class)), rep('igraph', 2), info = m)
  }
})

test_that('Evaluation of precomputed results ', {
  for(m in dcMethods()) {
    nets = dcPipeline(sim102, dc.func = m, precomputed = TRUE)
    res = dcEvaluate(sim102, dclist = nets, perf.method = 'f.measure', combine = TRUE)
    expect_equal(length(res), 1, info = m)
    res = dcEvaluate(sim102, dclist = nets, perf.method = 'f.measure', combine = FALSE)
    expect_equal(length(res), length(nets), info = m)
    expect_equal(names(res), names(nets), info = m)
  }
})
