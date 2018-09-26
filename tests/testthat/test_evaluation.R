library(dcevalr)
library(igraph)
library(stringr)
library(EBcoexpress)

context('Simulation accessors ')

data("sim102")

test_that('Result result of pipeline ', {
  for(m in c('zscore', 'diffcoex', 'dicer', 'ebcoexpress')) {
    res = dcPipeline(sim102, dc.func = m)
    expect_equal(length(res), 2, info = m)
    expect_equal(as.numeric(lapply(res, function(x) length(V(x)))), rep(nrow(sim102$data), 2), info = m)
    expect_equal(as.character(lapply(res, class)), rep('igraph', 2), info = m)
  }
})

test_that('Result result of pipeline ', {
  for(m in c('zscore', 'diffcoex', 'dicer', 'ebcoexpress')) {
    nets = dcPipeline(sim102, dc.func = m)
    res = dcEvaluate(sim102, dclist = nets, perf.method = 'f.measure', combine = TRUE)
    expect_equal(length(res), 1, info = m)
    res = dcEvaluate(sim102, dclist = nets, perf.method = 'f.measure', combine = FALSE)
    expect_equal(length(res), length(nets), info = m)
    expect_equal(names(res), names(nets), info = m)
  }
})

