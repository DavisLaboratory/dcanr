library(dcevalr)
library(igraph)
library(stringr)

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

