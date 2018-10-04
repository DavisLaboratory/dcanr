library(dcevalr)
library(igraph)
library(stringr)
library(Matrix)

context('Simulation accessors ')

data("sim102")

test_that('Network retreival works ', {
  for (type in c('direct', 'influence', 'association')) {
    #networks retrieved
    nets = lapply(getConditionNames(sim102), function(x)
      getTrueNetwork(sim102, x, truth.type = type))
    nets = lapply(nets, graph_from_adjacency_matrix, mode = 'undirected')
    netlens = lapply(nets, function(x)
      length(E(x)))

    #true networks
    infnet = sim102$infnet
    truesize = sum(get.edge.attribute(infnet, str_to_title(type)))
    expect_equal(sum(unlist(netlens)), truesize)
  }
})

test_that('Edges are the same when using the full dataset ', {
  for (type in c('direct', 'influence', 'association')) {
    #networks retrieved
    nets = lapply(getConditionNames(sim102), function(x)
      sum(getTrueNetwork(sim102, x, truth.type = type)))
    nets_full = lapply(getConditionNames(sim102), function(x)
      sum(getTrueNetwork(sim102, x, TRUE, truth.type = type)))

    expect_equal(sum(unlist(nets)), sum(unlist(nets_full)))
  }
})

mat2vec <- function(m) {
  m = pmax(m, t(m))
  v = m[upper.tri(m)]
  attr(v, 'feature.names') = rownames(m) #store names to enable reconstruction
  attr(v, 'mat.attrs') = attributes(m) #store names to enable reconstruction

  return(v)
}

test_that('Retrieved network vector is same as stored network ', {
  vecnames = c('truthvecDct', 'truthvecInf', 'truthvec')
  names(vecnames) = c('direct', 'influence', 'association')
  for (type in names(vecnames)) {
    #networks retrieved
    nets = lapply(getConditionNames(sim102), function(x){
      mat = getTrueNetwork(sim102, x, truth.type = type)
      mat = mat[sort(rownames(mat)), sort(rownames(mat))]
      return(mat2vec(mat))
    })
    nets = unlist(nets)

    expect_equal(length(nets), length(as.matrix(sim102$scores)[, vecnames[type]]))
    expect_equal(nets, as.matrix(sim102$scores)[, vecnames[type]])
  }
})
