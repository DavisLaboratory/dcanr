#' @include multtest_adjust.R performance_metrics.R statistical_tests.R
#' @importFrom igraph E V E<- V<-
NULL

#' @title Run a DC pipeline on a simulation
#' @description Run a differential co-expression pipeline on data from a
#'   simulation experiment. A default pipeline can be used which consists of
#'   methods in the package or custom pipelines can be provided.
#'
#' @inheritParams getSimData
#' @param dc.func a function or character. Character represents one of the
#'   method names from \code{dcMethods} which is run with the default settings.
#'   A function can be used to provide custom processing pipelines (see details)
#' @param precomputed a logical, indicating whether the precomputed inference
#'   should be used or a new one computed (default FALSE)
#' @param continuous a logical, indicating whether binary or continuous
#'   conditions should be used (default FALSE). No methods implemented currently
#'   use continuous conditions. This is to allow custom methods that require
#'   continuous conditions
#' @param ... additional parameters to \code{dc.func}
#'
#' @details If \code{dc.func} is a character, the existing methods in the
#'   package will be run with their default parameters. The pipeline is as such:
#'   dcScore -> dcTest -> dcAdjust -> dcNetwork, resulting in a igraph object.
#'   Parameters to the independent processing steps can also be provided to this
#'   function as shown in the examples.
#'
#'   If \code{precomputed} is TRUE while \code{dc.func} is a character,
#'   pre-computed results will be used. These can then be evaluated using
#'   \code{dcEvaluate}.
#'
#'   Custom pipelines need to be coded into a function which can then be
#'   provided instead of a character. Functions must have the following
#'   structure:
#'
#'   \code{ function(emat, condition, ...) }
#'
#'   They must return either an igraph object or an adjacency matrix stored in a
#'   base R 'matrix' or the S4 'Matrix' class, containing all genes in the
#'   expression matrix 'emat'. See examples for how the in-built functions are
#'   combined into a pipeline.
#'
#' @return a list of igraphs, representing the differential network for each
#'   independent condition (knock-out).
#'
#' @seealso \code{\link[igraph]{plot.igraph}}, \code{\link{dcScore}},
#'   \code{\link{dcTest}}, \code{\link{dcAdjust}}, \code{\link{dcNetwork}},
#'   \code{\link{dcMethods}}
#'
#' @examples
#' data(sim102)
#'
#' #run a standard pipeline
#' resStd <- dcPipeline(sim102, dc.func = 'zscore')
#'
#' #run a standard pipeline and specify params
#' resParam <- dcPipeline(sim102, dc.func = 'zscore', cor.method = 'pearson')
#'
#' #retrieve pre-computed results
#' resPrecomputed <- dcPipeline(sim102, dc.func = 'zscore', precomputed = TRUE)
#'
#' #run a custom pipeline
#' analysisInbuilt <- function(emat, condition, dc.method = 'zscore', ...) {
#'   #compute scores
#'   score = dcScore(emat, condition, dc.method, ...)
#'   #perform statistical test
#'   pvals = dcTest(score, emat, condition, ...)
#'   #adjust tests for multiple testing
#'   adjp = dcAdjust(pvals, ...)
#'   #threshold and generate network
#'   dcnet = dcNetwork(score, adjp, ...)
#'
#'   return(dcnet)
#' }
#' resCustom <- dcPipeline(sim102, dc.func = analysisInbuilt)
#'
#' plot(resCustom[[1]])
#'
#' @export
dcPipeline <- function(simulation, dc.func = 'zscore', precomputed = FALSE, continuous = FALSE, ...) {
  if (is.character(dc.func)) {
    stopifnot(dc.func %in% dcMethods())
    dc.method = dc.func

    if (precomputed) {
      return(retrieveResult(simulation, dc.method))
    }

    dc.func <- function(emat, condition, ...) {
      return(analysisInbuilt(emat, condition, dc.method, ...))
    }
  }
  stopifnot(is.function(dc.func))

  #perform inference for each condition
  infnets = lapply(getConditionNames(simulation), function (cond) {
    #get data from simulation
    simdata = getSimData(simulation, cond.name = cond, full = FALSE)
    emat = simdata$emat
    if (continuous) {
      condition = simdata$condition_c
    } else {
      condition = simdata$condition
    }

    #evaluate function
    dcnet = suppressWarnings(dc.func(emat = emat, condition = condition, ...))
    #validate result and convert to matrix
    dcnet = validateResult(emat, dcnet)
    dcnet = annotateig(simulation, dcnet)

    return(dcnet)
  })
  names(infnets) = getConditionNames(simulation)

  return(infnets)
}

#' @title Evaluate performance of DC methods on simulations
#' @description Quantify the performance of a differential co-expression
#'   pipeline on simulated data.
#'
#' @inheritParams getSimData
#' @inheritParams performanceMeasure
#' @param dclist a list of igraphs, produced using \code{dcPipeline}
#' @param combine a logical, indicating whether differential networks from
#'   independent knock-outs should be treated as a single inference or
#'   independent inferences (defaults to \code{TRUE})
#' @param ... additional parameters to be passed on to the performance metric
#'   method (see \code{performanceMeasure})
#'
#' @return a numeric, representing the performance metric. A single value if
#'   \code{combine = TRUE} and a named vector otherwise.
#'
#' @seealso \code{\link{dcPipeline}}, \code{\link{performanceMeasure}},
#'   \code{\link{perfMethods}}
#'
#' @examples
#' data(sim102)
#'
#' #run a standard pipeline
#' resStd <- dcPipeline(sim102, dc.func = 'zscore')
#' dcEvaluate(sim102, resStd)
#' dcEvaluate(sim102, resStd, combine = FALSE)
#'
#' @export
dcEvaluate <-
  function(simulation,
           dclist,
           truth.type = c('association', 'influence', 'direct'),
           perf.method = 'f.measure',
           combine = TRUE,
           ...) {
  truth.type = match.arg(truth.type)

  #retrieve the truth vectors
  truthvecs = lapply(names(dclist), function(cond.name) {
    adjmat = getTrueNetwork(simulation, cond.name, truth.type = truth.type)
    genes = rownames(getSimData(simulation, cond.name)$emat)
    genes = sort(genes)
    adjmat = adjmat[genes, genes]

    return(mat2vec(adjmat))
  })

  #convert adjacency matrices to vectors
  infvecs = lapply(names(dclist), function(cond.name) {
    #convert igraphs to adjacency matrices, subset and sort rows and cols
    dcnet = dclist[[cond.name]]
    adjmat = igraph::as_adj(dcnet, sparse = FALSE)
    genes = rownames(getSimData(simulation, cond.name)$emat)
    genes = sort(genes)
    adjmat = adjmat[genes, genes]

    return(mat2vec(adjmat))
  })

  #combine the independent vectors to a single vector
  if (combine) {
    infvecs = list(unlist(infvecs))
    truthvecs = list(unlist(truthvecs))
  }

  #compute performance metrics
  metrics = mapply(performanceMeasure,
                   infvecs,
                   truthvecs,
                   MoreArgs = c(list('perf.method' = perf.method), list(...)))

  #name the results if not combined
  if (!combine) {
    names(metrics) = names(dclist)
  }

  return(metrics)
}

validateResult <- function(emat, res) {
  #check for correct output format
  if (!class(res) %in% c('matrix', 'igraph', 'Matrix')) {
    stop('result of the function must be either a matrix or igraph object')
  }

  #if matrix, check that it is binary
  if (class(res) %in% 'matrix') {
    if (!all(table(res) %in% 0:1)) {
      stop('result of the function must be a binary matrix')
    }

    #helps with reordering vertices and clearing all attributes
    res = igraph::graph_from_adjacency_matrix(res)
  }

  #if Matrix, check that it is binary
  if (class(res) %in% 'Matrix') {
    if (!all(Matrix::as.array(res) %in% 0:1)) {
      stop('result of the function must be a binary matrix')
    }

    #helps with reordering vertices and clearing all attributes
    res = igraph::graph_from_adjacency_matrix(res)
  }

  #if igraph, check all nodes in data are present
  if (!setequal(V(res)$name, rownames(emat))) {
    stop('result of the function must contain all genes in the data matrix')
  }

  #convert result to adjacency matrix
  adjmat = igraph::as_adjacency_matrix(res, sparse = FALSE)
  adjmat = adjmat[rownames(emat), rownames(emat)]
  ig = igraph::graph_from_adjacency_matrix(adjmat, mode = 'undirected')

  return(ig)
}

#add attributes to the graph to enable plotting
annotateig <- function(simulation, ig) {
  #merge inferred graph and true graph
  truenet = igraph::as.undirected(simulation$infnet, mode = 'collapse', edge.attr.comb = 'first')
  infnet = igraph::intersection(truenet, ig, byname = TRUE)
  ig = igraph::difference(ig, infnet, byname = TRUE)
  infnet = igraph::union(infnet, ig, byname = TRUE)

  #add plotting attributes for inferred graph
  E(infnet)$color[is.na(E(infnet)$color)] = '#88888898'
  E(infnet)$width[is.na(E(infnet)$width)] = E(infnet)$width[!is.na(E(infnet)$width)][1]
  E(infnet)$arrow.size = 0

  #add layout as attributes
  V(infnet)$x = simulation$netlayout[, 1]
  V(infnet)$y = simulation$netlayout[, 2]

  return(infnet)
}

analysisInbuilt <- function(emat, condition, dc.method = 'zscore', ...) {
  #compute scores
  score = dcScore(emat, condition, dc.method, ...)
  #perform statistical test
  pvals = dcTest(score, emat, condition, ...)
  #adjust tests for multiple testing
  adjp = dcAdjust(pvals, ...)
  #threshold and generate network
  dcnet = dcNetwork(score, adjp, ...)

  return(dcnet)
}

retrieveResult <- function(simulation, dc.method = 'zscore') {
  scores = Matrix::as.matrix(simulation$scores)
  scores = split(scores[, dc.method], scores[, 'bimid'])
  names(scores) = getConditionNames(simulation)

  #convert vectors to matrices and then to igraphs
  dcnets = lapply(getConditionNames(simulation), function(b) {
    gnames = rownames(getSimData(simulation, b)$emat)

    #add attributes to transfer
    attr(scores[[b]], 'feature.names') = sort(gnames)
    attr(scores[[b]], 'mat.attrs') = list(
      'dim' = rep(length(gnames), 2),
      'dimnames' = list(sort(gnames), sort(gnames))
    )
    scmat = vec2mat(scores[[b]])
    scmat = scmat[gnames, gnames]
    diag(scmat) = 0

    #convert to igraphs
    ig = igraph::graph_from_adjacency_matrix(scmat, mode = 'undirected', weighted = 'score')
    ig = annotateig(simulation, ig)
    return(ig)
  })
  names(dcnets) = getConditionNames(simulation)

  return(dcnets)
}
