#' @include multtest_adjust.R
NULL

#' @title Generate a differential network from a DC analysis
#' @description Threshold the results from a differential co-expression analysis
#'   and create a differential network.
#'
#' @inheritParams dcTest
#' @param dcpvals a matrix, raw or adjusted p-values resulting from
#'   \code{dcTest} or \code{dcAdjust} respectively
#' @param thresh a numeric, threshold to apply. If \code{NULL}, defaults to 0.05
#'   for methods that generate a p-value, 0.9 for posterior probabilities from
#'   EBcoexpress and 0.05 on the absolute score from DiffCoEx
#'
#' @name dcNetwork
#' @return an igraph object, representing the differential network. Scores are
#'   added as edge attributes with the name 'score'
#' @seealso \code{\link{dcScore}}, \code{\link{dcTest}}, \code{\link{dcAdjust}}
#'
#' @examples
#' #create data
#' x <- matrix(rnorm(120), 4, 30)
#' cond <- rep(1:2, 15)
#'
#' #perform analysis
#' zscores <- dcScore(x, cond)
#' pvals <- dcTest(zscores, emat = x, condition = cond)
#' pvals <- dcAdjust(pvals, p.adjust, method = 'fdr')
#' ig <- dcNetwork(zscores, pvals, 0.5)
#'
#' \dontrun{
#' igraph::plot.igraph(ig)
#' }
#'
#' @export
dcNetwork <- function(dcscores, dcpvals, thresh = NULL) {
  if (!all(c('dc.test', 'dc.method') %in% names(attributes(dcpvals)))) {
    stop('Please ensure dcpvals has not been modified')
  }
  if (!all(c('dc.method', 'dc.method') %in% names(attributes(dcscores)))) {
    stop('Please ensure dcscores has not been modified')
  }

  #default thresh if not provided
  dc.method = attr(dcpvals, 'dc.method')
  if (is.null(thresh)) {
    warning('default thresholds being selected')
    thresh = methodmap[dc.method, 'default_thresh']
  }

  #perform thresholding
  proc_annot = attributes(dcpvals) [-(1:2)]
  if (dc.method %in% c('ebcoexpress', 'diffcoex')) {
    dcpvals = abs(dcpvals) > thresh
  } else {
    dcpvals = dcpvals < thresh
  }

  #return results as a igraph
  diag(dcpvals) = FALSE
  dcscores[!dcpvals] = 0
  ig = igraph::graph_from_adjacency_matrix(dcscores, mode = 'undirected', weighted = 'score')

  return(ig)
}
