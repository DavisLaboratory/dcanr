#' @include multtest_adjust.R
NULL

#' @title Generate a differential network from a DC analysis
#' @description Threshold the results from a differential co-expression analysis
#'   and create a differential network.
#'
#' @inheritParams dcTest
#' @param dcpvals a matrix or NULL, raw or adjusted p-values resulting from
#'   \code{dcTest} or \code{dcAdjust} respectively. Should be left NULL only if
#'   method is EBcoexpress or DiffCoEx
#' @param thresh a numeric, threshold to apply. If \code{NULL}, defaults to 0.1
#'   for methods that generate a p-value, 0.9 for posterior probabilities from
#'   EBcoexpress and 0.1 on the absolute score from DiffCoEx
#'
#' @details No extra arguments required for this function. The ellipsis are used
#'   to allow flexibility in pipelines.
#'
#' @name dcNetwork
#' @return an igraph object, representing the differential network. Scores are
#'   added as edge attributes with the name 'score'
#' @seealso \code{\link{dcScore}}, \code{\link{dcTest}}, \code{\link{dcAdjust}}
#'
#' @examples
#' #create data
#' set.seed(360)
#' x <- matrix(rnorm(120), 4, 30)
#' cond <- rep(1:2, 15)
#'
#' #perform analysis - z-score
#' zscores <- dcScore(x, cond)
#' pvals <- dcTest(zscores, emat = x, condition = cond)
#' pvals <- dcAdjust(pvals, p.adjust, method = 'fdr')
#' ig <- dcNetwork(zscores, pvals, 0.1)
#'
#' #perform analysis - DiffCoEx
#' dcscores <- dcScore(x, cond, dc.method = 'diffcoex')
#' ig <- dcNetwork(dcscores, thresh = 0.001)
#'
#' #plot the resulting differential co-expression network
#' igraph::plot.igraph(ig)
#'
#' @export
dcNetwork <- function(dcscores, dcpvals = NULL, thresh = NULL, ...) {
  if (!all(c('dc.method', 'dc.method') %in% names(attributes(dcscores)))) {
    stop('Please ensure dcscores has not been modified')
  }
  #if statistical tests are performed for method, check pval matrix
  if (methodmap[attr(dcscores, 'dc.method'), 'adjust']) {
    stopifnot(!is.null(dcpvals))
    if (!all(c('dc.test', 'dc.method') %in% names(attributes(dcpvals)))) {
      stop('Please ensure dcpvals has not been modified')
    }
  } else {
    dcpvals = dcscores
  }
  #default thresh if not provided
  dc.method = attr(dcpvals, 'dc.method')
  if (is.null(thresh)) {
    warning('default thresholds being selected')
    thresh = methodmap[dc.method, 'default_thresh']
  }

  #perform thresholding
  proc_annot = attributes(dcpvals)
  proc_annot = proc_annot[!grepl('dim', names(proc_annot))]
  if (dc.method %in% c('ebcoexpress', 'diffcoex', 'ldgm')) {
    dcpvals = abs(dcpvals) > thresh
  } else {
    dcpvals = dcpvals < thresh
  }

  #return results as a igraph
  diag(dcpvals) = FALSE
  dcscores[!dcpvals] = 0
  ig = igraph::graph_from_adjacency_matrix(dcscores, mode = 'undirected', weighted = 'score')

  #visualisation properties
  #node sizes based on degree
  V(ig)$size = log(igraph::degree(ig, mode = 'out') + 3) * 2
  V(ig)$color = '#F7F7F7B3'
  if (length(E(ig)) > 0){
    minmax = stats::quantile(abs(E(ig)$score), 0.75, na.rm = TRUE)
    E(ig)$color = colfunc(-minmax, minmax, RColorBrewer::brewer.pal(11, 'PRGn'))(E(ig)$score)
    E(ig)$color = stringr::str_replace(E(ig)$color, 'FF$', 'B3')
  }
  attributes(ig) = c(attributes(ig), proc_annot)

  return(ig)
}

colfunc <- function(start, end, colvec) {
  col = circlize::colorRamp2(seq(start, end, length.out = length(colvec)), colvec)
  return(col)
}
