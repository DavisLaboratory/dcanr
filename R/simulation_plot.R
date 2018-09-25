#' @include statistical_tests.R
NULL

#' @title Plot source and true differential networks from simulations
#' @description Plots either the source network or the true differential network
#'   for all KOs performed in the simulation. KO nodes are coloured with their
#'   resulting differential networks coloured accordingly.
#'
#' @inheritParams getSimData
#' @param what a character, indicating which network to plot, 'source'
#'   (default), 'direct', 'influence' or 'association'
#' @param ... additional parameters to \code{plot.igraph}
#'
#' @details The direct, influence and association networks represent different
#'   levels of true differential networks. The direct network contains
#'   differential regulatory interactions present in the original network. The
#'   influence network includes upstream interactions and the association
#'   network includes non-causative differential interactions.
#'
#' @return a plot of the network
#'
#' @seealso \code{\link[igraph]{plot.igraph}}
#'
#' @examples
#' data(sim102)
#' plotSimNetwork(sim102)
#' plotSimNetwork(sim102, what = 'direct')
#' plotSimNetwork(sim102, what = 'influence')
#' plotSimNetwork(sim102, what = 'association')
#'
#' @export
plotSimNetwork <- function(simulation, what = c('source', 'direct', 'influence', 'association')[1], ...) {
  lyt = simulation$netlayout
  if (what == 'source') {
    igraph::plot.igraph(simulation$staticnet, layout = lyt, ...)
  } else {
    what = stringr::str_to_title(what)
    infnet = simulation$infnet
    E(infnet)$color[!igraph::get.edge.attribute(infnet, what)] = NA
    igraph::plot.igraph(infnet, layout = lyt, ...)
  }
}
