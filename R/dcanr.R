#' @details
#' There are three categories of functions available
#' \enumerate{
#'   \item Differential co-expression methods (DC) - These functions are used to
#'   perform a differential co-expression analysis on experimental data with binary
#'   conditions.
#'   \item Functions to evaluate DC methods - These functions are used to evaluate
#'   methods implemented in the package and novel methods on simulated data.
#'   Expression data is simulated for 2 conditions, wild-type and knock-down of
#'   given genes.
#'   \item By-products of implementations
#' }
#'
#' @section Differential co-expression methods (DC):
#' \itemize{
#'    \item{\code{\link{dcMethods}}}
#'    \item{\code{\link{dcScore}}}
#'    \item{\code{\link{dcTest}}}
#'    \item{\code{\link{dcAdjust}}}
#'    \item{\code{\link{dcNetwork}}}
#' }
#'
#' @section Functions to evaluate DC methods:
#' Accessors of simulated data:
#' \itemize{
#'    \item{\code{\link{getConditionNames}}}
#'    \item{\code{\link{getSimData}}}
#'    \item{\code{\link{getTrueNetwork}}}
#'    \item{\code{\link{plotSimNetwork}}}
#' }
#'
#' Functions for evaluating inference methods
#' \itemize{
#'    \item{\code{\link{dcPipeline}}}
#'    \item{\code{\link{dcEvaluate}}}
#' }
#'
#' @section By-products of implementations:
#' These are functions used in the package but have further uses in general:
#' \itemize{
#'    \item{\code{\link{cor.pairs}}} - a faster implementation of pairwise
#'    correlation computation
#'    \item{\code{\link{mi.ap}}} - pairwise computation of mutual information
#'    MI with data discretisation performed using adaptive partitioning
#'    \item{\code{\link{perfMethods}}} - available performance metrics
#'    \item{\code{\link{performanceMeasure}}} - performance measures of prediction
#'    algorithms. Predictions have to be binary
#' }
#'
#' @keywords internal classif graphs
#'
"_PACKAGE"
