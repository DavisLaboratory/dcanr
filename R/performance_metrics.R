#' Performance metrics to evaluate classification
#'
#' @description Quantify the performance of a classification algorithm.
#'   Predictions and truth both have to be binary.
#'
#' @param pred a logical or numeric, where 0 and FALSE represent control, and, 1
#'   and TRUE represent cases
#' @param obs a logical or numeric, where 0 and FALSE represent control, and, 1
#'   and TRUE represent cases
#' @param perf.method a character, specifying the method to use. Available
#'   methods can be accessed using \code{perfMethods}
#' @param ... additional parameters to methods. see details
#'
#' @details The F-measure requires the beta parameter which can be specified
#'   using \code{f.beta} which defaults to 1 thereby computing the F1-measure.
#'
#' @return a numeric, representing the performance
#' @seealso \code{\link{perfMethods}}
#'
#' @examples
#' pred <- sample(0:1, 100, replace = TRUE, prob = c(0.75, 0.25))
#' obs <- sample(0:1, 100, replace = TRUE, prob = c(0.75, 0.25))
#'
#' #compute the F1 and F2 scores
#' f1 <- performanceMeasure(pred, obs)
#' f2 <- performanceMeasure(pred, obs, f.beta = 2)
#'
#' @export
performanceMeasure <- function(pred, obs, perf.method = 'f.measure', ...) {
  stopifnot(all(c(pred, obs) %in% c(0:1)) || is.logical(pred) || is.logical(obs))
  stopifnot(perf.method %in% perfMethods())

  #compute contingency matrix
  contVec = contingencyVec(pred, obs)
  #compute performance measure
  measure = do.call(perf.method, c(list(quote(contVec)), list(...)))

  return(measure)
}

#' @title Get names of performance metric methods
#'
#' @description Returns a list of performance metrics
#' @return names of methods implemented
#' @export
#'
#' @examples
#' perfMethods()
perfMethods <- function() {
  return(sort(c('f.measure', 'accuracy', 'precision', 'recall', 'MCC', 'AC')))
}

#create a contingency vector
contingencyVec <- function(pred, obs) {
  P = sum(obs)
  N = sum(1 - obs)
  TP = sum(obs & pred)
  TN = sum(!obs & !pred)
  FP = sum(!obs & pred)
  FN = sum(obs & !pred)

  return(c('TP' = TP, 'TN' = TN, 'FP' = FP, 'FN' = FN, 'P' = P, 'N' = N))
}

#F-score
f.measure <- function(contVec, f.beta = 1) {
  f.beta = f.beta ^ 2
  Fsc = with(as.list(contVec), (1 + f.beta) * TP / ((1 + f.beta) * TP + f.beta * FN + FP))
  return(Fsc)
}

#accuracy
accuracy <- function(contVec) {
  return(with(as.list(contVec), (TP + TN)/(P + N)))
}

#precision
precision <- function(contVec) {
  return(with(as.list(contVec), TP/(TP + FP)))
}

#recall
recall <- function(contVec) {
  return(with(as.list(contVec), TP/P))
}

#Mathew's correlation
MCC <- function(contVec) {
  mcc = with(as.list(contVec), (TP*TN - FP*FN)/sqrt((TP+FP) * (TP+FN) * (TN+FP) * (TN+FN)))
  return(mcc)
}

#Approximate correlation
AC <- function(contVec) {
  acp = with(as.list(contVec), TP / (TP + FN) + TP / (TP + FP) + TN / (TN + FP) + TN / (TN + FN))
  acp = acp * 0.25
  ac = 2 * (acp - 0.5)
  return(ac)
}
