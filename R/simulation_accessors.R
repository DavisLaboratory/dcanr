#' @title Get data and conditions from a given knock-out (KO)
#' @description Retrieves the simulated expression matrix and sample
#'   classification for a specific knock-out experiment.
#'
#' @param simulation a list, storing data and results generated from simulations
#' @param cond.name a character, indicating the knock-out to use to derive
#'   conditions. Multiple knock-outs (KOs) are performed per simulation. If
#'   \code{NULL}, the first KO is chosen
#' @param full a logical, indicating whether genes associated with the condition
#'   should be excluded. Defaults to \code{FALSE} and is recommended
#'
#' @details Genes discarded when \code{full} is \code{FALSE} are those that are
#'   solely dependent on the condition. These genes are discarded from the
#'   analysis to focus on those that are differentially co-expressed, not
#'   co-ordinately co-expressed.
#'
#'   The names of all genes knocked-out can be retrieved using
#'   \code{getConditionNames}
#'
#' @return a list, containing \code{emat}, a matrix representing the expression
#'   data and \code{condition}, a numeric containing the classification of
#'   samples
#' @seealso \code{\link{dcScore}}
#'
#' @examples
#' data(sim102)
#' KOs <- getConditionNames(sim102)
#' simdata <- getSimulatedData(sim102, KOs[2])
#' cond <- simdata$condition
#' emat <- simdata$emat
#' zscores <- dcScore(emat, cond)
#'
#' @describeIn getSimData get the expression matrix and sample
#'   classification
#'
#' @export
getSimData <- function(simulation, cond.name = NULL, full = FALSE) {
  #subset the condition matrix with the condition of interest
  emat = simulation$data
  condmat = attr(emat, 'classf')

  #checks for cond.name
  if (is.null(cond.name)) {
    cond.name = rownames(condmat)[1]
  }
  stopifnot(cond.name %in% getConditionNames(simulation))
  cond = condmat[cond.name, ]

  #filter out genes directly dependent on the condition
  if (!full) {
    noncond = attr(simulation$triplets, 'condcoex')[cond.name, ]
    noncond = names(noncond)[noncond == 0]

    emat = emat[noncond, ]
  } else {
    warning('directly regulated genes should be discarded for differential co-expression analysis')
  }

  return(list('emat' = emat, 'condition' = cond))
}

#' @describeIn getSimData get names of the conditions (KOs)
#' @export
getConditionNames <- function(simulation) {
  return(rownames(attr(simulation$data, 'classf')))
}
