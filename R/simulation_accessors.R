#'@title Get data and conditions from a given knock-down (KD)
#'@description Retrieves the simulated expression matrix and sample
#'  classification for a specific knock-down experiment.
#'
#'@param simulation a list, storing data and results generated from simulations
#'@param cond.name a character, indicating the knock-down to use to derive
#'  conditions. Multiple knock-downs (KDs) are performed per simulation. If
#'  \code{NULL}, the first KD is chosen
#'@param truth.type a character, specifying which level of the true network to
#'  retrieve: 'association' (default), 'influence' or 'direct'
#'@param full a logical, indicating whether genes associated with the condition
#'  should be excluded. Defaults to \code{FALSE} and is recommended
#'
#'@details Genes discarded when \code{full} is \code{FALSE} are those that are
#'  solely dependent on the condition. These genes are discarded from the
#'  analysis to focus on those that are differentially co-expressed, not
#'  coordinately co-expressed.
#'
#'  The names of all genes knocked-out can be retrieved using
#'  \code{getConditionNames}.
#'
#'  The direct, influence and association networks represent different levels of
#'  true differential networks. The direct network contains differential
#'  regulatory interactions present in the original network. The influence
#'  network includes upstream interactions and the association network includes
#'  non-causative differential interactions.
#'
#'@return a list, containing \code{emat}, a matrix representing the expression
#'  data, \code{condition}, a numeric containing the classification of samples,
#'  and , \code{condition_c}, a numeric containing the expression levels of the
#'  KD gene (continuous condition) for \code{getSimData}; the names of all genes
#'  that are KD for \code{getConditionNames}; and an adjacency matrix for
#'  \code{getTrueNetwork}.
#'@seealso \code{\link{dcScore}}
#'
#' @examples
#' data(sim102)
#' KDs <- getConditionNames(sim102)
#'
#' #get simulated data
#' simdata <- getSimData(sim102, KDs[2])
#' cond <- simdata$condition
#' emat <- simdata$emat
#' zscores <- dcScore(emat, cond)
#'
#' #get the true network to evaluate against
#' truenet <- getTrueNetwork(sim102, KDs[2], truth.type = 'association')
#'
#'@describeIn getSimData get the expression matrix and sample classification
#'
#'@export
getSimData <- function(simulation, cond.name = NULL, full = FALSE) {
  #subset the condition matrix with the condition of interest
  emat = simulation$data
  condmat = attr(emat, 'classf')

  #checks for cond.name
  if (is.null(cond.name)) {
    cond.name = getConditionNames(simulation)[1]
  }
  stopifnot(cond.name %in% getConditionNames(simulation))
  cond = condmat[cond.name, ]

  #filter out genes directly dependent on the condition
  cond_c = emat[cond.name, ]
  if (!full) {
    noncond = attr(simulation$triplets, 'condcoex')[cond.name, ]
    noncond = names(noncond)[noncond == 0]
    noncond = setdiff(noncond, cond.name)

    emat = emat[noncond, ]
  } else {
    warning('directly regulated genes should be discarded for differential co-expression analysis')
  }

  return(list('emat' = emat, 'condition' = cond, 'condition_c' = cond_c))
}

#' @describeIn getSimData get names of the conditions (KDs)
#' @export
getConditionNames <- function(simulation) {
  return(rownames(attr(simulation$data, 'classf')))
}

#' @describeIn getSimData get the true differential network
#' @export
getTrueNetwork <- function(simulation, cond.name = NULL, truth.type = c('association', 'influence', 'direct')[1], full = FALSE) {
  truth.type = stringr::str_to_title(truth.type)
  if (is.null(cond.name)) {
    cond.name = getConditionNames(simulation)[1]
  }
  genes = rownames(getSimData(simulation, cond.name, full)$emat)

  #initialise adjacency matrix
  adjmat = matrix(0, length(genes), length(genes))
  rownames(adjmat) = colnames(adjmat) = genes

  #populate the matrix
  diffpairs = simulation$triplets
  diffpairs = diffpairs[diffpairs$cond %in% cond.name, ]
  diffpairs = diffpairs[diffpairs[, truth.type],  c('TF', 'Target')]
  adjmat[as.matrix(diffpairs)] = adjmat[as.matrix(diffpairs)[, 2:1], drop = FALSE] = 1

  return(adjmat)
}



