#' @include inference_methods.R
#' @importFrom foreach foreach %:% %dopar%
#' @importFrom graphics par
#' @importFrom stats anova binomial glm p.adjust pnorm qchisq var
#' @import methods
NULL

#' @title Compute scores from differential association analysis
#' @description Implementations and wrappers for existing implementations for
#'   methods inferring differential associations/co-expression. This method
#'   requires a matrix of expression and a binary condition to compute the
#'   differential association scores for all pairs of features (genes).
#'   Applications are not limited to analysis of gene expression data and may be
#'   used for differential associations in general.
#'
#' @param emat a matrix, Matrix, data.frame, ExpressionSet, SummarizedExperiment
#'   or DGEList
#' @param condition a numeric, (with 1's and 2's representing a binary
#'   condition), a factor with 2 levels or a character representing 2 conditions
#' @param dc.method a character, representing the method to use. Use
#'   \code{dcMethods()} to get a list of methods
#' @param ... see details
#'
#' @details When using data from sequencing experiments, make sure appropriate
#'   filtering for low counts and data transformation has been performed. Not
#'   doing so will affect estimation of correlation coefficients which most
#'   methods rely on.
#'
#'   Additional method specific parameters can be supplied to the function.
#'   \code{cor.method} can be set to either 'pearson' (default) or 'spearman' to
#'   determine the method to use for estimating correlations. These are the two
#'   measures currently supported in the package. We recommend using the
#'   'spearman' correlation when dealing with sequencing data.
#'
#'   The beta parameter in the DiffCoEx method can be specified using
#'   \code{diffcoex.beta} (defaults to 6). This enable soft thresholding of
#'   correlations similar to WGCNA.
#'
#'   EBcoexpress specific parameters include \code{ebcoexpress.seed}
#'   representing the randomisation seed for the EM algorithm,
#'   \code{ebcoexpress.useBWMC} (defaults to \code{TRUE}) representing whether
#'   to use the bi-weight mid-correlation coefficient or not and
#'   \code{ebcoexpress.plot} which plots the diagnostic plots if set to
#'   \code{TRUE} (defaults to \code{FALSE}).
#'
#' @name dcScore
#' @return a matrix, of scores/statistics representing differential
#'   associations; p-values will be returned if FTGI is used and posterior
#'   probabilities if EBcoexpress is used.
#' @seealso \code{\link{dcMethods}}
#'
#' @examples
#' x <- matrix(rnorm(60), 2, 30)
#' cond <- rep(1:2, 15)
#' dcScore(x, cond) #defaults to zscore
#' dcScore(x, cond, dc.method = 'diffcoex')
#'
#' @exportMethod dcScore
setGeneric(
  name = 'dcScore',
  def = function(emat,
                 condition,
                 dc.method,
                 ...) {
    standardGeneric('dcScore')
  }
)

#' @rdname dcScore
setMethod(
  f = 'dcScore',
  signature = c('matrix', 'ANY', 'ANY'),
  definition = function(emat,
                        condition,
                        dc.method = 'zscore',
                        ...) {
    stopifnot(missing(dc.method) || dc.method %in% dcMethods())
    stopifnot(is.vector(condition) | is.factor(condition))
    stopifnot(length(condition) == ncol(emat))
    stopifnot(length(levels(as.factor(condition))) == 2)

    #default method
    if (missing(dc.method)) {
      dc.method = 'zscore'
    }

    #names of features
    if (is.null(rownames(emat))) {
      warning('rownames missing, assigning default rownames')
      rownames(emat) = 1:nrow(emat)
    }

    #names of observations
    if (is.null(colnames(emat))) {
      warning('colnames missing, assigning default colnames')
      rownames(emat) = 1:nrow(emat)
    }

    #convert conditions to numeric
    condition = as.numeric(as.factor(condition))

    #run method
    # scmat = eval(parse(text = methodmap[dc.method]))(emat, condition, ...)
    scmat = do.call(methodmap[dc.method, 'scoref'], list(quote(emat), quote(condition), ...))
    diag(scmat) = NA #no self loops: none of the methods are capable

    #add method name to attributes
    attr(scmat, 'dc.method') = dc.method

    return(scmat)
  }
)

#' @rdname dcScore
setMethod(
  f = 'dcScore',
  signature = c('Matrix', 'ANY', 'ANY'),
  definition = function(emat,
                        condition,
                        dc.method = 'zscore',
                        ...) {
    #default method
    if (missing(dc.method)) {
      dc.method = 'zscore'
    }

    emat = Matrix::as.matrix(emat)
    scmat = callGeneric(emat, condition, dc.method, ...)

    return(scmat)
  }
)

#' @rdname dcScore
setMethod(
  f = 'dcScore',
  signature = c('data.frame', 'ANY', 'ANY'),
  definition = function(emat,
                        condition,
                        dc.method = 'zscore',
                        ...) {
    stopifnot(all(sapply(emat, class) %in% 'numeric'))

    #default method
    if (missing(dc.method)) {
      dc.method = 'zscore'
    }

    emat = as.matrix(emat)
    scmat = callGeneric(emat, condition, dc.method, ...)

    return(scmat)
  }
)

#' @rdname dcScore
setMethod(
  f = 'dcScore',
  signature = c('ExpressionSet', 'ANY', 'ANY'),
  definition = function(emat,
                        condition,
                        dc.method = 'zscore',
                        ...) {
    if (!requireNamespace("Biobase", quietly = TRUE)){
      stop('\'Biobase\' needed for this function to work. Please install it.', call. = FALSE)
    }

    #default method
    if (missing(dc.method)) {
      dc.method = 'zscore'
    }

    emat = Biobase::exprs(emat)
    scmat = callGeneric(emat, condition, dc.method, ...)

    return(scmat)
  }
)

#' @rdname dcScore
setMethod(
  f = 'dcScore',
  signature = c('SummarizedExperiment', 'ANY', 'ANY'),
  definition = function(emat,
                        condition,
                        dc.method = 'zscore',
                        ...) {
    if (!requireNamespace("SummarizedExperiment", quietly = TRUE)){
      stop('\'SummarizedExperiment\' needed for this function to work. Please install it.', call. = FALSE)
    }

    #default method
    if (missing(dc.method)) {
      dc.method = 'zscore'
    }

    emat = SummarizedExperiment::assay(emat)
    scmat = callGeneric(emat, condition, dc.method, ...)

    return(scmat)
  }
)

#' @rdname dcScore
setMethod(
  f = 'dcScore',
  signature = c('DGEList', 'ANY', 'ANY'),
  definition = function(emat,
                        condition,
                        dc.method = 'zscore',
                        ...) {
    #default method
    if (missing(dc.method)) {
      dc.method = 'zscore'
    }

    emat = emat$counts
    scmat = methods::callGeneric(emat, condition, dc.method, ...)

    return(scmat)
  }
)
