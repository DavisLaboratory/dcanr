methodmap = data.frame(
  'scoref' = c(
    'z.score',
    'ggm.score',
    'magic.score',
    'ftgi.score',
    'diffcoex.score',
    'ebcoexpress.score',
    'dicer.score',
    'ecf.score',
    'ent.score',
    'mindy.score'
  ),
  'testf' = c('z.test', c('perm.test', 'no.test')[c(1, 1, 2, 2, 2, 1, 1, 1, 1)]),
  'default_thresh' = c(rep(0.05, 5), 0.95, rep(0.05, 4)),
  row.names = c(
    'zscore',
    'ggm-based',
    'magic',
    'ftgi',
    'diffcoex',
    'ebcoexpress',
    'dicer',
    'ecf',
    'entropy',
    'mindy'
  ),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

#' Get names of differential co-expression methods
#'
#' @description Returns a list of differential co-expression methods
#' @return names of methods implemented
#' @export
#'
#' @examples
#' dcMethods()
dcMethods <- function() {
  return(sort(rownames(methodmap)))
}

#' Fast pairwise correlation estimation
#'
#' @description Fast estimation of pairwise correlation coefficients.
#'
#' @param emat a numeric matrix
#' @param cor.method a character, specifying the method to use for estimation.
#'   Possible values are 'pearson' (default) and 'spearman'
#'
#' @return a numeric matrix with estimated correlation coefficients
#' @export
#'
#' @examples
#' x <- matrix(rnorm(200), 100, 2)
#' cor.pairs(x)
#' cor.pairs(x, cor.method = 'spearman')
cor.pairs <- function(emat, cor.method = 'pearson') {
  stopifnot(cor.method %in% c('pearson', 'spearman'))

  if (cor.method %in% 'spearman') {
    #spearman is pearson on ranked data
    emat = apply(emat, 2, rank)
  }

  design = matrix(1, nrow(emat), 1)
  QR = qr(design)
  E = qr.qty(QR, emat)
  s2 = colMeans(E[-1, ]^2)
  U = t(t(E[-1, ])/sqrt(s2))
  c = crossprod(U)/nrow(U)
  diag(c) = 1  # to correct for numerical inaccuracies

  return(c)
}

z.score <- function(emat, condition, cor.method = 'pearson', ...) {
  expr1 = emat[, condition == 1, drop = FALSE]
  expr2 = emat[, condition == 2, drop = FALSE]

  # apply the Fisher transformation
  const = 1
  if (cor.method %in% 'spearman'){
    const = 1.06
  }

  r1 = cor.pairs(t(expr1), cor.method)
  r2 = cor.pairs(t(expr2), cor.method)
  z1 = atanh(r1)
  z2 = atanh(r2)
  z = (z1 - z2)/sqrt(const/(sum(condition == 1) - 3) + const/(sum(condition == 2) - 3))
  z[is.nan(z)] = 0
  z[is.infinite(z)] = 0

  #add run parameters as attributes
  attributes(z) = c(attributes(z),
                    'cor.method' = cor.method,
                    'call' = match.call())

  return(z)
}

# paper: Quantifying differential gene connectivity between disease states for objective identification of disease-relevant genes
ggm.score <- function(emat, condition, ...) {
  if (!requireNamespace("GeneNet", quietly = TRUE)){
    stop('\'GeneNet\' needed for this function to work. Please install it.', call. = FALSE)
  }

  # apply GGM method from genenet and estimate partial correlations
  expr1 = emat[, condition == 1, drop = FALSE]
  expr2 = emat[, condition == 2, drop = FALSE]
  rmat1 = GeneNet::ggm.estimate.pcor(t(expr1))
  rmat2 = GeneNet::ggm.estimate.pcor(t(expr2))

  # get posterior probabilities
  df1 = GeneNet::network.test.edges(rmat1, plot = FALSE, verbose = FALSE)
  pmat1 = matrix(0, ncol(rmat1), ncol(rmat1))
  pmat1[apply(df1[, 2:3], 2, as.numeric)] = df1$prob
  pmat1[apply(df1[, 3:2], 2, as.numeric)] = df1$prob

  df2 = GeneNet::network.test.edges(rmat2, plot = FALSE, verbose = FALSE)
  pmat2 = matrix(0, ncol(rmat2), ncol(rmat2))
  pmat2[apply(df2[, 2:3], 2, as.numeric)] = df2$prob
  pmat2[apply(df2[, 3:2], 2, as.numeric)] = df2$prob
  colnames(pmat1) = rownames(pmat1) = colnames(pmat2) = rownames(pmat2) = colnames(rmat1)

  # assign zero-values a non-zero value
  nz = min(min(pmat1[pmat1 != 0]), min(pmat2[pmat2 != 0]))
  pmat1[pmat1 == 0] = nz
  pmat2[pmat2 == 0] = nz

  # compute posterior odds
  logOR = log((pmat1/(1 - pmat1))/(pmat1/(1 - pmat2)))

  #add run parameters as attributes
  attributes(logOR) = c(attributes(logOR),
                        'call' = match.call())

  return(logOR)
}

magic.score <- function(emat, condition, cor.method = 'pearson', ...) {
  expr1 = emat[, condition == 1, drop = FALSE]
  expr2 = emat[, condition == 2, drop = FALSE]

  const = 1
  if (cor.method %in% 'spearman'){
    const = 1.06
  }

  r1 = cor.pairs(t(expr1), cor.method)
  r2 = cor.pairs(t(expr2), cor.method)
  z1 = atanh(r1)/sqrt(const/(sum(condition == 1) - 3))
  z2 = atanh(r2)/sqrt(const/(sum(condition == 2) - 3))

  # compute magic score
  tgtsize = round(ncol(emat)/2) - 3
  expz1 = exp(const/sqrt(tgtsize) * 2 * z1)
  expz2 = exp(const/sqrt(tgtsize) * 2 * z2)
  IM = abs((expz1 - 1)/(expz1 + 1)) - abs((expz2 - 1)/(expz2 + 1))
  IM[is.nan(IM)] = 0

  #add run parameters as attributes
  attributes(IM) = c(attributes(IM),
                     'cor.method' = cor.method,
                     'call' = match.call())

  return(IM)
}

ftgi.score <- function(emat, condition, ...) {
  message('Results will be a matrix of p-values, not scores')

  score <- foreach(i = 1:nrow(emat), .combine = cbind) %:%
    foreach(j = 1:nrow(emat), .combine = rbind) %dopar% {
      e1 = as.numeric(emat[i, ])
      e2 = as.numeric(emat[j, ])
      m1 = glm(as.factor(condition) ~ e1 + e2, family = binomial(link = "logit"))
      m2 = glm(as.factor(condition) ~ e1 * e2, family = binomial(link = "logit"))
      sc = anova(m1, m2, test = "Chisq")[2, 5]
      return(sc)
    }
  rownames(score) = colnames(score) = rownames(emat)
  score[is.na(score)] = 1

  #add run parameters as attributes
  attributes(score) = c(attributes(score),
                        'call' = match.call())

  return(score)
}

diffcoex.score <- function(emat, condition, cor.method = 'pearson', diffcoex.beta = 6, ...) {
  expr1 = emat[, condition == 1, drop = FALSE]
  expr2 = emat[, condition == 2, drop = FALSE]

  # calculate correlations
  r1 = cor.pairs(t(expr1), cor.method)
  r2 = cor.pairs(t(expr2), cor.method)
  D = sqrt(0.5 * abs(sign(r1) * r1^2 - sign(r2) * r2^2))
  D = D^diffcoex.beta
  T.ovlap = D %*% D + ncol(D) * D  #calc topological ovlap

  mins = matrix(rep(rowSums(D), ncol(D)), nrow = ncol(D))
  mins = pmin(mins, matrix(rep(colSums(D), each = ncol(D)), nrow = ncol(D)))
  T.ovlap = 1 - (T.ovlap/(mins + 1 - D))

  diag(T.ovlap) = 1

  #add run parameters as attributes
  attributes(T.ovlap) = c(
    attributes(T.ovlap),
    'cor.method' = cor.method,
    'diffcoex.beta' = diffcoex.beta,
    'call' = match.call()
  )

  return(1 - T.ovlap)
}

ebcoexpress.score <- function(emat, condition, ebcoexpress.seed = sample.int(1e6, 1), ebcoexpress.useBWMC = TRUE, ebcoexpress.plot = FALSE, ...) {
  if (!requireNamespace("EBcoexpress", quietly = TRUE)){
    stop('\'EBcoexpress\' needed for this function to work. Please install it.', call. = FALSE)
  }
  if (!requireNamespace("EBarrays", quietly = TRUE)){
    stop('\'EBarrays\' needed for this function to work. Please install it.', call. = FALSE)
  }

  message('Results will be a matrix of posterior probabilities')

  #map functions required by EBcoexpress
  Mclust = mclust::Mclust

  scoremat = matrix(rep(0, nrow(emat)^2), nrow = nrow(emat))
  colnames(scoremat) = rownames(scoremat) = rownames(emat)

  set.seed(ebcoexpress.seed)
  pat = EBarrays::ebPatterns(c("1,1", "1,2"))
  D = NULL
  tryCatch(expr = {
    D = EBcoexpress::makeMyD(emat, condition, useBWMC = ebcoexpress.useBWMC)
  }, error = function(e) {
    warning(e)
  })
  if (is.null(D)) {
    return(scoremat)
  }

  tryCatch(expr = {
    initHP = EBcoexpress::initializeHP(D, condition)
  }, error = function(e) {
    stop(e, '\nApplication of this method will require loading the EBcoexpress library')
  })

  result1 = NULL
  tryCatch(expr = {
    oout = EBcoexpress::ebCoexpressOneStep(D, condition, pat, initHP)
    result1 = oout$POSTPROBS
  }, error = function(e) {
    warning(e)
  })

  if (is.null(result1)) {
    return(scoremat)
  }
  ppbDC1 = 1 - result1[, 1]

  # diagnostic plots if required
  if (ebcoexpress.plot) {
    par(mfrow = c(1, 2))
    EBcoexpress::priorDiagnostic(D, condition, oout, 1)
    EBcoexpress::priorDiagnostic(D, condition, oout, 2)
    par(mfrow = c(1, 1))
  }

  # convert to matrix
  corpairs = plyr::ldply(stringr::str_split(names(ppbDC1), "~"))
  corpairs["prob"] = ppbDC1
  for (i in 1:nrow(corpairs)) {
    g1 = as.character(corpairs[i, 1])
    g2 = as.character(corpairs[i, 2])
    scoremat[g1, g2] = scoremat[g2, g1] = corpairs[i, 3]
  }

  #add run parameters as attributes
  attributes(scoremat) = c(
    attributes(scoremat),
    'ebcoexpress.seed' = ebcoexpress.seed,
    'ebcoexpress.useBWMC' = ebcoexpress.useBWMC,
    'call' = match.call()
  )

  return(scoremat)
}

dicer.score <- function(emat, condition, cor.method = 'pearson', ...) {
  expr1 = emat[, condition == 1, drop = FALSE]
  expr2 = emat[, condition == 2, drop = FALSE]

  # apply the Fisher transformation
  r1 = cor.pairs(t(expr1), cor.method)
  r2 = cor.pairs(t(expr2), cor.method)
  mu1 = mean(r1)
  mu2 = mean(r2)
  var1 = var(as.numeric(r1))
  var2 = var(as.numeric(r2))
  tscore = ((r1 - r2) - (mu1 - mu2))/sqrt(var1 + var2)
  diag(tscore) = 0

  #add run parameters as attributes
  attributes(tscore) = c(attributes(tscore),
                         'cor.method' = cor.method,
                         'call' = match.call())

  return(tscore)
}

ecf.score <- function(emat, condition, ...) {
  if (!requireNamespace("COSINE", quietly = TRUE)){
    stop('\'COSINE\' needed for this function to work. Please install it.', call. = FALSE)
  }

  expr1 = emat[, condition == 1, drop = FALSE]
  expr2 = emat[, condition == 2, drop = FALSE]

  # apply the Fisher transformation
  ecfscore = COSINE::diff_gen(t(expr1), t(expr2))[[2]]
  colnames(ecfscore) = rownames(ecfscore) = rownames(emat)

  #add run parameters as attributes
  attributes(ecfscore) = c(attributes(ecfscore),
                           'call' = match.call())

  return(ecfscore)
}

ent.score <- function(emat, condition, cor.method = "pearson", ...) {
  expr1 = emat[, condition == 1, drop = FALSE]
  expr2 = emat[, condition == 2, drop = FALSE]

  # calculate correlations
  r1 = cor.pairs(t(expr1), cor.method)
  r2 = cor.pairs(t(expr2), cor.method)
  rall = cor.pairs(t(emat), cor.method)
  I1 = -((1 + abs(r1))/2 * log((1 + abs(r1))/2) + (1 - abs(r1))/2 * log((1 - abs(r1))/2))
  I2 = -((1 + abs(r2))/2 * log((1 + abs(r2))/2) + (1 - abs(r2))/2 * log((1 - abs(r2))/2))
  Iall = -((1 + abs(rall))/2 * log((1 + abs(rall))/2) + (1 - abs(rall))/2 * log((1 - abs(rall))/2))
  ent = (I1 + I2)/2 - Iall

  diag(ent) = 0

  #add run parameters as attributes
  attributes(ent) = c(attributes(ent),
                      'cor.method' = cor.method,
                      'call' = match.call())

  return(ent)
}

mi.ap.single <- function(x, y) {
  x = rank(x, ties.method = "min")
  y = rank(y, ties.method = "min")
  sz = length(x)
  testgrid = expand.grid(0, 0, sz, sz, length(x), length(x), length(x))
  finalgrid = c()
  cnames = c("sx", "sy", "ex", "ey", "count", "cx", "cy")
  colnames(testgrid) = cnames
  chithresh = qchisq(0.95, 3)

  while (!is.null(testgrid)) {
    newtestgrid = c()
    for (i in 1:nrow(testgrid)) {
      # partition grid
      ptsx = seq(testgrid$sx[i], testgrid$ex[i], length.out = 3)[-3]
      ptsy = seq(testgrid$sy[i], testgrid$ey[i], length.out = 3)[-3]
      gridsz = ptsx[2] - ptsx[1]
      df = expand.grid(ptsx, ptsy)
      df = cbind(df, df + gridsz, 0, 0, 0)
      colnames(df) = cnames

      # grid counts and chi-squared test
      ntot = testgrid[i, "count"]
      df$count = apply(df, 1, function(t) {
        sum(x > t["sx"] & x <= t["ex"] & y > t["sy"] & y <= t["ey"])
      })
      df$cx = apply(df, 1, function(t) {
        sum(x > t["sx"] & x <= t["ex"])
      })
      df$cy = apply(df, 1, function(t) {
        sum(y > t["sy"] & y <= t["ey"])
      })

      chival = sum((ntot/4 - df$count)^2/(ntot/4))

      # store/partition further
      if (chival <= chithresh) {
        finalgrid = rbind(finalgrid, testgrid[i, ])
      } else {
        newtestgrid = rbind(newtestgrid, df)
      }

      # remove approve any grid with only 1 element
      is1 = newtestgrid$count == 1
      finalgrid = rbind(finalgrid, newtestgrid[is1, ])
      newtestgrid = newtestgrid[!is1 & newtestgrid$count != 0, ]
    }

    testgrid = newtestgrid
  }

  prob = finalgrid[, 5:7]/length(x)
  mi = sum(prob[, 1] * log(prob[, 1]/(prob[, 2] * prob[, 3])))

  return(mi)
}

#' Mutual information using adaptive partitioning
#'
#' @description Computes the mutual information between all pairs of variables
#'   in the matrix (along the columns). Variables are discretised using the
#'   adaptive partitioning algorithm
#'
#' @param mat a numeric matrix
#'
#' @return matrix of pairwise mutual information estimates
#' @export
#'
#' @examples
#' x <- matrix(rnorm(200), 100, 2)
#' mi.ap(x)
mi.ap <- function(mat) {
  if (is.null(colnames(mat))){
    colnames(mat) = 1:ncol(mat)
  }

  gpairs = expand.grid(colnames(mat), colnames(mat))
  gpairs = unique(t(apply(gpairs, 1, sort)))
  gpairs = as.data.frame(gpairs)

  # calculate MIs for gene pairs - try parallelizing later
  gpairs$mi = apply(gpairs, 1, function(x) {
    mi.ap.single(mat[, x[1]], mat[, x[2]])
  })

  # convert result to matrix
  mimat = reshape2::acast(gpairs, V1 ~ V2, value.var = "mi")
  mimat = mimat[colnames(mat), colnames(mat)]
  mimat[is.na(mimat)] = -1
  mimat = pmax(mimat, t(mimat))

  return(mimat)
}

mindy.score <- function(emat, condition, ...) {
  expr1 = emat[, condition == 1, drop = FALSE]
  expr2 = emat[, condition == 2, drop = FALSE]

  m1 = mi.ap(t(expr1))
  m2 = mi.ap(t(expr2))
  diag(m1) = diag(m2) = 0

  scoremat = m1 - m2
  scoremat[is.na(scoremat)] = 0

  #add run parameters as attributes
  attributes(scoremat) = c(attributes(scoremat),
                           'call' = match.call())

  return(scoremat)
}
