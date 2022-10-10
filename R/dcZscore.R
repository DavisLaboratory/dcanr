#' DC analysis using the z-score method
#'
#' This function packs the entire DC analysis pipeline using the z-score method.
#' It simplifies the implementation of the analysis and increases the
#' flexibility of the analysis (not just limited to all pairwise comparisons).
#'
#' @param from a character vector, with the names of nodes from which
#'   comparisons need to be performed.
#' @param to a character vector, with the names of nodes to which comparisons
#'   need to be performed.
#' @param fdrthresh a numeric, specifying the FDR cutoff to apply to the
#'   inferred network.
#' @param cor.method a character, either 'spearman' (default) or 'pearson'
#'   specifying the correlation computation method to use.
#'
#' @return an igraph object, containing the differential coexpression network.
#' @export
#'
#' @inheritParams dcScore
#'
#' @examples
#' x <- matrix(rnorm(60), 10, 30)
#' rownames(x) = 1:10
#' cond <- rep(1:2, 15)
#' dcZscore(x, cond)
#' dcZscore(x, cond, to = 1:2)
#'
dcZscore <-
  function(emat,
           condition,
           from = NULL,
           to = NULL,
           fdrthresh = 0.1,
           cor.method = c('spearman', 'pearson')) {

  #check params
  cor.method = match.arg(cor.method)
  condition = as.factor(condition)
  stopifnot(length(levels(condition)) == 2)
  stopifnot(fdrthresh <=1 | fdrthresh >= 0)
  if (is.null(from)) {
    from = rownames(emat)
  } else{
    stopifnot(all(from %in% rownames(emat)))
  }
  if (is.null(to)) {
    to = rownames(emat)
  } else{
    stopifnot(all(to %in% rownames(emat)))
  }

  #separate group expression
  expr1 = emat[, !is.na(condition) & as.numeric(condition) == 1, drop = FALSE]
  expr2 = emat[, !is.na(condition) & as.numeric(condition) == 2, drop = FALSE]

  # apply the Fisher transformation
  const = 1
  if (cor.method %in% 'spearman'){
    const = 1.06
  }

  #compute z-scores
  r1 = cor(t(expr1[from, , drop = FALSE]), t(expr1[to, , drop = FALSE]), method = cor.method)
  r2 = cor(t(expr2[from, , drop = FALSE]), t(expr2[to, , drop = FALSE]), method = cor.method)
  z1 = atanh(r1)
  z2 = atanh(r2)
  z = (z1 - z2)/sqrt(const/(sum(as.numeric(condition) == 1) - 3) + const/(sum(as.numeric(condition) == 2) - 3))
  z[is.nan(z)] = 0
  z[is.infinite(z)] = 0

  #compute p-vals
  pvals = z.test(z)

  #adjust p-vals
  fdr = reshape2::melt(pvals)
  fdr[, 1:2] = sapply(fdr[, 1:2], as.character)
  fdr$min = apply(fdr[, 1:2], 1, min)
  fdr$max = apply(fdr[, 1:2], 1, max)
  fdr = fdr[, -(1:2)]
  fdr = unique(fdr)
  fdr$fdr = stats::p.adjust(fdr$value, method = 'fdr')
  fdr = fdr[fdr$fdr < fdrthresh, -1]
  colnames(fdr) = c('from', 'to', 'fdr')

  #annotate df
  fdr = annotateZdf(fdr, z, 'z.score')
  fdr = annotateZdf(fdr, r1, paste0('cor.', levels(condition)[1]))
  fdr = annotateZdf(fdr, r2, paste0('cor.', levels(condition)[2]))

  #create network
  dcnet = igraph::graph_from_data_frame(fdr, directed = FALSE)

  #visualisation properties
  #node sizes based on degree
  V(dcnet)$size = log(igraph::degree(dcnet, mode = 'out') + 3) * 2
  V(dcnet)$color = '#F7F7F7B3'
  if (length(E(dcnet)) > 0){
    minmax = stats::quantile(abs(E(dcnet)$z.score), 0.75, na.rm = TRUE)
    E(dcnet)$color = colfunc(-minmax, minmax, RColorBrewer::brewer.pal(11, 'PRGn'))(E(dcnet)$z.score)
    E(dcnet)$color = stringr::str_replace(E(dcnet)$color, 'FF$', 'B3')
  }

  return(dcnet)
}

annotateZdf <- function(df, mat, value.name) {
  #initialise matrix
  allg = union(rownames(mat), colnames(mat))
  vals = matrix(NA, length(allg), length(allg))
  rownames(vals) = colnames(vals) = allg

  #make a complete matrix
  valdf = reshape2::melt(mat)
  valdf[, 1:2] = sapply(valdf[, 1:2], as.character)
  vals[as.matrix(valdf[, 1:2])] = valdf$value
  vals[as.matrix(valdf[, 2:1])] = valdf$value

  #annotate df
  df$value = vals[as.matrix(df[, 1:2])]
  colnames(df)[ncol(df)] = value.name

  return(df)
}

# if (cor.method %in% 'spearman') {
#   eranks = t(apply(emat, 1, rank))
# } else {
#   eranks = emat
# }
#
# # .Call(stats:::C_cov, x, NULL, 1, FALSE)
# # cov.wt is faster for pairwise cors than cov
#
# #compute pairwise weighted correlation
# r1 = cov.wt(t(emat), wt = wt1, cor = TRUE)$cor
# r2 = cov.wt(t(emat), wt = wt2, cor = TRUE)$cor
# #subset edges
# r1 = r1[from, to, drop = FALSE]
# r2 = r2[from, to, drop = FALSE]
