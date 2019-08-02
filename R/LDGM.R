# Soft thresholding
softThreshold = function(w, t, lambda){
  w = sign(w) * pmax(0, abs(w) - lambda * t);
  return(w)
}

L1GeneralCompositeGradientAccelerated <- function(gradFunc, w, lambda, params) {
  #process options
  maxiter = 500
  optTol = 1e-6
  L = c()

  p = length(w)

  # Compute Evaluate Function
  fg_list = gradFunc(w)
  f = fg_list[[1]]
  g = fg_list[[2]]
  funEvals = 1

  if (length(L) == 0) {
    alpha_max = 1
    alpha = 1
  }
  else{
    alpha_max = L
  }
  alpha = alpha_max

  a = 1
  z = w

  for (i in seq_len(maxiter)) {
    w_old = w
    w_new = softThreshold(w - alpha * g, alpha, lambda)

    # Compute Evaluate Function
    fg_list_new = gradFunc(w)
    f_new = fg_list_new[[1]]
    g_new = fg_list_new[[2]]
    funEvals = funEvals + 1

    phi_T = f_new + sum(lambda * (abs(w_new)))
    mL = f + t(g) %*% (w_new - w) + t(w_new - w) %*% (w_new - w) / (2 * alpha) + sum(lambda * (abs(w_new)))

    if (phi_T > mL) {
      alpha = alpha/2
      w_new = softThreshold(w - alpha * g, alpha, lambda)

      # Compute Evaluate Function
      fg_list_new = gradFunc(w_new)
      f_new = fg_list_new[[1]]
      g_new = fg_list_new[[2]]
      funEvals = funEvals + 1

      phi_T = f_new + sum(lambda * (abs(w_new)))
      mL = f + t(g) %*% (w_new - w) + t(w_new - w) %*% (w_new - w) / (2 * alpha) + sum(lambda * (abs(w_new)))
    }

    #Extrapolation step
    z_new = w_new
    a_new = (1 + sqrt(1 + 4 * a * a)) / 2
    w_new = z_new + ((a - 1) / a_new) * (z_new - z)
    # Compute Evaluate Function
    fg_list_new = gradFunc(w_new)
    f_new = fg_list_new[[1]]
    g_new = fg_list_new[[2]]
    funEvals = funEvals + 1

    a = a_new
    z = z_new
    w = w_new
    f = f_new
    g = g_new

    if (sum(abs(w - w_old)) < optTol)
      break
  }

  return(list(w, funEvals))
}

DGLoss <- function(w, sigma1, sigma2, b) {
  d1 = ncol(sigma1)
  d2 = ncol(sigma2)
  W = matrix(w, d1, d2)
  tmp = sigma1 %*% W %*% sigma2
  Qw = matrix(tmp, d1 * d2, 1)
  f = 1 / 2 * t(w) %*% Qw -  t(w) %*% b
  g = Qw - b

  return(list(f, g))
}

differential_graph <- function(sigma1, sigma2, lambda) {
  d = ncol(sigma1)
  b = sigma1 - sigma2
  b = as.vector(b)
  nVars = d * d
  w_init = rep(0, nVars)
  funObj = function(w) DGLoss(w, sigma1, sigma2, b);

  params = c()
  params.verbose = 0
  optres = L1GeneralCompositeGradientAccelerated(funObj, w_init, lambda, params)
  w = optres[[1]]

  Theta = matrix(w, d, d)
  Theta = pmax(Theta, t(Theta))
  rownames(Theta) = colnames(Theta) = rownames(sigma1)
  attr(Theta, 'funEvals') = optres[[2]]

  return(Theta)
}

find_lambda_max <- function(sigma1, sigma2) {
  lambda_init = 1
  delta = differential_graph(sigma1, sigma2, lambda_init)
  iter = 0

  if (sum(delta != 0) / 2 == 0) {
    while (sum(delta[upper.tri(delta)] != 0) == 0 & iter < 30) {
      lambda_init = lambda_init / 2
      delta = differential_graph(sigma1, sigma2, lambda_init)
      iter = iter + 1
    }
    lambda_max = lambda_init * 5
  } else {
    while (sum(delta != 0) / 2 != 0 & iter < 30) {
      lambda_init = lambda_init * 2
      delta = differential_graph(sigma1, sigma2, lambda_init)
      iter = iter + 1
    }
    lambda_max = lambda_init
  }

  return(lambda_max)
}

find_lambda_min <- function(sigma1, sigma2, lambda_max, n_target) {
  lambda_min = 1 / 1.2 * lambda_max
  tp = 0
  iter = 0

  while (tp < n_target & iter < 30) {
    lambda_min = 1 / 1.2 * lambda_min
    delta = differential_graph(sigma1, sigma2, lambda_min)
    tp = sum(delta != 0) / 2
    iter = iter + 1
  }

  return(lambda_min)
}

#binary search
search_lambda <- function(sigma1, sigma2, lambda_min, lambda_max, n_target, maxiter = 50) {
  tp = 0
  iter = 0

  a = lambda_min
  b = lambda_max
  fa = sum(differential_graph(sigma1, sigma2, a) != 0) / 2
  fb = sum(differential_graph(sigma1, sigma2, b) != 0) / 2
  res = matrix(NA, maxiter, 2)

  while (tp != n_target & iter < maxiter) {
    c = b - (b - a) / 2
    net_c = differential_graph(sigma1, sigma2, c)
    fc = sum(net_c != 0) / 2
    if (fa >= n_target & n_target >= fc) {
      b = c
      fb = fc
    } else {
      a = c
      fa = fc
    }
    # print(c(c, fc))
    iter = iter + 1
    res[iter, ] = c(c, fc)
  }

  res = res[!is.na(res[, 1]), ]
  res[, 2] = abs(res[, 2] - n_target)
  c = res[which.min(res[, 2]), 1]

  return(c)
}

ldgm.score <- function(emat, condition, lambda = NA, n_target = NA, search.iter = 50, binary = FALSE) {
  expr1 = emat[, condition == 1, drop = FALSE]
  expr2 = emat[, condition == 2, drop = FALSE]

  #compute latent correlation matrices in each condition
  kendall1 = cor(t(expr1), method = 'kendall')
  kendall2 = cor(t(expr2), method = 'kendall')
  lcor1 = sin(pi/2 * kendall1) #latent correlation estimate
  lcor2 = sin(pi/2 * kendall2) #latent correlation estimate
  diag(lcor1) = diag(lcor2) = 1

  if (is.na(lambda) & is.na(n_target)) {
    stop('Need to specify either lambda or num of edges in true dc network')
  }

  if (!is.na(lambda)) {
    #place holder to prioritise lambda choice over n_target
  } else if (!is.na(n_target)) {
    #search for lambda
    lambda_max = find_lambda_max(lcor1, lcor2)
    lambda_min = find_lambda_min(lcor1, lcor2, lambda_max, n_target)
    lambda = search_lambda(lcor1, lcor2, lambda_min, lambda_max, n_target, search.iter)
  }

  delta = differential_graph(lcor1, lcor2, lambda)

  if (binary)
    delta = delta != 0
  diag(delta) = NA

  #add run parameters as attributes
  attributes(delta) = c(attributes(delta),
                    'lambda' = lambda,
                    'call' = match.call())

  return(delta)
}
