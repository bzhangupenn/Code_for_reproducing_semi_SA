###############################################################
###############################################################
# Some preliminaries
library(rootSolve)
library(matlib)
library(numDeriv)

expit <- function(x){
  return(1/(1+exp(-x)))
}

# Specify (possibly misspecify) a distribution on U
# discretize U \in [0, 1]
u_vec = seq(0, 1, 0.1)

# Approximate the integral I_{ij} given $beta, kappa, X_i$

p_z <- function(z, x, u, kappa, lambda){
  if (z == 1)
    return(expit(sum(kappa*x) + lambda*u))
  else
    return(1 - expit(sum(kappa*x) + lambda*u))
}

# p(Y | X, Z, U) under a logistic regression model without intercept
p_y <- function(y, x, z, u, beta, delta){
  if (y == 1)
    return(expit(sum(beta*c(x,z)) + delta*u))
  else
    return(1 - expit(sum(beta*c(x,z)) + delta*u))
}


###############################################################################
###############################################################################
# Approximate the integrals using quadrature

# Z and Y are both binary. We only need to sum 4 terms.
compute_I_ij<- function(beta, kappa, lambda, delta, u_i, u_j, u_j_prob, X_i, u_vec, u_prob, func){
  sum = 0
  for (y in 0:1){
    for (z in 0:1){
      integrand = func(y, beta, kappa, lambda, delta, u_j, u_j_prob, X_i, z, u_vec, u_prob)
      sum = sum + integrand*p_y(y, X_i, z, u_i, beta, delta)*p_z(z, X_i, u_i, kappa, lambda)
      #cat(sum,'\n')
    }
  }
  return(sum)
}

# Integrand of the integral I_{ij}
func <- function(y, beta, kappa, lambda, delta, u_j, u_j_prob, X_i, z, u_vec, u_prob){
  num = p_z(z, X_i, u_j, kappa, lambda)*p_y(y, X_i, z, u_j, beta, delta)*u_j_prob
  den = 0
  for (i in 1:length(u_vec)){
    den = den + p_z(z, X_i, u_vec[i], kappa, lambda)*
          p_y(y, X_i, z, u_vec[i], beta, delta)*u_prob[i]
  }
  return(num/den)
}

compute_b_i<- function(beta, kappa, lambda, delta, u_i, X_i, u_vec, u_prob, func_2){
  sum = 0
  for (y in 0:1){
    for (z in 0:1){
      integrand = func_2(y, beta, kappa, lambda, delta, u_i, X_i, z, u_vec, u_prob)
      sum = sum + integrand*p_y(y, X_i, z, u_i, beta, delta)*p_z(z, X_i, u_i, kappa, lambda)
    }
  }
  return(sum)
}

# Integrand of the integral b_i
func_2 <- function(y, beta, kappa, lambda, delta, u_i, X_i, z, u_vec, u_prob){
  num = 0
  for (i in 1:length(u_vec)){
    num = num + score(beta, kappa, lambda, delta, X_i, u_vec[i], z, y)*
      p_z(z, X_i, u_vec[i], kappa, lambda)*p_y(y, X_i, z, u_vec[i], beta, delta)*u_prob[i]
  }
  den = 0
  for (i in 1:length(u_vec)){
    den = den + p_z(z, X_i, u_vec[i], kappa, lambda)*p_y(y, X_i, z, u_vec[i], beta, delta)*
          u_prob[i]
  }
  return(num/den)
}

# Compute the score w.r.t beta and kappa
score <- function(beta, kappa, lambda, delta, x, u, z, y){
  score_beta = (y == 1)*c(x,z) - expit(sum(beta*c(x,z)) + delta*u) * c(x,z)
  score_kappa = (z == 1)*x - expit(sum(kappa*x) + lambda*u) * x
  return(c(score_beta, score_kappa))
}


##########################################################################################
# Compute the projection

# Solve for a(X_i): I^(-1)*b, where b is m by q, with b_i the ith row
# return a m by q (10 by 6) matrix
obtain_I <- function(beta, kappa, lambda, delta, X_i, u_vec, u_prob){
  I = matrix(0, nrow = length(u_vec), ncol = length(u_vec))
  for (i in 1:length(u_vec)){
    for (j in 1:length(u_vec)){
      I[i,j] = compute_I_ij(beta, kappa, lambda, delta, u_vec[i], u_vec[j], u_prob[j], X_i, u_vec, u_prob, func)
    }
  }
  return(I)
}

obtain_b <- function(beta, kappa, lambda, delta, X_i, u_vec, u_prob){
  b = NULL
  for (i in 1:length(u_vec)){
    b_i = compute_b_i(beta, kappa, lambda, delta, u_vec[i], X_i, u_vec, u_prob, func_2)
    b = rbind(b, b_i)
  }
  return(b)
}

solve_a <- function(beta, kappa, lambda, delta, X_i, u_vec, u_prob){
  I = obtain_I(beta, kappa, lambda, delta, X_i, u_vec, u_prob)
  b = obtain_b(beta, kappa, lambda, delta, X_i, u_vec, u_prob)
  #print(I)
  #print(b)
  #cat('X_i', X_i, '\n')
  #alpha = 0.1
  alpha = 0.01
  k = length(u_vec)
  W = diag(c(0.5, rep(1,k-2),0.5), nrow = k)
  h = 1/(k-1)
  #design_matrix = h*I %*% W
  design_matrix = I
  return(inv(t(design_matrix) %*% design_matrix + diag(alpha, nrow = k)) %*% t(design_matrix) %*% b)
}

efficient_score <- function(beta, kappa, lambda, delta, Y, Z, X_i, u_vec, u_prob){
  a = solve_a(beta, kappa, lambda, delta, X_i, u_vec, u_prob)
  num = 0
  for (i in 1:length(u_vec)){
    sc = score(beta, kappa, lambda, delta, X_i, u_vec[i], Z, Y)
    a_i = a[i,]
    #cat('i', i, u_vec[i], '\n')
    #cat('score', sc, '\n')
    #cat('a_i', a_i, '\n')
    temp_num = (sc - a_i)*p_z(Z, X_i, u_vec[i], kappa, lambda)*p_y(Y, X_i, Z, u_vec[i], beta, delta)*u_prob[i]
    #cat('temp_num', temp_num, '\n')
    num = num + temp_num
  }
  #cat(num,'\n')
  den = 0
  for (i in 1:length(u_vec)){
    den = den + p_z(Z, X_i, u_vec[i], kappa, lambda)*p_y(Y, X_i, Z, u_vec[i], beta, delta)*u_prob[i]
  }
  return(num/den)
}

# For solving Jacobian
efficient_score_jac <- function(param, q, lambda, delta, Y, Z, X_i, u_vec, u_prob){
  beta = param[1 : (q+1)]
  kappa = param[(q+2):(2*q+1)]
  a = solve_a(beta, kappa, lambda, delta, X_i, u_vec, u_prob)
  num = 0
  for (i in 1:length(u_vec)){
    sc = score(beta, kappa, lambda, delta, X_i, u_vec[i], Z, Y)
    a_i = a[i,]
    temp_num = (sc - a_i)*p_z(Z, X_i, u_vec[i], kappa, lambda)*p_y(Y, X_i, Z, u_vec[i], beta, delta)*u_prob[i]
    num = num + temp_num
  }
  #cat(num,'\n')
  den = 0
  for (i in 1:length(u_vec)){
    den = den + p_z(Z, X_i, u_vec[i], kappa, lambda)*p_y(Y, X_i, Z, u_vec[i], beta, delta)*u_prob[i]
  }
  return(num/den)
}

# param is c(beta, kappa)
# 4-dim in this toy example
estimating_eq <- function(param, lambda = 1, delta = 1, u_vec = seq(0,1,1), u_prob, dataset = data1){
  num_row = dim(dataset)[1]
  EE = 0
  for (i in 1:num_row){
    #cat(i, '\n')
    beta = param[1:3]
    kappa = param[4:5]
    Y = dataset[i, 4]
    Z = dataset[i, 3]
    X_i = as.matrix(dataset[i, 1:2])
    #cat(beta, kappa, Y, Z, X_i, '\n')
    e_s = efficient_score(beta, kappa, lambda, delta, Y, Z, X_i, u_vec, u_prob)
    #cat(e_s, '\n')
    EE = EE + e_s
  }
  return(EE)
}


# apply to real data
# q is number of X, followed by treatment Z, and finally by response Y
# Throughout, X includes the intercept and (beta, kappa) include the intercept
estimating_eq_2 <- function(param, q, lambda = 1, delta = 1, u_vec = seq(0,1,1), u_prob, dataset = dt_real){
  num_row = dim(dataset)[1]
  EE = 0
  for (i in 1:num_row){
    #cat(i, '\n')
    beta = param[1:(q+1)]
    kappa = param[(q+2):(2*q+1)]
    Y = dataset[i, q+2]
    Z = dataset[i, q+1]
    X_i = as.matrix(dataset[i, 1:q])
    #cat(beta, kappa, Y, Z, X_i, '\n')
    e_s = efficient_score(beta, kappa, lambda, delta, Y, Z, X_i, u_vec, u_prob)
    #cat(e_s, '\n')
    EE = EE + e_s
  }
  return(EE)
}

# Construct confidence band
# Assume lambda = delta
# lambda and delta mean supremum
mul_boot_once <- function(beta, kappa, q, lambda_sup, delta_sup, u_vec, u_prob, dataset = dt_real){
  num_row = dim(dataset)[1]
  eps = rnorm(num_row, 0, 1)
  Z_vec = NULL
  lambda_vec = seq(0, lambda_sup, 0.2)
  for (lambda in lambda_vec){
    param_0 = multiroot(estimating_eq_2, c(beta,kappa), q = q, lambda = lambda, delta = delta,
                        u_vec = u_vec, u_prob = u_prob, dataset = dataset)$root
    EE = 0
    for (i in 1:num_row){
      #cat(i, '\n')
      beta = param_0[1:(q+1)]
      kappa = param_0[(q+2):(2*q+1)]
      Y = dataset[i, q+2]
      Z = dataset[i, q+1]
      X_i = as.matrix(dataset[i, 1:q])
      #cat(beta, kappa, Y, Z, X_i, '\n')
      e_s = efficient_score(beta, kappa, lambda, delta, Y, Z, X_i, u_vec, u_prob)
      #cat(e_s, '\n')
      EE = EE + e_s[q+1]*eps[i]
    }
    est = param_0[q+1]
    var_beta = sandwich_variance(param_0 = param_0, q = q, lambda = lambda, delta = delta,
                                 u_vec = u_vec, u_prob = u_prob, dataset = dataset)
    sd_beta = sqrt(var_beta[q+1, q+1])
    Z_temp = abs(EE/sd_beta)/sqrt(num_row)
    cat(lambda, Z_temp, '\n')
    Z_vec = c(Z_vec, Z_temp)
  }
  return(max(Z_vec))
}

# Compute the robust sandwich estimator of the variance
sandwich_variance <- function(param_0, q, lambda, delta, u_vec, u_prob, dataset){
  beta_0 = param_0[1:(q+1)]
  kappa_0 = param_0[(q+2):(2*q+1)]
  A_m = 0
  B_m = 0
  m = dim(dataset)[1]
  for (i in 1:m){
    X_i = as.matrix(dataset[i,1:q])
    Z = dataset[i, (q+1)]
    Y = dataset[i, (q+2)]
    phi_i = efficient_score(beta_0, kappa_0, lambda, delta, Y, Z, X_i, u_vec, u_prob)
    B_i = phi_i %*% t(phi_i)
    A_i = -jacobian(efficient_score_jac, c(beta_0, kappa_0), q = q, lambda = lambda,
                   delta = delta, Y = Y, Z = Z, X_i = X_i, u_vec = u_vec, u_prob = u_prob)
    A_m = A_m + A_i
    B_m = B_m + B_i
  }
  A_m = A_m/m
  B_m = B_m/m
  var_cov_mat = (solve(A_m)%*%B_m%*%(t(solve(A_m))))/m
  return(var_cov_mat)
}


# Newton's method for solving the EE
newton <- function(estimating_eq, start_value, lambda, delta, u_vec, dataset){
  x_n = start_value
  flag = TRUE
  while (flag){
    cat(x_n, '\n')
    cat('Computing Jacobian', '\n')
    jac = jacobian(estimating_eq, x_n, method = 'simple', lambda = lambda, delta = delta, u_vec = u_vec, dataset = dataset)
    #jac = jac(estimating_eq, x_n, eps = 1e-4, lambda = lambda, delta = delta, u_vec = u_vec, dataset = dataset)
    cat('Jacobian computed', '\n')
    jac_inv = solve(jac)
    print(jac_inv)
    f_x_n = estimating_eq(x_n, lambda = lambda, delta = delta, u_vec = u_vec, dataset = dataset)
    cat('EE', f_x_n, '\n')
    delta_x = jac_inv %*% f_x_n
    cat('delta x', delta_x, '\n')
    x_new = x_n - delta_x
    cat('X_n updated', '\n')
    if (sum((x_new - x_n)^2) <= 1e-20) flag = FALSE
    x_n = x_new
  }
  return(x_n)
}

# Compute Jacobian numerically
jac <- function(f, x, eps, lambda, delta, u_vec, dataset){
  res = NULL
  f_0 = f(x, lambda = lambda, delta = delta, u_vec = u_vec, dataset = dataset)
  for (i in 1:length(x)){
    x_new = x
    x_new[i] = x_new[i] + eps
    column_i = (f(x_new, lambda = lambda, delta = delta, u_vec = u_vec, dataset = dataset) - f_0)/eps
    res = cbind(res, column_i)
  }
  return(res)
}


# Try solving the EE with different starting point
# Try at most 10/20 times
find_zero <- function(estimating_eq_2, start_point, q, lambda, delta, u_vec, u_prob, dataset, i_max = 20){
  i = 0
  solution_ready = FALSE

  while (i <= i_max & ! solution_ready){
    cat(i, '\n')
    solution = tryCatch(multiroot(estimating_eq_2, start_point, q = q, lambda = lambda, delta = delta, u_vec = u_vec, u_prob = u_prob, dataset = dataset),
                        error = function(e){return('It fails.')}, warning = function(e){return('diag zero')})
    if (is.list(solution)){
      if (!is.nan(solution$f.root[1])) solution_ready = TRUE
      else {
        start_point = start_point + rnorm(length(start_point), 0, 1.5)
        i = i + 1
        solution_ready = FALSE
      }
    }
    else{
      start_point = start_point + rnorm(length(start_point), 0, 2)
      i = i + 1
      solution_ready = FALSE
    }
  }
  if (solution_ready) return(c(solution$root, solution$f.root))
  else return('It fails')
}

