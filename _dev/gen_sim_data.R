gen_sim_data <- function(n = 1000) {
  #Generating data similar to Austin (2009) for demonstrating treatment effect estimation
  gen_X <- function(n) {
    X <- matrix(rnorm(9 * n), nrow = n, ncol = 9)
    X[,5] <- as.numeric(X[,5] < .5)
    X
  }

  #~20% treated
  gen_A <- function(X) {
    LP_A <- - 1.2 + log(2)*X[,1] - log(1.5)*X[,2] + log(2)*X[,4] - log(2.4)*X[,5] + log(2)*X[,7] - log(1.5)*X[,8]
    P_A <- plogis(LP_A)
    rbinom(nrow(X), 1, P_A)
  }

  # Continuous outcome
  gen_Y_C <- function(A, X) {
    2*A + 2*X[,1] + 2*X[,2] + 2*X[,3] + 1*X[,4] + 2*X[,5] + 1*X[,6] + rnorm(length(A), 0, 5)
  }
  #Conditional:
  #  MD: 2
  #Marginal:
  #  MD: 2

  # Binary outcome
  gen_Y_B <- function(A, X) {
    LP_B <- -2 + log(2.4)*A + log(2)*X[,1] + log(2)*X[,2] + log(2)*X[,3] + log(1.5)*X[,4] + log(2.4)*X[,5] + log(1.5)*X[,6]
    P_B <- plogis(LP_B)
    rbinom(length(A), 1, P_B)
  }
  #Conditional:
  #  OR:   2.4
  #  logOR: .875
  #Marginal:
  #  RD:    .144
  #  RR:   1.54
  #  logRR: .433
  #  OR:   1.92
  #  logOR  .655

  # Survival outcome
  gen_Y_S <- function(A, X) {
    LP_S <- -2 + log(2.4)*A + log(2)*X[,1] + log(2)*X[,2] + log(2)*X[,3] + log(1.5)*X[,4] + log(2.4)*X[,5] + log(1.5)*X[,6]
    sqrt(-log(runif(length(A)))*2e4*exp(-LP_S))
  }
  #Conditional:
  #  HR:   2.4
  #  logHR: .875
  #Marginal:
  #  HR:   1.57
  #  logHR: .452

  X <- gen_X(n)
  A <- gen_A(X)

  Y_C <- gen_Y_C(A, X)
  Y_B <- gen_Y_B(A, X)
  Y_S <- gen_Y_S(A, X)

  d <- data.frame(A, X, Y_C, Y_B, Y_S)
  simdata <<- d
  simform <<- A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9

  invisible(d)
}