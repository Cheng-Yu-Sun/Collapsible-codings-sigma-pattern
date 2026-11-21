#####################################################################
# Common setup and simulation routines for Section 5.4.
# The four simulation scripts source this file automatically.
#####################################################################


# Install these packages if needed:
# install.packages(c("doSNOW", "foreach", "hetGP", "LHD", "TestFunctions"))
library(doSNOW)
library(foreach)
library(parallel)
library(hetGP)
library(LHD)
library(TestFunctions)

# Run one simulation replicate.
sim.fun <- function(i, N, m, test.f, n_active, is.log = TRUE,
                    N.test = 10^4) {
  # Classes (i) and (ii): randomly expand each SOA to an LHD.
  L1 <- (apply(D_Sigma, 2, rank, ties.method = "random") - 0.5) / N
  L2 <- (apply(D_SP, 2, rank, ties.method = "random") - 0.5) / N

  # Classes (iii)-(vi): use the ith saved benchmark design.
  L3 <- Uniform_LHD[[i]]
  L4 <- Random_LHD[[i]]
  L5 <- Maximin_LHD[[i]]
  L6 <- MaxPro_LHD[[i]]

  # Apply the same column permutation to all six designs.
  perm <- sample.int(m)
  L1 <- L1[, perm]
  L2 <- L2[, perm]
  L3 <- L3[, perm]
  L4 <- L4[, perm]
  L5 <- L5[, perm]
  L6 <- L6[, perm]
  list.design <- list(L1, L2, L3, L4, L5, L6)

  # Use a common test LHD on (0,1)^m within this replicate.
  d.test <- (rLHD(N.test, m) - 0.5) / N.test
  y.test <- apply(d.test[, seq_len(n_active), drop = FALSE], 1, test.f)
  if (is.log) y.test <- log(y.test)

  NRMSE <- numeric(length(list.design))
  for (l in seq_along(list.design)) {
    trainX <- list.design[[l]]
    trainY <- apply(trainX[, seq_len(n_active), drop = FALSE], 1, test.f)
    if (is.log) trainY <- log(trainY)

    fit <- mleHomGP(X = trainX, Z = trainY, covtype = "Matern5_2")
    pre <- predict(object = fit, x = d.test)
    NRMSE[l] <- sqrt(mean((pre$mean - y.test)^2) /
                      mean((y.test - mean(y.test))^2))
  }
  NRMSE
}

# Run n_rep independent simulation replicates in parallel.
run_sim <- function(N, m, test.f, n_active, is.log = TRUE,
                    n_rep = 300, seed = 12345, N.test = 10^4,
                    n_cores = max(1, detectCores() - 2)) {
  cl <- makeCluster(n_cores)
  registerDoSNOW(cl)
  clusterSetRNGStream(cl, seed)
  on.exit({
    stopCluster(cl)
    registerDoSEQ()
  }, add = TRUE)

  NRMSE <- foreach(
    i = seq_len(n_rep),
    .combine = "rbind",
    .init = matrix(numeric(0), ncol = 6),
    .packages = c("LHD", "hetGP", "TestFunctions"),
    .export = c("sim.fun", "D_Sigma", "D_SP", "Uniform_LHD", "Random_LHD",
                "Maximin_LHD", "MaxPro_LHD")
  ) %dopar% {
    # The same replicate has the same seed on any worker.
    set.seed(seed + i)
    sim.fun(i, N, m, test.f, n_active, is.log, N.test)
  }

  colnames(NRMSE) <- c("Sigma", "SP", "UD", "Random", "Maximin", "MaxPro")
  NRMSE
}
