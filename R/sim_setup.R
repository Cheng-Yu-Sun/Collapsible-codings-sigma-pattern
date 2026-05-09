#####################################################################
# Common setup and simulation routines for Section 5.4 of the paper.
#####################################################################

# If the script is run from a subfolder, move the working directory
# back to the repository root. We allow at most two upward moves.
n_up <- 0

while (
  !(dir.exists("R") && dir.exists("data")) &&
  dirname(getwd()) != getwd() &&
  n_up < 2
) {
  setwd("..")
  n_up <- n_up + 1
}

if (!(dir.exists("R") && dir.exists("data"))) {
  stop(
    "Could not locate the repository root. ",
    "Please set the working directory to the repository root ",
    "or to a subfolder within it."
  )
}

# Install required packages if needed:
# install.packages(c("doSNOW", "foreach", "SLHD", "MaxPro", "LHD", "TestFunctions", "hetGP"))

library(doSNOW)
library(foreach)
library(hetGP)
library(LHD)
library(MaxPro)
library(parallel)
library(SLHD)
library(TestFunctions)

D_Sigma <- readRDS("data/D_Sigma.rds")
D_SP <- readRDS("data/D_SP.rds")

# The fixed designs D_Sigma and D_SP are 32 x 15 four-level arrays.
# They can be reproduced by running "R/construct_D_Sigma.R" and "R/construct_D_SP.R".
# Repeated executions of these scripts may return different arrays,
# but the resulting designs have the same Sigma/SP patterns.

sim.fun <- function(test.f, n_active, is.log = TRUE) {
  N <- 32
  N.test <- 1000
  m <- 15
  
  # L1 and L2: LHDs obtained from D_Sigma and D_SP
  L1 <- (apply(D_Sigma, 2, rank, ties.method = "random") - 0.5) / N
  L2 <- (apply(D_SP, 2, rank, ties.method = "random") - 0.5) / N
  
  # L3-L5: benchmark designs
  L3 <- (rLHD(N, m) - 0.5) / N
  L4 <- maximinSLHD(1, N, m)$StandDesign
  L5 <- MaxProLHD(N, m)$Design
  
  # Randomly permute the columns
  perm <- sample.int(m)
  L1 <- L1[, perm]
  L2 <- L2[, perm]
  L3 <- L3[, perm]
  L4 <- L4[, perm]
  L5 <- L5[, perm]
  
  list.design <- list(L1, L2, L3, L4, L5)
  
  # Test set on (0,1)^15
  d.test <- (rLHD(N.test, m) - 0.5) / N.test
  y.test <- apply(d.test[, 1:n_active], 1, test.f)
  if (is.log) y.test <- log(y.test)
  
  NRMSE <- numeric(length(list.design))
  
  for (l in seq_along(list.design)) {
    trainX <- list.design[[l]]
    trainY <- apply(trainX[, 1:n_active], 1, test.f)
    if (is.log) trainY <- log(trainY)
    
    fit <- mleHomGP(X = trainX, Z = trainY, covtype = "Matern5_2")
    pre <- predict(object = fit, x = d.test)
    
    NRMSE[l] <- sqrt(
      mean((pre$mean - y.test)^2) /
        mean((mean(trainY) - y.test)^2)
    )
  }
  
  NRMSE
}

run_sim <- function(test.f, n_active, is.log = TRUE,
                    n_rep = 300, seed = 12345) {
  cpu.cores <- detectCores()
  cl <- makeCluster(max(1, cpu.cores - 2))
  registerDoSNOW(cl)
  clusterSetRNGStream(cl, seed)
  
  on.exit(stopCluster(cl), add = TRUE)
  
  NRMSE <- foreach(
    i = 1:n_rep,
    .combine = "rbind",
    .packages = c("SLHD", "MaxPro", "LHD", "hetGP", "TestFunctions"),
    .export = c("sim.fun", "D_Sigma", "D_SP")
  ) %dopar% sim.fun(test.f, n_active, is.log)
  
  colnames(NRMSE) <- c("L1", "L2", "L3", "L4", "L5")
  NRMSE
}

plot_nrmse <- function(NRMSE, main_title) {
  boxplot(
    NRMSE,
    names = c("L1(Sigma)", "L2(SP)", "L3", "L4", "L5"),
    main = main_title
  )
}