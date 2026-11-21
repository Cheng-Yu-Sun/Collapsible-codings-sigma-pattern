#####################################################################
# Self-defined test function for z = 4: lambda = 3 and 30.
# Run All fits the GPs and saves the NRMSE matrices.
# Figure2.R combines the z = 4 and z = 6 results.
#####################################################################

# Set the working directory to the repository root or a subfolder.
n_up <- 0
while (!(dir.exists("R") && dir.exists("design_matrices")) &&
       dirname(getwd()) != getwd() && n_up < 3) {
  setwd("..")
  n_up <- n_up + 1
}
if (!(dir.exists("R") && dir.exists("design_matrices"))) {
  stop("Open GitHub_Sigma.Rproj, or set the working directory to the repository root.")
}

source("R/Section 5/sim_setup.R")

n_rep <- 300
seed <- 12345
N.test <- 10^4
n_cores <- max(1, detectCores() - 2)
dir.create("detailed_outputs/simulations", recursive = TRUE, showWarnings = FALSE)

# Fixed inputs; see README for their search and generation scripts.
D_Sigma <- readRDS("design_matrices/SOA_z4_Sigma.rds")
D_SP <- readRDS("design_matrices/SOA_z4_SP.rds")
Uniform_LHD <- readRDS("design_matrices/class_iii/Uniform_LHD_z4.rds")
Random_LHD <- readRDS("design_matrices/class_iv/Random_LHD_z4.rds")
Maximin_LHD <- readRDS("design_matrices/class_v/Maximin_LHD_z4.rds")
MaxPro_LHD <- readRDS("design_matrices/class_vi/MaxPro_LHD_z4.rds")

make_f_test <- function(lambda) {
  force(lambda)
  function(x) lambda * sin(pi * x[1] * x[2]) + sum(x[3:10])
}

for (lambda in c(3, 30)) {
  message("Self-defined function: z = 4, lambda = ", lambda)
  NRMSE <- run_sim(
    N = nrow(D_Sigma), m = ncol(D_Sigma),
    test.f = make_f_test(lambda), n_active = 10, is.log = FALSE,
    n_rep = n_rep, seed = seed, N.test = N.test, n_cores = n_cores
  )
  saveRDS(NRMSE, paste0("detailed_outputs/simulations/NRMSE_z4_lambda",
                        lambda, ".rds"))
  print(colMeans(NRMSE))
}
