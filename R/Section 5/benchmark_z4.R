#####################################################################
# Benchmark test functions for z = 4.
# Run All fits the GPs and saves the NRMSE matrices.
# Table3_benchmarks.R combines both design settings.
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

test_functions <- list(
  OTL_Circuit = TestFunctions::OTL_Circuit,
  Borehole = TestFunctions::borehole,
  Steel_Column = TestFunctions::steelcolumnstress,
  Wing_Weight = TestFunctions::wingweight
)
n_active <- c(6, 8, 9, 10)

for (k in seq_along(test_functions)) {
  name <- names(test_functions)[k]
  message("Benchmark function: z = 4, ", name)
  NRMSE <- run_sim(
    N = nrow(D_Sigma), m = ncol(D_Sigma),
    test.f = test_functions[[k]], n_active = n_active[k], is.log = TRUE,
    n_rep = n_rep, seed = seed, N.test = N.test, n_cores = n_cores
  )
  saveRDS(NRMSE, paste0("detailed_outputs/simulations/NRMSE_z4_", name, ".rds"))
  print(colMeans(NRMSE))
}
