#####################################################################
# Class (iii): Uniform LHDs for z = 4 and z = 6.
# Run All replaces the two saved lists in design_matrices/class_iii/.
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

# Install these packages if needed:
# install.packages(c("doSNOW", "foreach", "SFDesign"))
library(doSNOW)
library(foreach)
library(parallel)

n_rep <- 300
seed <- 12345
n_cores <- max(1, detectCores() - 2)

cl <- makeCluster(n_cores)
registerDoSNOW(cl)
clusterSetRNGStream(cl, seed)

for (z in c(4, 6)) {
  N <- if (z == 4) 32 else 72
  m <- if (z == 4) 15 else 13
  message("Generating uniform LHDs for z = ", z)

  Uniform_LHD <- foreach(i = seq_len(n_rep)) %dopar% {
    # A seed per replicate makes the draws independent of worker scheduling.
    set.seed(seed + i)
    SFDesign::uniformLHD(n = N, p = m)$design
  }

  saveRDS(Uniform_LHD,
          paste0("design_matrices/class_iii/Uniform_LHD_z", z, ".rds"))
}

stopCluster(cl)
registerDoSEQ()
