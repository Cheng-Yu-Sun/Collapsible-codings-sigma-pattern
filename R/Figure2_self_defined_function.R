############################################################
# Simulation code for the self-defined test function
# in Section 5.4 of the paper.
#
# This script reproduces Figure 2(a) and Figure 2(b),
# corresponding to lambda = 30 and lambda = 3, respectively.
############################################################

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

source("R/sim_setup.R")

if (!dir.exists("results")) {
  dir.create("results")
}

if (!dir.exists("outputs")) {
  dir.create("outputs")
}

if (!dir.exists("outputs/figures")) {
  dir.create("outputs/figures")
}

make_f_test <- function(lambda) {
  force(lambda)
  
  function(x) {
    lambda * sin(pi * x[1] * x[2]) + sum(x[3:length(x)])
  }
}

# -----------------------------
# Figure 2(a): lambda = 30
# -----------------------------

NRMSE_Figure2a <- run_sim(
  test.f = make_f_test(30),
  n_active = 10,
  is.log = FALSE,
  seed = 12345
)

saveRDS(
  NRMSE_Figure2a,
  "results/NRMSE_Figure2a_self_defined_lambda30.rds"
)

pdf("outputs/figures/Figure2a_self_defined_lambda30.pdf", width = 6, height = 4)
plot_nrmse(NRMSE_Figure2a, "Figure 2(a): lambda = 30")
dev.off()


# -----------------------------
# Figure 2(b): lambda = 3
# -----------------------------

NRMSE_Figure2b <- run_sim(
  test.f = make_f_test(3),
  n_active = 10,
  is.log = FALSE,
  seed = 12345
)

saveRDS(
  NRMSE_Figure2b,
  "results/NRMSE_Figure2b_self_defined_lambda3.rds"
)

pdf("outputs/figures/Figure2b_self_defined_lambda3.pdf", width = 6, height = 4)
plot_nrmse(NRMSE_Figure2b, "Figure 2(b): lambda = 3")
dev.off()


# Display the plots in RStudio for visual checking.

plot_nrmse(NRMSE_Figure2a, "Figure 2(a): lambda = 30")
plot_nrmse(NRMSE_Figure2b, "Figure 2(b): lambda = 3")