############################################################
# Simulation code for the additional test functions
# in Section 5.4 of the paper.
#
# This script reproduces Figure 3(a)--Figure 3(d),
# corresponding to the OTL Circuit, Borehole, Steel Column,
# and Wing Weight functions, respectively.
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

# Test functions used in Figure 3.
Circuit6 <- function(x) TestFunctions::OTL_Circuit(x)
Borehole8 <- function(x) TestFunctions::borehole(x)
Steel9 <- function(x) TestFunctions::steelcolumnstress(x)
Wing10 <- function(x) TestFunctions::wingweight(x)


# -----------------------------
# Figure 3(a): OTL Circuit
# -----------------------------

NRMSE_Figure3a <- run_sim(
  test.f = Circuit6,
  n_active = 6,
  is.log = TRUE,
  seed = 12345
)

saveRDS(
  NRMSE_Figure3a,
  "results/NRMSE_Figure3a_OTL_Circuit.rds"
)

pdf("outputs/figures/Figure3a_OTL_Circuit.pdf", width = 6, height = 4)
plot_nrmse(NRMSE_Figure3a, "Figure 3(a): OTL Circuit")
dev.off()


# -----------------------------
# Figure 3(b): Borehole
# -----------------------------

NRMSE_Figure3b <- run_sim(
  test.f = Borehole8,
  n_active = 8,
  is.log = TRUE,
  seed = 12345
)

saveRDS(
  NRMSE_Figure3b,
  "results/NRMSE_Figure3b_Borehole.rds"
)

pdf("outputs/figures/Figure3b_Borehole.pdf", width = 6, height = 4)
plot_nrmse(NRMSE_Figure3b, "Figure 3(b): Borehole")
dev.off()


# -----------------------------
# Figure 3(c): Steel Column
# -----------------------------

NRMSE_Figure3c <- run_sim(
  test.f = Steel9,
  n_active = 9,
  is.log = TRUE,
  seed = 12345
)

saveRDS(
  NRMSE_Figure3c,
  "results/NRMSE_Figure3c_Steel_Column.rds"
)

pdf("outputs/figures/Figure3c_Steel_Column.pdf", width = 6, height = 4)
plot_nrmse(NRMSE_Figure3c, "Figure 3(c): Steel Column")
dev.off()


# -----------------------------
# Figure 3(d): Wing Weight
# -----------------------------

NRMSE_Figure3d <- run_sim(
  test.f = Wing10,
  n_active = 10,
  is.log = TRUE,
  seed = 12345
)

saveRDS(
  NRMSE_Figure3d,
  "results/NRMSE_Figure3d_Wing_Weight.rds"
)

pdf("outputs/figures/Figure3d_Wing_Weight.pdf", width = 6, height = 4)
plot_nrmse(NRMSE_Figure3d, "Figure 3(d): Wing Weight")
dev.off()


# Display the plots in RStudio for visual checking.

plot_nrmse(NRMSE_Figure3a, "Figure 3(a): OTL Circuit")
plot_nrmse(NRMSE_Figure3b, "Figure 3(b): Borehole")
plot_nrmse(NRMSE_Figure3c, "Figure 3(c): Steel Column")
plot_nrmse(NRMSE_Figure3d, "Figure 3(d): Wing Weight")