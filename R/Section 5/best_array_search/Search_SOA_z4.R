#####################################################################
# Search for the z=4 SOAs underlying classes (i) and (ii).
# D = 2A + B, where A ranges over the 144 regular 32-run designs
# and B ranges over all binary arrays satisfying the SOA constraints.
#
# Run this script to perform both searches. Common functions are in
# soa_z4_functions.R. The selected arrays are saved in design_matrices/ and
# the search scores are saved in detailed_outputs/best_array_search/.
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
# install.packages(c("lpSolve", "doSNOW", "foreach"))

library(lpSolve)
library(doSNOW)
library(foreach)
library(parallel)

functions_file <- normalizePath(
  "R/Section 5/best_array_search/soa_z4_functions.R"
)
source(functions_file)

n_starts <- 100 # Number of random initial B arrays for each A.
seed <- 41001
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
n_cores <- max(1, detectCores() - 2)

catalog <- read_catalog("design_matrices/regular_32run_catalog.txt")
X <- binary_columns()

# nu33 depends only on A.
nu33 <- vapply(seq_len(nrow(catalog)), function(i) {
  coarse_score(X[, catalog[i, ]])
}, numeric(1))


#####################################################################
# SP: minimize nu33 over A, then search B to minimize nu42.
#####################################################################

sp_ids <- which(nu33 == min(nu33))
sp <- vector("list", length(sp_ids))

sp_time <- system.time({
  for (k in seq_along(sp_ids)) {
    i <- sp_ids[k]
    message("SP: searching A ", i)
    A <- X[, catalog[i, ]]
    sp[[k]] <- multistart_B(A, n_starts, seed + i)
  }
})[["elapsed"]]

sp_scores <- vapply(sp, function(x) x$score, numeric(1))
sp_position <- order(sp_scores, sp_ids)[1]
sp_id <- sp_ids[sp_position]
sp_best <- sp[[sp_position]]


#####################################################################
# Sigma: search B for every A in parallel.
# Then minimize nu42 followed by nu33.
#####################################################################

cl <- makeCluster(n_cores)
registerDoSNOW(cl)

# Load the common functions and set the RNG method on each worker.
clusterExport(cl, "functions_file")
invisible(clusterEvalQ(cl, {
  source(functions_file)
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
}))

# Show progress whenever the search for an A is completed.
progress <- function(n) {
  message("Sigma: completed ", n, "/", nrow(catalog), " choices of A")
}

sigma_time <- system.time({
  sigma <- foreach(
    i = seq_len(nrow(catalog)),
    .packages = "lpSolve",
    .options.snow = list(progress = progress)
  ) %dopar% {
    A <- X[, catalog[i, ]]
    # The same A uses the same random starts in both searches.
    multistart_B(A, n_starts, seed + i)
  }
})[["elapsed"]]

stopCluster(cl)
registerDoSEQ()

nu42 <- vapply(sigma, function(x) x$score, numeric(1))
sigma_id <- order(nu42, nu33, seq_along(nu33))[1]
sigma_best <- sigma[[sigma_id]]


#####################################################################
# Save the selected arrays and search scores.
#####################################################################

D_Sigma <- sigma_best$design
D_SP <- sp_best$design
colnames(D_Sigma) <- colnames(D_SP) <- paste0("x", seq_len(ncol(D_Sigma)))

# These replace the saved z=4 SOAs used by the GP simulations.
saveRDS(D_Sigma, "design_matrices/SOA_z4_Sigma.rds")
saveRDS(D_SP, "design_matrices/SOA_z4_SP.rds")

# For this family, nu33 = 16*sigma33 and nu42 = 4*sigma42.
selected_scores <- data.frame(
  criterion = c("Sigma", "SP"),
  candidate = c(sigma_id, sp_id),
  nu33 = nu33[c(sigma_id, sp_id)],
  nu42 = c(sigma_best$score, sp_best$score),
  sigma33 = nu33[c(sigma_id, sp_id)] / 16,
  sigma42 = c(sigma_best$score, sp_best$score) / 4,
  starts_per_A = n_starts,
  seconds = c(sigma_time, sp_time)
)

# One row per A in the Sigma search.
candidate_scores <- data.frame(
  candidate = seq_along(nu33),
  nu33 = nu33,
  best_nu42 = nu42,
  sigma33 = nu33 / 16,
  best_sigma42 = nu42 / 4,
  starts = n_starts
)

dir.create(
  "detailed_outputs/best_array_search",
  recursive = TRUE,
  showWarnings = FALSE
)

write.csv(
  selected_scores,
  "detailed_outputs/best_array_search/SOA_z4_selected_scores.csv",
  row.names = FALSE
)

write.csv(
  candidate_scores,
  "detailed_outputs/best_array_search/SOA_z4_candidate_scores.csv",
  row.names = FALSE
)

print(selected_scores, row.names = FALSE)
