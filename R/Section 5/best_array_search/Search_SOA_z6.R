#####################################################################
# Search for the z=6 SOAs underlying classes (i) and (ii).
# D = 2A + B, where A = rbind(A_star, A_star[, p]) and A_star is
# a fixed OA(36,13,3,2). B satisfies the SOA constraints.
#
# Run this script to perform both searches. Common functions are in
# soa_z6_functions.R. The selected arrays are saved in design_matrices/ and
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
  "R/Section 5/best_array_search/soa_z6_functions.R"
)
source(functions_file)

n_sp <- 1000   # Number of random permutations screened for SP.
n_sigma <- 100 # Number of random permutations for Sigma, plus identity.
n_starts <- 50 # Number of random initial B arrays for each selected P.
sp_seed <- 62002
sigma_seed <- 62001
binary_seed <- 63001
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
n_cores <- max(1, detectCores() - 2)

# A_star is taken from N. J. A. Sloane's library of orthogonal arrays:
# https://neilsloane.com/oadir/oa.36.13.3.2.txt
# The CSV preserves the original row and column order and levels 0, 1, 2.
A_star <- as.matrix(read.csv("design_matrices/OA_36_13_3_2.csv"))
m <- ncol(A_star)
N <- 2 * nrow(A_star)


#####################################################################
# SP: screen random P by nu33, then search B for every minimizer.
#####################################################################

set.seed(sp_seed)
sp_permutations <- replicate(n_sp, sample.int(m), simplify = FALSE)

# nu33 depends only on A, hence only on P for this fixed A_star.
sp_nu33 <- vapply(sp_permutations, function(p) {
  coarse_score(rbind(A_star, A_star[, p]))
}, numeric(1))

sp_ids <- which(abs(sp_nu33 - min(sp_nu33)) < 1e-10)
sp <- vector("list", length(sp_ids))
message("SP: retained ", length(sp_ids), " of ", n_sp, " permutations")

sp_time <- system.time({
  for (k in seq_along(sp_ids)) {
    i <- sp_ids[k]
    message("SP: searching P ", i)
    A <- rbind(A_star, A_star[, sp_permutations[[i]]])
    sp[[k]] <- multistart_B(A, n_starts, binary_seed + i)
  }
})[["elapsed"]]

sp_scores <- vapply(sp, function(x) x$score, numeric(1))
sp_position <- order(round(sp_scores, 10), sp_ids)[1]
sp_id <- sp_ids[sp_position]
sp_best <- sp[[sp_position]]


#####################################################################
# Sigma: identity plus random P, with the B searches in parallel.
# Then minimize nu42 followed by nu33.
#####################################################################

set.seed(sigma_seed)
sigma_permutations <- c(
  list(seq_len(m)),
  replicate(n_sigma, sample.int(m), simplify = FALSE)
)

sigma_nu33 <- vapply(sigma_permutations, function(p) {
  coarse_score(rbind(A_star, A_star[, p]))
}, numeric(1))

cl <- makeCluster(n_cores)
registerDoSNOW(cl)

# Load the common functions and set the RNG method on each worker.
clusterExport(cl, "functions_file")
invisible(clusterEvalQ(cl, {
  source(functions_file)
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
}))

# Show progress whenever the search for a P is completed.
progress <- function(n) {
  message("Sigma: completed ", n, "/", length(sigma_permutations),
          " permutations")
}

sigma_time <- system.time({
  sigma <- foreach(
    i = seq_along(sigma_permutations),
    .packages = "lpSolve",
    .options.snow = list(progress = progress)
  ) %dopar% {
    A <- rbind(A_star, A_star[, sigma_permutations[[i]]])
    # A separate seed for each P makes the starts independent of scheduling.
    multistart_B(A, n_starts, binary_seed + n_sp + i)
  }
})[["elapsed"]]

stopCluster(cl)
registerDoSEQ()

sigma_nu42 <- vapply(sigma, function(x) x$score, numeric(1))
sigma_id <- order(round(sigma_nu42, 10), round(sigma_nu33, 10),
                  seq_along(sigma_nu33))[1]
sigma_best <- sigma[[sigma_id]]


#####################################################################
# Save the selected arrays and search scores.
#####################################################################

D_Sigma <- sigma_best$design
D_SP <- sp_best$design
colnames(D_Sigma) <- colnames(D_SP) <- paste0("x", seq_len(m))

# These replace the saved z=6 SOAs used by the GP simulations.
saveRDS(D_Sigma, "design_matrices/SOA_z6_Sigma.rds")
saveRDS(D_SP, "design_matrices/SOA_z6_SP.rds")

# nu33 = (72/27)^2*sigma33 and nu42 = (72/36)^2*sigma42.
# The elapsed times below cover the B searches.
selected_scores <- data.frame(
  criterion = c("Sigma", "SP"),
  permutation = c(sigma_id, sp_id),
  nu33 = c(sigma_nu33[sigma_id], sp_nu33[sp_id]),
  nu42 = c(sigma_best$score, sp_best$score),
  sigma33 = c(sigma_nu33[sigma_id], sp_nu33[sp_id]) * (27 / N)^2,
  sigma42 = c(sigma_best$score, sp_best$score) * (36 / N)^2,
  starts_per_P = n_starts,
  seconds = c(sigma_time, sp_time)
)

# One row per sampled P, including identity for Sigma (permutation 1).
# NA means the SP permutation was screened out before the B search.
sp_nu42 <- rep(NA_real_, length(sp_permutations))
sp_nu42[sp_ids] <- sp_scores

candidate_scores <- data.frame(
  criterion = c(rep("Sigma", length(sigma_permutations)),
                rep("SP", length(sp_permutations))),
  permutation = c(seq_along(sigma_permutations), seq_along(sp_permutations)),
  nu33 = c(sigma_nu33, sp_nu33),
  best_nu42 = c(sigma_nu42, sp_nu42),
  sigma33 = c(sigma_nu33, sp_nu33) * (27 / N)^2,
  best_sigma42 = c(sigma_nu42, sp_nu42) * (36 / N)^2,
  starts = c(rep(n_starts, length(sigma_permutations)),
             ifelse(seq_along(sp_permutations) %in% sp_ids, n_starts, 0))
)

# p1,...,p13 record the column order of A_star in the second block.
permutation_columns <- do.call(rbind, c(sigma_permutations, sp_permutations))
colnames(permutation_columns) <- paste0("p", seq_len(m))
candidate_scores <- cbind(candidate_scores, permutation_columns)

dir.create(
  "detailed_outputs/best_array_search",
  recursive = TRUE,
  showWarnings = FALSE
)

write.csv(
  selected_scores,
  "detailed_outputs/best_array_search/SOA_z6_selected_scores.csv",
  row.names = FALSE
)

write.csv(
  candidate_scores,
  "detailed_outputs/best_array_search/SOA_z6_candidate_scores.csv",
  row.names = FALSE
)

print(selected_scores, row.names = FALSE)
