#####################################################################
# Exhaustive search for Table 1:
# Sigma-optimal 4^{m-(m-2)} regular designs, m = 3, 4, 5.
#
# This script uses GF(4) and searches over all relevant regular
# 16-run designs. The final results are saved as both .rds and .csv.
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

source("R/regular_design_functions.R")

#####################################################################
# Setup for GF(4)
#####################################################################

AT <- cbind(
  c(0, 1, 2, 3),
  c(1, 0, 3, 2),
  c(2, 3, 0, 1),
  c(3, 2, 1, 0)
) # addition table of GF(4)

MT <- cbind(
  0,
  0:3,
  c(0, 2, 3, 1),
  c(0, 3, 1, 2)
) # multiplication table of GF(4)

SAT <- t(full(4, 2)[c(2, 5:8), 2:1])
SAT.f <- t(full(4, 2)[, 2:1])


#####################################################################
# Search function for a fixed m
#####################################################################

search_table1_for_m <- function(m) {
  C1 <- combn(3, m - 2) + 2
  C2 <- combn(m, 3)
  FF <- full(3, m - 2) + 1
  
  R.store <- list()
  G.list <- list()
  
  for (i in 1:ncol(C1)) {
    for (j in 1:nrow(FF)) {
      G <- cbind(SAT[, C1[, i]])
      
      if (m == 3) {
        G <- mul(G, FF[j])
      }
      
      if (m != 3) {
        G <- Mm(G, diag(FF[j, ]))
      }
      
      G <- cbind(diag(2), G)
      g <- t(apply(G, 2, one)) %*% 4^(0:1)
      
      R <- c(0, 0, 0)
      
      for (k in 1:ncol(C2)) {
        w <- defW(red(G[, C2[, k]]))
        r <- n_dis(w)
        R[r] <- R[r] + 1
      }
      
      R.store[[length(R.store) + 1]] <- R
      G.list[[length(G.list) + 1]] <- t(G) %*% 4^(0:1)
    }
  }
  
  R.store <- sapply(R.store, identity)
  G.list <- sapply(G.list, identity)
  
  # Remove duplicate designs.
  U <- Unique(G.list)
  G.list <- U[[1]]
  R.store <- R.store[, U[[2]], drop = FALSE]
  
  # Sequentially minimize n1, n2, and n3.
  idx <- which(R.store[1, ] == min(R.store[1, ]))
  G.list <- G.list[, idx, drop = FALSE]
  R.store <- R.store[, idx, drop = FALSE]
  
  idx <- which(R.store[2, ] == min(R.store[2, ]))
  G.list <- G.list[, idx, drop = FALSE]
  R.store <- R.store[, idx, drop = FALSE]
  
  idx <- which(R.store[3, ] == min(R.store[3, ]))
  G.list <- G.list[, idx, drop = FALSE]
  R.store <- R.store[, idx, drop = FALSE]
  
  # Use the first optimal design found.
  generator <- sort(as.vector(G.list[, 1]))
  rank_counts <- as.vector(R.store[, 1])
  
  data.frame(
    m = m,
    generating_matrix = paste(generator, collapse = ", "),
    n1 = rank_counts[1],
    n2 = rank_counts[2],
    n3 = rank_counts[3],
    search_type = "exhaustive"
  )
}


#####################################################################
# Run the searches for m = 3, 4, 5
#####################################################################

Table1_results <- do.call(
  rbind,
  lapply(3:5, search_table1_for_m)
)

print(Table1_results)


#####################################################################
# Save results
#####################################################################

if (!dir.exists("outputs")) {
  dir.create("outputs")
}

if (!dir.exists("outputs/tables")) {
  dir.create("outputs/tables")
}

saveRDS(
  Table1_results,
  file = "outputs/tables/Table1_regular_4_k2.rds"
)

write.csv(
  Table1_results,
  file = "outputs/tables/Table1_regular_4_k2.csv",
  row.names = FALSE
)