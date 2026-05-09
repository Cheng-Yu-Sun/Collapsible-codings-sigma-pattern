#####################################################################
# Search for Table 3:
# Sigma-optimal / nearly Sigma-optimal 9^{m-(m-2)} regular designs,
# m = 3, ..., 10.
#
# This script uses GF(9). The search has two parts:
#   1. Exhaustive search for m = 3, 4, 5.
#   2. Partial iterative search for m = 6, ..., 10.
#
# For m = 3, 4, 5, candidates are selected by sequentially minimizing
# (n1, n2, n3), and then the Sigma-pattern is computed to select the
# Sigma-optimal design.
#
# For m = 6, ..., 10, exhaustive search is computationally intensive.
# We therefore use the partial iterative procedure described in the
# Supplement. These designs are reported as nearly Sigma-optimal.
#
# The final Table 3 results are saved as both .rds and .csv.
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
# Setup for GF(9)
#####################################################################

AT <- matrix(
  c(
    0, 1, 2, 3, 4, 5, 6, 7, 8,
    1, 2, 0, 4, 5, 3, 7, 8, 6,
    2, 0, 1, 5, 3, 4, 8, 6, 7,
    3, 4, 5, 6, 7, 8, 0, 1, 2,
    4, 5, 3, 7, 8, 6, 1, 2, 0,
    5, 3, 4, 8, 6, 7, 2, 0, 1,
    6, 7, 8, 0, 1, 2, 3, 4, 5,
    7, 8, 6, 1, 2, 0, 4, 5, 3,
    8, 6, 7, 2, 0, 1, 5, 3, 4
  ),
  9
) # addition table of GF(9)

MT <- matrix(
  c(
    0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 1, 2, 3, 4, 5, 6, 7, 8,
    0, 2, 1, 6, 8, 7, 3, 5, 4,
    0, 3, 6, 7, 1, 4, 5, 8, 2,
    0, 4, 8, 1, 5, 6, 2, 3, 7,
    0, 5, 7, 4, 6, 2, 8, 1, 3,
    0, 6, 3, 5, 2, 8, 7, 4, 1,
    0, 7, 5, 8, 3, 1, 4, 2, 6,
    0, 8, 4, 2, 7, 3, 1, 6, 5
  ),
  9
) # multiplication table of GF(9)

# Generating matrix of the saturated 81-run design.
SAT <- t(full(9, 2)[c(2, 10, 11:18), 2:1])

# Includes equivalent columns and the null vector.
SAT.f <- t(full(9, 2)[, 2:1])

# All additive subgroups of size 3 in GF(9), excluding the zero-only case.
# This is used to compute the rank of a defining word.
H <- MT[1:3, c(2, 4, 5, 6)]


#####################################################################
# Output folders
#####################################################################

if (!dir.exists("outputs")) {
  dir.create("outputs")
}

if (!dir.exists("outputs/tables")) {
  dir.create("outputs/tables")
}

if (!dir.exists("outputs/Table3_intermediate")) {
  dir.create("outputs/Table3_intermediate")
}


#####################################################################
# Helper function: word rank over GF(9)
#####################################################################

word_rank_GF9 <- function(w) {
  # The rank is the number of additive subgroups H_theta that contain
  # at least one nonzero entry of the word.
  r <- apply(matrix(H %in% w, 3), 2, sum)
  return(sum(r != 0))
}


#####################################################################
# Step 1: exhaustive search for m = 3, 4, 5
#####################################################################

cat("Starting exhaustive search for m = 3, 4, 5...\n\n")

for (m in 3:5) {
  cat("Searching", m, "-factor designs...\n")
  tic <- proc.time()
  
  C1 <- combn(8, m - 2) + 2
  C2 <- combn(m, 3)
  FF <- full(8, m - 2) + 1
  
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
      g <- t(apply(G, 2, one)) %*% 9^(0:1)
      
      R <- c(0, 0, 0)
      
      for (k in 1:ncol(C2)) {
        w <- defW(red(G[, C2[, k]]))
        r <- word_rank_GF9(w)
        R[r] <- R[r] + 1
      }
      
      R.store[[length(R.store) + 1]] <- R
      G.list[[length(G.list) + 1]] <- t(G) %*% 9^(0:1)
    }
  }
  
  R.store <- sapply(R.store, identity)
  G.list <- sapply(G.list, identity)
  
  # Sequentially minimize n1 and then n2.
  # For fixed m, n3 is determined after n1 and n2 are fixed because
  # every three-factor projection has either no word or a word of rank 1, 2, or 3.
  idx <- which(R.store[1, ] == min(R.store[1, ]))
  G.list <- G.list[, idx, drop = FALSE]
  R.store <- R.store[, idx, drop = FALSE]
  
  idx <- which(R.store[2, ] == min(R.store[2, ]))
  G.list <- G.list[, idx, drop = FALSE]
  R.store <- R.store[, idx, drop = FALSE]
  
  # Remove duplicate designs.
  U <- Unique(G.list)
  G.list <- U[[1]]
  R.store <- R.store[, U[[2]], drop = FALSE]
  
  saveRDS(
    G.list,
    paste0("outputs/Table3_intermediate/G_list_", m, ".rds")
  )
  
  saveRDS(
    R.store,
    paste0("outputs/Table3_intermediate/rank_counts_", m, ".rds")
  )
  
  cat(
    m,
    "-factor exhaustive search finished in",
    round((proc.time() - tic)[3], 1),
    "seconds.\n\n"
  )
}


#####################################################################
# Step 2: partial iterative search for m = 6, ..., 10
#####################################################################

cat("Starting partial iterative search for m = 6, ..., 10...\n\n")

F92 <- t(full(9, 2)[, 2:1])

G.list.old <- readRDS("outputs/Table3_intermediate/G_list_5.rds")

# The stored 5-factor designs have the current minimum n1 and n2.
# In the original computation, these values are initialized as below.
m1 <- 0
m2 <- 2

for (l in 5:9) {
  # The old designs have l columns. We add one column to obtain l+1 columns.
  cat("Extending to", l + 1, "factors by partial search...\n")
  tic <- proc.time()
  
  G.list <- c()
  N1 <- c()
  N2 <- c()
  C2 <- combn(l, 2)
  
  for (i in 1:ncol(G.list.old)) {
    for (j in 1:80) {
      g <- c(G.list.old[, i], j)
      
      # Exclude duplicated columns.
      if (n_dis(g) != l + 1) next
      
      G <- F92[, g + 1]
      
      # Exclude equivalent columns.
      g.one <- t(apply(G, 2, one)) %*% 9^(0:1)
      if (n_dis(g.one) != l + 1) next
      
      n1 <- 0
      n2 <- 0
      
      # Only new three-factor projections need to be checked:
      # each old pair of columns plus the newly added column.
      for (k in 1:ncol(C2)) {
        w <- defW(red(G[, c(C2[, k], l + 1)]))
        r <- word_rank_GF9(w)
        
        if (r == 1) n1 <- n1 + 1
        if (r == 2) n2 <- n2 + 1
      }
      
      G.list <- cbind(G.list, t(G) %*% 9^(0:1))
      N1[length(N1) + 1] <- n1 + m1
      N2[length(N2) + 1] <- n2 + m2
    }
  }
  
  # Sequentially minimize n1 and then n2.
  idx <- which(N1 == min(N1))
  G.list <- G.list[, idx, drop = FALSE]
  N2 <- N2[idx]
  
  idx <- which(N2 == min(N2))
  G.list <- G.list[, idx, drop = FALSE]
  
  U <- Unique(G.list)
  G.list.unique <- U[[1]]
  
  saveRDS(
    G.list.unique,
    paste0("outputs/Table3_intermediate/G_list_", l + 1, ".rds")
  )
  
  m1 <- min(N1)
  m2 <- min(N2)
  n3 <- choose(l + 1, 3) - m1 - m2
  
  saveRDS(
    c(n1 = m1, n2 = m2, n3 = n3),
    paste0("outputs/Table3_intermediate/rank_counts_", l + 1, ".rds")
  )
  
  cat(
    l + 1,
    "-factor partial search finished in",
    round((proc.time() - tic)[3], 1),
    "seconds; (n1, n2, n3) = (",
    m1, ", ", m2, ", ", n3, ").\n\n",
    sep = ""
  )
  
  G.list.old <- G.list.unique
}


#####################################################################
# Step 3: compute Sigma-patterns for m = 3, 4, 5
#####################################################################

cat("Computing Sigma-patterns for m = 3, 4, 5...\n\n")

for (l in 3:5) {
  G.list <- readRDS(paste0("outputs/Table3_intermediate/G_list_", l, ".rds"))
  
  # WLP(d, 3, 2) returns the dimension-by-weight table.
  # The Sigma-pattern is obtained by column-major comparison,
  # implemented by comparing t(Sigma_table).
  Sigma_table <- array(0, c(l, 2 * l, ncol(G.list)))
  
  for (i in 1:ncol(G.list)) {
    G <- SAT.f[, G.list[, i] + 1]
    d <- Mm(t(SAT.f), G)
    Sigma_table[, , i] <- WLP(d, 3, 2)
    
    if (i %% 5 == 0 || i == ncol(G.list)) {
      cat(paste0(
        l, "-factor case: ",
        round(i / ncol(G.list) * 100, 2),
        "% completed.\n"
      ))
    }
  }
  
  saveRDS(
    Sigma_table,
    paste0("outputs/Table3_intermediate/Sigma_table_", l, ".rds")
  )
}


#####################################################################
# Step 4: select final designs for Table 3
#####################################################################

Table3_results <- data.frame(
  m = integer(),
  generating_matrix = character(),
  n1 = integer(),
  n2 = integer(),
  n3 = integer(),
  search_type = character()
)

# m = 3, 4, 5: exhaustive search + Sigma-pattern comparison.
for (l in 3:5) {
  Sigma_table <- readRDS(
    paste0("outputs/Table3_intermediate/Sigma_table_", l, ".rds")
  )
  
  G.list <- readRDS(
    paste0("outputs/Table3_intermediate/G_list_", l, ".rds")
  )
  
  R.store <- readRDS(
    paste0("outputs/Table3_intermediate/rank_counts_", l, ".rds")
  )
  
  best <- 1
  
  if (dim(Sigma_table)[3] >= 2) {
    for (i in 2:dim(Sigma_table)[3]) {
      if (compare(t(Sigma_table[, , best]), t(Sigma_table[, , i])) == 1) {
        best <- i
      }
    }
  }
  
  generator <- sort(as.vector(G.list[, best]))
  rank_counts <- as.vector(R.store[, best])
  
  Table3_results <- rbind(
    Table3_results,
    data.frame(
      m = l,
      generating_matrix = paste(generator, collapse = ", "),
      n1 = rank_counts[1],
      n2 = rank_counts[2],
      n3 = rank_counts[3],
      search_type = "exhaustive search with complete Sigma-pattern comparison"
    )
  )
}

# m = 6, ..., 10: partial / nearly optimal search.
for (l in 6:10) {
  G.list <- readRDS(
    paste0("outputs/Table3_intermediate/G_list_", l, ".rds")
  )
  
  rank_counts <- readRDS(
    paste0("outputs/Table3_intermediate/rank_counts_", l, ".rds")
  )
  
  generator <- sort(as.vector(G.list[, 1]))
  
  Table3_results <- rbind(
    Table3_results,
    data.frame(
      m = l,
      generating_matrix = paste(generator, collapse = ", "),
      n1 = rank_counts["n1"],
      n2 = rank_counts["n2"],
      n3 = rank_counts["n3"],
      search_type = "partial iterative search; reported as nearly Sigma-optimal"
    )
  )
}

print(Table3_results)


#####################################################################
# Save final Table 3 results
#####################################################################

saveRDS(
  Table3_results,
  file = "outputs/tables/Table3_regular_9level_k2.rds"
)

write.csv(
  Table3_results,
  file = "outputs/tables/Table3_regular_9level_k2.csv",
  row.names = FALSE
)