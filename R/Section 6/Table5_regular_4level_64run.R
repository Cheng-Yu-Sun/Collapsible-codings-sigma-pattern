#####################################################################
# Search for Table 5:
# Sigma-optimal 4^{m-(m-3)} regular designs, m = 7, ..., 16.
#
# This script uses GF(4). The search has three main steps:
#   1. Find all 4^{7-4} designs with n1 = n2 = 0.
#   2. Sequentially extend them to m = 8, ..., 16 while preserving
#      n1 = n2 = 0.
#   3. Among designs with minimum n3, compute the Sigma-pattern and
#      select the Sigma-optimal design.
#
# Detailed results are saved as RDS; the paper table is saved as CSV.
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

source("R/Section 6/regular_design_functions.R")

dir.create("detailed_outputs/regular_designs", recursive = TRUE, showWarnings = FALSE)
dir.create("detailed_outputs/regular_designs/regular_4level_64run_intermediate", recursive = TRUE, showWarnings = FALSE)


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

# Generating matrix of the saturated 64-run design.
SAT <- t(full(4, 3)[c(2, 5, 17, 6:8, 18:32), 3:1])

# Includes equivalent columns and the null vector.
SAT.f <- t(full(4, 3)[, 3:1])


#####################################################################
# Step 1: find all 4^{7-4} designs with n1 = n2 = 0
#####################################################################

C1 <- combn(18, 4) + 3
C2 <- combn(7, 3)
F34 <- full(3, 4) + 1

N3 <- c()
G.list <- c()

cat("Searching 7-factor designs with n1 = n2 = 0...\n")

tic <- proc.time()

for (i in 1:ncol(C1)) {
  for (j in 1:3^4) {
    G <- SAT[, C1[, i]]
    G <- Mm(G, diag(F34[j, ]))
    G <- cbind(diag(3), G)
    
    n3 <- 0
    fail <- FALSE
    
    for (k in 1:ncol(C2)) {
      w <- defW(red(G[, C2[, k]]))
      
      # If the word has rank 1 or 2, then n1 > 0 or n2 > 0.
      if (!all(w == 0) && !all((1:3) %in% w)) {
        fail <- TRUE
        break
      }
      
      # In GF(4), any nonzero acceptable length-three word here has rank 3.
      if (!all(w == 0)) {
        n3 <- n3 + 1
      }
    }
    
    if (!fail) {
      G.list <- cbind(G.list, t(G) %*% 4^(0:2))
      N3 <- c(N3, n3)
    }
  }
  
  if (i %% 100 == 0 || i == ncol(C1)) {
    pct <- round(i / ncol(C1) * 100)
    
    if (is.null(G.list)) {
      cat(paste0(pct, "% completed; no 64x7 design found yet.\n"))
    } else {
      cat(paste0(
        pct, "% completed; ",
        ncol(G.list),
        " 64x7 designs with n1 = n2 = 0 found.\n"
      ))
    }
  }
}

cat("Step 1 finished in", round((proc.time() - tic)[3], 1), "seconds.\n\n")

# Remove duplicates up to column permutations / coordinate permutations.
U <- Unique(G.list)
G.list <- U[[1]]
N3 <- N3[U[[2]]]

saveRDS(G.list, "detailed_outputs/regular_designs/regular_4level_64run_intermediate/G_list_7.rds")
saveRDS(N3, "detailed_outputs/regular_designs/regular_4level_64run_intermediate/N3_7.rds")


#####################################################################
# Step 2: sequentially extend to m = 8, ..., 16
#####################################################################

G.list.old <- G.list
N3.old <- N3

for (l in 7:15) {
  # The old designs have l columns. We add one column to obtain l+1 columns.
  G.list <- c()
  N3 <- c()
  C2 <- combn(l, 2)
  
  cat("Extending to", l + 1, "factors...\n")
  tic <- proc.time()
  
  for (i in 1:ncol(G.list.old)) {
    for (j in 1:63) {
      g <- c(G.list.old[, i], j)
      
      # Exclude duplicated columns.
      if (n_dis(g) != l + 1) next
      
      G <- SAT.f[, g + 1]
      
      # Exclude equivalent columns.
      g.one <- t(apply(G, 2, one)) %*% 4^(0:2)
      if (n_dis(g.one) != l + 1) next
      
      n3 <- N3.old[i]
      fail <- FALSE
      
      # Only new three-factor projections need to be checked:
      # each old pair of columns plus the newly added column.
      for (k in 1:ncol(C2)) {
        w <- defW(red(G[, c(C2[, k], l + 1)]))
        
        if (!all(w == 0) && !all((1:3) %in% w)) {
          fail <- TRUE
          break
        }
        
        if (!all(w == 0)) {
          n3 <- n3 + 1
        }
      }
      
      if (!fail) {
        G.list <- cbind(G.list, t(G) %*% 4^(0:2))
        N3 <- c(N3, n3)
      }
    }
  }
  
  U <- Unique(G.list)
  G.list.unique <- U[[1]]
  N3.unique <- N3[U[[2]]]
  
  saveRDS(
    G.list.unique,
    paste0("detailed_outputs/regular_designs/regular_4level_64run_intermediate/G_list_", l + 1, ".rds")
  )
  
  saveRDS(
    N3.unique,
    paste0("detailed_outputs/regular_designs/regular_4level_64run_intermediate/N3_", l + 1, ".rds")
  )
  
  cat(
    l + 1,
    "-factor case finished in",
    round((proc.time() - tic)[3], 1),
    "seconds.\n\n"
  )
  
  G.list.old <- G.list.unique
  N3.old <- N3.unique
}


#####################################################################
# Step 3: compute Sigma-patterns for candidates with minimum n3
#####################################################################

cat("Computing Sigma-patterns for minimum-n3 candidates...\n")

for (l in 7:16) {
  G.list <- readRDS(paste0("detailed_outputs/regular_designs/regular_4level_64run_intermediate/G_list_", l, ".rds"))
  N3 <- readRDS(paste0("detailed_outputs/regular_designs/regular_4level_64run_intermediate/N3_", l, ".rds"))
  
  idx_min <- which(N3 == min(N3))
  G.list.min <- G.list[, idx_min, drop = FALSE]
  N3.min <- N3[idx_min]
  
  saveRDS(
    G.list.min,
    paste0("detailed_outputs/regular_designs/regular_4level_64run_intermediate/G_list_min_", l, ".rds")
  )
  
  # WLP(d, 2, 2) returns the dimension-by-weight table.
  # The Sigma-pattern is obtained by column-major comparison,
  # implemented by comparing t(Sigma_table).
  Sigma_table <- array(0, c(l, 2 * l, length(N3.min)))
  
  for (i in 1:length(N3.min)) {
    G <- SAT.f[, G.list.min[, i] + 1]
    d <- Mm(t(SAT.f), G)
    Sigma_table[, , i] <- WLP(d, 2, 2)
    
    if (i %% 5 == 0 || i == length(N3.min)) {
      cat(paste0(
        l, "-factor case: ",
        round(i / length(N3.min) * 100, 2),
        "% completed.\n"
      ))
    }
  }
  
  saveRDS(
    Sigma_table,
    paste0("detailed_outputs/regular_designs/regular_4level_64run_intermediate/Sigma_table_", l, ".rds")
  )
}


#####################################################################
# Step 4: select the Sigma-optimal design for each m
#####################################################################

results <- data.frame(
  m = integer(),
  generating_matrix = character(),
  n1 = integer(),
  n2 = integer(),
  n3 = integer(),
  search_type = character()
)

for (l in 7:16) {
  Sigma_table <- readRDS(
    paste0("detailed_outputs/regular_designs/regular_4level_64run_intermediate/Sigma_table_", l, ".rds")
  )
  
  G.list.min <- readRDS(
    paste0("detailed_outputs/regular_designs/regular_4level_64run_intermediate/G_list_min_", l, ".rds")
  )
  
  N3 <- readRDS(
    paste0("detailed_outputs/regular_designs/regular_4level_64run_intermediate/N3_", l, ".rds")
  )
  
  N3.min <- min(N3)
  
  best <- 1
  
  if (dim(Sigma_table)[3] >= 2) {
    for (i in 2:dim(Sigma_table)[3]) {
      if (compare(t(Sigma_table[, , best]), t(Sigma_table[, , i])) == 1) {
        best <- i
      }
    }
  }
  
  generator <- sort(as.vector(G.list.min[, best]))
  
  results <- rbind(
    results,
    data.frame(
      m = l,
      generating_matrix = paste(generator, collapse = ", "),
      n1 = 0,
      n2 = 0,
      n3 = N3.min,
      search_type = "iterative search with complete Sigma-pattern comparison"
    )
  )
}

print(results)


#####################################################################
# Save final Table 5 results
#####################################################################

saveRDS(
  results,
  file = "detailed_outputs/regular_designs/regular_4level_64run.rds"
)


#####################################################################
# Paper table: this section can also use the saved results alone.
#####################################################################

results <- readRDS("detailed_outputs/regular_designs/regular_4level_64run.rds")
dir.create("Tables", showWarnings = FALSE)
write.csv(results, "Tables/Table5.csv", row.names = FALSE)
message("Saved Tables/Table5.csv")
