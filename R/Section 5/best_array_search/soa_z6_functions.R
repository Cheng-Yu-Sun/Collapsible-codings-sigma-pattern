#####################################################################
# Common functions for Search_SOA_z6.R.
# A has 72 rows and 13 columns, with levels 0, 1, 2; B has levels 0, 1.
# The constraints ensure B is an OA(72,13,2,1) and 2A+B is an SOA2+.
#####################################################################


#####################################################################
# Calculate nu33 and nu42.
#####################################################################

coarse_score <- function(A) {
  triples <- combn(ncol(A), 3)
  sum(apply(triples, 2, function(j) {
    cells <- 9 * A[, j[1]] + 3 * A[, j[2]] + A[, j[3]] + 1
    mean((tabulate(cells, nbins = 27) - nrow(A) / 27)^2)
  }))
}

fine_score <- function(D) {
  pairs <- combn(ncol(D), 2)
  sum(apply(pairs, 2, function(j) {
    cells <- 6 * D[, j[1]] + D[, j[2]] + 1
    mean((tabulate(cells, nbins = 36) - nrow(D) / 36)^2)
  }))
}


#####################################################################
# Admissibility constraints and random initial B.
#####################################################################

balance_constraints <- function(A, i) {
  # Each eight-row (A_i, A_j) cell must have four zeros and four ones
  # in B_i, for every j != i. These constraints also balance B_i.
  do.call(rbind, lapply(setdiff(seq_len(ncol(A)), i), function(j) {
    cells <- 3 * A[, i] + A[, j]
    t(vapply(0:8, function(k) as.numeric(cells == k), numeric(nrow(A))))
  }))
}

solve_binary <- function(objective, constraints, directions, rhs,
                         binary = seq_along(objective)) {
  fit <- lpSolve::lp(
    direction = "min", objective.in = objective,
    const.mat = constraints, const.dir = directions, const.rhs = rhs,
    binary.vec = binary
  )
  if (fit$status != 0) {
    stop("Integer programming failed; solver status: ", fit$status)
  }
  fit$solution
}

random_B <- function(balance) {
  # Random objectives give feasible starting columns, not uniform draws
  # from the set of all admissible columns.
  N <- ncol(balance[[1]])
  B <- matrix(0L, N, length(balance))
  for (i in seq_along(balance)) {
    M <- balance[[i]]
    solution <- solve_binary(rnorm(N), M, rep("=", nrow(M)),
                             rep(4, nrow(M)))
    B[, i] <- as.integer(round(solution))
  }
  B
}


#####################################################################
# Update one B column by mixed-integer linear programming.
#####################################################################

update_column <- function(A, D, i, balance) {
  # Each (A_i=a, D_j=d) group has four rows. If k entries of B_i are
  # one, its contribution to nu42 is (k-2)^2/18.
  groups <- do.call(rbind, lapply(setdiff(seq_len(ncol(D)), i), function(j) {
    cells <- 6 * A[, i] + D[, j]
    t(vapply(0:17, function(k) as.numeric(cells == k), numeric(nrow(A))))
  }))

  b <- integer(nrow(A))

  # The three A_i strata can be optimized separately (24 binary entries each).
  for (a in 0:2) {
    rows <- which(A[, i] == a)
    M <- balance[rowSums(balance[, rows, drop = FALSE]) > 0, rows, drop = FALSE]
    G <- groups[rowSums(groups[, rows, drop = FALSE]) > 0, rows, drop = FALSE]
    ng <- nrow(G)

    # For integer k in 0:4, (k-2)^2 = max(2-k, k-2, 4-3*k, 3*k-8).
    # One continuous auxiliary variable per group represents this square.
    constraints <- rbind(
      cbind(M, matrix(0, nrow(M), ng)),
      cbind(G, diag(ng)),
      cbind(-G, diag(ng)),
      cbind(3 * G, diag(ng)),
      cbind(-3 * G, diag(ng))
    )
    directions <- c(rep("=", nrow(M)), rep(">=", 4 * ng))
    rhs <- c(rep(4, nrow(M)), rep(2, ng), rep(-2, ng),
             rep(4, ng), rep(-8, ng))
    objective <- c(rep(0, length(rows)), rep(1, ng))

    solution <- solve_binary(objective, constraints, directions, rhs,
                             binary = seq_along(rows))
    b[rows] <- as.integer(round(solution[seq_along(rows)]))
  }
  b
}


#####################################################################
# Coordinate updates and multiple random starts for a fixed A.
#####################################################################

search_B <- function(A, B, balance) {
  D <- 2 * A + B
  score <- fine_score(D)

  repeat {
    improved <- FALSE
    for (i in sample.int(ncol(A))) {
      candidate <- D
      candidate[, i] <- 2 * A[, i] + update_column(A, D, i, balance[[i]])
      candidate_score <- fine_score(candidate)

      # Accept only strict improvements, allowing for rounding error.
      if (candidate_score < score - 1e-10) {
        D <- candidate
        score <- candidate_score
        improved <- TRUE
      }
    }
    if (!improved) break
  }

  list(design = D, score = score)
}

multistart_B <- function(A, n_starts, seed) {
  set.seed(seed)
  balance <- lapply(seq_len(ncol(A)), function(i) balance_constraints(A, i))
  best_score <- Inf

  for (s in seq_len(n_starts)) {
    B <- random_B(balance)
    result <- search_B(A, B, balance)

    # Retain the first array attaining the best score found.
    if (result$score < best_score - 1e-10) {
      best <- result
      best_score <- result$score
    }
    if (s %% 10 == 0) message("  Completed ", s, "/", n_starts, " starts")
  }

  best
}
