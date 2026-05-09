#####################################################################
# Common functions for the regular-design searches in Tables 1--3.
#
# This file contains:
#   1. finite-field arithmetic based on addition/multiplication tables;
#   2. utilities for generating full factorial designs;
#   3. functions for defining words and word-rank calculations;
#   4. functions for computing WLP/SP/SFP/Sigma-related patterns.
#
# The scripts Table1_regular_4level_16run.R, Table2_regular_4level_64run.R,
# and Table3_regular_9level_81run.R source this file.
#####################################################################


#####################################################################
# Full factorial design
#####################################################################

full <- function(l, n) {
  # Return an l^n full factorial design with levels 0, ..., l-1.
  if (n == 1) return(cbind(1:l - 1))
  
  a <- as.integer(gl(l, l^n / l, l^n))
  
  if (n > 1) {
    for (i in 2:n) {
      a <- cbind(a, as.integer(gl(l, l^(n - i), l^n)))
    }
  }
  
  d <- a - 1
  colnames(d) <- 1:ncol(d)
  return(d)
}


#####################################################################
# Finite-field arithmetic
#
# These functions use the global objects AT and MT, which should be
# defined in the calling script.
#
# AT: addition table of GF(z)
# MT: multiplication table of GF(z)
#####################################################################

add <- function(x, y) {
  # Compute x + y over GF(z).
  # y may be a scalar.
  if (length(y) == 1) y <- rep(y, length(x))
  return(AT[x * ncol(AT) + y + 1])
}

mul <- function(x, y) {
  # Compute x * y over GF(z).
  # x may be a scalar.
  if (length(x) == 1) x <- rep(x, length(y))
  return(MT[x * ncol(MT) + y + 1])
}

SUM <- function(x) {
  # Sum all elements in a vector over GF(z).
  z <- 0
  for (i in 1:length(x)) {
    z <- AT[z + 1, x[i] + 1]
  }
  return(z)
}

PROD <- function(x) {
  # Multiply all elements in a vector over GF(z).
  z <- 1
  for (i in 1:length(x)) {
    z <- MT[z + 1, x[i] + 1]
  }
  return(z)
}

Mm <- function(X, Y) {
  # Matrix multiplication X %*% Y over GF(z).
  M <- matrix(0, nrow(X), ncol(Y))
  
  for (i in 1:nrow(X)) {
    for (j in 1:ncol(Y)) {
      M[i, j] <- SUM(mul(X[i, ], Y[, j]))
    }
  }
  
  return(M)
}

rec <- function(x) {
  # Pointwise reciprocal over GF(z), with 0^{-1} defined as 0.
  c(0, (which(MT == 1) - 1) %% nrow(MT))[x + 1]
}


#####################################################################
# Row reduction, defining words, and canonical representatives
#####################################################################

red <- function(X) {
  # Return a reduced form of matrix X over GF(z).
  Y <- rbind(X)
  mini <- min(nrow(Y), ncol(Y))
  
  if (mini == 1) return(X)
  
  neg1 <- which(AT[, 2] == 0) - 1  # additive inverse of 1
  
  for (i in 1:mini) {
    if (Y[i, i] == 0 && i < nrow(Y)) {
      for (k in (i + 1):nrow(Y)) {
        if (Y[k, i] == 0) next
        Y[c(i, k), ] <- Y[c(k, i), ]
        break
      }
    }
    
    if (Y[i, i] == 0) next
    
    for (j in (1:nrow(Y))[-i]) {
      r <- PROD(c(rec(Y[i, i]), Y[j, i], neg1))
      if (r == 0) next
      Y[j, ] <- add(mul(r, Y[i, ]), Y[j, ])
    }
  }
  
  return(Y)
}

defW <- function(G) {
  # Given a reduced three-column generating matrix G,
  # return the corresponding defining word.
  #
  # If no length-three word exists, the function returns c(0, 0, 0).
  
  if (nrow(G) == 3) {
    if (!all(G[3, ] == 0)) return(c(0, 0, 0))
  }
  
  for (i in 1:2) {
    G[i, ] <- mul(rec(G[i, i]), G[i, ])
  }
  
  neg1 <- which(AT[, 2] == 0) - 1
  return(c(G[1:2, 3], neg1))
}

one <- function(x) {
  # Normalize a nonzero vector x by multiplying a scalar on the left
  # so that the last nonzero entry becomes 1.
  r <- rev(x[which(x != 0)])[1]
  return(mul(rec(r), x))
}

n_dis <- function(x) {
  # Number of distinct values in x.
  return(length(unique(x)))
}


#####################################################################
# Removing duplicate designs
#
# This function removes duplicate generating matrices after accounting
# for column permutations. It also checks coordinate permutations for
# k = 2 or k = 3, where k is the number of independent columns.
#
# The function uses the global object SAT.f, which should be defined
# in the calling script.
#####################################################################

Unique <- function(X) {
  # Remove duplicate designs from X, where each column of X encodes
  # one generating matrix.
  
  z <- ncol(SAT.f)^(1 / nrow(SAT.f))
  
  if (ncol(X) == 1) return(list(X, 1))
  
  X.u <- list(X[, 1])
  index <- 1
  
  for (i in 2:ncol(X)) {
    IN <- FALSE
    
    for (j in 1:length(X.u)) {
      if (setequal(X[, i], X.u[[j]])) {
        IN <- TRUE
        break
      }
      
      if (nrow(SAT.f) == 2) {
        if (setequal(
          z^(0:1) %*% SAT.f[, X[, i] + 1],
          z^(1:0) %*% SAT.f[, X.u[[j]] + 1]
        )) {
          IN <- TRUE
          break
        }
      }
      
      if (nrow(SAT.f) == 3) {
        if (setequal(
          z^(0:2) %*% SAT.f[, X[, i] + 1],
          z^(c(0, 2, 1)) %*% SAT.f[, X.u[[j]] + 1]
        )) {
          IN <- TRUE
          break
        }
        
        if (setequal(
          z^(0:2) %*% SAT.f[, X[, i] + 1],
          z^(c(1, 0, 2)) %*% SAT.f[, X.u[[j]] + 1]
        )) {
          IN <- TRUE
          break
        }
        
        if (setequal(
          z^(0:2) %*% SAT.f[, X[, i] + 1],
          z^(c(1, 2, 0)) %*% SAT.f[, X.u[[j]] + 1]
        )) {
          IN <- TRUE
          break
        }
        
        if (setequal(
          z^(0:2) %*% SAT.f[, X[, i] + 1],
          z^(c(2, 0, 1)) %*% SAT.f[, X.u[[j]] + 1]
        )) {
          IN <- TRUE
          break
        }
        
        if (setequal(
          z^(0:2) %*% SAT.f[, X[, i] + 1],
          z^(c(2, 1, 0)) %*% SAT.f[, X.u[[j]] + 1]
        )) {
          IN <- TRUE
          break
        }
      }
    }
    
    if (!IN) {
      X.u[[length(X.u) + 1]] <- X[, i]
      index[length(index) + 1] <- i
    }
  }
  
  return(list(sapply(X.u, identity), index))
}


#####################################################################
# Comparing two patterns sequentially
#####################################################################

compare <- function(a, b) {
  # Sequentially compare two vectors a and b.
  #
  # Return:
  #   0   if a is preferred to b;
  #   1   if b is preferred to a;
  #   0.5 if they are identical.
  
  for (i in 1:length(a)) {
    if (a[i] < b[i]) return(0)
    if (b[i] < a[i]) return(1)
  }
  
  return(0.5)
}


#####################################################################
# Functions for calculating GWLP, beta-WLP, SFP, SP, and Sigma-pattern
#####################################################################

WLP <- function(d, s, p, type = "SP", r = 3) {
  # Return the wordlength pattern of a specified type.
  #
  # type = "SP"   : stratification pattern / dimension-by-weight table
  # type = "SFP"  : space-filling pattern
  # type = "GWLP" : generalized wordlength pattern
  # type = "beta" : beta wordlength pattern
  #
  # References:
  #   Tang and Xu (2021)
  #   Tian and Xu (2024)
  
  if (type == "SP") {
    E <- Ed(d, s, p)
    V <- Vd(ncol(d) * (p * ncol(d) + 1))
    ans <- round(matrix(Re(V %*% (E - 1)), ncol(d)), r)
    rownames(ans) <- 1:ncol(d)
    colnames(ans) <- 0:(p * ncol(d))
    return(ans[, -1])
  }
  
  if (type == "SFP") {
    E <- Ed(d, s, p, type = "SFP")
    V <- Vd(p * ncol(d))
    return(round(as.vector(Re(V %*% (E - 1))), r))
  }
  
  if (type == "GWLP") {
    E <- Ed(d, s, p, "GWLP")
    V <- Vd(ncol(d))
    return(round(as.vector(Re(V %*% (E - 1))), r))
  }
  
  if (type == "beta") {
    E <- Ed(d, s, p, "beta")
    V <- Vd((s - 1) * ncol(d))
    return(round(as.vector(Re(V %*% (E - 1))), r))
  }
}

HA <- function(s, p) {
  # Construct the real-valued full-factorial-based coding matrix.
  
  h <- contr.poly(s) * s^0.5
  h <- cbind(1, h)
  
  if (p == 1) return(h)
  
  ff <- full(s, p) + 1
  H <- c()
  
  for (i in 1:nrow(ff)) {
    int <- h[, ff[i, 1]]
    
    for (j in 2:ncol(ff)) {
      int <- h[, ff[i, j]] %x% int
    }
    
    H <- cbind(H, int)
  }
  
  colnames(H) <- 1:ncol(H) - 1
  return(H)
}

Vd <- function(K, inv = TRUE) {
  # Fourier matrix used in the wordlength-pattern calculation.
  
  w <- exp(2 * pi * 1i / K)
  V <- matrix(0, K, K)
  
  for (i in 1:nrow(V)) {
    for (j in 1:ncol(V)) {
      V[i, j] <- w^(i * j)
    }
  }
  
  if (inv) return(Conj(t(V)) / K)
  return(V)
}

Ed <- function(d, s, p, type = "STP") {
  # Helper function for WLP().
  
  z <- s^p
  H <- HA(s, p)
  
  a <- c(0, rep(1, z - 1))
  b <- ceiling(log(1:z, s))
  
  if (type == "GWLP") K <- ncol(d)
  if (type == "STP") K <- ncol(d) * (p * ncol(d) + 1)
  if (type == "SFP") K <- p * ncol(d)
  if (type == "beta") K <- (s - 1) * ncol(d)
  
  w <- exp(2 * pi * 1i / K)^(1:K)
  R <- array(0, c(z, z, K))
  
  for (k in 1:K) {
    for (i in 1:z) {
      for (j in 1:z) {
        if (type == "GWLP") {
          R[i, j, k] <- sum(H[i, ] * H[j, ] * w[k]^a)
        }
        
        if (type == "STP") {
          R[i, j, k] <- sum(H[i, ] * H[j, ] * w[k]^(b * ncol(d) + a))
        }
        
        if (type == "SFP") {
          R[i, j, k] <- sum(H[i, ] * H[j, ] * w[k]^b)
        }
        
        if (type == "beta") {
          R[i, j, k] <- sum(H[i, ] * H[j, ] * w[k]^(1:s - 1))
        }
      }
    }
  }
  
  E <- rep(0, K)
  
  for (k in 1:K) {
    for (i in 1:nrow(d)) {
      for (j in 1:nrow(d)) {
        E[k] <- E[k] + prod(diag(R[, , k][d[i, ] + 1, d[j, ] + 1]))
      }
    }
    
    E[k] <- E[k] / nrow(d)^2
  }
  
  return(E)
}