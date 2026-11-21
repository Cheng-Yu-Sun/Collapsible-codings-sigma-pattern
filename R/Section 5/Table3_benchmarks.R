#####################################################################
# Table 3: GP prediction for the four benchmark functions.
# Run All performs both simulations and saves Tables/Table3.csv.
# With saved NRMSE matrices, run only the table section below.
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


source("R/Section 5/benchmark_z4.R")
source("R/Section 5/benchmark_z6.R")

#####################################################################
# Table section: summarize the saved NRMSE matrices.
#####################################################################

function_names <- c("OTL_Circuit", "Borehole", "Steel_Column", "Wing_Weight")
function_labels <- c("OTL Circuit", "Borehole", "Steel Column", "Wing Weight")
methods <- c("Sigma", "SP", "UD", "Random", "Maximin", "MaxPro")
summary_rows <- list()
table_rows <- list()

for (z in c(4, 6)) {
  for (k in seq_along(function_names)) {
    NRMSE <- readRDS(paste0("detailed_outputs/simulations/NRMSE_z", z, "_",
                            function_names[k], ".rds"))
    mu <- colMeans(NRMSE)
    s <- apply(NRMSE, 2, sd)
    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      z = z, test_function = function_labels[k], method = methods,
      mean = mu, sd = s, n_rep = nrow(NRMSE)
    )

    # Cells show mean NRMSE (standard deviation in units of 0.001).
    cells <- sprintf("%.3f (%.1f)", mu, 1000 * s)
    table_rows[[length(table_rows) + 1]] <- data.frame(
      z = z, test_function = function_labels[k],
      as.list(setNames(cells, methods)), check.names = FALSE
    )
  }
}

summary_table <- do.call(rbind, summary_rows)
write.csv(summary_table, "detailed_outputs/simulations/benchmark_summary.csv",
          row.names = FALSE)
Table3 <- do.call(rbind, table_rows)
rownames(Table3) <- NULL

dir.create("Tables", showWarnings = FALSE)
write.csv(Table3, "Tables/Table3.csv", row.names = FALSE)
print(Table3, row.names = FALSE)
message("Saved Tables/Table3.csv")
