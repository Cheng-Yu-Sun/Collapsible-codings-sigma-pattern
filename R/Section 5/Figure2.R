#####################################################################
# Figure 2: self-defined test function, z = 4 and z = 6.
# Run All performs both simulations and saves one complete Figure2.pdf.
# With saved NRMSE matrices, run only the plotting section below.
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


# Recalculate the two design settings.
source("R/Section 5/self_defined_z4.R")
source("R/Section 5/self_defined_z6.R")

#####################################################################
# Plotting section: read the saved results and draw the complete figure.
#####################################################################

NRMSE <- list(
  z4_lambda3 = readRDS("detailed_outputs/simulations/NRMSE_z4_lambda3.rds"),
  z4_lambda30 = readRDS("detailed_outputs/simulations/NRMSE_z4_lambda30.rds"),
  z6_lambda3 = readRDS("detailed_outputs/simulations/NRMSE_z6_lambda3.rds"),
  z6_lambda30 = readRDS("detailed_outputs/simulations/NRMSE_z6_lambda30.rds")
)

# Leave NULL for automatic limits, or supply c(lower, upper).
ylim_z4 <- NULL
ylim_z6 <- NULL
colors <- c("#80A6BF", "#E6A160", "#9FBF97", "#E5E5E5", "#B9ACCC", "#CEB891")
methods <- c("Sigma", "SP", "UD", "Random", "Maximin", "MaxPro")
classes <- c("(i)", "(ii)", "(iii)", "(iv)", "(v)", "(vi)")

draw_figure2 <- function() {
  layout(matrix(c(1, 1, 2, 3), nrow = 2, byrow = TRUE),
         heights = c(0.15, 0.85))
  par(mar = c(0, 0, 0, 0))
  plot.new()
  legend("center", legend = paste(classes, methods), fill = colors, border = "#444444",
         horiz = TRUE, bty = "n", cex = 0.95)

  for (z in c(4, 6)) {
    values <- c(as.list(as.data.frame(NRMSE[[paste0("z", z, "_lambda3")]])),
                as.list(as.data.frame(NRMSE[[paste0("z", z, "_lambda30")]])))
    limits <- if (z == 4) ylim_z4 else ylim_z6
    if (is.null(limits)) {
      limits <- range(vapply(values, function(x) boxplot.stats(x)$stats, numeric(5)))
      limits <- limits + c(-1, 1) * 0.04 * diff(limits)
    }

    # Each panel has its own scale; the right panel uses a right-hand axis.
    par(mar = if (z == 4) c(4.4, 4.0, 2.1, 0.8) else c(4.4, 0.8, 2.1, 4.0),
        mgp = c(2, 0.6, 0), las = 1)
    boxplot(values, at = c(1:6, 9:14), col = rep(colors, 2),
            border = "#444444", boxwex = 0.72, outline = FALSE,
            axes = FALSE, xlim = c(0.1, 14.9), ylim = limits,
            yaxs = "i", main = paste("z =", z), font.main = 1)
    axis(if (z == 4) 2 else 4)
    # Class labels identify the boxes in both color and grayscale.
    axis(1, at = c(1:6, 9:14), labels = rep(classes, 2),
         tick = FALSE, cex.axis = 0.75)
    axis(1, at = c(3.5, 11.5), labels = expression(lambda == 3, lambda == 30),
         tick = FALSE, line = 2)
    abline(v = 7.5, col = "#DDDDDD", lty = 3)
    box()
  }
}

dir.create("Figures", showWarnings = FALSE)
pdf("Figures/Figure2.pdf", width = 10.5, height = 4.8, bg = "white")
draw_figure2()
dev.off()
message("Saved Figures/Figure2.pdf")

# Display the same complete figure in RStudio.
if (interactive()) draw_figure2()
