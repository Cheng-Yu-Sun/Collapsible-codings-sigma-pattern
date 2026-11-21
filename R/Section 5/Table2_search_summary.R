#####################################################################
# Table 2: SOA search settings, elapsed times, and selected scores.
# Run All reads the saved search results; it does not rerun the searches.
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


selected_z4 <- read.csv("detailed_outputs/best_array_search/SOA_z4_selected_scores.csv")
selected_z6 <- read.csv("detailed_outputs/best_array_search/SOA_z6_selected_scores.csv")
candidates_z4 <- read.csv("detailed_outputs/best_array_search/SOA_z4_candidate_scores.csv")
candidates_z6 <- read.csv("detailed_outputs/best_array_search/SOA_z6_candidate_scores.csv")

criteria <- c("Sigma", "SP")
selected_z4 <- selected_z4[match(criteria, selected_z4$criterion), ]
selected_z6 <- selected_z6[match(criteria, selected_z6$criterion), ]

considered <- c(rep(nrow(candidates_z4), 2),
                sum(candidates_z6$criterion == "Sigma"),
                sum(candidates_z6$criterion == "SP"))
retained <- c(nrow(candidates_z4),
              sum(abs(candidates_z4$nu33 - min(candidates_z4$nu33)) < 1e-10),
              sum(candidates_z6$criterion == "Sigma" & candidates_z6$starts > 0),
              sum(candidates_z6$criterion == "SP" & candidates_z6$starts > 0))
starts <- c(selected_z4$starts_per_A, selected_z6$starts_per_P)

Table2 <- data.frame(
  z = c(4, 4, 6, 6), criterion = rep(criteria, 2),
  considered = considered, retained = retained,
  starts_per_A = starts, search_runs = retained * starts,
  seconds = round(c(selected_z4$seconds, selected_z6$seconds), 2),
  sigma33 = round(c(selected_z4$sigma33, selected_z6$sigma33), 2),
  sigma42 = round(c(selected_z4$sigma42, selected_z6$sigma42), 2)
)

dir.create("Tables", showWarnings = FALSE)
write.csv(Table2, "Tables/Table2.csv", row.names = FALSE)
print(Table2, row.names = FALSE)
message("Saved Tables/Table2.csv")
