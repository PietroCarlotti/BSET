rm(list = ls())
library(BSET)
source("Parast_et_al_2024_simulations.R")
graphics.off()

args <- commandArgs(trailingOnly = TRUE)
batch_id <- as.numeric(args[1])

batch_size <- 50

start_time <- Sys.time()

results <- Parast_et_al_2024_simulations(
  seed = 123 + batch_id,
  n_simulations = batch_size,
  parallel = FALSE
)

end_time <- Sys.time()

cat("Batch:", batch_id, "\n")
cat("Total simulation time:", end_time - start_time, "\n")

write.csv(
  results,
  paste0("Parast_et_al_2024_batch_", batch_id, ".csv"),
  row.names = FALSE
)
