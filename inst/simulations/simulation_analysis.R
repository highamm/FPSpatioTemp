library(tidyverse)

files <- list.files("inst/simulations/raw_sims", full.names = TRUE)

sims_df <- map(files, read_csv, col_names = FALSE)
bind_rows(sims_df) |> print(width = Inf) %>%
  summarise(sum(X15), mean(X15))
