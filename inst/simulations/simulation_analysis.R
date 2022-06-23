library(tidyverse)

files <- list.files(here("inst/simulations/raw_sims"), full.names = TRUE)

sims_readin <- map(files, read_csv, col_names = FALSE)
sims_df <- bind_rows(sims_readin) 

names(sims_df) <- c("pred", "se", "lb", "ub", "truetotal",
                    "parms.yo",
                    "parms.sigma_parsil_spat",
                    "parms.range",
                    "parms.sigma_nugget_spat",
                    "parms.sigma_parsil_time",
                    "parms.rho",
                    "parms.sigma_nugget_time",
                    "parms.sigma_nugget_spacetime",
                    "parms.sigma_parsil_spacetime",
                    "conf_ind")
sims_df                     
