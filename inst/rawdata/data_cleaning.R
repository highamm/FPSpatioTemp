## script for cleaning the raw data sets

library(tidyverse)
library(here)

load(file = "data/raw_moose.rda")
all_moose <- raw_moose |>
  mutate(stratfact = factor(str_to_upper(Stratname)))

all_moose <- all_moose |>
  mutate(totalmoosena = if_else(Counted == 1,
                                true = totalmoose,
                                false = NA_real_)) |>
  mutate(prop.cov.1 =
           if_else(is.na(prop.cov.1), true = 0, false = prop.cov.1)) |>
  mutate(prop.cov.6_7 =
           if_else(is.na(prop.cov.6_7), true = 0, false = prop.cov.6_7)) |>
  dplyr::select(totalmoosena, everything())

all_moose |> filter(Surveyyear == 2020) |> select(stratfact)

## change lat/lon to TM
library(sptotal)
xcoords <- LLtoTM(cm = mean(all_moose$centrlon), lat = all_moose$centrlat, lon = all_moose$centrlon)$xy[, "x"]
ycoords <- LLtoTM(cm = mean(all_moose$centrlon), lat = all_moose$centrlat, lon = all_moose$centrlon)$xy[, "y"]
coords <- data.frame(xcoords, ycoords)

moose_tm <- as.data.frame(bind_cols(all_moose, coords))

## all observations in the data set are in the sampling frame right now
moose_tm <- moose_tm |> mutate(samp_frame = 1) |>
  as_tibble()

## drop extra columns for simplicity
moose_small <- moose_tm |>
  select(totalmoosena, SurveyID, SurveyName, Surveyyear, ID, 
         Stratname, Counted, AreaMi, ELEV_MEAN, stratfact,
         xcoords, ycoords, samp_frame)
moose_full <- moose_small |> complete(Surveyyear, nesting(ID),
                                       fill = list(samp_frame = 0))

## add coordinates
moose_full <- moose_full |> group_by(ID) |> 
  fill(xcoords, .direction = "downup") |>
  fill(ycoords, .direction = "downup") |>
  ungroup()

moose_full <- moose_full |> arrange(Surveyyear) |>
  group_by(ID) |> 
  fill(stratfact, .direction = "down") |>
  ungroup()

## finally, add 2016 NAs in (no survey was done in 2016 because of
## lack of snowfall)
## assume 2016 has same sampling frame as 2015

moose_2016 <- moose_full |> filter(Surveyyear == 2015)
moose_2016$totalmoosena <- NA
moose_2016$Surveyyear <- 2016

## moose_complete <- bind_rows(moose_full, moose_2016)

## output the full data set



load(file = "data/moose_complete.rda")
test <- moose_complete |> group_by(ID) |>
  count()
table(test$n)

## add 2021 data to moose_complete
library(tidyverse)
moose_21 <- read_csv("inst/rawdata/moose_21.csv")

moose_21$totalmoose
moose_21$Counted
library(sptotal)

moose_21 <- moose_21 |> mutate(totalmoosena = if_else(Counted == TRUE,
                                           true = totalmoose,
                                           false = NA_real_)) |>
  select(totalmoosena, everything()) |>
  mutate(ELEV_MEAN = NA) |>
  mutate(Stratname = str_to_upper(Stratname)) |>
  mutate(stratfact = factor(Stratname)) |>
  mutate(xcoords = LLtoTM(cm = mean(moose_21$centrlon), lat = moose_21$centrlat, lon = moose_21$centrlon)$xy[, "x"]) |>
  mutate(ycoords = LLtoTM(cm = mean(moose_21$centrlon), lat = moose_21$centrlat, lon = moose_21$centrlon)$xy[, "y"])

moose_21 <- moose_21 |> mutate(samp_frame = 1) |> select(Surveyyear, ID, totalmoosena, SurveyID, SurveyName,
                   Stratname, Counted, AreaMi,
                   ELEV_MEAN, stratfact, xcoords, ycoords, samp_frame) 
## next step: expand to include all coordinates but assign samp_frame to be 1 for those in the actual sampling frame for 2021

moose_21_outframe <- moose_complete |> filter(Surveyyear == 2020) |> anti_join(moose_21,
                                                          by = "ID") |>
  mutate(Surveyyear = 2021)

moose_21_comp <- bind_rows(moose_21, moose_21_outframe)
moose_through_21 <- bind_rows(moose_complete, moose_21_comp)


moose_22 <- read_csv("inst/rawdata/2022_20E_GSPE_survey_final_results.csv")
moose_22$totalmoose
moose_22$Counted

moose_22 <- moose_22 |> mutate(totalmoosena = if_else(Counted == TRUE,
                                                      true = totalmoose,
                                                      false = NA_real_)) |>
  select(totalmoosena, everything()) |>
  mutate(ELEV_MEAN = NA) |>
  mutate(Stratname = str_to_upper(Stratname)) |>
  mutate(stratfact = factor(Stratname)) |>
  mutate(xcoords = LLtoTM(cm = mean(moose_22$centrlon), lat = moose_22$centrlat, lon = moose_22$centrlon)$xy[, "x"]) |>
  mutate(ycoords = LLtoTM(cm = mean(moose_22$centrlon), lat = moose_22$centrlat, lon = moose_22$centrlon)$xy[, "y"])


moose_22 <- moose_22 |> mutate(samp_frame = 1) |> select(Surveyyear, ID, totalmoosena, SurveyID, SurveyName,
                                                         Stratname, Counted, AreaMi,
                                                         ELEV_MEAN, stratfact, xcoords, ycoords, samp_frame)

moose_22_outframe <- moose_complete |> filter(Surveyyear == 2020) |> anti_join(moose_22,
                                                                               by = "ID") |>
  mutate(Surveyyear = 2022)

moose_22_comp <- bind_rows(moose_22, moose_22_outframe)
moose_through_22 <- bind_rows(moose_through_21, moose_22_comp)
moose_complete <- moose_through_22
# save(moose_complete, file = "data/moose_complete.rda")
# load(file = "data/moose_complete.rda")
# moose_complete |> filter(Surveyyear == 2022)
moose_complete |> filter(Surveyyear == 2022)
slmfit(totalmoosena ~ stratfact, data = moose_complete |> filter(Surveyyear == 2022 & samp_frame == 1),
       xcoordcol = "xcoords", ycoordcol = "ycoords") |>
  predict()
