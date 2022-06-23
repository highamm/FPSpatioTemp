## script for cleaning the raw data sets

library(tidyverse)
library(here)
# moose_upto2019 <- read_csv(here("inst", "rawdata",
#                                 "All_20E_Surveys_1998_2019_with_covariate_data_final.csv"))
# ## issues reading in searchmin, datecounted, and comments
# ## none of these are used for now so that shouldn't be an issue
# 
# moose_2020 <- read_csv(here("inst", "rawdata",
#                             "2020_Taylor_Corridor_survey.csv"))
all_moose <- bind_rows(moose_upto2019, moose_2020)
raw_moose <- all_moose |> select(-c(YBULL_SF, YBULL_GTSF, BULL_30_40, 
                        BULL_41_50, BULL_GT_50, COW_W_0, COW_W_1,
                        COW_W_2, COW_W_3, CALF, UNKNOWN, TotalAdults,
                        TotalCows, TotalBulls, TotalCalves)) |>
  ungroup()

##save(raw_moose, file = "data/raw_moose.rda")

load(file = "data/raw_moose.rda")
all_moose <- raw_moose |>
  mutate(stratfact = factor(str_to_upper(Stratname)))

## change some 0's that should be NA to NA and vice versa
all_moose <- all_moose |>
  mutate(totalmoosena = if_else(Counted == 1,
                                true = totalmoose, false = NA_real_)) |>
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

## A super-complete data set would have all sites across all years
## with an indicator for whether the site was in the sampling frame
## for that year.

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

moose_full |> filter(Surveyyear == 2014, Stratname == "HIGH") |> View()
## add coordinates
moose_full <- moose_full |> group_by(ID) |> 
  fill(xcoords, .direction = "downup") |>
  fill(ycoords, .direction = "downup") |>
  ungroup()

## if using sites outside the sampling frame for that year, need a method
## to classify them into a strata (this may not come up though)
## 
## for now, give sites the same strata as the most recent previous year

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

moose_complete <- bind_rows(moose_full, moose_2016)

## output the full data set

## save(file = "data/moose_complete.rda", moose_complete)


load(file = "data/moose_complete.rda")
test <- moose_complete |> group_by(ID) |>
  count()
table(test$n)

moose_test <- moose_complete |>
  dplyr::filter(Surveyyear >= 2014 & Surveyyear <= 2018) |>
  dplyr::filter(stratfact == "HIGH")
test2 <- moose_test |> group_by(ID) |> 
  count()
table(test2$n)

## make some useful subsets

## 2014 - 2018
moose_14_18 <- moose_complete |>
  filter(Surveyyear >= 2014 & Surveyyear <= 2018)

## only keep rows where the site was in the sampling frame at least
## once across the years.
moose_14_18 <- moose_14_18 |> group_by(ID) |>
  filter((sum(samp_frame) != 0))

load(file = "data/moose_14_18.rda")
moose_14_18
## 2014 - 2020

moose_14_18 <- moose_14_18 |> ungroup()
## save(moose_14_18, file = "data/moose_14_18.rda")


moose_14_20 <- moose_complete |>
  filter(Surveyyear >= 2014 & Surveyyear <= 2020)
moose_14_20 <- moose_14_20 |> group_by(ID) |>
  filter((sum(samp_frame) != 0))

moose_14_20 <- moose_14_20 |> ungroup()
## save(moose_14_20, file = "data/moose_14_20.rda")


## 2004 - 2012
moose_04_12 <- moose_complete |>
  filter(Surveyyear >= 2004 & Surveyyear <= 2012)
moose_04_12 <- moose_04_12 |> group_by(ID) |>
  filter((sum(samp_frame) != 0))

moose_04_12 <- moose_04_12 |> ungroup()
## save(moose_04_12, file = "data/moose_04_12.rda")
