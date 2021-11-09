## phase_test

# script to analyze 7-day averages of 2019 vs 2020

##### load and process data

# load libraries
library(sf)
library(tidyr)
library(dplyr)
library(lme4)
library(lmerTest)
library(spdep)
library(tmap)

# rm(list = ls())
source("lmer_sf.R")

## Load data
dat2 <- readRDS("./data/data_for_COVIDmodel.rds")

# ## Subset for uncomment for testing
# library(lubridate)
# dat2$TIME1 <- ymd(dat2$TIME1)
# dat2 <- dat2 %>%
#   filter(month(dat2$TIME1) %in% c(4, 5) & year(dat2$TIME1) == 2020)
# dat2$TAVG_90[sample(1:nrow(dat2), 100)] <- 1 ## Meaningless dummy
# dat2$snow_dummy[sample(1:nrow(dat2), 100)] <- 1 ## Meaningless dummy

## Subset to create the graph network
dat3 <- distinct(dat2, SIGNAL, .keep_all = TRUE)
coords <- st_coordinates(dat3)
## Gabriel graph
nb <- graph2nb(relativeneigh(coords), sym = TRUE, row.names = dat3$SIGNAL)

## Model formula
myform1_re <- "lpedest ~ phase*(ln_popden_000_qtmi+ln_empden_000_qtmi + 
             per_com_qtmi + per4wy_qtmi + stops_qtmi+ schools_qtmi+ worship_qtmi +
             ln_park_acre_qtmi + income_000_qtmi +  ln_hhsize_qtmi +  major_road + SLC) + 
             weekend + TAVG + TAVG_90 + PRCP + snow_dummy + recall + (1 | SIGNAL)"

## Run the spatial filter
final_model <- lmer_sf(myform1_re, dat2, 
                       idcol = "SIGNAL", 
                       nb = nb,
                       alpha = 0.05, ## Stopping threshold (MI p-val)
                       evlimit = 200, ## Limit number of EVs to test
                       parallel = FALSE, ## Change to use multiple cores
                       ncores = 1, ## Sets number of cores if parallel = T
                       verbose = TRUE)

save(final_model, file = "covid_sf_phase.RData")
