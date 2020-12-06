# script to analyze 7-day averages of 2019 vs 2020

##### load and process data

# load libraries
library("sf")
library("tidyr")
library("dplyr")
library(lme4)
library("lmerTest")
library(spdep)

#multilevel modeling (test)
dat2 <- readRDS("data_for_COVIDmodel.rds")
names(dat2)
myform1 <- "log(pedest+1) ~ ln_popden_000_hami+ln_empden_000_qtmi + ln_hhsize_qtmi+ income_000_hami + avgveh_hami + 
            per_res_qtmi + per_com_qtmi + intden_hami + per4wy_hami + schools_qtmi + worship_hami + stops_qtmi +
             ln_park_acre_hami + major_road + SLC + TAVG + PRCP + snow_dummy + Count + phase + recall + days_since_COVID"
lm<-lm(myform1,dat2);summary(lm)

myform2 <- formula(paste0(myform1, "+ (1|SIGNAL)"))
system.time(mlm<-lmer(myform2,dat2))
summary(mlm)


#multilevel spatial filtering

#converting a data frame to a sf object
dat2 <- as_Spatial(dat2)
class(dat2)

#build spatial weight matrix #http://rstudio-pubs-static.s3.amazonaws.com/5009_aa49bd54e40f41f19352af3d9e00b036.html
coords<-coordinates(dat2)
nb = dnearneigh(coords, d1 = 0, d2 = 0.75 * 0.7654856)
##the maximum distance between two neighboring intersections is 0.7654856
#0.7654856*0.75 = 0.57

lw = nb2listw(nb, style='W',zero.policy=TRUE)

#inversed distance
dsts = nbdists(nb, coords); invd = function(x) {1/(x/1000)}; idw = lapply(dsts, invd)
lw_idwB = nb2listw(nb, glist = idw, style = "B")


# end