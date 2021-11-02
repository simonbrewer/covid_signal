# script to analyze 7-day averages of 2019 vs 2020

##### load and process data

# load libraries
library(sf)
library(tidyr)
library(dplyr)
library(INLA)
library(lmerTest)
library(spdep)
library(ggregplot)

#multilevel modeling (test)
load("newdata.RData")

dat2$lpedest <- log(dat2$pedest+1)

myform1 <- lpedest ~ ln_popden_000_hami + ln_empden_000_qtmi + 
  ln_hhsize_qtmi+ income_000_hami + avgveh_hami + 
  per_res_qtmi + per_com_qtmi + intden_hami + per4wy_hami + 
  schools_qtmi + worship_hami + stops_qtmi + 
  ln_park_acre_hami + major_road + 
  SLC + TAVG + PRCP + snow_dummy + 
  Count + phase + recall + days_since_COVID

fit <- lm(myform1, dat2)
summary(fit)

## Cross sectional

# dat4 <- dat2 %>%
#   filter(TIME2 == "2020-04-01")

#### --------------------------------------------------------------------- ####
dat4 <- dat2 %>%
  filter(TIME2 < "2020-04-01")

mod1 <- inla(myform1,
             data = dat4,
             control.inla = list(strategy = "gaussian", int.strategy = "eb"),
             control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
             control.predictor = list(compute = TRUE))
summary(mod1)

Efxplot(mod1)

#### --------------------------------------------------------------------- ####
## IID in SIGNAL
myform2 <- lpedest ~ ln_popden_000_hami + ln_empden_000_qtmi + 
  ln_hhsize_qtmi+ income_000_hami + avgveh_hami + 
  per_res_qtmi + per_com_qtmi + intden_hami + per4wy_hami + 
  schools_qtmi + worship_hami + stops_qtmi + 
  ln_park_acre_hami + major_road + 
  SLC + TAVG + PRCP + snow_dummy + 
  Count + phase + recall + days_since_COVID + 
  f(SIGNAL, model = "iid", hyper = prec.prior)

# myform2 <- lpedest ~ 1 + 
#   f(SIGNAL, model = "iid", hyper = prec.prior)

prec.prior <- list(prec = list(param = c(0.001, 0.001)))
system.time(
  mod2 <- inla(myform2,
               data = dat2,
               control.inla = list(strategy = "gaussian", int.strategy = "eb"),
               control.compute = list(dic = TRUE, waic = TRUE),
               control.predictor = list(compute = TRUE))
)
summary(mod2)
Efxplot(mod2)

#### --------------------------------------------------------------------- ####
## Graph
signal.W <- nb2mat(nb)

## IID in SIGNAL
myform3 <- lpedest ~ ln_popden_000_hami + ln_empden_000_qtmi + 
  ln_hhsize_qtmi+ income_000_hami + avgveh_hami + 
  per_res_qtmi + per_com_qtmi + intden_hami + per4wy_hami + 
  schools_qtmi + worship_hami + stops_qtmi + 
  ln_park_acre_hami + major_road + 
  SLC + TAVG + PRCP + snow_dummy + 
  Count + phase + recall + days_since_COVID + 
  f(SIGNAL2, model = "bym2", hyper = bym2.prior, graph = signal.W)

# myform2 <- lpedest ~ 1 + 
#   f(SIGNAL, model = "iid", hyper = prec.prior)

bym2.prior <- list(
  prec = list(
    prior = "pc.prec",
    param = c(0.5 / 0.31, 0.01)),
  phi = list(
    prior = "pc",
    param = c(0.5, 2 / 3))
)

system.time(
  mod3 <- inla(myform3,
               data = dat2,
               control.inla = list(strategy = "gaussian", int.strategy = "eb"),
               control.compute = list(dic = TRUE, waic = TRUE),
               control.predictor = list(compute = TRUE))
)
summary(mod3)
Efxplot(mod3)

#### --------------------------------------------------------------------- ####

## IID in SIGNAL
myform4 <- lpedest ~ ln_popden_000_hami + ln_empden_000_qtmi + 
  ln_hhsize_qtmi+ income_000_hami + avgveh_hami + 
  per_res_qtmi + per_com_qtmi + intden_hami + per4wy_hami + 
  schools_qtmi + worship_hami + stops_qtmi + 
  ln_park_acre_hami + major_road + 
  SLC + TAVG + PRCP + snow_dummy + 
  Count + phase + recall + 
  f(SIGNAL2, model = "bym2", hyper = bym2.prior, graph = signal.W) +
  f(days_since_COVID, model = "rw1", hyper = rw1_prior)

# myform2 <- lpedest ~ 1 + 
#   f(SIGNAL, model = "iid", hyper = prec.prior)

bym2.prior <- list(
  prec = list(
    prior = "pc.prec",
    param = c(0.5 / 0.31, 0.01)),
  phi = list(
    prior = "pc",
    param = c(0.5, 2 / 3))
)

rw1_prior <- list(prec = list(param = c(0.001, 0.001)))

system.time(
  mod4 <- inla(myform4,
               data = dat2,
               control.inla = list(strategy = "gaussian", int.strategy = "eb"),
               control.compute = list(dic = TRUE, waic = TRUE),
               control.predictor = list(compute = TRUE))
)
summary(mod4)
Efxplot(mod4)