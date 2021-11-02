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

#### --------------------------------------------------------------------- ####
## Graph
signal.W <- nb2mat(nb)

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

library(ggpubr)
time.re <- mod4$summary.random$days_since_COVID
ggline(time.re, x = "ID", y = "0.5quant", xlab = "Year",
       numeric.x.axis = TRUE)

p1 <- ggplot(time.re, aes(x = ID)) + 
  geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`), fill = 'gray70') +
  geom_line(aes(y = mean)) + 
  scale_x_continuous("Days since COVID") + 
  scale_y_continuous("Time random effect") + theme_bw()

pdf("covid_inla_days_since_coivd_re.pdf")
print(p1)
dev.off()

mod4$marginals.random$SIGNAL2
dat2$yhat <- exp(mod4$summary.fitted.values$mean)
dat2$cilo <- exp(mod4$summary.fitted.values$`0.025quant`)
dat2$cihi <- exp(mod4$summary.fitted.values$`0.975quant`)

library(tmap)
dat_map <- dat2 %>% 
  filter(days_since_COVID == 0)

tm_shape(dat_map) + tm_symbols("yhat")
save(mod4, file = "covid_inla_mod4.RData")


mysignal = 1001
dat_ts <- dat2 %>% 
  filter(SIGNAL == mysignal)

p1 <- ggplot(dat_ts, aes(x = days_since_COVID)) + 
  geom_ribbon(aes(ymin = cilo, ymax = cihi), fill = 'gray70') +
  geom_line(aes(y = yhat)) + 
  scale_x_continuous("Days since COVID") + 
  scale_y_continuous("yhat") + theme_bw() + ggtitle(paste("SIGNAL", mysignal))
print(p1)