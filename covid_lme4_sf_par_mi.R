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
library(foreach)

#multilevel modeling (test)
load("newdata.RData")

dat2$lpedest <- log(dat2$pedest+1)

## Subset for testing
dat2 <- dat2 %>%
  filter(days_since_COVID > 20 & days_since_COVID < 40)
dat3 <- dat2 %>% distinct(SIGNAL)
coords <- st_coordinates(dat3)
nb <- graph2nb(relativeneigh(coords), sym = TRUE, row.names = dat3$SIGNAL)

myform1 <- lpedest ~ ln_popden_000_hami + ln_empden_000_qtmi + 
  ln_hhsize_qtmi+ income_000_hami + avgveh_hami + 
  per_res_qtmi + per_com_qtmi + intden_hami + per4wy_hami + 
  schools_qtmi + worship_hami + stops_qtmi + 
  ln_park_acre_hami + major_road + 
  SLC + TAVG + PRCP + snow_dummy + 
  Count + phase + recall + days_since_COVID

sig.lw = nb2listw(nb, style = "B")

# generates a series of hypothetical potential spatial patterns 
n <- length(nb)
C <- listw2mat(sig.lw)
M <- diag(1,n)-1/n
MCM <- M%*%C%*%M
eig <- eigen(MCM, symmetric = TRUE)
E <- eig$vector
E <- E[, which(eig$values > 0.01)]
E.df <- as.data.frame(E)
colnames(E.df) <- paste0("E", 1:ncol(E))
E.df$SIGNAL <- attr(nb, "region.id")

dat2 <- merge(dat2, E.df, by = "SIGNAL")

## Map to show MEs
myME <- "E10"
mydsc <- 20
# dat2 %>% 
#   filter(days_since_COVID > 19 & days_since_COVID < 20) %>%
#   tm_shape() + tm_symbols(col = myME, size = 0.25, palette = "PRGn", alpha = 0.5)

#multilevel modeling (test)
myform1 <- "log(pedest+1) ~ ln_popden_000_hami+ln_empden_000_qtmi + ln_hhsize_qtmi+ income_000_hami + avgveh_hami + 
            per_res_qtmi + per_com_qtmi + intden_hami + per4wy_hami + schools_qtmi + worship_hami + stops_qtmi +
             ln_park_acre_hami + major_road + SLC + TAVG + PRCP + snow_dummy + Count + phase + recall + days_since_COVID"
myform1 <- "log(pedest+1) ~ ln_popden_000_hami+ln_empden_000_qtmi + ln_hhsize_qtmi+ income_000_hami + avgveh_hami + 
            per_res_qtmi + per_com_qtmi + intden_hami + per4wy_hami + schools_qtmi + worship_hami + stops_qtmi +
             ln_park_acre_hami + major_road + SLC + TAVG + PRCP + snow_dummy + Count + recall + days_since_COVID"
lm<-lm(myform1,dat2);summary(lm)

myform2 <- formula(paste0(myform1, "+ (1|SIGNAL)"))
system.time(mlm<-lmer(myform2, dat2))
AIC(mlm)
summary(mlm)

## Moran's test on base model
## ranef gets second level residuals
dat3$resid_lev2 <- ranef(mlm)[[1]][,1]
tm_shape(dat3) + tm_symbols(col = "resid_lev2", style = "fisher")

target_MI <- moran.test(dat3$resid_lev2, sig.lw)
target_z <- target_MI$statistic
target_p <- target_MI$p.value
# var<-as.vector(c("E1", "E2", "E3", "E4", "E5"))
# 
# system.time(
#   mlm.sf <- forward.lmer(mlm,var)
# )

threshold_p <- 0.05
model_basis <- mlm
dat3$resid_lev2 <- ranef(model_basis)[[1]][,1]
target_MI <- moran.test(dat3$resid_lev2, sig.lw)
target_z <- target_MI$statistic
target_p <- target_MI$p.value

Evar <- as.vector(paste0("E", seq(1,100)))

cl <- parallel::makeCluster(4)
doParallel::registerDoParallel(cl)

## Outer loop
## iterate over all candidate MEs - maybe replace with a condition?
while ((target_p < threshold_p) & (length(Evar) > 0))  {
  # Iteratively updating the model with addition of one block of variable(s)
  # Also: extracting the loglikelihood of each estimated model
  
  candidate_p <- target_p
  list_p <- foreach(j=1:length(Evar), 
                      .combine = 'c',
                      .packages = c("lme4", "spdep"),
                      .export = "dat2") %dopar% ## Make this parallel
    {
      #print(paste(i,j, "AIC", target_AIC))
      candidate_model <- update(model_basis, as.formula(paste(". ~ . + ", Evar[j])))
      dat3$resid_lev2 <- ranef(candidate_model)[[1]][,1]
      moran.test(dat3$resid_lev2, sig.lw)$p.value
    }
  candidate_j <- which.max(list_p)
  candidate_E <- Evar[candidate_j]
  candidate_p <- max(list_p)
  
  model_basis <- update(model_basis, as.formula(paste(". ~ . + ", candidate_E)))
  dat3$resid_lev2 <- ranef(model_basis)[[1]][,1]
  target_MI <- moran.test(dat3$resid_lev2, sig.lw)
  target_z <- target_MI$statistic
  target_p <- target_MI$p.value
  Evar <- Evar[-candidate_j]
  
  print(paste("Selected:", candidate_E, "z:", target_z, "p-val:", target_p))
    
}

parallel::stopCluster(cl)
