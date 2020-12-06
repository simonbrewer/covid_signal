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

#multilevel modeling (test)
load("newdata.RData")

dat2$lpedest <- log(dat2$pedest+1)

dat2 <- dat2 %>% 
  filter(days_since_COVID > 20 & days_since_COVID < 40)

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
dat2 %>% 
  filter(days_since_COVID > 19 & days_since_COVID < 100) %>%
  tm_shape() + tm_symbols(col = myME, size = 0.25, palette = "PRGn", alpha = 0.5)

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

# var<-as.vector(c("E1", "E2", "E3", "E4", "E5"))
# 
# system.time(
#   mlm.sf <- forward.lmer(mlm,var)
# )

threshold_AIC <- 5
model_basis <- mlm
model_null <- update(model_basis, . ~ -. + 1 + (1|SIGNAL))
target_AIC <- AIC(model_basis)
improve_AIC <- AIC(model_null) - target_AIC

Evar <- as.vector(paste0("E", seq(1,5)))

## Outer loop
## iterate over all candidate MEs - maybe replace with a condition?
while (improve_AIC > threshold_AIC)  {
  # Iteratively updating the model with addition of one block of variable(s)
  # Also: extracting the loglikelihood of each estimated model
  
  candidate_AIC <- target_AIC
  for(j in 1:length(Evar)) ## Make this parallel
  {
    #print(paste(i,j, "AIC", target_AIC))
    candidate_model <- update(model_basis, as.formula(paste(". ~ . + ", Evar[j])))
    #print(paste(candidate_AIC, AIC(candidate_model)))
    if (AIC(candidate_model) < candidate_AIC) {
      candidate_AIC <- AIC(candidate_model)
      candidate_E <- Evar[j]
      candidate_j <- j
    }
    
  }
  improve_AIC <- target_AIC - candidate_AIC
  if (improve_AIC > threshold_AIC) {
    model_basis <- update(model_basis, as.formula(paste(". ~ . + ", candidate_E)))
    target_AIC <- AIC(model_basis)
    Evar <- Evar[-candidate_j]
    
    print(paste("Selected", candidate_E, "AIC", candidate_AIC, "Improvement", improve_AIC))
    
  }
}



stop()
## Code modified from Rense Nieuwenhuis `forward.lmer()` 

threshold_AIC <- 5
model_basis <- mlm
target_AIC <- AIC(model_basis)
models <- list()

Evar <- as.vector(paste0("E", seq(1,10)))
niter <- 5
## Outer loop
## iterate over all candidate MEs - maybe replace with a condition?
for (i in 1:niter) {
  # Iteratively updating the model with addition of one block of variable(s)
  # Also: extracting the loglikelihood of each estimated model
  for(j in 1:length(Evar))
  {
    print(paste(i,j))
    models[[j]] <- update(model_basis, as.formula(paste(". ~ . + ", Evar[j])))
  }
  
  step_AIC <- unlist(lapply(models, AIC))
  diff_AIC <- target_AIC - min(step_AIC) 
  sel_E <- which.min(step_AIC)
  
  if (diff_AIC > threshold_AIC) {
    model_basis <- models[[sel_E]]
    target_AIC <- AIC(model_basis)
    
    Evar <- Evar[-sel_E]
  } else {
    stop()
  }
  stop()
  
}
