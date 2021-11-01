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

rm(list = ls())

## Load data
dat2 <- readRDS("./data/data_for_COVIDmodel.rds")
head(dat2)
dim(dat2)
summary(dat2$SIGNAL)
summary(dat2$sigid)

## Subset to create the graph network
dat3 <- distinct(dat2, SIGNAL, .keep_all = TRUE)
coords <- st_coordinates(dat3)

## Gabriel graph
nb <- graph2nb(relativeneigh(coords), sym = TRUE, row.names = dat3$SIGNAL)

## Spatial weight matrix
sig.lw = nb2listw(nb, style = "B")

## Generate set of Eigenvectors
n <- length(nb)
## Weight matrix (as matrix)
C <- listw2mat(sig.lw)
## M - diagonal
M <- diag(1,n)-1/n
## Griffith's MCM matrix
MCM <- M%*%C%*%M
## Eigen decomposition
eig <- eigen(MCM, symmetric = TRUE)
## Extract vectors and sub out positive vectors
E <- eig$vector
E <- E[, which(eig$values > 0.01)]
## Make up as data frame for modeling
E.df <- as.data.frame(E)
colnames(E.df) <- paste0("E", 1:ncol(E))
E.df$SIGNAL <- attr(nb, "region.id")

## Merge full set of eigenvectors back with data
dat2 <- merge(dat2, E.df, by = "SIGNAL")

## ----------------------------------------------------------------------------
## Model

## Set up model formula
myform1 <- "lpedest ~ phase*(ln_popden_000_qtmi+ln_empden_000_qtmi + 
            per_res_qtmi + per_com_qtmi + per4wy_qtmi + stops_qtmi+ schools_qtmi+ worship_qtmi +
             ln_park_acre_qtmi + income_000_qtmi +  ln_hhsize_qtmi + avgveh_qtmi + major_road + SLC) + weekend + TAVG + TAVG_90 + PRCP + snow_dummy + recall"
myform1_re <- "lpedest ~ phase*(ln_popden_000_qtmi+ln_empden_000_qtmi + 
             per_com_qtmi + per4wy_qtmi + stops_qtmi+ schools_qtmi+ worship_qtmi +
             ln_park_acre_qtmi + income_000_qtmi +  ln_hhsize_qtmi +  major_road + SLC) + weekend + 
TAVG + TAVG_90 + PRCP + snow_dummy + recall"

## Simple linear model for VIF check
lm_re<-lm(myform1_re,dat2);summary(lm_re)
car::vif(lm_re)

## MLM formula
myform2 <- formula(paste0(myform1_re, "+ (1|SIGNAL)"))
mlm<-lmer(myform2, dat2)
# AIC(mlm)
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

threshold_p <- 0.07
model_basis <- mlm
dat3$resid_lev2 <- ranef(model_basis)[[1]][,1]
target_MI <- moran.test(dat3$resid_lev2, sig.lw)
target_z <- target_MI$statistic
target_p <- target_MI$p.value

Evar <- as.vector(paste0("E", seq(1,200))) ## Increase
# 
cl <- parallel::makeCluster(8)
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

summary(model_basis)

stop()


#############
#### COUNT model ########
#####################
start_time <- Sys.time()
myform_count <- "lpedest ~ Count*(ln_popden_000_qtmi+ln_empden_000_qtmi + 
             per_com_qtmi + per4wy_qtmi + stops_qtmi+ schools_qtmi+ worship_qtmi +
             ln_park_acre_qtmi + income_000_qtmi +  ln_hhsize_qtmi +  major_road + SLC) + weekend + 
TAVG + TAVG_90 + PRCP + snow_dummy + recall - Count"

# summary(dat2[,1:30])
lm<-lm(myform_count,dat2)
car::vif(lm)
# summary(lm)

myform2 <- formula(paste0(myform_count, "+ (1|SIGNAL)"))
system.time(mlm<-lmer(myform2, dat2))
# AIC(mlm)
# summary(mlm)

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

Evar <- as.vector(paste0("E", seq(1,200))) ## Increase
# 
# cl <- parallel::makeCluster(8)
# doParallel::registerDoParallel(cl)

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

# parallel::stopCluster(cl)
end_time <- Sys.time()

end_time - start_time
summary(model_basis)



