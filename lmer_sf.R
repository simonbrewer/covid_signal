lmer_sf <- function(formula, data, idcol,
                    nb, pthresh = 0.05, 
                    evlimit = NULL, parallel = FALSE,
                    verbose = TRUE) {
  require(lme4)
  require(spdep)
  require(foreach)
  
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
  if (is.null(evlimit)) {
    nE <- ncol(E)
  } else {
    nE <- evlimit
  }
  ## Make up as data frame for modeling
  E.df <- as.data.frame(E[, 1:nE])
  colnames(E.df) <- paste0("E", 1:nE)
  E.df[idcol] <- attr(nb, "region.id")
  ## Create list of EVs for selection
  Evar <- as.vector(paste0("E", seq(1,nE))) 
  print(Evar)
  
  ## Merge full set of eigenvectors back with data
  data <- merge(data, E.df, by = idcol)
  
  ## Model basis (non-filtered)
  model_basis <- lmer(formula, data)
  
  ## Moran's test on base model
  ## Use ranef to get second level residuals
  resid_lev2 <- ranef(model_basis)[[1]][,1]
  
  ## Estimate Moran's I and store values + p-val
  target_MI <- moran.test(resid_lev2, sig.lw)
  target_z <- target_MI$statistic
  target_p <- target_MI$p.value
  
  if (verbose)
    print(paste("Basis:", "z:", target_z, "p-val:", target_p))
  
  ## Set up for parallel run (adjust core number as needed)
  if (parallel) {
    cl <- parallel::makeCluster(8)
    doParallel::registerDoParallel(cl)
  }
  
  ## Outer loop
  while ((target_p < pthresh) & (length(Evar) > 0))  {
    ## Inner loop
    candidate_p <- target_p 
    if (parallel) {
      list_p <- foreach(j=1:length(Evar), 
                        .combine = 'c',
                        .packages = c("lme4", "spdep"),
                        .export = "dat2") %dopar% ## Parallel
        {
          #print(paste(i,j, "AIC", target_AIC))
          candidate_model <- update(model_basis, as.formula(paste(". ~ . + ", Evar[j])))
          resid_lev2 <- ranef(candidate_model)[[1]][,1]
          moran.test(resid_lev2, sig.lw)$p.value
        }
    } else {
      list_p <- foreach(j=1:length(Evar), 
                        .combine = 'c',
                        .packages = c("lme4", "spdep"),
                        .export = "dat2") %do% ## Sequential
        {
          #print(paste(i,j, "AIC", target_AIC))
          candidate_model <- update(model_basis, as.formula(paste(". ~ . + ", Evar[j])))
          resid_lev2 <- ranef(candidate_model)[[1]][,1]
          moran.test(resid_lev2, sig.lw)$p.value
        }
    }
    candidate_j <- which.max(list_p)
    candidate_E <- Evar[candidate_j]
    candidate_p <- max(list_p)
    
    ## Update model basis with selected EV
    model_basis <- update(model_basis, as.formula(paste(". ~ . + ", candidate_E)))
    ## Get updated residuals
    resid_lev2 <- ranef(model_basis)[[1]][,1]
    ## Calculate Moran's I on new residuals
    target_MI <- moran.test(resid_lev2, sig.lw)
    target_z <- target_MI$statistic
    target_p <- target_MI$p.value
    ## Remove selected EV from list
    Evar <- Evar[-candidate_j]
    
    ## Some output
    if (verbose)
      print(paste("Selected:", candidate_E, "z:", target_z, "p-val:", target_p))
    
  }
  
  if (parallel) 
    parallel::stopCluster(cl)

  return(list(model_basis = model_basis,
              mi = target_MI))
  
}