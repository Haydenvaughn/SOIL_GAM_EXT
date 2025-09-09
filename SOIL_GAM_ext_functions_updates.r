library(mboost)
library(mgcv)
library(tidyverse)
library(MASS)

# Boosting GAM and Extraction Function
boost_extract = function(iters, y, X){
  candidate_mat = matrix(0, nrow = iters,ncol = ncol(X))
  data_comb = data.frame(y = y, X)
  gam.fit = gamboost(y ~ ., 
                     data = data_comb,
                     baselearner = "bbs",
                     control = boost_control(mstop = iters))
  print(gam.fit$xselect())
  for (i in seq_along(gam.fit$xselect())) {
      candidate_mat[i:nrow(candidate_mat),gam.fit$xselect()[i]] = 1
  }
  return(candidate_mat)
}

# ARM weighting for GAM ext
arm_gam <- function(candidate_mat, L, y, X) {
  
  # Storage for final weights
  weights.mat <- matrix(0, nrow = L, ncol = nrow(candidate_mat))
  
  # Combine y and X into one data frame
  data_comb <- data.frame(y = y, X)
  #print(head(data_comb))
  
  #initial weights for candidate models
  psi <- 1/nrow(candidate_mat)
  
  # Loop over repetitions
  for (l in 1:L) {
    
    # Train/test split
    train.index <- sample(1:nrow(data_comb), trunc(nrow(data_comb)/2))
    test.index  <- setdiff(1:nrow(data_comb), train.index)
    test.y      <- y[test.index]
    # Initialize storage for this repetition
    count.vec   <- rep(0, nrow(candidate_mat))
    unique.row  <- rep(0, nrow(candidate_mat))
    sig.hat     <- rep(NA, nrow(candidate_mat))
    df.resid    <- rep(NA, nrow(candidate_mat))
    fitted.vals <- matrix(NA, nrow = nrow(candidate_mat), ncol = length(train.index))
    preds       <- matrix(NA, nrow = nrow(candidate_mat), ncol = length(test.index))
    edf         <- rep(NA, nrow(candidate_mat))
    C           <- rep(NA, nrow(candidate_mat))
    num.weight  <- rep(NA, nrow(candidate_mat))
    weight      <- rep(NA, nrow(candidate_mat))
    # Loop over each candidate row
    for (row_idx in 1:nrow(candidate_mat)) {
      
      # Check for duplicates among earlier rows
      if (row_idx > 1 && all(candidate_mat[row_idx, ] == candidate_mat[row_idx-1, ])) {
            sig.hat[row_idx]     <- sig.hat[row_idx-1]
            df.resid[row_idx]    <- df.resid[row_idx-1]
            fitted.vals[row_idx, ] <- fitted.vals[row_idx-1, ]
            preds[row_idx, ]     <-  preds[row_idx-1,]
            edf[row_idx]         <-  edf[row_idx-1]
            C[row_idx]           <-  C[row_idx-1]
      }
      
      else {
        # Extract the subset of variables for this candidate
        selected_vars <- which(candidate_mat[row_idx, ] == 1)
        data.sub  <- data_comb[, c(1, selected_vars + 1)] # +1 for y in col1
    
        
        train.sub <- data.sub[train.index, ]
        test.sub  <- data.sub[test.index, ]
        #print(head(test.sub))
        
          
        covariates <- setdiff(names(train.sub), "y")
        smooth_terms <- paste0("s(", covariates, ', bs = "ps")')
        formula_str  <- paste("y ~", paste(smooth_terms, collapse = " + "))
          
        gam_fit <- gam(as.formula(formula_str), data = train.sub)
        
        df.resid[row_idx]     <- gam_fit$df.residual 
        fitted.vals[row_idx, ] <- gam_fit$fitted.values
        sig.hat[row_idx] <-(sum((train.sub$y-fitted.vals[row_idx, ])^2))/
                            df.resid[row_idx]
        preds[row_idx,]  <- predict.gam(gam_fit, newdata = test.sub)
        edf[row_idx]     <- sum(gam_fit$edf)
        
        
      }
      C[row_idx]       <- edf[row_idx] * log(exp(1)*max(edf,na.rm = TRUE)/edf[row_idx]) + 
                            2 * log(edf[row_idx] + 2)
      
      num.weight[row_idx]  <- (exp(-psi * C[row_idx]) * 
                          sig.hat[row_idx]^(-nrow(candidate_mat)) * 
                          prod(exp( -sig.hat[row_idx]^(-2) * 
                          (test.y - preds[row_idx,])^2 / 2 )))
                          
    }
    # calculate weights and set them in lth vector in weights mat
    weights.mat[l, ] <- num.weight/sum(num.weight)
    
  }
  
  final_weights <- colMeans(weights.mat)
  return(final_weights)
  
}

# SOIL_gam
soil_gam = function(candidate_mat,arm_weights){
  soil.vals = rep(NA, ncol(candidate_mat))
  for (j in 1:ncol(candidate_mat)){
    soil.vals[j] = sum(candidate_mat[,j]*arm_weights)
  }
  return(soil.vals)
}


# Creating a combined function that executes all three functions
soil.gam = function(boosting_iters, y, X, num_splits) {
  candidate.mat = boost_extract(boosting_iters, y, X)
  arm.weights = arm_gam(candidate.mat, num_splits, y, X)
  soil.scores = soil_gam(candidate.mat,arm.weights)
  return(soil.scores)
}

# Revising arm_gam to use mboost instead of mgcv
arm_gam_mboost_old <- function(candidate_mat, L, y, X, nu, mstop) {
  
  # Storage for final weights
  weights.mat <- matrix(0, nrow = L, ncol = nrow(candidate_mat))
  
  # Combine y and X into one data frame
  data_comb <- data.frame(y = y, X)
  #print(head(data_comb))
  
  #initial weights for candidate models
  psi <- 1/nrow(candidate_mat)
  
  # Loop over repetitions
  for (l in 1:L) {
    
    # Train/test split
    train.index <- sample(1:nrow(data_comb), trunc(nrow(data_comb)/2))
    test.index  <- setdiff(1:nrow(data_comb), train.index)
    test.y      <- y[test.index]
    # Initialize storage for this repetition
    count.vec   <- rep(0, nrow(candidate_mat))
    unique.row  <- rep(0, nrow(candidate_mat))
    sig.hat     <- rep(NA, nrow(candidate_mat))
    df.resid    <- rep(NA, nrow(candidate_mat))
    fitted.vals <- matrix(NA, nrow = nrow(candidate_mat), ncol = length(train.index))
    preds       <- matrix(NA, nrow = nrow(candidate_mat), ncol = length(test.index))
    edf         <- rep(NA, nrow(candidate_mat))
    C           <- rep(NA, nrow(candidate_mat))
    num.weight  <- rep(NA, nrow(candidate_mat))
    # Loop over each candidate row
    for (row_idx in 1:nrow(candidate_mat)) {
      
      # Check for duplicates among earlier rows
      if (row_idx > 1 && all(candidate_mat[row_idx, ] == candidate_mat[row_idx-1, ])) {
            sig.hat[row_idx]     <- sig.hat[row_idx-1]
            df.resid[row_idx]    <- df.resid[row_idx-1]
            fitted.vals[row_idx, ] <- fitted.vals[row_idx-1, ]
            preds[row_idx, ]     <-  preds[row_idx-1,]
            edf[row_idx]         <-  edf[row_idx-1]
            C[row_idx]           <-  C[row_idx-1]
      }
      
      else {
        # Extract the subset of variables for this candidate
        selected_vars <- which(candidate_mat[row_idx, ] == 1)
        data.sub  <- data_comb[, c(1, selected_vars + 1)] # +1 for y in col1
    
        
        train.sub <- data.sub[train.index, ]
        test.sub  <- data.sub[test.index, ]
        #print(head(test.sub))
        
          
        covariates <- setdiff(names(train.sub), "y")
        smooth_terms <- paste0("bbs(", covariates, ', degree = 2)')
        formula_str  <- paste("y ~", paste(smooth_terms, collapse = " + "))
          
        gam_fit <- gamboost(as.formula(formula_str), 
                            data = train.sub,
                            control = boost_control(mstop = mstop, nu = nu))
        
        #ifelse(length(unique(gam_fit$xselect())) == length(covariates),
               #print("yay"), 
               #print("nay"))
        
        df.resid[row_idx]     <- nrow(train.sub) - sum(hatvalues(gam_fit)) 
        fitted.vals[row_idx, ] <- gam_fit$fitted()
        sig.hat[row_idx] <-(sum((train.sub$y-fitted.vals[row_idx, ])^2))/
                            df.resid[row_idx]
        preds[row_idx,]  <- predict(gam_fit, newdata = test.sub)
        edf[row_idx]     <- sum(hatvalues(gam_fit))
        
        
      }
      C[row_idx]       <- edf[row_idx] * log(exp(1)*max(edf,na.rm = TRUE)/edf[row_idx]) + 
                            2 * log(edf[row_idx] + 2)
      
      num.weight[row_idx]  <- exp(-psi * C[row_idx] - (nrow(data_comb)) * sig.hat[row_idx] - sum(sig.hat[row_idx]^(-2) * (test.y - preds[row_idx,])^2 / 2 ))
                          
    }
    #print("DF Residual")
    #print(df.resid)
    #print("Fitted Vals")
    #print(fitted.vals)
    #print("Sigma Hats")
    #print(sig.hat)
    #print("Predicted Vals")
    #print(preds)
    #print("Predicted Vals Last Row")
    #print(preds[250,])
    #print(test.y)
    #print("Effective Degrees of Freedom")
    #print(edf)
    #print("AIC vals")
    #print(C)
    #print("Numerator Weights")
    #print(num.weight)
    # calculate weights and set them in lth vector in weights mat
    weights.mat[l, ] <- num.weight/sum(num.weight)
    print(l)
  }
  #print("Weights Matrix")
  #print(weights.mat)
  final_weights <- colMeans(weights.mat)
  return(final_weights)
}

# arm_gam_mboost with progress bar
arm_gam_mboost <- function(candidate_mat, L, y, X, nu, mstop) {
  
  # Storage for final weights
  weights.mat <- matrix(0, nrow = L, ncol = nrow(candidate_mat))
  
  # Combine y and X into one data frame
  data_comb <- data.frame(y = y, X)
  
  # initial weights for candidate models
  psi <- 1/nrow(candidate_mat)
  
  # set up progress bar
  pb <- txtProgressBar(min = 0, max = L, style = 3)
  
  # Loop over repetitions
  for (l in 1:L) {
    
    # Train/test split
    train.index <- sample(1:nrow(data_comb), trunc(nrow(data_comb)/2))
    test.index  <- setdiff(1:nrow(data_comb), train.index)
    test.y      <- y[test.index]
    
    # Initialize storage for this repetition
    sig.hat     <- rep(NA, nrow(candidate_mat))
    df.resid    <- rep(NA, nrow(candidate_mat))
    fitted.vals <- matrix(NA, nrow = nrow(candidate_mat), ncol = length(train.index))
    preds       <- matrix(NA, nrow = nrow(candidate_mat), ncol = length(test.index))
    edf         <- rep(NA, nrow(candidate_mat))
    C           <- rep(NA, nrow(candidate_mat))
    num.weight  <- rep(NA, nrow(candidate_mat))
    
    # Loop over each candidate row
    for (row_idx in 1:nrow(candidate_mat)) {
      
      if (row_idx > 1 && all(candidate_mat[row_idx, ] == candidate_mat[row_idx-1, ])) {
        sig.hat[row_idx]      <- sig.hat[row_idx-1]
        df.resid[row_idx]     <- df.resid[row_idx-1]
        fitted.vals[row_idx,] <- fitted.vals[row_idx-1,]
        preds[row_idx, ]      <- preds[row_idx-1,]
        edf[row_idx]          <- edf[row_idx-1]
        C[row_idx]            <- C[row_idx-1]
      } else {
        # Extract the subset of variables for this candidate
        selected_vars <- which(candidate_mat[row_idx, ] == 1)
        data.sub  <- data_comb[, c(1, selected_vars + 1)] 
        
        train.sub <- data.sub[train.index, ]
        test.sub  <- data.sub[test.index, ]
        
        covariates <- setdiff(names(train.sub), "y")
        smooth_terms <- paste0("bbs(", covariates, ", degree = 2)")
        formula_str  <- paste("y ~", paste(smooth_terms, collapse = " + "))
        
        gam_fit <- gamboost(as.formula(formula_str), 
                            data = train.sub,
                            control = boost_control(mstop = mstop, nu = nu))
        if (length(unique(gam_fit$xselect())) != length(covariates)){
          print(l)
          print(row_idx)
          print(covariates)
          print(unique(gam_fit$xselect()))
        }
        
        df.resid[row_idx]      <- nrow(train.sub) - sum(hatvalues(gam_fit)) 
        fitted.vals[row_idx, ] <- gam_fit$fitted()
        sig.hat[row_idx]       <- sum((train.sub$y - fitted.vals[row_idx, ])^2) / df.resid[row_idx]
        preds[row_idx, ]       <- predict(gam_fit, newdata = test.sub)
        edf[row_idx]           <- sum(hatvalues(gam_fit))
      }
      
      C[row_idx] <- edf[row_idx] * log(exp(1) * max(edf, na.rm = TRUE)/edf[row_idx]) + 
        2 * log(edf[row_idx] + 2)
      
      num.weight[row_idx] <- exp(-psi * C[row_idx] -
                                   (nrow(data_comb)) * sig.hat[row_idx] -
                                   sum(sig.hat[row_idx]^(-2) * 
                                   (test.y - preds[row_idx,])^2 / 2 ))
    }
    #print("DF Residual")
    #print(df.resid)
    #print("Fitted Vals")
    #print(fitted.vals)
    #print("Sigma Hats")
    #print(sig.hat)
    #print("Predicted Vals")
    #print(preds)
    #print("Predicted Vals Last Row")
    #print(preds[250,])
    #print(test.y)
    #print("Effective Degrees of Freedom")
    #print(edf)
    #print("AIC vals")
    #print(C)
    #print("Numerator Weights")
    #print(num.weight)
    # calculate weights and set them in lth vector in weights mat
    weights.mat[l, ] <- num.weight/sum(num.weight)
    
    # update progress bar
    setTxtProgressBar(pb, l)
  }
  
  close(pb)  # close progress bar
  
  final_weights <- colMeans(weights.mat)
  return(final_weights)
}

# Creating a combined function that executes all three functions
soil.gam.mboost = function(boosting_iters, y, X, num_splits, nu, mstop) {
  candidate.mat = boost_extract(boosting_iters, y, X)
  arm.weights = arm_gam_mboost(candidate.mat, num_splits, y, X, nu, mstop)
  soil.scores = soil_gam(candidate.mat,arm.weights)
  return(soil.scores)
}
