##############################
# MATCHING ESTIMATORS		     #
##############################

library(Matching)

matching_estimators <- function(data) {
  
  R <- length(data)
  
  # output matrix
  
  out <- matrix(NA_real_, nrow = R, ncol = 8)
  colnames(out) <- c(
    "mat_lin_pgs_t",  "mat_lin_pgs_f",
    "mat_nn_pgs_t",   "mat_nn_pgs_f",
    "mat_lin_fpgs_t", "mat_lin_fpgs_f",
    "mat_nn_fpgs_t",  "mat_nn_fpgs_f"
  )
  
  for (i in seq_len(R)) {
    
    datat <- as.data.frame(data[[i]])
    
    
    if (ncol(datat) != 10) {
      stop(paste("Replicate", i, "does not have 10 columns."))
    }
    
    # ---- force numeric vectors ----
    Y  <- as.numeric(datat$y)
    Tr <- as.numeric(datat$tr)
    
    # numeric matrix 
    makeX <- function(df) {
      X <- unname(as.matrix(df))
      storage.mode(X) <- "double"
      X
    }
    
    # -------- PGS (1-dimensional) --------
    x_pgs_t    <- makeX(datat["mu0"])
    x_pgs_f    <- makeX(datat["mu0f"])
    x_pgs_nl_t <- makeX(datat["murf0"])
    x_pgs_nl_f <- makeX(datat["murf0_f"])
    
    # -------- FPGS (2-dimensional) --------
    x_fpgs_t    <- makeX(datat[c("mu0","mu1")])
    x_fpgs_f    <- makeX(datat[c("mu0f","mu1f")])
    x_fpgs_nl_t <- makeX(datat[c("murf0","murf1")])
    x_fpgs_nl_f <- makeX(datat[c("murf0_f","murf1_f")])
    
    # -------- Matching estimators--------
    
    out[i,1] <- Match(Y = Y, Tr = Tr, X = x_pgs_t,    estimand = "ATE", M = 1)$est
    out[i,2] <- Match(Y = Y, Tr = Tr, X = x_pgs_f,    estimand = "ATE", M = 1)$est
    out[i,3] <- Match(Y = Y, Tr = Tr, X = x_pgs_nl_t, estimand = "ATE", M = 1)$est
    out[i,4] <- Match(Y = Y, Tr = Tr, X = x_pgs_nl_f, estimand = "ATE", M = 1)$est
    
    out[i,5] <- Match(Y = Y, Tr = Tr, X = x_fpgs_t,    estimand = "ATE", M = 1)$est
    out[i,6] <- Match(Y = Y, Tr = Tr, X = x_fpgs_f,    estimand = "ATE", M = 1)$est
    out[i,7] <- Match(Y = Y, Tr = Tr, X = x_fpgs_nl_t, estimand = "ATE", M = 1)$est
    out[i,8] <- Match(Y = Y, Tr = Tr, X = x_fpgs_nl_f, estimand = "ATE", M = 1)$est
  }
  
  return(out)
}
