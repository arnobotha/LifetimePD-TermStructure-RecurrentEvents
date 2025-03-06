# ============================== SURVIVAL FUNCTIONS ==============================
# Defining bespoke customs relating to various generic aspects of survival analysis
# --------------------------------------------------------------------------------
# PROJECT TITLE: Default Survival Modelling
# SCRIPT AUTHOR(S): Marcel Muller, Dr Arno Botha, Bernard Scheepers

# VERSION: 1.0 (July-2023)
# DESCRIPTION: 
# This script defines various functions specific to survival modelling
# that are used elsewhere in this project or, indeed, used across other projects.
# Functions are grouped thematically.
# ================================================================================




# ----------------- 1. Functions related to the modelling process------------------

# --- Function to ensure that the column is indeed part of the given dataset.
# Input:  [cols], Vector containing the column names to be checked if they exist in the dataset.
#         [dataset], Dataset to check on if the columns exist.
# Output: print on whether the columns are in the dataset or which columns are not in the dataset.
colCheck <- function(cols, dataset) {
  # Check if all columns exist in the dataset
  missing_cols <- cols[!(cols %in% colnames(dataset))]
  
  if (length(missing_cols) == 0) { # All columns are in the dataset.
    paste0("SAFE: All columns are in the dataset")
  } else { # List columns not in the dataset.
    warning("WARNING: Some columns are not in the dataset.")
    paste0("Columns not in the dataset: ", paste(missing_cols, collapse = ", "))
  }
}



# --- Function to add and/or remove certain values from a vector
# Input:  [mat], (2 x n) Matrix on which changes will be made. The first column contains the variable names (vars) and the second their respective variable type (vartypes)
#         [Remove], Vector containing the entries to be removed.
# Output:        [m], Updated matrix.
vecChange <- function(mat,Remove=FALSE,Add=FALSE){
  m <- mat
  if (all(Remove != FALSE)) {m <- mat[!(mat$vars %in% Remove)]} # Remove values in [Remove]
  if (all(Add != FALSE)) {m <- rbind(m,Add[!(Add$vars %in% mat$vars)])} # Add values in [Add]
  return(m)
}



# --- Function to detect significant correlations (abs(cor) > 0.6) between vectors.
# Input:  [data], Dataset containing the variables in varlist to dertermine correlations.
#         [varlist], list of variable names to determine correlations.
#         [corrThresh], the absolute correlation threshold above which correlations are deemed significant.
#         [method], the method by which correlations are determined.
#         [Remove], Vector containing the entries to be removed.
# Output: <graph>  Upper half of correlation matrix.
#         <print> Text indicating the variable pairs with high correlation.
corrAnalysis <- function(data, varlist, corrThresh = 0.6, method = 'spearman') {
  # Compute the correlation matrix
  corrMat <- as.data.table(data) %>% subset(select = varlist) %>% cor(method = method)
  
  # Visualize the correlation matrix
  corrplot(corrMat, type = 'upper', addCoef.col = 'black', tl.col = 'black', diag=FALSE,
           tl.srt = 45)
  
  # Find correlation coordinates exceeding the threshold
  corrCoordinates <- which(abs(corrMat) > corrThresh & abs(corrMat) < 1 & upper.tri(corrMat), arr.ind = TRUE)
  
  if(nrow(corrCoordinates) != 0){
    # Create a data table with correlation pairs
    corrProbs <- data.table(x = rownames(corrMat)[corrCoordinates[, 1]], y = colnames(corrMat)[corrCoordinates[, 2]])
    
    # Print the identified correlations
    for (i in 1:nrow(corrProbs)) {
      cat("Absolute correlations of ",percent(corrMat[corrProbs[i, x], corrProbs[i, y]]),
          " found for ", corrProbs[i, x], " and ", corrProbs[i, y],"\n")
    }
  }else{
    cat("No significant correlations were detected")
  }
}



# --- Function to return the appropriate formula object based on the time definition.
#         [TimeDef], Time definition on which the cox ph models are based on.
#         [var], Single variable name
TimeDef_Form <- function(TimeDef="TFD", vars){
  # Create formula based on time definition of the dataset.
  if(TimeDef=="TFD"){# Formula for time to first default time defintion (containing only the fist performance spell).
    formula <- as.formula(paste0("Surv(Start,End,Default_Ind) ~ ",
                                 paste(vars,collapse=" + ")))
  }else if(TimeDef=="PWP_ST"){# Formula for Prentice-Williams-Peterson Spell time definition (containing only the fist performance spell).
    formula <- as.formula(paste0("Surv(Start,End,Default_Ind) ~ PerfSpell_Grp + ",
                                 paste(vars,collapse=" + ")))
  }
}



# --- Function to fit a given formula within a Cox regression model towards extracting Harrell's C-statistic and related quantities
calc_HarrellC <- function(formula, data_train, data_valid, variable="", it=NA, logPath="") {
  model <- coxph(formula,id=get(fldSpellID), data = data_train) # Fit Cox model
  if (!is.na(it)) {
    cat(paste0("\n\t ", it,") Single-factor survival model built. "),
        file=paste0(logPath,"HarrelsC_log.txt"), append=T)
  }
  c <- concordance(model, newdata=data_valid) # Calculate concordance of the model based on the validation set.
  conc <- as.numeric(c[1])# Extract concordance
  sd <- sqrt(c$var)# Extract concordance variability
  lr_stat <- round(2 * (model$loglik[2] - model$loglik[1]),0)# Extract LRT from the model's log-likelihood
  # Return results as a data.table
  return(data.table(Variable = variable, Concordance = conc, SD = sd, LR_Statistic = lr_stat))
}



# --- Function to extract the concordances (Harrell's C) from single-factor models
# Input:  [data_train], Training dataset on which the models are built on.
#         [data_valid], Validation dataset on which the models' concordance is validated on.
#         [variables], Vector containing a list of the variables for the models.
#         [TimeDef], Time definition on which the cox ph models are based on.
# Output: [matResults]  Table containing the concordance, se(concordance) and log ratio of the singular models.
concTable <- function(data_train, data_valid, variables, fldSpellID="PerfSpell_Key",
                      TimeDef="TFD", numThreads=6, genPath) {
  # - Testing conditions
  # data_valid <- datCredit_valid_TFD; TimeDef="TFD"' numThreads=6
  
  # - Iterate across loan space using a multi-threaded setup
  ptm <- proc.time() #IGNORE: for computation time calculation
  cl.port <- makeCluster(round(numThreads)); registerDoParallel(cl.port) # multi-threading setup
  cat("New Job: Estimating B-statistic (1-KS) for each variable as a single-factor survival model ..",
      file=paste0(genPath,"HarrelsC_log.txt"), append=F)
  
  results <- foreach(j=1:length(variables), .combine='rbind', .verbose=F, .inorder=T,
                                 .packages=c('data.table', 'survival'), .export=c('calc_HarrellC')) %dopar%
    { # ----------------- Start of Inner Loop -----------------
      # - Testing conditions
      # j <- 1
      calc_HarrellC(formula=TimeDef_Form(TimeDef,variables[j]), variable=variables[j],
                    data_train=data_train, data_valid=data_valid, it=j, logPath=genPath)
    } # ----------------- End of Inner Loop -----------------
  stopCluster(cl.port); proc.time() - ptm  
  
  # Sort by concordance in descending order.
  setorder(results, -Concordance)
  
  # Return resulting table.
  return(results)
}



# --- Function to calcualte the complement of the KS test statistic "B-statistic"
# Inputs: [formula]: Cox regression formula object; [data_train]: training data
#         [fldSpellID]: Field name of spell-level ID; [vEvents]: spell-level vector of event indicators
#         [seedVal]: Seed value for random number generation; [it]: optional iteration parameter for logging purposes;
#         [logPath]: Optional path for log file for logging purposes
# Outputs: b-statistic (single value)
calcBStat <- function(formula, data_train, fldSpellID="PerfSpell_Key", vEvents, seedVal, it=NA, logPath=NA) {
  # Fit Model
  model <- coxph(formula, data = data_train, id=get(fldSpellID))
  
  if (!is.na(it)) {
    cat(paste0("\n\t ", it,") Single-factor survival model built. "),
        file=paste0(logPath,"BStat_log.txt"), append=T)
  }
  
  # Calculate Cox-Snell (adjusted) residuals
  vCS <- calc_CoxSnell_Adj(model, vIDs=data_train[[fldSpellID]], vEvents=vEvents)
  # Initialize a unit exponential distribution
  set.seed(seedVal, kind = "Mersenne-Twister")
  vExp <- rexp(length(vCS),1)
  # Perform the two-sample Kolmogorov-Smirnov test of distribution equality
  #   H_0: vCS and vExp originates from the same distribution
  #   NOTE: We only desire the KS test statistic in measuring distributional dissimilarity
  #   Then, we subtract this from 1 in creating a coherent statistic; greater is better
  bStat <- 1 - round(suppressWarnings(ks.test(vCS,vExp))$statistic,4)
  return(bStat)
}



# --- Function to extract the B-statistic from a range of models built on a list of variables based on a time definition.
# Input:  [data_train], Training dataset on which the models are built on.
#         [seedVal], Seed value to ensure results are reproducible.
#         [numIt], Number of simulations to calculate the B-statistic.
#         [TimeDef], Time definition on which the cox ph models are based on.
# Output: [Results]  Table containing the concordance, se(concordance) and log ratio of the singular models.
csTable <- function(data_train, variables, TimeDef="TFD", seedVal=1, numIt=5, 
                    fldSpellID="PerfSpell_Key", fldLstRowInd="PerfSpell_Exit_Ind", fldEventInd="Default_Ind",
                    numThreads=6, genPath=NA){
  
  # - Testing conditions
  # data_train <- datCredit_train_TFD; variables<-vars2; TimeDef<-"TFD"; seedVal<-1; numIt<-5; 
  # fldLstRowInd="PerfSpell_Exit_Ind";  fldSpellID="PerfSpell_Key"; fldEventInd="Default_Ind"; numThreads=6
  
  # - Initialize results
  results <- data.frame(Variable = variables, B_Statistic = NA_real_)
  
  # - Data preparation
  # Subset last row per performing spell for Goodness-of-Fit (GoF) purposes
  datLstRow <- copy(data_train[get(fldLstRowInd)==1,])
  vLstRow_Events <- datLstRow[, get(fldEventInd)]
  
  # - Simulate null distribution if seedVal is not NA
  if (!is.na(seedVal)) {
    
    
    # - Iterate across loan space using a multi-threaded setup
    ptm <- proc.time() #IGNORE: for computation time calculation
    cl.port <- makeCluster(round(numThreads)); registerDoParallel(cl.port) # multi-threading setup
    cat("New Job: Estimating B-statistic (1-KS) for each variable as a single-factor survival model ..",
        file=paste0(genPath,"BStat_log.txt"), append=F)
    
    results$B_Statistic <- foreach(j=1:length(variables), .combine='rbind', .verbose=F, .inorder=T,
                      .packages=c('data.table', 'survival'), .export=c('calc_CoxSnell_Adj', 'calcBStat', 'TimeDef_Form')) %dopar%
      
      { # ----------------- Start of Inner Loop -----------------
        # - Testing conditions
        # var <- variables[1]
        calcBStat(formula=TimeDef_Form(TimeDef,variables[j]), data_train=data_train, fldSpellID=fldSpellID, vEvents=vLstRow_Events,
                  seedVal=seedVal, it=j, logPath=genPath)
        
      } # ----------------- End of Inner Loop -----------------
    stopCluster(cl.port); proc.time() - ptm
    
    # Sort results by B statistic in descending order
    results <- results[order(-results$B_Statistic, na.last = TRUE), ]
    
    # Return results and range of B statistics
    return(list(Results = results, Range = diff(range(results$B_Statistic, na.rm = TRUE))))
    
  } else {
    # Perform iterative B calculation when seedVal is NA
    # Initialize Results matrix to contain the number of interations
    matResults <- matrix(NA, nrow = length(variables), ncol = numIt,
                         dimnames = list(variables,
                                         paste0("Iteration_", 1:numIt))) %>%
                          as.data.table()
    
    # - Iterate across loan space using a multi-threaded setup
    ptm <- proc.time() #IGNORE: for computation time calculation
    cl.port <- makeCluster(round(numThreads)); registerDoParallel(cl.port) # multi-threading setup
    cat("New Job: Estimating B-statistics (1-KS) ..",
        file=paste0(genPath,"BStat_log.txt"), append=F)
    
    for (it in seq_len(numIt)) {
      
      cat(paste0("\n Estimating B-statistic (1-KS) for each variable as a single-factor survival model for iteration ", it, " .."),
          file=paste0(genPath,"BStat_log.txt"), append=T)
      
      matResults[, it] <- foreach(j=1:length(variables), .combine='rbind', .verbose=F, .inorder=T,
                                     .packages=c('data.table', 'survival'), .export=c('calc_CoxSnell_Adj', 'calcBStat', 'TimeDef_Form')) %dopar%
        
        { # ----------------- Start of Inner Loop -----------------
          # - Testing conditions
          # var <- variables[1]
          calcBStat(formula=TimeDef_Form(TimeDef,variables[j]), data_train=data_train, fldSpellID=fldSpellID, vEvents=vLstRow_Events,
                    seedVal=seedVal*it, it=j, logPath=genPath)
          
        } # ----------------- End of Inner Loop -----------------
    }
    stopCluster(cl.port); proc.time() - ptm
    
    # Compute additional statistics for the results matrix
    colRanges <- matResults[, lapply(.SD, function(x) diff(range(x, na.rm = TRUE)))] # Calculate the Range of B-statistic value for each iteration
    matResults[, Average := rowMeans(.SD, na.rm = TRUE), .SDcols = patterns("^Iteration_")]# Calculate the average B-statistic for each variable
    matResults[, Variable := variables] 
    matResults <- matResults %>% relocate(Variable, .before=Iteration_1)
    
    #matResults <- cbind(Variables = c(variables,"Range"),matResults)# Add a column to the Results matrix to cross reference the variables with their respective B-statistics
    setorder(matResults,-Average)# Arrange matrix according to average
    
    # Return matrix of B statistics
    return(list(Results=matResults, IterationRanges=colRanges))
  }
}



# --- Function to calculate various survival-related quantities for a given loan history
# Input:    [datGiven]: given loan history; [coxGiven]: fitted cox PH model; [it]: current iteration index; 
#           [numKeys]: total keys; [genPath]: path of reporting log
# Output:   Survival probability, cumulative hazard, hazard, and event probability
survQuants <- function(datGiven, coxGiven, it=1, numKeys, genPath="") {
  # datGiven <- subset(test,PerfSpell_Key == vSpellKeys[2]); coxGiven <- cox_TFD
  # it=1; numKeys <- numSpellKeys
  
  # - Compute individual survival curve from fitted Cox model
  survFit_pred <- survfit(coxGiven, centered=F, newdata=datGiven, id=PerfSpell_Key)
  
  cat("\n\t", it, "of", numKeys, "| Estimation completed for spell key:", unique(datGiven$PerfSpell_Key),
      file=paste0(genPath,"survQuants_log.txt"), append=T)
  
  datSurv <- data.table(PerfSpell_Key = unique(datGiven$PerfSpell_Key), End=datGiven$End, # composite key
                        CHaz=survFit_pred$cumhaz, #RiskSetSize=survFit_pred$n.risk,
                        #NumEvents=survFit_pred$n.event, NumCensored=survFit_pred$n.censor,
                        Survival=round(survFit_pred$surv,digits=15))
  # plot(survFit_pred)
  
  # - Approximate baseline hazard h_0(t) from cumulative baseline hazard H_0(t)
  datSurv[, Hazard := c(datSurv$CHaz[1]/datSurv$End[1], diff(datSurv$CHaz) / diff(datSurv$End))]
  #plot(datSurv[, Time], datSurv[, Hazard_base], type="b") # mean survival probability over time
  datSurv[, EventProb := Hazard * shift(Survival,n=1,type="lag",fill=1)] # f(t|X) = S(t-1|X) . h(t|X)
  
  # - Render survival predictions: risk scores exp(\beta . x) and calculate hazard and survival probability
  #datSurv[, RiskScore := predict(coxGiven, newdata=datGiven, type="risk")]
  #datSurv[, Hazard := Hazard_base * RiskScore] # by definition of Cox regression model
  
  return(datSurv)
}


# --- Function to calculate various survival-related quantities for a given loan history
# Input:    [datGiven]: given loan history; [coxGiven]: fitted cox PH model; [it]: current iteration index; 
#           [numKeys]: total keys; [genPath]: path of reporting log
# Output:   Survival probability, cumulative hazard, hazard, and event probability
survQuants.data <- function(datGiven, vars, beta, vStartTimes, vStopTimes, vStatus, id, centered = F, datBaselineHaz,
                            it=1, numKeys, genPath="") {
  # datGiven <- subset(test,PerfSpell_Key == vSpellKeys[2]); coxGiven <- cox_TFD
  # it=1; numKeys <- numSpellKeys; centered <- T
  # vars <- vecVars_TFD; beta <- cox_TFD$coefficients
  #X <- as.matrix(datGiven[, mget(vars)])
  
  # Compute linear predictors and risk scores, depending on centering
  if (centered==F) {
    lp <- (X %*% beta)
  } else {
    lp_sample <- colMeans(X) %*% beta # assuming only numeric covariates; 0 for categorical
    lp_indiv <- (X %*% beta)
    lp <- lp_indiv - as.numeric(lp_sample)
  }
  risk_scores <- exp(lp)  # exp(Xβ)
  
  # --- Calculate baseline cumulative hazard | Breslow estimate
  
  ### AB: Draft further from here, by removing the dependence on coxGiven.
  # X <- as.matrix(datCredit_train_TFD[, mget(vars)])
  datHaz <- baseHaz.coxph.data(X=as.matrix(datCredit_train_TFD[, mget(vars)]), beta=beta, vStartTimes=datCredit_train_TFD$Start, vStopTimes=datCredit_train_TFD$End, 
                               vStatus=datCredit_train_TFD$Default_Ind, id=datCredit_train_TFD$PerfSpell_Key, centered=F, ties="Breslow")
  datHaz[1:6,]
  
  # - Compute individual survival curve from fitted Cox model
  survFit_pred <- survfit(coxGiven, centered=F, newdata=datGiven, id=PerfSpell_Key)
  
  cat("\n\t", it, "of", numKeys, "| Estimation completed for spell key:", unique(datGiven$PerfSpell_Key),
      file=paste0(genPath,"survQuants_log.txt"), append=T)
  
  (datSurv <- data.table(PerfSpell_Key = unique(datGiven$PerfSpell_Key), End=datGiven$End, # composite key
                        CHaz=survFit_pred$cumhaz, RiskSetSize=survFit_pred$n.risk,
                        NumEvents=survFit_pred$n.event, NumCensored=survFit_pred$n.censor,
                        Survival=round(survFit_pred$surv,digits=15)))
  # plot(survFit_pred)
  
  # - Approximate baseline hazard h_0(t) from cumulative baseline hazard H_0(t)
  datSurv[, Hazard := c(datSurv$CHaz[1]/datSurv$End[1], diff(datSurv$CHaz) / diff(datSurv$End))]
  #plot(datSurv[, Time], datSurv[, Hazard_base], type="b") # mean survival probability over time
  datSurv[, EventProb := Hazard * shift(Survival,n=1,type="lag",fill=1)] # f(t|X) = S(t-1|X) . h(t|X)
  
  # - Render survival predictions: risk scores exp(\beta . x) and calculate hazard and survival probability
  #datSurv[, RiskScore := predict(coxGiven, newdata=datGiven, type="risk")]
  #datSurv[, Hazard := Hazard_base * RiskScore] # by definition of Cox regression model
  
  return(datSurv)
}


baseHaz.coxph.data <- function(X, beta, vStartTimes, vStopTimes, vStatus, id, centered = F, ties="Breslow",
                               vStrata=NULL) {
  # - Testing conditions
  # X=as.matrix(datCredit_train_TFD[, mget(vars)]); beta <- cox_TFD$coefficients
  # X=as.matrix(datGiven[, mget(vars)])
  # vStartTimes=datCredit_train_TFD$Start; vStopTimes=datCredit_train_TFD$End
  # vStatus=datCredit_train_TFD$Default_Ind; id=datCredit_train_TFD$PerfSpell_Key; centered=F; ties="Breslow"; vStrata=NULL
  
  # Compute linear predictors and risk scores, depending on centering
  if (centered==F) {
    lp <- (X %*% beta)
  } else {
    lp_sample <- colMeans(X) %*% beta # assuming only numeric covariates; 0 for categorical
    lp_indiv <- (X %*% beta)
    lp <- lp_indiv - as.numeric(lp_sample)
  }
  
  risk_scores <- exp(lp)  # exp(Xβ)
  
  lp2 <- predict(cox_TFD, type="lp", centered=T)
  head(X)
  head(lp); head(lp2)
  all.equal(as.numeric(risk_scores), exp(lp2))
  risk_scores <- as.matrix(exp(lp2))
  
  # Compute mean linear predictor for centering
  mean_lp <- mean(lp)
  
  # Identify stratification (if applicable)
  if (all(!is.null(vStrata))) {
    strata_levels <- unique(vStrata)
    strata_indices <- rep(strata_levels, length(vStrata))
  } else {
    strata_indices <- rep("Main", length(vStartTimes))
  }
  
  # Initialize storage for baseline hazards
  baseline_hazard_list <- list()
  
  # Loop over strata to compute baseline hazard separately
  for (stratum in unique(strata_indices)) {
    # stratum <- strata_indices[1]
    
    # Filter data for the current stratum
    stratum_mask <- strata_indices == stratum
    vStartTimes_s <- vStartTimes[stratum_mask]
    vStopTimes_s <- vStopTimes[stratum_mask]
    status_s <- vStatus[stratum_mask]
    risk_scores_s <- risk_scores[stratum_mask]
    
    # Further subset for individuals with a single event (vStatus == 1)
    event_mask <- (status_s == 1)
    
    if (sum(event_mask) == 0) {
      next  # Skip this stratum since no events occurred
    }
    
    # Unique failure times for event status
    #event_times <- sort(unique(vStopTimes_s[event_mask]))
    event_times <- sort(unique(vStopTimes_s)) # corresponds with basehaz(), and therefore includes censoring times
    
    # Initialize cumulative baseline hazard & related quantities
    H0 <- numeric(length(event_times))
    riskSetSize <- numeric(length(event_times))
    numEvents <- numeric(length(event_times))
    
    # Compute hazard at each failure time using Breslow's method
    for (i in seq_along(event_times)) {
      # i <- 2
      t_i <- event_times[i]
      
      # Identify individuals who had an event at t_i
      vEventIDs <- id[vStopTimes == t_i & status_s == 1]
      d_i <- length(vEventIDs) # number of events at t_i
      numEvents[i] <- d_i
      
      # Risk set: sum of exp(Xβ) for individuals still at risk at t_i
      risk_set_ind <- vStartTimes_s < t_i & vStopTimes_s >= t_i & !(id %in% vEventIDs & vStopTimes < t_i)
      (riskSetSize[i] <- sum(risk_set_ind)) # same as sum(vStartTimes_s < t_i) when not handling ids
      risk_set <- sum(risk_scores_s[vStartTimes_s < t_i])
      (risk_set <- sum(risk_scores_s[risk_set_ind]))
      
      # Calculate the Efron adjustment for ties
      if (d_i > 1 & ties=="Efron") {
        # Adjust for ties at the event time using Efron estimator
        tied_contribution <- sum(risk_scores_s[vStartTimes_s < t_i & vStopTimes_s == t_i]) / risk_set
        risk_set <- risk_set + tied_contribution
      }
      
      # Efron estimator: \delta(t) = d_i / risk_set
      H0[i] <- (ifelse(i > 1, H0[i-1], 0)) + (d_i / risk_set)
    }
    
    # Apply centering correction if requested
    if (centered) {
      H0 <- H0 * exp(-mean_lp)
    }
    
    # Store results
    baseline_hazard_list[[stratum]] <- data.frame(
      time = event_times,
      CHaz = H0,
      RiskSetSize = riskSetSize,
      NumEvents = numEvents,
      strata = stratum
    )
  }
  
  # Combine all strata into a single dataframe
  baseline_hazard_df <- do.call(rbind, baseline_hazard_list)
  
  return(baseline_hazard_df)
}



# HW BS: Comment better the spline function
### AB: Marked for deletion, rather use lm(output ~ ns(input, df=3), data=dat) from the splines package
# -- Cubic spline function
spline_estimation <- function(times, hazard, nknots, degree) {
  n <- length(times)  # Number of observations
  
  # Calculate quantiles for knot placement
  qtiles <- (1:nknots + 1) / (nknots + 2)  # Generate nknots evenly spaced quantiles
  knots <- quantile(times, probs = qtiles, na.rm = TRUE)
  
  # Construct the T matrix efficiently
  matT <- cbind(
    sapply(0:degree, function(j) times^j),                           # Polynomial terms
    sapply(knots, function(k) pmax(0, (times - k)^degree))           # Truncated power basis
  )
  
  # Ordinary Least Squares: Compute the parameter estimates
  coef <- solvet(crossprod(matT), tol=1e-15) %*% crossprod(matT, hazard)  # More efficient than t(T) %*% T and t(T) %*% y
  
  # Polynomial terms: ∑_{j=0}^{d} α_j * t^j
  poly_terms <- sapply(0:degree, function(j) coef[j + 1] * times^j)
  
  # Truncated power basis terms: ∑_{p=1}^{r} α_{p+d} * (t - K_p)^d_+
  truncated_terms <- sapply(1:nknots, function(p) {
    coef[1 + degree + p] * pmax(0, (times - knots[p])^degree)
  })
  
  # Combine polynomial and truncated power terms
  y <- rowSums(poly_terms) + rowSums(truncated_terms)
  y <- pmax(y,0)
  return(y)
}







### AB: Given the work of Bernard, I'm no longer sure of the utility of the below.

# --- function to compute the (unscaled) Schoenfeld residuals for a Cox PH model
#   1) Schoenfeld residuals are computed for each specified variable in the training dataset
#   2) Tests are conducted for significance of the correlation of the Schoenfeld residuals against time
#   3) Graphs are created of the Schoenfeld residuals against time
# Input:  cph - Cox PH model to be assessed
#         dat_train - Dataset used to train the Cox PH model
#         var - Name of the variables for which the Schoenfeld residuals are to be computed
#         id - Name of the column uniquely identifying each row of the dat_train
#         time - Name of the column identifying the associated time in dat_train
#         status - Name of the column identifying the target variable in dat_train
#         verbose - Indicator variable for supressing graphs created by the function
#         max_time - The maximum time for which the graph should be plotted
# Output: A list containing the following:
#         data - A dataset containing the id, time, raw variable's value, and the associated Schoenfeld residual
#         CorTest - A list of the results from test(s) for significance of the correlation of the Schoenfeld residuals against time
#         plots - A list of graphs of the Schoenfeld residuals against time
cph_schoen <- function(cph, var=NULL, dat_train, id, time, status, verbose=T, max_time=NULL){
  # UNIT TEST (VARIABLE INITIALISATION)
  # cph <- cph_Default_PH_test; var <- c("Principal_wins"); id <- "PerfSpell_Key"; time <- "TimeInPerfSpell"; status <- "DefaultStatus1"; max_time <- 240; verbose <- F
  
  # Copying the training dataset to ensure no contamination and ensuring that the dataset is of the correct class (this step usually takes a considerable amount of time for large datasets)
  if (any(class(dat_train) %in% "tbl_df")){
    dat_train2 <- as.data.table(dat_train)
  } else {
    dat_train2 <- copy(dat_train)
  }
  if (any(class(dat_train2) %in% "grouped_df")){
    dat_train2 <- ungroup(dat_train2)
  }
  
  cph_sum <- summary(cph) # Getting a summary of the model
  row_names <- rownames(cph_sum$coefficients) # Getting the names of all the variables in the Cox model (includes levels of categorical variables)
  
  # Object for graphs
  gplots <- list()
  
  # Object for correlation test
  cor_test <- list()
  
  # Getting all the variable names in the model if no names are specified
  if (is.null(var)){
    var <- unlist(strsplit(toString(cph_sum$call$formula[[3]]), '[,+ ]+'))
    var <- var[var!=""]
  }
  
  # Initialising the dataset to be returned
  dat_return <- data.table(ID = numeric(),
                           Time = numeric(),
                           Var_Val = numeric(),
                           Sch_Res = numeric(),
                           Var_Name = as.character(),
                           Var_Name_Base = as.character())
  
  # Initialising the temporary dataset used in the loop for all variables
  col_names <- c(id, time, status, var) # Getting the names of the columns from the training dataset to be subsetted
  dat_temp <- dat_train2[, ..col_names] # Subsetting from the main dataset
  colnames(dat_temp)[1:3] <- c("ID", "Time", "Status") # Renaming the columns for conveinience
  dat_temp <- cbind(dat_temp, predict(cph, dat_train2, type="risk")); colnames(dat_temp)[length(colnames(dat_temp))] <- "RiskScore" # Getting the risk score of each observation
  dat_temp <- merge(dat_temp, dat_temp[, list(SumScore=sum(RiskScore)), by=list(Time)], by="Time", all.x=T) # Getting the total risk score at each time point and merging it back into the dataset
  
  # Computing the Schoenfeld residuals for each selected variable and appending the dataset which to return
  k <- 0 # Counting variable for lists in for loop below
  for (i in 1:length(var)){
    # i <- 1
    var_type <- class(dat_temp[[var[i]]]) # The class of the selected variable
    var_name <- var[i]
    
    col_names2 <- c("ID", "Time", "Status", "RiskScore", "SumScore", var_name)
    dat_temp2 <- dat_temp[, ..col_names2] # Subsetting from the temporary dataset to only include the i'th variable's values
    colnames(dat_temp2)[6] <- "Var_Val"
    dat_temp2[, Var_Name := var_name]; dat_temp2[, Var_Name_Base := var_name] # [Var_Name] is the name of the variable as in the training dataset; [Var_Name_Base] is the level of the categorical variable (equal to Var_Name for numeric variables)
    
    # Computing the Schoenfeld residuals for a numeric variable
    if (var_type == "numeric"){
      # Computing the risk-weighted value of the chosen variable
      dat_temp2[, RW_Val := Var_Val*RiskScore/SumScore]
      # Computing the expected value of the numerical variable at each time point and merging it back into the dataset
      dat_temp2 <- merge(dat_temp2, dat_temp2[, list(Exp_Val=sum(RW_Val)), by=list(Time)], by="Time", all.x = TRUE)
      
      # Subsetting to only include defaulted accounts || Change naming accordingly as structure of data changes
      dat_temp3 <- dat_temp2[Status==1,]; rm(dat_temp2)
      # Computing the Schoenfeld residuals
      dat_temp3[, Sch_Res := Var_Val - Exp_Val] # -447886.85507
      
      # Appending the temporary dataset to the main dataset that is to be returned
      dat_return <- rbind(dat_return, dat_temp3[, list(ID, Time, Var_Val, Sch_Res, Var_Name, Var_Name_Base)])
      
      # Updating the counter variable
      k <- k + 1
      
      # Plotting if verbose = True
      if (verbose==F){
        if (is.null(max_time)){
          max_time <- max(dat_return$Time)
        } # if

        g <- ggplot(dat_temp3, aes(x=Time, y=Sch_Res)) +
             theme_minimal() + xlab(expression(Default~Time~italic(tau[d]))) + ylab(expression(italic(hat(Sc)(x[k~i~tau[d]])))) +
             theme(text=element_text(family=chosenFont),legend.position = "bottom",
                   strip.background=element_rect(colour="grey", fill="#D3D3D3")) +
             facet_wrap(~ Var_Name, strip.position = "right") +
             geom_point(shape=1, col="#3F702F") + geom_smooth(aes(linetype="solid"), col="#043927", fill="#043927") +
             scale_x_continuous(limits = c(NA, max_time)) +
             scale_linetype_manual(name = "Loess-smoother", values = "solid", labels = NULL)
        
        g_name <- paste0("SchRes_",var[i])
        
        gplots[[g_name]] <- g; names(gplots[k])  
      } # if

      # Correlation Test (Formal)
      cor_test[[k]] <- cor.test(dat_temp3$Time, dat_temp3$Sch_Res, method = "spearman"); names(cor_test)[k] <- var[i] # Spearman rank correlation is used for robustness
      
      # Clean up
      rm(dat_temp3)
      
    } else if (var_type %in% c("character", "factor")){
      levels_n <- length(grep(var[i], row_names)) # Computing the number of levels (-1) of the categorical variable | k-1 levels of the categorical variable
      levels <-  substring(row_names[grep(var[i], row_names)], nchar(var[i])+1) # Getting the levels of the categorical variable (-1)
      
      # Computing the Schoenfeld residuals for each level of the categorical variable
      for (j in 1:levels_n){
        # j <- 1
        dat_temp3 <- copy(dat_temp2) # Copying the temporary dataset to increase efficiency with multiple levels of a categorical variable
        dat_temp3[, Var_Name := levels[j]] # Setting [Var_Name] to level j of the categorical variable
        
        # Computing the risk-weighted value of the chosen variable
        dat_temp3[, RW_Val := as.numeric(Var_Val==levels[j])*RiskScore/SumScore]
        # Computing the expected value of the chosen variable at each time point and merging it back into the dataset
        dat_temp3 <- merge(dat_temp3, dat_temp3[, list(Exp_Val=sum(RW_Val)), by=list(Time)], by="Time", all.x = TRUE)
        # Subsetting to only include defaulted accounts
        dat_temp4 <- dat_temp3[Status==1,]; rm(dat_temp3)
        # Computing the Schoenfeld residuals
        dat_temp4[, Sch_Res := as.numeric(Var_Val==levels[j]) - Exp_Val]
        
        # Appending the temporary dataset to the main dataset
        dat_return <- rbind(dat_return, dat_temp4[, list(ID, Time, Var_Val, Sch_Res, Var_Name, Var_Name_Base)])
        
        # Updating the counter variable
        k <- k+1
        
        # Plotting if verbose = True
        if (verbose==F){
          if (is.null(max_time)){
            max_time <- max(dat_return$Time)
          } # if
          
          # Modification for faceting to plotting dataset
          dat_temp4[, fac := paste0(Var_Name_Base, " ~ ", Var_Name)]
          
          g <- ggplot(dat_temp4, aes(x=Time, y=Sch_Res)) +
               theme_minimal() + xlab(expression(Default~Time~italic(tau[d]))) + ylab(expression(italic(hat(Sc)(x[k~i~tau[d]])))) +
               theme(text=element_text(family=chosenFont),legend.position = "bottom",
                     strip.background=element_rect(colour="grey", fill="#D3D3D3")) +
               facet_wrap(~ fac, strip.position = "right") +
               geom_point(shape=1, col="#3F702F") + geom_smooth(aes(linetype="solid"), col="#043927", fill="#043927") +
               scale_x_continuous(limits = c(NA, max_time)) +
               scale_linetype_manual(name = "Loess-smoother", values = "solid", labels = NULL)
          g_name <- paste0("SchRes_",var[i], levels[j])  
          gplots[[g_name]] <- g; names(gplots[k])    
        } # if
        
        # Correlation Test
        cor_test[[k]] <- cor.test(dat_temp4$Time, dat_temp4$Sch_Res, method = "spearman"); names(cor_test)[k] <- paste(var[i], levels[j]) # Spearman rank correlation is used for robustness
      }
    } # else if
    
  } # for
  
  dat_return_wider <- dat_return %>% pivot_wider(names_from = c(Var_Name_Base, Var_Name), values_from = c(Var_Val, Sch_Res)) %>% setDT()
  
  # Small correction in naming
  names(dat_return_wider)[grep("Var_Val*", names(dat_return_wider))[names(dat_return_wider)[grep("Var_Val*", names(dat_return_wider))] == paste0("Var_Val", "_", var, "_", var)]] <- paste0("Var_Val_", var[sapply(dat_train2[,..var], is.numeric)])
  names(dat_return_wider)[grep("Sch_Res*", names(dat_return_wider))[names(dat_return_wider)[grep("Sch_Res*", names(dat_return_wider))] == paste0("Sch_Res", "_", var, "_", var)]] <- paste0("Sch_Res_", var[sapply(dat_train2[,..var], is.numeric)])
  
  # Clean up
  suppressWarnings(rm(dat_train2, dat_temp, dat_temp2, dat_temp3, dat_temp4, var, var_type, var_name, col_names, g, g_name, k, levels, levels_n))

  return(list(data = dat_return_wider, CorTest = cor_test, plots = gplots))
  
  # rm(dat_temp, dat_temp2, dat_return, id, var, time, status, max_time, verbose)
} # function

# --- UNIT TEST (Breslow Approximation)
# - Setup
# method <- "breslow"
# cph <- coxph(Surv(time=TimeInPerfSpell-1, time2=TimeInPerfSpell, event=DefaultStatus1==1) ~
#              Principal + Instalment + LN_TPE + slc_pmnt_method,
#              id= PerfSpell_Key, data = dat_train %>% group_by(PerfSpell_Key, TimeInPerfSpell),
#              ties = method)
# cph_bres_scaled_schoenfeld <- residuals(cph, type="schoenfeld") # | Use residuals() for classical Schoenfeld residuals and not Scaled; Schoefeld residuals computed for each level

# - Function Execution
# return <- cph_schoen(cph=cph, id = "PerfSpell_Key", dat_train = dat_train, time = "TimeInPerfSpell", status = "DefaultStatus1", verbose = F, max_time = 240)

# - Comparison of Schoenfeld residuals for a numerical variable
# a <- cph_bres_scaled_schoenfeld[,1]
# b <- return$data[, Sch_Res_Principal]
# ab <- data.table(Sch_Res_Principal_Residuals = a,
#                  Sch_Res_Principal_CustFunc = b)
# all.equal(ab$Sch_Res_Principal_Residuals, ab$Sch_Res_Principal_CustFunc)
### RESULTS:~ TRUE
# - Plot comparison for numeric variable
# par(mfcol=c(1,2)); plot(x=rownames(cph_bres_scaled_schoenfeld), y=cph_bres_scaled_schoenfeld[,1], xlim=c(0,240)); plot(x=return$data$Time, y=return$data$Sch_Res_Principal, xlim=c(0,240))

# - Comparison of Schoenfeld residuals for a categorical variable
# c <- cph_bres_scaled_schoenfeld[,4]
# d <- as.numeric(return$data$`Sch_Res_slc_pmnt_method_Debit Order other bank`)
# cd <- data.table(Sch_Res_slc_pmnt_method_Debit_Order_Other_Bank_Residuals = c,
#                  Sch_Res_slc_pmnt_method_Debit_Order_Other_Bank_CustFunc = d)
# all.equal(cd$Sch_Res_slc_pmnt_method_Debit_Order_Other_Bank_Residuals,cd$Sch_Res_slc_pmnt_method_Debit_Order_Other_Bank_CustFunc)
### RESULTS:~ TRUE
# - Plot comparison for categorical variable
# par(mfcol=c(1,2)); plot(x=rownames(cph_bres_scaled_schoenfeld), y=cph_bres_scaled_schoenfeld[,8], xlim=c(0,240)); plot(x=return$data$Time, y=return$data$Sch_Res_slc_pmnt_method_Suspense, xlim=c(0,240))
