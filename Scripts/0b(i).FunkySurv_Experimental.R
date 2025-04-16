### The following are experimental constructs, though not used (yet).


calc_LPs <- function(X, beta, centered=F, coxGiven=NULL) {
  # - Testing conditions
  # X=X_train; coxGiven=NULL
  
  if (centered==F) {
    lp <- (X %*% beta)
  } else {
    lp_sample <- colMeans(X) %*% beta # assuming only numeric covariates; 0 for categorical
    lp_indiv <- (X %*% beta)
    lp <- lp_indiv - as.numeric(lp_sample)
  }
  
  # internal cox-model based output
  if (!is.null(coxGiven)) {
    lp2 <- predict(coxGiven, type="lp", centered=F)
    
    # Comparisons
    # head(lp); head(coxGiven$linear.predictors)
    # head(coxGiven$linear.predictors - mean(coxGiven$linear.predictors))
    #all.equal(as.numeric(lp), coxGiven$linear.predictors - mean(coxGiven$linear.predictors)) ### TRUE
    lp <- lp2
  }
  
  return(as.numeric(lp))
}


TESTSPACE <- function() {
  
  
  # Select one case from validation set
  test <- subset(datCredit_valid_TFD, LoanID %in% unique(datCredit_valid_TFD[PerfSpell_Age > 5 & PerfSpell_Num > 2,LoanID])[1],
                 select=c("LoanID", "Date", vecVars_TFD, "PerfSpell_Key", "PerfSpell_Num","PerfSpell_Counter","Start", "End", "Default_Ind"))
  
  # - Estimate survival curves
  survFit_test <- survfit(cox_TFD, centered=F, newdata=test, id=PerfSpell_Key)
  
  # - Create survival data object for graphing and comparison purposes | New data
  (datSurv_test_indiv <- data.table(End=survFit_test$time,
                                    CHaz=survFit_test$cumhaz, RiskSetSize=survFit_test$n.risk,
                                    NumEvents=survFit_test$n.event, NumCensored=survFit_test$n.censor,
                                    Survival=round(survFit_test$surv,digits=15)))
  plot(datSurv_test_indiv$Survival, type="b") # contains survival curves of multiple spells, hence why it jumps around
  
  # - Create survival data object for graphing and comparison purposes | No new data
  # NOTE: This creates the baseline hazard function since all covariates are set to 0
  # NOTE: This is useful as a comparison against our own baseline cumulative hazard function, which
  # I suppose must be estimated from the original training set
  survFit_test2 <- survfit(cox_TFD, centered=F)
  # - Create survival data object for graphing and comparison purposes
  (datSurv_test_train <- data.table(End=survFit_test2$time,
                                    CHaz=survFit_test2$cumhaz, RiskSetSize=survFit_test2$n.risk,
                                    NumEvents=survFit_test2$n.event, NumCensored=survFit_test2$n.censor,
                                    Survival=round(survFit_test2$surv,digits=15)))
  plot(datSurv_test_train$Survival, type="b")
  plot(datSurv_test_train$CHaz, type="b")
  
  # Obtain some hidden detail from the fitted CoxPH object
  coxDetail <- coxph.detail(cox_TFD) # Not useful in our context; riskmat=T fails
  
  # - Approximate baseline hazard h_0(t) from cumulative baseline hazard H_0(t)
  # datSurv[, Hazard := c(datSurv$CHaz[1]/datSurv$End[1], diff(datSurv$CHaz) / diff(datSurv$End))]
  #plot(datSurv[, Time], datSurv[, Hazard_base], type="b") # mean survival probability over time
  # datSurv[, EventProb := Hazard * shift(Survival,n=1,type="lag",fill=1)] # f(t|X) = S(t-1|X) . h(t|X)
  
  # - Render survival predictions: risk scores exp(\beta . x) and calculate hazard and survival probability
  #datSurv[, RiskScore := predict(coxGiven, newdata=datGiven, type="risk")]
  #datSurv[, Hazard := Hazard_base * RiskScore] # by definition of Cox regression model
}



survQuants.data <- function(datGiven_train, datGiven_score, vars, beta, fldID="PerfSpell_Key", 
                            fldStart="Start", fldStop="End", fldEvent="Default_Ind",
                            centered = F, ties="Breslow", coxGiven=NULL) {
  # - Testing conditions
  #test <- subset(datCredit_valid_TFD, LoanID %in% unique(datCredit_valid_TFD[PerfSpell_Age > 5 & PerfSpell_Num > 2,LoanID])[1],
                 #select=c("LoanID", "Date", vecVars_TFD, "PerfSpell_Key", "PerfSpell_Num","PerfSpell_Counter","Start", "End", "Default_Ind"))
  # vars <- vecVars_TFD; fldID<-"PerfSpell_Key"; fldStart<-"Start"; fldStop<-"End"; fldEvent<-"Default_Ind"
  # datGiven_train <- datCredit_train_TFD; datGiven_score <- test
  # centered <- F; beta <- cox_TFD$coefficients; coxGiven <- NULL; ties<-"Breslow";
  
  
  # --- Initialising data structures
  
  # - Set specific order since Breslow-estimator is sensitive to ordering
  setDT(datGiven_train, key=c(fldID, fldStop))
  
  # - Coalesce input variables as matrices
  # Use training set for calculating risk scores in the Breslow-estimator in determining baseline cumulative hazard H_0(t),
  # while relegating the to-be-scored set for calculating risk scores in deriving f(t|x)
  X_train <- as.matrix(datGiven_train[, mget(vars)])
  
  # Transform given X to a numeric matrix, encoding nominal/character columns as numeric factors
  X_transform <- matrix(NA, nrow=NROW(X_train), ncol=1)
  jj <- 0 # counter variable of number of new columns added (relating to character vectors)
  for (i in 1:NCOL(X_train)) {
    # i <- 10
    if ( all(!is.na(suppressWarnings(as.numeric(X_train[,i])))) ) {# If true, then numeric
      X_transform[, i+jj] <- as.numeric(X_train[,i])
      suppressWarnings(X_transform <- matrix(X_transform, nrow=NROW(X_train), i+jj+1)) # Grow the dimensions of transformed X
    } else {# If false, then character
      sLvls <- NROW(unique(X_train[,i]))
      # Scan [beta] and position new columns in the correct order corresponding to the coefficients' order
      strtPoint <- jj # Pick up column-creation from last point (if applicable)
      for (j in strtPoint:(strtPoint+sLvls-1-1)) { # Iterate only over (p-1) levels of the nominal variable (given reference cell coding)
        # j <- 0
        sColName <- sub(colnames(X_train)[i], "", names(beta)[i+jj]) # column name extracted
        X_transform[,i+jj] <- ifelse(X_train[,i]==sColName,1,0 )
        suppressWarnings(X_transform <- matrix(X_transform, nrow=NROW(X_train), i+jj+1)) # Grow the dimensions of transformed X
        jj <- jj + 1
      }
      jj <- jj -1 # preventing subscript out of bounds error in next i-iteration
    }
  }
  # Remove added column
  suppressWarnings(X_transform <- matrix(X_transform, nrow=NROW(X_train), i+jj))
  colnames(X_transform) <- names(beta)
  
  # - Extract start and stop times
  vStartTimes_train <- datGiven_train[[fldStart]]
  vStopTimes_train <- datGiven_train[[fldStop]]
  
  # - Extract status and id fields
  vStatus_train <- datGiven_train[[fldEvent]]
  vID_train <- datGiven_train[[fldID]]
  vID_score <- datGiven_score[[fldID]]
  
  
  # --- Estimate cumulative baseline hazard and survival probabilities over unique event times
  
  # Obtain unique failure times for event status
  #vEventTimes <- sort(unique(vStopTimes_train[event_mask]))
  vEventTimes <- sort(unique(vStopTimes_train)) # corresponds with basehaz(), and therefore includes censoring times
  
  # Initialize cumulative baseline hazard & related quantities
  H0 <- numeric(length(vEventTimes))
  vSurv <- numeric(length(vEventTimes))
  riskSetSize <- numeric(length(vEventTimes))
  numEvents <- numeric(length(vEventTimes))
  
  # - Compute linear predictors and risk scores, depending on centering
  vLP_train <- calc_LPs(X=X_transform, beta=beta, centered=centered, coxGiven=coxGiven)
  vRskScores_train <- exp(vLP_train)  # exp(\beta.X)
  
  # Compute cumulative baseline hazard at each failure time using Breslow's method
  for (i in seq_along(vEventTimes)) {
    # i <- 2
    t_i <- vEventTimes[i]
    
    # Identify individuals who had an event at t_i
    vEventIDs <- vID_train[vStopTimes_train == t_i & vStatus_train == 1]
    d_i <- length(vEventIDs) # number of events at t_i
    numEvents[i] <- d_i
    
    # Risk set: sum of exp(Xβ) for individuals still at risk at t_i
    vRskSet_ind <- vStartTimes_train < t_i & vStopTimes_train >= t_i # & !(vID_train %in% vEventIDs & vStopTimes_train < t_i)
    (riskSetSize[i] <- sum(vRskSet_ind)) # same as sum(vStatus_train < t_i) when not handling ids
    (riskSet <- sum(vRskScores_train[vRskSet_ind]) )
    
    # Calculate the Efron adjustment for ties by adjusting the risk set
    if (d_i > 1 & ties=="Efron") {
      # Adjust for ties at the event time using Efron estimator
      tied_contribution <- sum(vRskScores_train[vStartTimes_train < t_i & vStopTimes_train == t_i]) / riskSet
      risk_set <- risk_set + tied_contribution
    }
    
    # Breslow estimator
    (H0[i] <- (ifelse(i > 1, H0[i-1], 0)) + (d_i / riskSet) )
    (vSurv[i] <- ifelse(i== 1, 1, exp(-H0[i])))
  }
  
  # Apply centering correction if requested
  if (centered) {
    #H0 <- H0 * exp(-mean(vLP_train))
  }
  
  # Store results
  datBaseHaz <- data.table(Time = vEventTimes, CHaz = H0, RiskSetSize = riskSetSize, 
                           NumEvents = numEvents, Survival = vSurv)
  # plot(datBaseHaz$Survival)
  # lines(datSurv_test_indiv2$Survival, col="red") # from within TESTSPACE()
  # plot(datBaseHaz$CHaz)
  
  
  # --- Iterate across spell keys and calculate survival-related quantities using survQuants()
  
  # - Preliminaries
  vSpellKeys <- unique(vID_score)
  numSpellKeys <- length(vSpellKeys)
  
  
  # - Iteration logic
  ptm <- proc.time() #IGNORE: for computation time calculation
  #cl.port <- makeCluster(round(6)); registerDoParallel(cl.port) # multi-threading setup
  cat("New Job: Estimating various survival quantities at the loan-period level for a given dataset ..",
      file=paste0(genPath,"survQuants_log.txt"), append=F)
  
  datSurvSet <- foreach(i=1:numSpellKeys, .combine='rbind', .verbose=F, .inorder=T,
                        .packages=c('data.table')) %do%
    { # ----------------- Start of Inner Loop -----------------
      # - Testing conditions
      # i <- 3
      vSubset <- datGiven_score[get(fldID)==vSpellKeys[i],]
      # Get start and stop times from subset, from which point the survival & cumulative hazard shall be calculated
      t <- vSubset[[fldStop]][1]
      t.end <- vSubset[[fldStop]][NROW(vSubset)]
      # Retrieve cumulative baseline hazard given start and stop times
      vBaseCHaz <- datBaseHaz[Time>=t & Time <= t.end, CHaz]
      #vBaseCHaz <- datSurv_test_indiv2[End>=t & End <= t.end, CHaz] # still gives near-flat survival curve
      
      # Calculate risk scores
      X_score <- as.matrix(vSubset[, mget(vars)])
      vLP_score <- calc_LPs(X=X_score, beta=beta, centered=centered)
      vRskScores_score <- exp(vLP_score)  # exp(\beta.X)
      # Calculate survival-related quantities given subject's risk score
      vCHaz_score <- vRskScores_score * vBaseCHaz # cumulative hazard H(t|x)
      vSurv_score <- exp(-vBaseCHaz*vRskScores_score) # survival probability S(t|x)=exp(_H_0(t))^ { exp(\beta . x)}
      vHaz_score <- c(vCHaz_score[1], diff(vCHaz_score)) # hazard h(t|x)
      vEventProb_score <-  vHaz_score * shift(data.table(vSurv_score), n=1, type="lag", fill=1)[[1]] # f(t) = h(t) . S(t-1)
      datResult <- data.table(Key=vSubset[[fldID]], Time=vSubset[[fldStop]], 
                              Survival=as.numeric(vSurv_score), CHaz=as.numeric(vCHaz_score), EventProb=vEventProb_score)
      return(datResult)
    } # ----------------- End of Inner Loop -----------------
  #stopCluster(cl.port); 
  proc.time() - ptm; 
  
  # plot(datSurvSet$Survival, type="b")
  
  return(datSurvSet)
}
