# ============================================ INPUT SPACE =========================================
# Divide data into thematic groups and perform data analysis on them to compile an input space for the TTFD model
# ------------------------------------------------------------------------------------------------------
# PROJECT TITLE: Default survival modelling
# SCRIPT AUTHOR(S): Bernard Scheepers
# ------------------------------------------------------------------------------------------------------
# -- Script dependencies:
#   - 0.Setup.R
#   - 1.Data_Import.R
#   - 2a.Data_Prepare_Credit_Basic.R
#   - 2b.Data_Prepare_Credit_Advanced.R
#   - 2c.Data_Prepare_Credit_Advanced2.R
#   - 2d.Data_Enrich.R
#   - 2f.Data_Fusion1.R
#   - 3b.Data_Enrich2.R

# -- Inputs:
#   - datCredit_train_TFD | Prepared from script 3b
#   - datCredit_valid_TFD | Prepared from script 3b

#
# -- Outputs:
#   - Input_Space
# ------------------------------------------------------------------------------------------------------

# ------ 1. Preliminaries
# - Confirm prepared datasets are loaded into memory
if (!exists('datCredit_train_TFD')) unpack.ffdf(paste0(genPath,"creditdata_train_TFD"), tempPath);gc()
if (!exists('datCredit_valid_TFD')) unpack.ffdf(paste0(genPath,"creditdata_valid_TFD"), tempPath);gc()

### HOMEWORK: Create [Removed], [slc_acct_roll_ever_24_imputed_med_f]
datCredit_train_TFD[,Removed := ifelse(Date==PerfSpell_Max_Date,T,F)]
datCredit_train_TFD[, slc_acct_arr_dir_3_Change_Ind := ifelse(slc_acct_arr_dir_3 != "SAME", 1,0)]
datCredit_train_TFD[,TimeInDelinqState_Lag_1 := shift(TimeInDelinqState,fill=0),by=LoanID]

datCredit_valid_TFD[,Removed := ifelse(Date==PerfSpell_Max_Date,T,F)]
datCredit_valid_TFD[, slc_acct_arr_dir_3_Change_Ind := ifelse(slc_acct_arr_dir_3 != "SAME", 1,0)]
datCredit_valid_TFD[,TimeInDelinqState_Lag_1 := shift(TimeInDelinqState,fill=0),by=LoanID]

datCredit_train_TFD[,slc_acct_roll_ever_24_imputed_med_f := factor(slc_acct_roll_ever_24_imputed_med)]

# Define functions for the analysis
# Spell check
# Ensures that the variables contained in a vector is indeed in a column.
colCheck <- function(cols, dataset) {
  # Check if all columns exist in the dataset
  missing_cols <- cols[!(cols %in% colnames(dataset))]
  
  if (length(missing_cols) == 0) {
    paste0("SAFE: All columns are in the dataset")
  } else {
    warning("WARNING: Some columns are not in the dataset.")
    paste0("Columns not in the dataset: ", paste(missing_cols, collapse = ", "))
  }
}

# Remove variables
# Used to remove and/or add variable names to a vector.
vecChange <- function(Vector,Remove=FALSE,Add=FALSE){
  v <- Vector
  if (all(Remove != FALSE)) {v <- Vector[!(Vector$vars %in% Remove)]} # Remove values in [Remove]
  if (all(Add != FALSE)) {v <- rbind(v,Add[!(Add$vars %in% Vector$vars)])} # Add values in [Add]
  return(v)
}

# Create Correlation function
# Used to detect variables that are highly correlated (> 0.6)
corrAnalysis <- function(data, varlist, corrThresh = 0.6, method = 'spearman') {
  # Compute the correlation matrix
  corrMat <- as.data.table(data) %>% subset(select = varlist) %>% cor(method = method)
  
  # Visualize the correlation matrix
  corrplot(corrMat, type = 'upper', addCoef.col = 'black', tl.col = 'black', diag=FALSE,
           tl.srt = 45)
  
  # Find correlations exceeding the threshold
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

# Table concordance statistic of single variable cox ph models
# Used to compare predictive performance of variables
concTable <- function(data, variables) {
  # Use lapply to efficiently compute concordances for all univariate models
  results <- lapply(variables, function(var) {
    formula <- as.formula(paste0("Surv(Start, End, Default_Ind) ~ ", var))
    tryCatch({
      model <- coxph(formula,id=LoanID, data = data)# Fit Cox model
      conc <- as.numeric(concordance(model)[1])# Extract concordance
      sd <- sqrt(concordance(model)$var)# Extract concordance variability
      lr_stat <- round(2 * (model$loglik[2] - model$loglik[1]),0)# Extract LRT from the model's log-likelihood
      # Return results as a data.table
      data.table(Variable = var, Concordance = conc, SD = sd, LR_Statistic = lr_stat)
    }, warning = function(w) {
      cat("Warning: ", w$message, " for variable: ", var, "\n")
      stop()
    }, error = function(e) {
      cat("Error: ", e$message, " for variable: ", var, "\n")
      stop()
    })
  })
  
  # Combine all results into a single data.table
  results <- rbindlist(results)
  
  # Sort by concordance in descending order
  setorder(results, -Concordance)
  
  return(results)
}

# Table KS statistics of single variable cox ph models
# Used to compare the goodness of fit of variables
csTable <- function(data,variables,seedVal=1,numIt=5){
  
  # Simulate null distribution if seedVal is not NA
  # Initialize results
  results <- data.frame(Variable = variables, KS_Statistic = NA_real_)
  
  # Simulate null distribution if seedVal is not NA
  if (!is.na(seedVal)) {
    set.seed(seedVal, kind = "Mersenne-Twister")
    #null_distribution <- rexp(nrow(data))
    
    # Vectorized calculation for KS statistics
    results$KS_Statistic <- sapply(variables, function(var) {
      formula <- as.formula(paste0("Surv(Start, End, Default_Ind) ~ ", var))
      tryCatch({
        model <- coxph(formula, data = data, id=LoanID)  # Fit a univariate Cox model
        GoF_CoxSnell_KS(model, data, GraphInd=F)$Stat  # Calculate KS statistic
      }, warning = function(w) {
        cat("Warning: ", w$message, " for variable: ", var, "\n")
        NA
      }, error = function(e) {
        cat("Error: ", e$message, " for variable: ", var, "\n")
        NA
      })
    })
    
    # Sort results by KS statistic in descending order
    results <- results[order(-results$KS_Statistic, na.last = TRUE), ]
    
    # Return results and range of KS statistics
    return(list(Table = results, Range = diff(range(results$KS_Statistic, na.rm = TRUE))))
    
  } else {
    # Perform iterative KS calculation when seedVal is NA
    matResults <- matrix(NA, nrow = length(variables), ncol = numIt,
                         dimnames = list(variables,
                                         paste0("Iteration_", 1:numIt))) %>%
      as.data.table()
    
    for (it in seq_len(numIt)) {
      #null_distribution <- rexp(nrow(data))
      
      # Vectorized iteration for KS statistics
      matResults[, it] <- sapply(variables, function(var) {
        formula <- as.formula(paste0("Surv(Start, End, Default_Ind) ~ ", var))
        tryCatch({
          model <- coxph(formula, data = data)
          GoF_CoxSnell_KS(model, data, GraphInd = FALSE)$Stat
        }, warning = function(w) {
          cat("Warning: ", w$message, " for variable: ", var, " in iteration: ", it, "\n")
          NA
        }, error = function(e) {
          cat("Error: ", e$message, " for variable: ", var, " in iteration: ", it, "\n")
          NA
        })
      })
    }
    
    # Compute additional statistics for the results matrix
    colRanges <- matResults[, lapply(.SD, function(x) diff(range(x, na.rm = TRUE)))]
    matResults <- rbind(matResults, Range = colRanges, fill=T)
    matResults[, Average := rowMeans(.SD, na.rm = TRUE), .SDcols = patterns("^Iteration_")]
    
    matResults <- cbind(Variables = c(variables,"Range"),matResults)
    setorder(matResults,-Average)
    
    # Return matrix of KS statistics
    return(matResults)
  }
}
#============================================================================================
# ------ 1. Delinquency measures
varlist <- data.table(vars=c("g0_Delinq","g0_Delinq_fac","PerfSpell_g0_Delinq_Num", "Arrears" ,
                             "TimeInDelinqState","g0_Delinq_Any_Aggr_Prop","g0_Delinq_Ave",
                             "slc_acct_arr_dir_3", "slc_acct_roll_ever_24_imputed_med"),
                      vartypes=c("acc", "cat", "acc", "acc", "acc", "dte", "dte", "cat", "acc"))

#=========================================================================================



# ------ 1.1 Which time window length is the best in calculating Delinquency volatility?

# Initialize variables to be tested
vars <- c("g0_Delinq_SD_4", "g0_Delinq_SD_5", "g0_Delinq_SD_6", "g0_Delinq_SD_9", "g0_Delinq_SD_12")

csTable(datCredit_train_TFD,vars)
#         Variable     KS
# 1  g0_Delinq_SD_4 0.6188
# 2  g0_Delinq_SD_5 0.6124
# 5 g0_Delinq_SD_12 0.6040
# 4  g0_Delinq_SD_9 0.5999
# 3  g0_Delinq_SD_6 0.5918

### RESULTS:  [g0_Delinq_SD_4] fits the data the best, slightly better than [g0_Delinq_SD_5]

concTable(datCredit_valid_TFD,vars)
#             Variable Concordance        SD LR_Statistic
# 1:  g0_Delinq_SD_4   0.9803661 0.001401872        48597
# 2:  g0_Delinq_SD_5   0.9732740 0.001754351        53030
# 3:  g0_Delinq_SD_6   0.9537569 0.002316001        52873
# 4:  g0_Delinq_SD_9   0.9213853 0.002923814        50256
# 5: g0_Delinq_SD_12   0.8885501 0.003291758        42174

### RESULTS: As the SD period increase, there is a slight decrease in concordance.
### NOTE: Concordance is extremely high with low variability

### Conclusion: Larger window are less influenced by large changes therefore significant changes
###             are less pronounced. Include [g0_Delinq_SD_4] in the varlist.

varlist <- vecChange(varlist,Remove="PerfSpell_g0_Delinq_SD",Add=data.table(vars=c("g0_Delinq_SD_4"), vartypes=c("acc")))


# ------ 1.2 Should we add a lag version of [g0_Delinq_Any_Aggr_Prop]?

vars <- c("g0_Delinq_Any_Aggr_Prop","g0_Delinq_Any_Aggr_Prop_Lag_1",
          "g0_Delinq_Any_Aggr_Prop_Lag_2","g0_Delinq_Any_Aggr_Prop_Lag_3",
          "g0_Delinq_Any_Aggr_Prop_Lag_4","g0_Delinq_Any_Aggr_Prop_Lag_5",
          "g0_Delinq_Any_Aggr_Prop_Lag_6","g0_Delinq_Any_Aggr_Prop_Lag_9",
          "g0_Delinq_Any_Aggr_Prop_Lag_12")

# Compare goodness of fit of different variables
csTable(datCredit_train_TFD,vars,seedVal = NA)
#                           Variables Iteration_1 Iteration_2 Iteration_3 Iteration_4 Iteration_5 Average
# 1:  g0_Delinq_Any_Aggr_Prop_Lag_9      0.6511      0.6504      0.6519      0.6452      0.6459 0.64890
# 2:  g0_Delinq_Any_Aggr_Prop_Lag_5      0.6497      0.6461      0.6471      0.6517      0.6494 0.64880
# 3:  g0_Delinq_Any_Aggr_Prop_Lag_2      0.6459      0.6507      0.6491      0.6484      0.6481 0.64844
# 4:        g0_Delinq_Any_Aggr_Prop      0.6491      0.6492      0.6472      0.6463      0.6493 0.64822
# 5:  g0_Delinq_Any_Aggr_Prop_Lag_4      0.6482      0.6495      0.6447      0.6491      0.6462 0.64754
# 6:  g0_Delinq_Any_Aggr_Prop_Lag_6      0.6480      0.6448      0.6475      0.6461      0.6480 0.64688
# 7: g0_Delinq_Any_Aggr_Prop_Lag_12      0.6490      0.6450      0.6489      0.6430      0.6470 0.64658
# 8:  g0_Delinq_Any_Aggr_Prop_Lag_3      0.6475      0.6466      0.6449      0.6447      0.6490 0.64654
# 9:  g0_Delinq_Any_Aggr_Prop_Lag_1      0.6478      0.6448      0.6459      0.6462      0.6475 0.64644
# 10:                          Range      0.0052      0.0059      0.0072      0.0087      0.0035 0.00610

### RESULTS:  After 5 iterations, [g0_Delinq_Any_Aggr_Prop_Lag_9] seems to have the best fit, albeit with 
###           a small range diminishing the validity of results (sampling variability may be present)

# Compare concordance of different variables
concTable(datCredit_valid_TFD,vars)
#                           Variable Concordance        SD  LR_Statistic
# 1:  g0_Delinq_Any_Aggr_Prop_Lag_3   0.5420038 0.004106057          103
# 2:  g0_Delinq_Any_Aggr_Prop_Lag_2   0.5414531 0.004076323          101
# 3:  g0_Delinq_Any_Aggr_Prop_Lag_5   0.5406850 0.004160495           99
# 4:  g0_Delinq_Any_Aggr_Prop_Lag_1   0.5406070 0.004049427           93
# 5: g0_Delinq_Any_Aggr_Prop_Lag_12   0.5397417 0.004163617          100
# 6:  g0_Delinq_Any_Aggr_Prop_Lag_6   0.5392982 0.004160057           99
# 7:  g0_Delinq_Any_Aggr_Prop_Lag_9   0.5388705 0.004186601          101
# 8:  g0_Delinq_Any_Aggr_Prop_Lag_4   0.5384452 0.004112573           96
# 9:        g0_Delinq_Any_Aggr_Prop   0.5379244 0.004031059           81

### RESULTS: [g0_Delinq_Any_Aggr_Prop] has the lowest concordance.
###           Although all concordances are quite close to one another with 
###           a range of 0.004 and low SD's.

### CONCLUSION: Include [g0_Delinq_Any_Aggr_Prop_Lag_3] in the model, since it has
###             the best concordance (disregard csTable due to high variability)

varlist <- vecChange(varlist,Add=data.table(vars=c("g0_Delinq_Any_Aggr_Prop_Lag_3"), vartypes=c("acc")))



# ------ 1.2 Which variables are highly correlated?

# Correlation analysis
corrAnalysis(datCredit_train_TFD, varlist[vartypes!="cat"]$vars, corrThresh = 0.6, method = 'spearman') # Obtain correlation groups

### RESULTS:  1) [g0_Delinq_Any_Aggr_Prop] and [g0_Delinq_Ave] with a correlation of 1
###           2) [PerfSpell_g0_Delinq_Num] and [slc_acct_roll_ever_24_imputed_med]
### NOTE: Group 1) are also highly correlated with [g0_Delinq_Any_Aggr_Prop_Lag_3],
###       which is to be expected.
###           3) [g0_Delinq] and [Arrears]

### CONCLUSION: A single variable from each group must be retained while the rest are removed.



# ------ 1.2.1 Which variable should be kept from group 1) [g0_Delinq_Any_Aggr_Prop] and [g0_Delinq_Ave]

vars <- c("g0_Delinq_Any_Aggr_Prop", "g0_Delinq_Ave")

# Goodness of fit
csTable(datCredit_train_TFD,vars,seedVal = NA)
#                   Variables Iteration_1 Iteration_2 Iteration_3 Iteration_4 Iteration_5 Average
# 1:           g0_Delinq_Ave      0.6482      0.6434      0.6482      0.6493       0.646 0.64702
# 2: g0_Delinq_Any_Aggr_Prop      0.6436      0.6504      0.6475      0.6476       0.644 0.64662
# 3:                   Range      0.0046      0.0070      0.0007      0.0017       0.002 0.00320


### RESULTS: [g0-Delinq_Ave] seems to have the better goodness of fit over 5 iterations.

# Accuracy
concTable(datCredit_valid_TFD,vars)
#                     Variable Concordance          SD LR_Statistic
# 1:           g0_Delinq_Ave   0.5397435 0.004030550           90
# 2: g0_Delinq_Any_Aggr_Prop   0.5379244 0.004031059           81

### RESULTS: [g0-Delinq_Ave] seems to have a slightly better concordance with the concordances have low SD.

### CONCLUSION: [g0-Delinq_Ave] seems to outperform [g0_Delinq_Any_Aggr_Prop] and therefore is kept in the model
###             and [g0_Delinq_Any_Aggr_Prop_Lag_3] is removed along with [g0_Delinq_Any_Aggr_Prop] due to the high correlation
###             I expect similar results.

varlist <- vecChange(varlist,Remove=c("g0_Delinq_Any_Aggr_Prop", "g0_Delinq_Any_Aggr_Prop_Lag_3"))



# ------ 1.2.2 Which variable should be kept from group 2) [PerfSpell_g0_Delinq_Num] and [slc_acct_roll_ever_24_imputed_med]

vars <- c("PerfSpell_g0_Delinq_Num", "slc_acct_roll_ever_24_imputed_med")

# Goodness of fit
csTable(datCredit_train_TFD,vars,seedVal = NA)
#                             Variables Iteration_1 Iteration_2 Iteration_3 Iteration_4 Iteration_5 Average
# 1:           PerfSpell_g0_Delinq_Num      0.6458      0.6463      0.6457      0.6426      0.6443 0.64494
# 2: slc_acct_roll_ever_24_imputed_med      0.6287      0.6307      0.6280      0.6300      0.6304 0.62956
# 3:                             Range      0.0171      0.0156      0.0177      0.0126      0.0139 0.01538

### RESULTS: [PerfSpell_g0_Delinq_Num] seems to have the better goodness of fit over 5 iterations.

# Accuracy
concTable(datCredit_valid_TFD,vars)
#                               Variable Concordance        SD LR_Statistic
# 1:           PerfSpell_g0_Delinq_Num   0.9474644 0.0006456311         4596
# 2: slc_acct_roll_ever_24_imputed_med   0.8588785 0.0031795664        18196

### RESULTS: [PerfSpell_g0_Delinq_Num] seems to have a significant better concordance with the concordances having low SD's.

### CONCLUSION: Keep [PerfSpell_g0_Delinq_Num] in the model and remove [slc_acct_roll_ever_24_imputed_med]

varlist <- vecChange(varlist,Remove="slc_acct_roll_ever_24_imputed_med")



# ------ 1.3 Which version of [g0_Delinq] should be kept in the model?

# [g0_Delinq]
cox <- coxph(Surv(Start,End,Default_Ind) ~ g0_Delinq, id=LoanID, datCredit_train_TFD)
summary(cox);rm(cox)
### RESULTS: Beta is unstable with high variability.

# [g0_Delinq_fac]
cox <- coxph(Surv(Start,End,Default_Ind) ~ g0_Delinq_fac, id=LoanID, datCredit_train_TFD)
summary(cox);rm(cox)
### RESULTS: coef is unstable with high variability.

### INVESTIGATE: Is quasi-complete seperation present?
datCredit_train_TFD[g0_Delinq==3 & Default_Ind==0, .N]
### RESULTS: 0
datCredit_train_TFD[g0_Delinq!=3 & Default_Ind==1, .N]
### RESULTS: 0

### CONCLUSION: Quasi-complete separation is present for [g0_Delinq]=3

varlist <- vecChange(varlist,Remove=c("g0_Delinq","g0_Delinq_fac"))

### INVESTIGATE: Should an indicator version of [g0_Delinq] be included in the model.

# An indicator for when the value is greater than 0
datCredit_train_TFD[,g0_Delinq_Ind := ifelse(g0_Delinq > 0, 1, 0)]
cox <- coxph(Surv(Start,End,Default_Ind) ~ g0_Delinq_Ind, id=LoanID, datCredit_train_TFD)
summary(cox); rm(cox)
### RESULTS: coef is unstable with high variability.
datCredit_train_TFD[,g0_Delinq_Ind := NULL]

### INVESTIGATE: Should an lagged version of [g0_Delinq] be included in the model.

# A lag version for [g0_Delinq]
datCredit_train_TFD[,g0_Delinq_Lag_1 := shift(g0_Delinq,fill=0),by=LoanID]
cox <- coxph(Surv(Start,End,Default_Ind) ~ g0_Delinq_Lag_1, id=LoanID, datCredit_train_TFD)
### RESULTS: exp(coef) tends to Inf
datCredit_train_TFD[,g0_Delinq_Lag_1 := NULL]

# Arrears
cox <- coxph(Surv(Start,End,Default_Ind) ~ Arrears, id=LoanID, datCredit_train_TFD)
summary(cox); rm(cox)
### RESULTS: Obtain a stable model

### CONCLUSION: Unable to add [Delinq_0] to the model, since the various forms'
###             coef is unstable, however, [Arrears] compiles a seemingly stable
###             model, therefore it can serve as a proxy for g0_Delinq.



# ------ 1.4 What is the performance of current thematic variables in univariate models?

# Build thematic model based on remaining delinquency variables.
vars <- c("PerfSpell_g0_Delinq_Num","TimeInDelinqState","g0_Delinq_Ave",
          "slc_acct_arr_dir_3","g0_Delinq_SD_4")  

# Test Goodness of fit
csTable(datCredit_train_TFD,vars)
#                     Variable KS_Statistic
# 3           g0_Delinq_Ave       0.6477
# 1 PerfSpell_g0_Delinq_Num       0.6458
# 5          g0_Delinq_SD_4       0.6151
# 2       TimeInDelinqState           NA
# 4      slc_acct_arr_dir_3           NA

### RESULTS: Fits are close to on another.
### NOTE: [TimeInDelinqState] ran out of iterations and did not converge
### NOTE: [slc_acct_arr_dir_3] exp overflow due to covariates

# ------ 1.4.1 Why does [TimeInDelinqState] not converge?

cox <- coxph(Surv(Start,End,Default_Ind) ~ TimeInDelinqState, datCredit_train_TFD)
### RESULTS:  Beta tends to Inf. After some inspection on the data it relates to all
###           Defaulting events starting in a new delinquency state,
###           i.e. [TimeInDelinqState] = 1. Therefore quasi-complete separation
###           seems to be present.

### INVESTIGATE: WHETHER QUASI-COMPLETE SEPERATION IS PRESENT
datCredit_train_TFD[TimeInDelinqState!=1 & Default_Ind==1, .N]
### RESULTS: 0
datCredit_train_TFD[TimeInDelinqState==1 & Default_Ind==0, .N]
### RESULTS: 167718

### CONCLUSION: Quasi-complete separation seems to be present for [g0_Delinq]=3
###             and should therefore be removed.

# Test a lagged version of TimeInDelinqState
datCredit_train_TFD[,TimeInDelinqState_Lag_1 := shift(TimeInDelinqState,fill=0),
                    by=LoanID] # REMOVE
cox <- coxph(Surv(Start,End,Default_Ind) ~ TimeInDelinqState_Lag_1, id=LoanID,
             datCredit_train_TFD)
summary(cox)
# Concordance= 0.942  (se = 0.002 )
datCredit_valid_TFD[,TimeInDelinqState_Lag_1 := shift(TimeInDelinqState,fill=0),
                    by=LoanID] # REMOVE

### CONCLUSION: Replace old variable with new variable

varlist <- vecChange(varlist,Remove=c("TimeInDelinqState") ,
                     Add=data.table(vars=c("TimeInDelinqState_Lag_1"),
                                    vartypes=c("acc")))


# ------ 1.4.2 Why does [slc_acct_arr_dir_3] exp overflow?

describe(datCredit_train_TFD[,slc_acct_arr_dir_3])
# Value            CURING MISSING_DATA      ROLLING         SAME
# Proportion        0.015        0.126        0.022        0.837

describe(datCredit_train_TFD[Default_Ind==1,slc_acct_arr_dir_3])
# Value            CURING MISSING_DATA      ROLLING         SAME
# Proportion        0.007        0.160        0.783        0.049

### RESULTS:  Although the ROLLING level is not prevalent in datCredit_train_TFD,
###           we can clearly see that it is highly predictive of a default event
###           occurring.

# Create an indicator function
datCredit_train_TFD[, slc_acct_arr_dir_3_ROLLING_Ind := 
                      ifelse(slc_acct_arr_dir_3 == "ROLLING", 1,0)]
cox <- coxph(Surv(Start,End,Default_Ind) ~ slc_acct_arr_dir_3_ROLLING_Ind,
             datCredit_train_TFD)
### RESULTS: exp overflows

# Make the indicator categorical
cox <- coxph(Surv(Start,End,Default_Ind) ~ factor(slc_acct_arr_dir_3_ROLLING_Ind),
             datCredit_train_TFD)
### RESULTS: exp overflows
datCredit_train_TFD[,slc_acct_arr_dir_3_ROLLING_Ind := NULL]

# Create an indicator variable for a change in account
datCredit_train_TFD[, slc_acct_arr_dir_3_Change_Ind := 
                      ifelse(slc_acct_arr_dir_3 != "SAME", 1,0)] # REMOVE
cox <- coxph(Surv(Start,End,Default_Ind) ~ slc_acct_arr_dir_3_Change_Ind,
             datCredit_valid_TFD)
summary(cox)
# Concordance= 0.843
datCredit_valid_TFD[, slc_acct_arr_dir_3_Change_Ind :=
                      ifelse(slc_acct_arr_dir_3 != "SAME", 1,0)] # REMOVE

### CONCLUSION: Replace old variable with new variable

varlist <- vecChange(varlist,Remove=c("slc_acct_arr_dir_3") ,
                     Add=data.table(vars=c("slc_acct_arr_dir_3_Change_Ind"),
                                    vartypes=c("bin")))



# ------ 1.4.3 What is the performance of current thematic variables in univariate models?

vars <- c("PerfSpell_g0_Delinq_Num","g0_Delinq_Ave","g0_Delinq_SD_4",
          "TimeInDelinqState_Lag_1","slc_acct_arr_dir_3_Change_Ind", "Arrears")

# Goodness of fit
csTable(datCredit_train_TFD,vars)
#                         Variable KS_Statistic
# 2                 g0_Delinq_Ave       0.6477
# 6                       Arrears       0.6476
# 1       PerfSpell_g0_Delinq_Num       0.6458
# 5 slc_acct_arr_dir_3_Change_Ind       0.6446
# 4       TimeInDelinqState_Lag_1       0.6331
# 3                g0_Delinq_SD_4       0.6151

### RESULTS: [g0_Delinq_SD_4] seems to have a notable worse fit than the other variables.

# Accuracy
concTable(datCredit_valid_TFD,vars)
#                         Variable Concordance           SD LR_Statistic
# 1:                g0_Delinq_SD_4   0.9803661 0.0014018721        48597
# 2:                       Arrears   0.9573776 0.0016578273         2438
# 3:       PerfSpell_g0_Delinq_Num   0.9474644 0.0006456311         4596
# 4: slc_acct_arr_dir_3_Change_Ind   0.8433911 0.0019345216        18802
# 5:       TimeInDelinqState_Lag_1   0.7507071 0.0056323114        12757
# 6:                 g0_Delinq_Ave   0.5397435 0.0040305498           90

### RESULTS: [g0_Delinq_Ave] has significant less concordance than the other variables

### CONCLUSION: Leave all variables in the model (including [g0_Delinq_Ave], since it has the best fit for the data).



# ------ 1.5 What is the performance of current thematic cox ph model?

# Goodness of fit of coxDelinq model

# Build cox model based on all thematic variables
coxDelinq_train <- coxph(Surv(Start,End,Default_Ind) ~ g0_Delinq_SD_4 + g0_Delinq_Ave +
                           PerfSpell_g0_Delinq_Num + slc_acct_arr_dir_3_Change_Ind +
                           TimeInDelinqState_Lag_1 + Arrears, id=LoanID,
                         data=datCredit_train_TFD)

# Kolmogorov-Smirnof of coxDelinq
cs_ks_test(coxDelinq_train,datCredit_train_TFD,GraphInd = FALSE) # 0.6171

# Accuracy
coxDelinq_valid<- coxph(Surv(Start,End,Default_Ind) ~ g0_Delinq_SD_4 + g0_Delinq_Ave +
                          PerfSpell_g0_Delinq_Num + slc_acct_arr_dir_3_Change_Ind +
                          TimeInDelinqState_Lag_1, id=LoanID,
                        data=datCredit_valid_TFD)

# (0,3) (4,12) (13,24) (0,12) (0,36)
#timedROC(datCredit_valid_TFD, coxDelinq_valid, month_Start=0, month_End=36,
         # fld_ID="LoanID", fld_Event="Default_Ind",fld_StartTime="Start",
         # fld_EndTime="End", numDigits=0, Graph=FALSE)
# AUC: 0.9760028
#HW: combine ggplot objects into graph facets

### RESULTS: Graph is extremely accurate with a AUC of 97.6%

# # ------ 1.6 Compare results with that of an algorithm
# 
# # My final cox model
# myCox <- coxph(Surv(Start, End, Default_Ind) ~ g0_Delinq_Ave)
# 
# # Cox model with all the other variables.
# Fullcox <- coxph(Surv(Start, End, Default_Ind) ~
#                    PerfSpell_g0_Delinq_Num + g0_Delinq_Any_Aggr_Prop +
#                    g0_Delinq_Ave + g0_Delinq_SD_4 +  g0_Delinq_SD_5 + 
#                    g0_Delinq_Any_Aggr_Prop + g0_Delinq_Any_Aggr_Prop_Lag_1 + 
#                    g0_Delinq_Any_Aggr_Prop_Lag_2 + g0_Delinq_Any_Aggr_Prop_Lag_3 + 
#                    g0_Delinq_Any_Aggr_Prop_Lag_5, datCredit_train_TFD)
# 
# compCox <- stepAIC(Fullcox, direction = "backward", trace = TRUE)
# #                                 Df    AIC
# # <none>                              98216
# # - g0_Delinq_Any_Aggr_Prop_Lag_5  1  98218
# # - g0_Delinq_Any_Aggr_Prop_Lag_1  1  98226
# # - g0_Delinq_SD_4                 1  98235
# # - PerfSpell_g0_Delinq_Num        1  98601
# # - g0_Delinq_SD_5                 1 105649
#============================================================================================

ModelVar <- c("PerfSpell_g0_Delinq_Num", "Arrears", "g0_Delinq_Ave", "g0_Delinq_SD_4",
              "TimeInDelinqState_Lag_1", "slc_acct_arr_dir_3_Change_Ind")

#============================================================================================
# ------ 2. Engineered measures
varlist <- data.table(vars=c("slc_acct_pre_lim_perc_imputed_med",
                             "slc_acct_prepaid_perc_dir_12_imputed_med",
                             "slc_pmnt_method"),
                      vartypes=c("prc", "dec", "cat"))

#=========================================================================================

# ------ 2.1 Which variables should be removed due to high correlation?

# - Correlation analysis
corGroups <- corrAnalysis(datCredit_train_TFD, varlist[vartypes!="cat"]$vars, corrThresh = 0.6, method = 'spearman') # Obtain correlation groups

### RESULTS:  1) [slc_acct_pre_lim_perc_imputed_med] and [slc_acct_prepaid_perc_dir_12_imputed_med] with a correlation of 0.82

# Initialize variables to be tested
vars <- c("slc_acct_pre_lim_perc_imputed_med", "slc_acct_prepaid_perc_dir_12_imputed_med")

# Compare goodness of fit of different variables
csTable(datCredit_train_TFD,vars)
#                                 Variable     KS
# 1        slc_acct_pre_lim_perc_imputed_med 0.6485
# 2 slc_acct_prepaid_perc_dir_12_imputed_med 0.6444

### RESULTS:  [slc_acct_prepaid_perc_dir_12_imputed_med] ran out of iterations and did not converge

### INVESTIGATE: Why did slc_acct_prepaid_perc_dir_12_imputed_med] ran out of iterations and not converge?

hist(datCredit_train_TFD[,slc_acct_prepaid_perc_dir_12_imputed_med])
### RESULTS: A significant portion of the values are 0 with possible extreme values.

# Proportion of values being 0
datCredit_train_TFD[slc_acct_prepaid_perc_dir_12_imputed_med==0,.N]/datCredit_train_TFD[,.N]
# 0.734379

hist(datCredit_train_TFD[Default_Ind==1, slc_acct_prepaid_perc_dir_12_imputed_med])
### Similar distribution with most of the values being 0, but all values are below 3.5

# Proportion of values being 0 | Default_Ind == 1
datCredit_train_TFD[slc_acct_prepaid_perc_dir_12_imputed_med==0 & Default_Ind == 1,.N]/datCredit_train_TFD[Default_Ind == 1,.N]
# 0.9980394

# Proportion of values defaulted | slc_acct_prepaid_perc_dir_12_imputed_med == 0
datCredit_train_TFD[slc_acct_prepaid_perc_dir_12_imputed_med==0 & Default_Ind == 1,.N]/datCredit_train_TFD[slc_acct_prepaid_perc_dir_12_imputed_med==0,.N]
# 0.002991395

# RESULTS: (Default_Ind == 1) => (slc_acct_prepaid_perc_dir_12_imputed_med == 0)

### CONCLUSION: Quasi-complete separation seems to be present for [slc_acct_prepaid_perc_dir_12_imputed_med],
###             therefore remove it from the varlist

varlist <- vecChange(varlist,Remove="slc_acct_prepaid_perc_dir_12_imputed_med")



# ------ 2.2 What is the predictive power of the variables left in varlist?

concTable(datCredit_valid_TFD,vars)
### RESULTS: [slc_pmnt_method] Loglik converged before variable  1,2,3,4,5,6 ; beta may be infinite.

### INVESTIGATE: Why does [slc_pmnt_method] have an infinite beta.

describe(datCredit_valid_TFD$slc_pmnt_method)
# Value       "Debit Order" "Debit Order FNB account"  "Debit Order other bank"  "MISSING_DATA" "Salary" "Statement" "Suspense"
# Proportion  0.000                   0.492                   0.191                   0.129       0.060    0.124        0.002

### "Debit Order" has extremely few values which might be a problem

describe(datCredit_valid_TFD[Default_Ind==1,slc_pmnt_method])
# Value       "Debit Order FNB account"  "Debit Order other bank"  "MISSING_DATA" "Salary" "Statement" "Suspense"
# Proportion          0.160                       0.094                  0.144      0.044      0.438    0.120

### RESULTS: Statement makes up a significant proportion of defaulted cases compared to the population and "Debit Order" does not exists.

# ------ 2.2.1 Can [slc_pmnt_method] be replaced with [pmnt_method_grp]?

# Distribution of [slc_pmnt_method]
describe(datCredit_train_TFD[,slc_pmnt_method])
# Value       "Debit Order" "Debit Order FNB account"  "Debit Order other bank"  "MISSING_DATA" "Salary" "Statement" "Suspense"
# Proportion  0.000                   0.510                   0.191                   0.135      0.063    0.100        0.001

# Distribution of [pmnt_method_grp]
describe(datCredit_train_TFD[,pmnt_method_grp])
# Value       "Debit Order"    "MISSING_DATA" "Salary/Suspense" "Statement"
# Proportion      0.701           0.135           0.064           0.100

### CONCLUSION: [pmnt_method_grp] is a refined grouping of [slc_pmnt_method]

varlist <- vecChange(varlist,Remove="slc_pmnt_method",Add=data.table(vars=c("pmnt_method_grp"), vartypes=c("cat")))

# ------ 2.2.3 What is the predictive power of the variables left in varlist?

# Initialize variables to be tested
vars <- c("pmnt_method_grp", "slc_acct_pre_lim_perc_imputed_med")

csTable(datCredit_train_TFD,vars)
#                               Variable KS_Statistic
# 2 slc_acct_pre_lim_perc_imputed_med       0.6477
# 1                   pmnt_method_grp       0.6453

### RESULTS:  [slc_acct_pre_lim_perc_imputed_med] fits the better and there is an improvement for [pmnt_method_grp]

concTable(datCredit_valid_TFD,vars)
#                             Variable Concordance           SD LR_Statistic
# 1:                   pmnt_method_grp   0.7363669 0.0036080644         5777
# 2: slc_acct_pre_lim_perc_imputed_med   0.6360895 0.0008227257         2865

### RESULTS: [pmnt_method_grp]  has a much better concordance than [slc_acct_pre_lim_perc_imputed_med]

### Conclusion: All values should be kept in the model ([slc_pmnt_method] should be replaced with [pmnt_method_grp])

varlist <- vecChange(varlist,Remove="slc_pmnt_method",Add=data.table(vars=c("pmnt_method_grp"), vartypes=c("cat")))

# ------ 2.3 How predictive is a single model based on the thematic variables?

vars <- c("pmnt_method_grp", "slc_acct_pre_lim_perc_imputed_med")

# Build Cox model based on each variable
coxEngineered_train <- coxph(Surv(Start, End, Default_Ind) ~ pmnt_method_grp +
                               slc_acct_pre_lim_perc_imputed_med, id=LoanID,
                             data=datCredit_train_TFD)
summary(coxEngineered_train)

### RESULTS: [slc_acct_pre_lim_perc_imputed_med]  has an insignificant coef of -154

# Goodness of fit
GoF_CoxSnell_KS(coxEngineered_train,datCredit_train_TFD) # 0.6407

coxEngineered_valid <- coxph(Surv(Start, End, Default_Ind) ~ pmnt_method_grp +
                               slc_acct_pre_lim_perc_imputed_med, id=LoanID,
                             data=datCredit_valid_TFD)

summary(coxEngineered_valid)
# Concordance= 0.784  (se = 0.003 ) which improved from univariate [pmnt_method_grp] model.

# Time dependent AUC
timedROC(datCredit_valid_TFD, coxEngineered_valid, month_Start=0, month_End=36,
         fld_ID="LoanID", fld_Event="Default_Ind",fld_StartTime="Start",
         fld_EndTime="End", numDigits=0, Graph=FALSE) # 0.7040435

# House keeping
rm(coxEngineered_valid, coxEngineered_valid); gc()

#============================================================================================
# Save variables to the model
ModelVar <- vecChange(ModelVar,Add=data.table(vars=c("slc_acct_pre_lim_perc_imputed_med",
                                                     "pmnt_method_grp"),
                                              vartypes=c("prc","cat")))

#============================================================================================
# ------ 3. Interest Rate
varlist <- data.table(vars=c("InterestRate_Nom", "InterestRate_Margin"),
                      vartypes=c("prc","prc"))

#=========================================================================================

# ------ 3.1 How correlated are the two variables?

# - Correlation analysis
corrAnalysis(datCredit_train_TFD, varlist[vartypes!="cat"]$vars, corrThresh = 0.6, method = 'spearman') # Obtain correlation groups

### RESULTS: Correlation of 0.41, therefore no significant correlations.

# Obtain the KS-statistics for goodness of fit tests.
csTable(datCredit_train_TFD, varlist$vars, seedVal=NA)
#             Variables Iteration_1 Iteration_2 Iteration_3 Iteration_4 Iteration_5 Average
# 1:    InterestRate_Nom      0.6448      0.6499      0.6499      0.6503      0.6483 0.64864
# 2: InterestRate_Margin      0.6461      0.6481      0.6464      0.6474      0.6520 0.64800
# 3:               Range      0.0013      0.0018      0.0035      0.0029      0.0037 0.00264

### RESULTS:  [InterestRate_Nom] is a better fit than [InterestRate_Margin]

# Obtain the concordance for different variables.
concTable(datCredit_valid_TFD, varlist$vars)
#                 Variable Concordance LR_Statistic
# 1:    InterestRate_Nom   0.5497724          165
# 2: InterestRate_Margin   0.5452677          174

### RESULTS:  [InterestRate_Nom] has a slightly higher concordance than [InterestRate_Margin].

### CONCLUSION: Variables are quite similar in predictave power, however are dissimilar
###              in correlation, therefore both can be kept in the model.



# ------ 3.2 Will the aggregate portfolio mean for [InterestRate_Margin] perform better?

vars <- c("InterestRate_Margin", "InterestRate_Margin_Aggr_Med")

# Obtain the KS-statistics for goodness of fit tests.
csTable(datCredit_train_TFD, vars, seedVal=NA)
# Variables Iteration_1 Iteration_2 Iteration_3 Iteration_4 Iteration_5 Average
# 1: InterestRate_Margin_Aggr_Med      0.6468      0.6457      0.6485      0.6510      0.6492 0.64824
# 2:          InterestRate_Margin      0.6484      0.6455      0.6476      0.6469      0.6453 0.64674
# 3:                        Range      0.0016      0.0002      0.0009      0.0041      0.0039 0.00214

### RESULTS:  [InterestRate_Margin_Aggr_Med] both values have relatively good fits.

# Obtain the concordance for different variables.
concTable(datCredit_valid_TFD, vars)
#                        Variable Concordance          SD LR_Statistic
# 1: InterestRate_Margin_Aggr_Med   0.5506784 0.004195174           96
# 2:          InterestRate_Margin   0.5452677 0.004064693          174

### RESULTS:  [InterestRate_Margin_Aggr_Med] has a slightly higher concordance than [InterestRate_Margin],
###           although both are quite low.

### CONCLUSIONS: Replace InterestRate_Margin with [InterestRate_Margin_Aggr_Med]

varlist <- vecChange(varlist,Remove="InterestRate_Margin",Add=data.table(vars=c("InterestRate_Margin_Aggr_Med"), vartypes=c("prc")))



# ------ 3.2 Should we add lagging variables to the model?

vars <- c("InterestRate_Margin_Aggr_Med","InterestRate_Margin_Aggr_Med_1",
          "InterestRate_Margin_Aggr_Med_12","InterestRate_Margin_Aggr_Med_2",
          "InterestRate_Margin_Aggr_Med_3","InterestRate_Margin_Aggr_Med_4",
          "InterestRate_Margin_Aggr_Med_5","InterestRate_Margin_Aggr_Med_6",
          "InterestRate_Margin_Aggr_Med_9")

# Obtain the KS-statistics for goodness of fit tests.
csTable(datCredit_train_TFD, vars, seedVal=NA)
# Variables Iteration_1 Iteration_2 Iteration_3 Iteration_4 Iteration_5 Average
# 1: InterestRate_Margin_Aggr_Med_12      0.6449      0.6520      0.6510      0.6472      0.6511 0.64924
# 2:  InterestRate_Margin_Aggr_Med_3      0.6500      0.6484      0.6492      0.6491      0.6492 0.64918
# 3:    InterestRate_Margin_Aggr_Med      0.6479      0.6474      0.6485      0.6501      0.6497 0.64872
# 4:  InterestRate_Margin_Aggr_Med_9      0.6504      0.6476      0.6471      0.6482      0.6495 0.64856
# 5:  InterestRate_Margin_Aggr_Med_1      0.6477      0.6483      0.6469      0.6506      0.6480 0.64830
# 6:  InterestRate_Margin_Aggr_Med_5      0.6499      0.6455      0.6463      0.6459      0.6507 0.64766
# 7:  InterestRate_Margin_Aggr_Med_2      0.6461      0.6460      0.6453      0.6482      0.6490 0.64692
# 8:  InterestRate_Margin_Aggr_Med_6      0.6464      0.6457      0.6462      0.6475      0.6466 0.64648
# 9:  InterestRate_Margin_Aggr_Med_4      0.6481      0.6468      0.6434      0.6478      0.6448 0.64618
# 10:                           Range      0.0055      0.0065      0.0076      0.0047      0.0063 0.00612

### RESULTS:  Difficult to make conclusive results with variability in KS-statistics.

# Obtain the concordance for different variables.
concTable(datCredit_valid_TFD, vars)
#                           Variable Concordance          SD LR_Statistic
# 1:  InterestRate_Margin_Aggr_Med_4   0.5611257 0.004199905          104
# 2:  InterestRate_Margin_Aggr_Med_3   0.5603674 0.004199684          105
# 3:  InterestRate_Margin_Aggr_Med_2   0.5575290 0.004200634          102
# 4:  InterestRate_Margin_Aggr_Med_5   0.5566768 0.004176251           99
# 5:  InterestRate_Margin_Aggr_Med_6   0.5566370 0.004185777           98
# 6:  InterestRate_Margin_Aggr_Med_1   0.5555107 0.004193798           99
# 7:  InterestRate_Margin_Aggr_Med_9   0.5546448 0.004170861           85
# 8:    InterestRate_Margin_Aggr_Med   0.5506784 0.004195174           96
# 9: InterestRate_Margin_Aggr_Med_12   0.5481947 0.004155018           68

### RESULTS:  [InterestRate_Margin_Aggr_Med_4] seems to have the best concordance, 
###           with [InterestRate_Margin_Aggr_Med_3] falling close behind.

### CONCLUSION: Add [InterestRate_Margin_Aggr_Med_3] to the model since it seems
###             to have the second best goodness of fit and concordance.

varlist <- vecChange(varlist,Add=data.table(vars=c("InterestRate_Margin_Aggr_Med_3"),
                                            vartypes=c("prc")))



# ------ 3.3 What is the performance of current thematic variables in univariate models?

vars <- c("InterestRate_Nom", "InterestRate_Margin_Aggr_Med", "InterestRate_Margin_Aggr_Med_3")

# Goodness of fit
csTable(datCredit_train_TFD,vars)
#                           Variable KS_Statistic
# 1               InterestRate_Nom       0.6479
# 2   InterestRate_Margin_Aggr_Med       0.6477
# 3 InterestRate_Margin_Aggr_Med_3       0.6449

### RESULTS: Variables seem to be a good fit.

# Accuracy
concTable(datCredit_valid_TFD,vars)
#                         Variable Concordance          SD LR_Statistic
# 1: InterestRate_Margin_Aggr_Med_3   0.5603674 0.004199684          105
# 2:   InterestRate_Margin_Aggr_Med   0.5506784 0.004195174           96
# 3:               InterestRate_Nom   0.5497724 0.004486873          165

### RESULTS: Concordances are quite low, this may lead to not including the variables into the model

### CONCLUSION: Take note of low concordances, but evaluate the concordance of the full model.



# ------ 3.4 What is the performance of current thematic cox ph model?

# Build cox model based on all thematic variables
coxInterest_train <- coxph(Surv(Start,End,Default_Ind) ~ InterestRate_Nom +
                           InterestRate_Margin_Aggr_Med + InterestRate_Margin_Aggr_Med_3,
                         id=LoanID, data=datCredit_train_TFD)

# Kolmogorov-Smirnof of coxDelinq
cs_ks_test(coxInterest_train,datCredit_train_TFD,GraphInd = FALSE) # 0.646

# Accuracy
coxInterest_valid<- coxph(Surv(Start,End,Default_Ind) ~ InterestRate_Nom +
                            InterestRate_Margin_Aggr_Med + InterestRate_Margin_Aggr_Med_3, id=LoanID,
                        data=datCredit_valid_TFD)
summary(coxInterest_valid)
# Concordance= 0.571

timedROC(datCredit_valid_TFD, coxInterest_valid, month_Start=0, month_End=36,
         fld_ID="LoanID", fld_Event="Default_Ind",fld_StartTime="Start",
         fld_EndTime="End", numDigits=0, Graph=FALSE)
# AUC: 0.5214215
#===========================================================================================

### CONCLUSION: Low concordance (0.571) and AUC implies weak accuracy, therefore do  
###             not include variables in final model.

#============================================================================================
# ------ 4. General
varlist <- data.table(vars=c("Balance","Instalment",
                             "Principal","Undrawn_Amt",
                             "LN_TPE","Term"),
                      vartypes=c("int", "fin", "fin","int","cat", "int"))

#=========================================================================================
# ------ 4.1 Which variables are highly correlated in the group?

# - Correlation analysis
corrAnalysis(datCredit_train_TFD, varlist[vartypes!="cat"]$vars, corrThresh = 0.6, method = 'spearman') # Obtain correlation groups

### RESULTS: 1) [Balance], [Installment] and [Principal] are highly correlated.

# ------ 4.2 Which variable in the correlated group should be kept in the model?

# Initialize variable
vars <- c("Balance", "Instalment", "Principal")

# Goodness of fit test of variables
csTable(datCredit_train_TFD, vars, seedVal=NA)
#     Variables Iteration_1 Iteration_2 Iteration_3 Iteration_4 Iteration_5 Average
# 1: Instalment      0.6469      0.6492      0.6462      0.6506      0.6491 0.64840
# 2:    Balance      0.6485      0.6453      0.6463      0.6501      0.6482 0.64768
# 3:  Principal      0.6510      0.6434      0.6471      0.6472      0.6478 0.64730
# 4:      Range      0.0041      0.0058      0.0009      0.0034      0.0013 0.00310

### RESULTS:  From the 5 iterations, [Instalment] seems to have the best fit.

# Obtain the concordance for different variables.
concTable(datCredit_valid_TFD, vars)
#     Variable Concordance          SD LR_Statistic
# 1:  Principal   0.5827355 0.004098595          335
# 2: Instalment   0.5636894 0.004286208           84
# 3:    Balance   0.5544691 0.004145147          101

## RESULTS: [Principal] has the highest concordance

### CONCOLUSION:  Considering Principal has the best concordance, but the worse fit
###               complicates the decision. However, [InstalmentToBalance_Aggr_Prop]
###               contains information of both.

# ------ 4.2.1 Can [Balance] and [Instalment] be replaced with [InstalmentToBalance_Aggr_Prop]?

# Initialize variable
vars <- c("Balance", "Instalment", "InstalmentToBalance_Aggr_Prop")

# Goodness of fit test of variables
csTable(datCredit_train_TFD, vars, seedVal=NA)
#                         Variables Iteration_1 Iteration_2 Iteration_3 Iteration_4 Iteration_5 Average
# 1:                    Instalment      0.6478      0.6507      0.6512      0.6501      0.6495 0.64986
# 2: InstalmentToBalance_Aggr_Prop      0.6510      0.6501      0.6489      0.6497      0.6482 0.64958
# 3:                       Balance      0.6464      0.6495      0.6457      0.6483      0.6478 0.64754
# 4:                         Range      0.0046      0.0012      0.0055      0.0018      0.0017 0.00296

### RESULTS:  From the 5 iterations, [InstalmentToBalance_Aggr_Prop] consistently seems to outfit [Balance]

# Accuracy of variables
concTable(datCredit_valid_TFD, vars)
#                         Variable Concordance          SD LR_Statistic
# 1:                    Instalment   0.5636894 0.004286208           84
# 2:                       Balance   0.5544691 0.004145147          101
# 3: InstalmentToBalance_Aggr_Prop   0.5146318 0.004341557           18

### RESULTS:  [InstalmentToBalance_Aggr_Prop] has a worse concordance than [Instalment] and [Balance].

### CONCLUSION: [InstalmentToBalance_Aggr_Prop] cannot replace [Instalment] and [Balance].

# ------ 4.2.1 Can [Balance] be replaced with [ArrearsToBalance_Aggr_Prop]?

# Initialize variable
vars <- c("Balance", "ArrearsToBalance_Aggr_Prop")

# Goodness of fit test of variables
csTable(datCredit_train_TFD, vars, seedVal=NA)
#                     Variables Iteration_1 Iteration_2 Iteration_3 Iteration_4 Iteration_5 Average
# 1:                    Balance      0.6475      0.6519      0.6479      0.6482      0.6483 0.64876
# 2: ArrearsToBalance_Aggr_Prop      0.6457      0.6440      0.6480      0.6446      0.6470 0.64586
# 3:                      Range      0.0018      0.0079      0.0001      0.0036      0.0013 0.00294

### RESULTS:  From the 5 iterations, [Balance] seems to have the best fit.

# Test for accuracy
concTable(datCredit_valid_TFD, vars)
#                     Variable Concordance          SD LR_Statistic
# 1:                    Balance   0.5544691 0.004145147          101
# 2: ArrearsToBalance_Aggr_Prop   0.5367480 0.004126016          118

## RESULTS: [Balance] has the highest concordance.

### CONCOLUSION:  [ArrearsToBalance_Aggr_Prop] cannot reaplace [Balance]

### CONCLUSION: Remove [Instalment] and [Balance].

varlist <- vecChange(varlist,Remove=c("Instalment","Balance"))

# ------ 4.3 How does [Term] compare to [AgeToTerm_Aggr_Mean]?

# Initialize variables
vars <- c("Term", "AgeToTerm_Aggr_Mean")

# Goodness of fit
csTable(datCredit_train_TFD,vars,seedVal = NA)
#               Variables Iteration_1 Iteration_2 Iteration_3 Iteration_4 Iteration_5 Average
# 1: AgeToTerm_Aggr_Mean      0.6462      0.6514      0.6487      0.6471      0.6509 0.64886
# 2:                Term      0.6455      0.6458      0.6486      0.6455      0.6465 0.64638
# 3:               Range      0.0007      0.0056      0.0001      0.0016      0.0044 0.00248

### RESULTS: [AgeToTerm_Aggr_Mean] seems to concistently outfit [Term].

# Accuracy
concTable(datCredit_valid_TFD,vars)
#               Variable Concordance          SD LR_Statistic
# 1: AgeToTerm_Aggr_Mean   0.5091579 0.003670705           24
# 2:                Term   0.5066915 0.002188255            2

### RESULTS: Both have poor concordance.

### CONCLUSION: Remove [Term].

varlist <- vecChange(varlist,Remove=c("Term"))

# ------ 4.4 What is the predictive power of the variables in univariate models?
#============================================================================================
vars <- c("Principal", "Undrawn_Amt", "LN_TPE")

# Goodness of fit
csTable(datCredit_train_TFD,vars)
#       Variable KS_Statistic
# 1   Principal       0.6479
# 2 Undrawn_Amt       0.6477
# 3      LN_TPE       0.6449

### RESULTS: [LN_TPE] appears to have the weakest goodness of fit.

# Accuracy
concTable(datCredit_valid_TFD,vars)
#       Variable Concordance           SD LR_Statistic
# 1: Undrawn_Amt   0.6883511 0.0006726495         2126
# 2:   Principal   0.5827355 0.0040985952          335
# 3:      LN_TPE   0.5121296 0.0018064070           29

### RESULTS:  [Undrawn_Amt] appears to have significant more concordance than the other variables.
###           [LN_TPE] could possibly be removed.

# ------ 4.4 What is the predictive power of the variables still in a single model?

# Goodness of fit
coxGen_train <- coxph(as.formula(paste0("Surv(Start,End,Default_Ind) ~ ",
                                  paste(vars,collapse=" + "))), id=LoanID,
                      datCredit_train_TFD)


GoF_CoxSnell_KS(coxGen_train,datCredit_train_TFD,GraphInd = FALSE) # 0.6485

# Accuracy
coxGen_valid <- coxph(as.formula(paste0("Surv(Start,End,Default_Ind) ~ ",
                                        paste(vars,collapse=" + "))), id=LoanID,
                      datCredit_valid_TFD)

concordance(coxGen_valid) # 0.6882

# House keeping
rm(coxGen_train, coxGen_valid)


### CONCLUSION: Leave all variables in the model.

#============================================================================================









#============================================================================================
# ------ 5. Portfolio Level
varlist <- data.table(vars=c("CuringEvents_Aggr_Prop",
                             "NewLoans_Aggr_Prop"),
                      vartypes=c("dec", "dec", "dec", "dec", "dec"))

#============================================================================================





#============================================================================================
# ------ 6. Macro Economic
varlist <- data.table(vars=c("M_DTI_Growth","M_Emp_Growth","M_Inflation_Growth",
                             "M_RealGDP_Growth","M_RealIncome_Growth","M_Repo_Rate"),
                      vartypes=c("prc", "prc", "prc", "prc", "prc", "prc"))
#============================================================================================



# ------ 6.1 Which economic variables are highly correlated with one another?

# - Correlation analysis
corrAnalysis(datCredit_train_TFD, varlist$vars[varlist$vartypes != 'cat'],
             corrThresh = 0.6, method = 'spearman') # Obtain correlation groups

### RESULTS: [M_Emp_Growth], [M_RealGDP_Growth] and [M_RealIncome_Growth] are highly correlated

# ------ 6.2 Which correlated variable should remain in the model?

# Initialize variables
vars <- c("M_Emp_Growth", "M_RealGDP_Growth", "M_RealIncome_Growth")

# Obtain the KS-statistics for goodness of fit tests.
csTable(datCredit_train_TFD, vars, seedVal=NA)
#             Variables Iteration_1 Iteration_2 Iteration_3 Iteration_4 Iteration_5 Average
# 1:        M_Emp_Growth      0.6469      0.6484      0.6468      0.6498      0.6474 0.64786
# 2:    M_RealGDP_Growth      0.6473      0.6463      0.6516      0.6449      0.6489 0.64780
# 3: M_RealIncome_Growth      0.6445      0.6467      0.6448      0.6509      0.6518 0.64774
# 4:               Range      0.0028      0.0021      0.0068      0.0060      0.0044 0.00442

### RESULTS:  Non-conclusive results as rankings fluctuated too much over 5 iterations.

# Obtain the concordance for different variables.
concTable(datCredit_valid_TFD, vars)
#               Variable Concordance          SD LR_Statistic
# 1: M_RealIncome_Growth   0.5341913 0.003951199           68
# 2:        M_Emp_Growth   0.5267075 0.003999905           45
# 3:    M_RealGDP_Growth   0.5211591 0.003910739           18

### RESULTS:  [M_RealIncome_Growth] has the highest concordance, although the concordances
###           are not that high.

### CONCLUSION: Keep M_RealIncome_Growth in the model.

varlist = vecChange(varlist, Remove=c("M_Emp_Growth", "M_RealGDP_Growth"))



# ------ 6.3 What are the predictive powers of the current variables in the models?

# Initialize variables
vars <- c("M_DTI_Growth", "M_Inflation_Growth", "M_RealIncome_Growth", "M_Repo_Rate")

# Compare goodness of fit of different variables
csTable(datCredit_train_TFD,vars,seedVal = NA)
#               Variables Iteration_1 Iteration_2 Iteration_3 Iteration_4 Iteration_5 Average
# 1:  M_Inflation_Growth      0.6479      0.6488      0.6452      0.6466      0.6503 0.64776
# 2:         M_Repo_Rate      0.6450      0.6501      0.6433      0.6472      0.6525 0.64762
# 3: M_RealIncome_Growth      0.6448      0.6474      0.6470      0.6473      0.6471 0.64672
# 4:        M_DTI_Growth      0.6454      0.6441      0.6466      0.6484      0.6471 0.64632
# 5:               Range      0.0031      0.0060      0.0037      0.0018      0.0054 0.00400

### RESULTS:  After 5 iterations, there is too much variability in kS-statistice, therefore no conclusions can be made.

# Compare concordance of different variables
concTable(datCredit_valid_TFD,vars)
#               Variable Concordance          SD LR_Statistic
# 1:        M_DTI_Growth   0.5350355 0.003955393          117
# 2: M_RealIncome_Growth   0.5341913 0.003951199           68
# 3:  M_Inflation_Growth   0.5202318 0.004253188           55
# 4:         M_Repo_Rate   0.5180882 0.004234430           29

### RESULTS: [M_DTI_Growth] has the highest concordance, albeit not that high.
### CONCLUSION: The low predictive power may be enhanced with lags.


# ------ 6.4 Which lags for current variables be included in the models?

# ------ 6.4.1 Should lags for [M_DTI_Growth] be included in the models?
vars <- c("M_DTI_Growth_1","M_DTI_Growth_10","M_DTI_Growth_11","M_DTI_Growth_12",
          "M_DTI_Growth_2","M_DTI_Growth_3","M_DTI_Growth_4","M_DTI_Growth_5",
          "M_DTI_Growth_6","M_DTI_Growth_7","M_DTI_Growth_8","M_DTI_Growth_9",
          "M_DTI_Growth")

cox <- coxph(as.formula(paste0("Surv(Start,End,Default_Ind) ~ ",
                               paste(vars,collapse=" + "))), id=LoanID, datCredit_train_TFD)
summary(cox); rm(cox)
### RESULTS:  Only [M_DTI_Growth] had a significant p-value, but concordance did improve.
###           The third and sixth lags had significant coef.

# Build a model on these lags.
vars <- c("M_DTI_Growth_3", "M_DTI_Growth_6","M_DTI_Growth")

cox <- coxph(as.formula(paste0("Surv(Start,End,Default_Ind) ~ ",
                               paste(vars,collapse=" + "))), id=LoanID, datCredit_train_TFD)
summary(cox); rm(cox)

### RESULTS: lag 3 was not significant. Therefore remove the lag from the model

vars <- c("M_DTI_Growth_6","M_DTI_Growth")

cox <- coxph(as.formula(paste0("Surv(Start,End,Default_Ind) ~ ",
                               paste(vars,collapse=" + "))), id=LoanID, datCredit_train_TFD)
summary(cox); rm(cox)

### RESULTS: The concordance was improved and all the variables are significant.

### CONCLUSION: Add [M_DTI_Growth_6] to the model.

varlist <- vecChange(varlist,Add=data.table(vars="M_DTI_Growth_6", vartypes="prc"))



# ------ 6.4.2 Should lags for [M_DTI_Growth] be included in the models?
vars <- c("M_Inflation_Growth_1","M_Inflation_Growth_10","M_Inflation_Growth_11",
          "M_Inflation_Growth_12","M_Inflation_Growth_2","M_Inflation_Growth_3",
          "M_Inflation_Growth_4","M_Inflation_Growth_5","M_Inflation_Growth_6",
          "M_Inflation_Growth_7","M_Inflation_Growth_8","M_Inflation_Growth_9",
          "M_Inflation_Growth")

cox <- coxph(as.formula(paste0("Surv(Start,End,Default_Ind) ~ ",
                               paste(vars,collapse=" + "))), id=LoanID, datCredit_train_TFD)
summary(cox); rm(cox)
### RESULTS:  [M_Inflation_Growth], [M_Inflation_Growth_12], [M_Inflation_Growth_11],
###           [M_Inflation_Growth_3] and [M_Inflation_Growth_1] has significant p-values,
###           with an improved concordance.

# Build a model on these lags.
vars <- c("M_Inflation_Growth", "M_Inflation_Growth_1","M_Inflation_Growth_3",
          "M_Inflation_Growth_11", "M_Inflation_Growth_12")

cox <- coxph(as.formula(paste0("Surv(Start,End,Default_Ind) ~ ",
                               paste(vars,collapse=" + "))), id=LoanID, datCredit_train_TFD)
summary(cox); rm(cox)

### RESULTS: lag 11 is not significant. Therefore remove the lag from the model

# Build a model on these lags.
vars <- c("M_Inflation_Growth", "M_Inflation_Growth_1","M_Inflation_Growth_3",
          "M_Inflation_Growth_12")

cox <- coxph(as.formula(paste0("Surv(Start,End,Default_Ind) ~ ",
                               paste(vars,collapse=" + "))), id=LoanID, datCredit_train_TFD)
summary(cox); rm(cox)

### RESULTS: lag 1 is not as significant as the other variables. Therefore remove the lag from the model

# Build a model on these lags.
vars <- c("M_Inflation_Growth","M_Inflation_Growth_3","M_Inflation_Growth_12")

cox <- coxph(as.formula(paste0("Surv(Start,End,Default_Ind) ~ ",
                               paste(vars,collapse=" + "))), id=LoanID, datCredit_train_TFD)
summary(cox); rm(cox)

### RESULTS: The concordance was improved and all the variables are significant.

### CONCLUSION: Add [M_Inflation_Growth_3] and [M_Inflation_Growth_12] to the model.

varlist <- vecChange(varlist,Add=data.table(vars=c("M_Inflation_Growth_3",
                                                   "M_Inflation_Growth_12"),
                                            vartypes=c("prc", "prc")))


# ------ 6.4.3 Should lags for [M_Repo_Rate] be included in the models?

vars <- c("M_Repo_Rate_1","M_Repo_Rate_10","M_Repo_Rate_11","M_Repo_Rate_12",
          "M_Repo_Rate_2","M_Repo_Rate_3","M_Repo_Rate_4","M_Repo_Rate_5",
          "M_Repo_Rate_6","M_Repo_Rate_7","M_Repo_Rate_8","M_Repo_Rate_9",
          "M_Repo_Rate")

cox <- coxph(as.formula(paste0("Surv(Start,End,Default_Ind) ~ ",
                               paste(vars,collapse=" + "))), id=LoanID, datCredit_train_TFD)
summary(cox); rm(cox)
### RESULTS:  [M_Repo_Rate], [M_Repo_Rate_3], [M_Repo_Rate_12],
###           [M_Inflation_Growth_3] are the most significant values.

# Build a model on these lags.
vars <- c("M_Repo_Rate", "M_Repo_Rate_3", "M_Repo_Rate_12")

cox <- coxph(as.formula(paste0("Surv(Start,End,Default_Ind) ~ ",
                               paste(vars,collapse=" + "))), id=LoanID, datCredit_train_TFD)
summary(cox); rm(cox)
### RESULTS: All variables are significant with [M_Repo_Rate_3] being the most significant.         

### CONCLUSION: Add [M_Inflation_Growth_3] and [M_Inflation_Growth_12] to the model.

varlist <- vecChange(varlist,Add=data.table(vars=c("M_Repo_Rate_3",
                                                   "M_Repo_Rate_12"),
                                            vartypes=c("prc", "prc")))


# ------ 6.4.4 Should lags for [M_RealIncome_Growth] be included in the models?
vars <- c("M_RealIncome_Growth_1","M_RealIncome_Growth_10","M_RealIncome_Growth_11",
          "M_RealIncome_Growth_12","M_RealIncome_Growth_2","M_RealIncome_Growth_3",
          "M_RealIncome_Growth_4","M_RealIncome_Growth_5","M_RealIncome_Growth_6",
          "M_RealIncome_Growth_7","M_RealIncome_Growth_8","M_RealIncome_Growth_9",
          "M_RealIncome_Growth")

cox <- coxph(as.formula(paste0("Surv(Start,End,Default_Ind) ~ ",paste(vars,collapse=" + "))), id=LoanID, datCredit_train_TFD)
summary(cox); rm(cox)
### RESULTS:  [M_RealIncome_Growth], [M_RealIncome_Growth_12], [M_RealIncome_Growth_11],
###           has significant p-values below the 0.001 threshold,
###           with an improved concordance.

# Build a model on these lags.
vars <- c("M_RealIncome_Growth", "M_RealIncome_Growth_12","M_RealIncome_Growth_11")

cox <- coxph(as.formula(paste0("Surv(Start,End,Default_Ind) ~ ",paste(vars,collapse=" + "))), id=LoanID, datCredit_train_TFD)
summary(cox); rm(cox)
### RESULTS: All values have significant p-values.


### CONCLUSION: Add [M_RealIncome_Growth_12] and [M_RealIncome_Growth_11] to the model.

varlist <- vecChange(varlist,Add=data.table(vars=c("M_RealIncome_Growth_11",
                                                   "M_RealIncome_Growth_12"),
                                            vartypes=c("prc", "prc")))



# ------ 6.5 What is the predictive performance of the currrent thematic model?

# Initiialize macro-economic thematic variables
vars <- c("M_DTI_Growth", "M_DTI_Growth_6", "M_Inflation_Growth",
          "M_Inflation_Growth_3", "M_Inflation_Growth_12", "M_RealIncome_Growth",
          "M_RealIncome_Growth_11", "M_RealIncome_Growth_12", "M_Repo_Rate",
          "M_Repo_Rate_3", "M_Repo_Rate_12")

coxMacro_train <- coxph(as.formula(paste0("Surv(Start,End,Default_Ind) ~ ",
                                          paste(vars,collapse=" + "))), id=LoanID,
                        datCredit_train_TFD)
summary(coxMacro_train)

### RESULTS: [M_RealIncome_Growth...], [M_Inflation_Growth] and [M_DTI_Growth_6]
###         does not have significant p-values.

vars <- c("M_DTI_Growth", "M_Repo_Rate", "M_Repo_Rate_3",
          "M_Repo_Rate_12")

coxMacro_train <- coxph(as.formula(paste0("Surv(Start,End,Default_Ind) ~ ",
                                          paste(vars,collapse=" + "))), id=LoanID,
                        datCredit_train_TFD)
summary(coxMacro_train)

### RESULTS: [M_Repo_Rate_12] does not have such significant p-values as the other variables.

vars <- c("M_DTI_Growth", "M_Repo_Rate", "M_Repo_Rate_3")

coxMacro_train <- coxph(as.formula(paste0("Surv(Start,End,Default_Ind) ~ ",
                                          paste(vars,collapse=" + "))), id=LoanID,
                        datCredit_train_TFD)
summary(coxMacro_train);rm(coxMacro_train)

### RESULTS: All variables have significant p-values.

#============================================================================================
### CONCLUSION: [M_DTI_Growth], [M_Repo_Rate] and [M_Repo_Rate_3] should be included in the model variables.

modelVar <- c(modelVar,vars)

#============================================================================================
