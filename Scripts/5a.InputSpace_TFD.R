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
  
  # Create a data table with correlation pairs
  corrProbs <- data.table(x = rownames(corrMat)[corrCoordinates[, 1]], y = colnames(corrMat)[corrCoordinates[, 2]])
  
  # Print the identified correlations
  for (i in 1:nrow(corrProbs)) {
    cat("Absolute correlations of ",percent(corrMat[corrProbs[i, x], corrProbs[i, y]]),
      " found for ", corrProbs[i, x], " and ", corrProbs[i, y],"\n")
  }
  
  # Return the correlation table
  return(corrProbs)
}

# Table concordance statistic of single variable cox ph models
# Used to compare predictive performance of variables
concTable <- function(data, variables) {
  # Use lapply to efficiently compute concordances for all univariate models
  results <- lapply(variables, function(var) {
    formula <- as.formula(paste0("Surv(Start, End, Default_Ind) ~ ", var))
    
    # Fit Cox model
    model <- coxph(formula,id=LoanID, data = data)
    
    # Extract concordance
    conc <- as.numeric(concordance(model)[1])
    
    # Extract concordance variability
    sd <- sqrt(concordance(model)$var)
    
    # Extract LRT from the model's log-likelihood
    # `loglik` contains log-likelihoods for null (intercept-only) and full models
    lr_stat <- round(2 * (model$loglik[2] - model$loglik[1]),0)
    
    # Return results as a data.table
    return(data.table(Variable = var, Concordance = conc, SD = sd, LR_Statistic = lr_stat))
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
        model <- coxph(formula, data = data)  # Fit a univariate Cox model
        cs_ks_test(model, data, GraphInd = FALSE)$Stat  # Calculate KS statistic
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
          cs_ks_test(model, data, GraphInd = FALSE)$Stat
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

# ------ 1. Delinquency measures
varlist <- data.table(vars=c("g0_Delinq","g0_Delinq_fac","PerfSpell_g0_Delinq_Num",
                             "TimeInDelinqState","g0_Delinq_Any_Aggr_Prop","g0_Delinq_Ave",
                             "slc_acct_arr_dir_3", "slc_acct_roll_ever_24_imputed_med"),
                      vartypes=c("acc", "cat", "acc", "acc", "dte", "dte", "cat", "acc"))

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

### CONCLUSION: Unable to add [Delinq_0] to the model, since the various forms'
###             coef is unstable.



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
          "TimeInDelinqState_Lag_1","slc_acct_arr_dir_3_Change_Ind")

# Goodness of fit
csTable(datCredit_train_TFD,vars)
#                         Variable KS_Statistic
# 2                 g0_Delinq_Ave       0.6477
# 1       PerfSpell_g0_Delinq_Num       0.6458
# 5 slc_acct_arr_dir_3_Change_Ind       0.6446
# 4       TimeInDelinqState_Lag_1       0.6331
# 3                g0_Delinq_SD_4       0.6151

### RESULTS: [g0_Delinq_SD_4] seems to have a notacible worse fit than the other variabeles.

# Accuracy
concTable(datCredit_valid_TFD,vars)
#                         Variable Concordance           SD LR_Statistic
# 1:                g0_Delinq_SD_4   0.9803661 0.0014018721        48597
# 2:       PerfSpell_g0_Delinq_Num   0.9474644 0.0006456311         4596
# 3: slc_acct_arr_dir_3_Change_Ind   0.8433911 0.0019345216        18802
# 4:       TimeInDelinqState_Lag_1   0.7507071 0.0056323114        12757
# 5:                 g0_Delinq_Ave   0.5397435 0.0040305498           90

### RESULTS: [g0_Delinq_Ave] has significant less concordance than the other variables

### CONCLUSION: Leave all variables in the model (including [g0_Delinq_Ave], since it has the best fit for the data).

# ------ 1.5 What is the performance of current thematic cox ph model?

# Goodness of fit of coxDelinq model

# Build cox model based on all thematic variables
coxDelinq_train <- coxph(Surv(Start,End,Default_Ind) ~ g0_Delinq_SD_4 + g0_Delinq_Ave +
                     PerfSpell_g0_Delinq_Num + slc_acct_arr_dir_3_Change_Ind +
                       TimeInDelinqState_Lag_1, id=LoanID,
                   data=datCredit_train_TFD)

# Kolmogorov-Smirnof of coxDelinq
cs_ks_test(coxDelinq_train,datCredit_train_TFD,GraphInd = FALSE) # 0.617

# Accuracy
coxDelinq_valid<- coxph(Surv(Start,End,Default_Ind) ~ g0_Delinq_SD_4 + g0_Delinq_Ave +
                           PerfSpell_g0_Delinq_Num + slc_acct_arr_dir_3_Change_Ind +
                          TimeInDelinqState_Lag_1, id=LoanID,
                         data=datCredit_valid_TFD)

# (0,3) (4,12) (13,24) (0,12) (0,36)
tROC(datCredit_valid_TFD, coxDelinq_valid, month_Start=0, month_End=36,
        fld_ID="LoanID", fld_Event="Default_Ind",fld_StartTime="Start",
        fld_EndTime="End", numDigits=0, Graph=FALSE)
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





# ------ 2. Engineered measures
varlist <- data.table(vars=c("slc_acct_pre_lim_perc_imputed_med",
                             "slc_acct_prepaid_perc_dir_12_imputed_med",
                             "slc_pmnt_method"),
                      vartypes=c("cat", "dec", "dec", "dec", "fin", "cat"))

# Sanity Check
colCheck(varlist[,vars],datCredit_train_TFD)

# ------ 2.1 Check correlations of variables to ascertain variables that are similar.

# - Correlation analysis
corGroups <- corrAnalysis(datCredit_train_TFD, varlist[vartypes!="cat"]$vars, corrThresh = 0.6, method = 'spearman') # Obtain correlation groups

### RESULTS:  1) [slc_acct_pre_lim_perc_imputed_med] and [slc_acct_prepaid_perc_dir_12_imputed_med] with a correlation of 0.82
#             2) All [value_...] indicators of missing values.

### INVESTIGATE: Compare [slc_acct_pre_lim_perc_imputed_med] vs [slc_acct_prepaid_perc_dir_12_imputed_med]

# Initialize variables to be tested
vars <- c("slc_acct_pre_lim_perc_imputed_med", "slc_acct_prepaid_perc_dir_12_imputed_med")

# Compare goodness of fit of different variables
getKS(datCredit_train_TFD,vars)
#                                 Variable     KS
# 1        slc_acct_pre_lim_perc_imputed_med 0.6485
# 2 slc_acct_prepaid_perc_dir_12_imputed_med 0.6444

### RESULTS:  The values are quite close to one another with a varying range, therefore
###           no decisions can be made from the results.

# Compare concordance of different variables
getConcs(datCredit_valid_TFD,vars)
#                                     Variable Concordance LR_Statistic
# 1:        slc_acct_pre_lim_perc_imputed_med   0.6360895         2865
# 2: slc_acct_prepaid_perc_dir_12_imputed_med   0.5857421          120

### RESULTS: [slc_acct_pre_lim_perc_imputed_med] has a higher concordance than that of [slc_acct_prepaid_perc_dir_12_imputed_med]
### CONCLUSION: Remove [slc_acct_prepaid_perc_dir_12_imputed_med]

varlist <- vecChange(varlist,Remove="slc_acct_prepaid_perc_dir_12_imputed_med")

### RESULTS:  [value...] functions are all identical. Therefore a value must be chosen to
###           represent the rest, e.g. [value_ind_slc_past_due_amt]

varlist <- vecChange(varlist,Remove=c("value_ind_slc_acct_pre_lim_perc",
                                      "value_ind_slc_acct_prepaid_perc_dir_12",
                                      "value_ind_slc_acct_roll_ever_24",
                                      "value_ind_slc_past_due_amt"))

# ------ 2.2 Test for proportionality in hazard models

# [slc_acct_arr_dir_3]
cox <- coxph(Surv(Start,End, Default_Ind) ~ slc_acct_arr_dir_3, datCredit_train_TFD)
### RESULTS: Coefficients does not converge

### INVESTIGATE [slc_acct_arr_dir_3]
describe(datCredit_train_TFD[Default_Ind==1,]$slc_acct_arr_dir_3)
# Value            CURING MISSING_DATA      ROLLING         SAME
# Frequency            68         1470         7189          454
# Proportion        0.007        0.160        0.783        0.049

### RESULTS: Merge small categories.

datCredit_train_TFD[, slc_acct_arr_dir_3_grp :=
                      ifelse(slc_acct_arr_dir_3 == "CURING" |
                               slc_acct_arr_dir_3 == "SAME", "OTHER",
                            slc_acct_arr_dir_3)]

cox <- coxph(Surv(Start,End, Default_Ind) ~ slc_acct_arr_dir_3, datCredit_train_TFD)
### RESULTS: Coefficients does not converge

datCredit_train_TFD[, slc_acct_arr_dir_3_grp :=
                      ifelse(slc_acct_arr_dir_3 == "ROLLING" ,1,0)]

cox <- coxph(Surv(Start,End, Default_Ind) ~ slc_acct_arr_dir_3_grp, datCredit_train_TFD)
### RESULTS: Coefficients does not converge

### CONCLUSION: Remove slc_acct_arr_dir_3

varlist <- vecChange(varlist,Remove="slc_acct_arr_dir_3")
datCredit_train_TFD[, slc_acct_arr_dir_3_grp := NULL]

# [slc_acct_pre_lim_perc_imputed_med]
cox <- coxph(Surv(Start,End, Default_Ind) ~ slc_acct_pre_lim_perc_imputed_med, datCredit_train_TFD)
sfResiduals(cox,datCredit_train_TFD,"slc_acct_pre_lim_perc_imputed_med", c(100,0.2))
### RESULTS: Residuals tend to zero
### CONCLUSION: [slc_acct_pre_lim_perc_imputed_med] should be kept in the model

# [slc_acct_roll_ever_24_imputed_med]
cox <- coxph(Surv(Start,End, Default_Ind) ~ slc_acct_roll_ever_24_imputed_med, datCredit_train_TFD)
sfResiduals(cox,datCredit_train_TFD,"slc_acct_roll_ever_24_imputed_med", c(100,0.2))
### RESULTS: Residuals does not tend to zero, however the graphs almost act like categorical values.
cox <- coxph(Surv(Start,End, Default_Ind) ~ slc_acct_roll_ever_24_imputed_med_f, datCredit_train_TFD)
sfResiduals(cox,datCredit_train_TFD,"slc_acct_roll_ever_24_imputed_med_f")
### RESULTS: Improved Concordance and loess curves around 0
### Conclusion: Keep slc_acct_roll_ever_24_imputed_med_f

# [slc_past_due_amt_imputed_med]
cox <- coxph(Surv(Start,End, Default_Ind) ~ slc_past_due_amt_imputed_med, datCredit_train_TFD)
sfResiduals(cox,datCredit_train_TFD,"slc_past_due_amt_imputed_med", c(100,100000))
### RESULTS: Residuals does not tend to zero.
### Conclusion: Reject proportionality for [slc_past_due_amt_imputed_med]

varlist <- vecChange(varlist,Remove="slc_past_due_amt_imputed_med")

# [slc_pmnt_method]
cox <- coxph(Surv(Start,End, Default_Ind) ~ slc_pmnt_method, datCredit_train_TFD)
sfResiduals(cox,datCredit_train_TFD,"slc_pmnt_method", c(100,0.2))
### RESULTS: Residuals does not tend to zero.
### Conclusion: Reject proportionality for [slc_pmnt_method]

varlist <- vecChange(varlist,Remove="slc_pmnt_method")

modelVar <- c(modelVar,varlist$vars)

# ------ 3. Interest Rate
varlist <- data.table(vars=c("InterestRate_Nom", "InterestRate_Margin"),
                      vartypes=c("prc","prc"))

# ------ 3.1 Check correlations of variables to ascertain variables that are similar.

# - Correlation analysis
corGroups <- corrAnalysis(datCredit_train_TFD, varlist[vartypes!="cat"]$vars, corrThresh = 0.6, method = 'spearman') # Obtain correlation groups

### RESULTS: Correlation of 0.41, therefore no significant correlations.

# Obtain the KS-statistics for goodness of fit tests.
getKS(datCredit_train_TFD, varlist$vars)
#           Variable     KS
# 1    InterestRate_Nom 0.6497
# 2 InterestRate_Margin 0.6490

### RESULTS:  After a few iterations, it would seem that [InterestRate_Nom] is a
###           a better fit than [InterestRate_Margin], although with slim differences.

# Obtain the concordance for different variables.
getConcs(datCredit_valid_TFD, varlist$vars)
#                 Variable Concordance LR_Statistic
# 1:    InterestRate_Nom   0.5497724          165
# 2: InterestRate_Margin   0.5452677          174

### RESULTS:  [InterestRate_Nom] has a slightly higher concordance than [InterestRate_Margin].

### CONCLUSION: Variables are quite similar in predictave power, however are dissimilar
###              in correlation, therefore both can be kept in the model.

# ------ 3.2 Explore how [InterestRate_Margin] can be improved.

### INVESTIGATE: [InterestRate_Margin_Aggr_Med]

vars <- c("InterestRate_Margin", "InterestRate_Margin_Aggr_Med")

# - Correlation analysis
datCredit_train_TFD %>% subset(select = vars) %>% cor(method = 'spearman')

### RESULTS:  [InterestRate_Margin] and [InterestRate_Margin_Aggr_Med] does not have
###           a significant correlation.
### CONCLUSION: Both variables can be kept in the model.

# Obtain the KS-statistics for goodness of fit tests.
getKS(datCredit_train_TFD, vars)
#                     Variable     KS
# 1          InterestRate_Margin 0.6504
# 2 InterestRate_Margin_Aggr_Med 0.6482

### RESULTS:  Non-conclusive results as KS fluctuate too much.

# Obtain the concordance for different variables.
getConcs(datCredit_valid_TFD, vars)
#                         Variable Concordance LR_Statistic
# 1: InterestRate_Margin_Aggr_Med   0.5506784           96
# 2:          InterestRate_Margin   0.5452677          174

### RESULTS:  [InterestRate_Margin_Aggr_Med] has a slightly higher concordance than [InterestRate_Margin].

### CONCLUSIONS: Add [InterestRate_Margin_Aggr_Med] to the model.

varlist <- vecChange(varlist,Add=data.table(vars=c("InterestRate_Margin_Aggr_Med"), vartypes=c("prc")))

### INVESTIGATE: Lagging [InterestRate_Margin_Aggr_Med]

vars <- c("InterestRate_Margin_Aggr_Med","InterestRate_Margin_Aggr_Med_1",
          "InterestRate_Margin_Aggr_Med_12","InterestRate_Margin_Aggr_Med_2",
          "InterestRate_Margin_Aggr_Med_3","InterestRate_Margin_Aggr_Med_4",
          "InterestRate_Margin_Aggr_Med_5","InterestRate_Margin_Aggr_Med_6",
          "InterestRate_Margin_Aggr_Med_9")

# - Correlation analysis
corGroups <- corrAnalysis(datCredit_train_TFD, vars, corrThresh = 0.6, method = 'spearman') # Obtain correlation groups

### RESULTS: All are significantly correlated

# Obtain the KS-statistics for goodness of fit tests.
getKS(datCredit_train_TFD, vars)
#                     Variable     KS
# 1          InterestRate_Margin 0.6504
# 2 InterestRate_Margin_Aggr_Med 0.6482

### RESULTS:  Non-conclusive results as KS fluctuate too much.

# Obtain the concordance for different variables.
getConcs(datCredit_valid_TFD, vars)
#                             Variable Concordance LR_Statistic
# 1:  InterestRate_Margin_Aggr_Med_4   0.5611257          104
# 2:  InterestRate_Margin_Aggr_Med_3   0.5603674          105
# 3:  InterestRate_Margin_Aggr_Med_2   0.5575290          102
# 4:  InterestRate_Margin_Aggr_Med_5   0.5566768           99
# 5:  InterestRate_Margin_Aggr_Med_6   0.5566370           98
# 6:  InterestRate_Margin_Aggr_Med_1   0.5555107           99
# 7:  InterestRate_Margin_Aggr_Med_9   0.5546448           85
# 8:    InterestRate_Margin_Aggr_Med   0.5506784           96
# 9: InterestRate_Margin_Aggr_Med_12   0.5481947           68

### RESULTS:  [InterestRate_Margin_Aggr_Med_4] has the best concordance, significantly
###           better than [InterestRate_Margin_Aggr_Med].

### CONCLUSION: Replace [InterestRate_Margin_Aggr_Med] with [InterestRate_Margin_Aggr_Med_4]

varlist <- vecChange(varlist,Remove=c("InterestRate_Margin_Aggr_Med"),
                     Add=data.table(vars=c("InterestRate_Margin_Aggr_Med_4"),
                                    vartypes=c("prc")))

# ------ 1.5 Test for proportionality in hazard models

# [InterestRate_Nom]
cox <- coxph(Surv(Start,End, Default_Ind) ~ InterestRate_Nom, datCredit_train_TFD)
sfResiduals(cox,datCredit_train_TFD,"InterestRate_Nom", c(50,0.1))
### RESULTS: NO trend away from 0, therefore do not reject proportionality.

# [InterestRate_Margin]
cox <- coxph(Surv(Start,End, Default_Ind) ~ InterestRate_Margin, datCredit_train_TFD)
sfResiduals(cox,datCredit_train_TFD,"InterestRate_Margin", c(50,0.1))
### RESULTS: NO trend away from 0, therefore do not reject proportionality.

# [InterestRate_Margin_Aggr_Med_4]
cox <- coxph(Surv(Start,End, Default_Ind) ~ InterestRate_Margin_Aggr_Med_4, datCredit_train_TFD)
sfResiduals(cox,datCredit_train_TFD,"InterestRate_Margin_Aggr_Med_4", c(50,0.01))
### RESULTS: Small trend away from 0, therefore do not reject proportionality.

# Add to model variables
modelVar <- c(modelVar,varlist$vars)

# ------ 4. General
varlist <- data.table(vars=c("Arrears","Balance","Instalment",
                             "pmnt_method_grp","Principal","Undrawn_Amt",
                             "LN_TPE","Term"),
                      vartypes=c("int", "fin", "fin", "cat", "int","fin", "cat", "int"))

# ------ 4.1 Check correlations of variables to ascertain variables that are similar.

# - Correlation analysis
corGroups <- corrAnalysis(datCredit_train_TFD, varlist[vartypes!="cat"]$vars, corrThresh = 0.6, method = 'spearman') # Obtain correlation groups

### RESULTS: [Balance], [Instalment] and [Principal] are highly corrolated.

### INVESTIGATE: Determine which variable to remove.

vars <- c("Balance", "Instalment", "Principal")

# Obtain the KS-statistics for goodness of fit tests.
getKS(datCredit_train_TFD, vars)

### RESULTS:  Non-conclusive results as KS fluctuate too much.

# Obtain the concordance for different variables.
getConcs(datCredit_valid_TFD, vars)
#     Variable Concordance LR_Statistic
# 1:  Principal   0.5827355          335
# 2: Instalment   0.5636894           84
# 3:    Balance   0.5544691          101

## RESULTS: [Principal] has the highest concordance

### CONCOLUSION: Keep [Principal] in the varlist

varlist <- vecChange(varlist,Remove=data.table(vars=c("Instalment","Balance")))

# ------ 4.2 Incorporate a timing component to Age and Balance with reference to term.

varlist <- vecChange(varlist,Add=data.table(vars=c("AgeToTerm","BalanceToTerm"), vartypes=c("dec","dec")))

# - Correlation analysis
corGroups <- corrAnalysis(datCredit_train_TFD, varlist$vars[varlist$vartypes != 'cat'], corrThresh = 0.6, method = 'spearman') # Obtain correlation groups

### RESULTS: [Principal] and [BalanceToTerm] are highly correlated

### INVESTIGATE: Remove highly correlated variables

vars <- c("Principal", "BalanceToTerm")

# Obtain the KS-statistics for goodness of fit tests.
getKS(datCredit_train_TFD, vars)
#       Variable     KS
# 2 BalanceToTerm 0.6480
# 1     Principal 0.6466

### RESULTS:  Non-conclusive results as KS fluctuate too much.

# Obtain the concordance for different variables.
getConcs(datCredit_valid_TFD, vars)
#         Variable Concordance LR_Statistic
# 1:     Principal   0.5827355          335
# 2: BalanceToTerm   0.5561134           99

### RESULTS:  [Principal] has a higher concordance than that of
###           [BalanceToTerm].

### CONCLUSION: Remove BalanceToTerm

varlist <- vecChange(varlist,Remove="BalanceToTerm")

# ------ 4.3 Test for proportionality in hazard models

# [Arrears]
cox <- coxph(Surv(Start,End, Default_Ind) ~ Arrears, datCredit_train_TFD)
sfResiduals(cox,datCredit_train_TFD,"Arrears", c(50,-50000))
### RESULTS: Clear trend away from 0, therefore reject proportionality.

# [pmnt_method_grp]
cox <- coxph(Surv(Start,End, Default_Ind) ~ pmnt_method_grp, datCredit_train_TFD)
sfResiduals(cox,datCredit_train_TFD,"pmnt_method_grp")
### RESULTS: Clear trends away from 0, therefore reject proportionality.

# [Principal]
cox <- coxph(Surv(Start,End, Default_Ind) ~ Principal, datCredit_train_TFD)
sfResiduals(cox,datCredit_train_TFD,"Principal", c(50,5000000))
### RESULTS:  Small trend away from 0, may be influenced by extreme values,
###           but do not reject proportionality.

# [Undrawn_Amt]
cox <- coxph(Surv(Start,End, Default_Ind) ~ Undrawn_Amt, datCredit_train_TFD)
sfResiduals(cox,datCredit_train_TFD,"Undrawn_Amt", c(50,1250000))
### RESULTS: Small trend away from 0, may be influenced by extreme values,
###           but do not reject proportionality.

# [LN_TPE]
cox <- coxph(Surv(Start,End, Default_Ind) ~ LN_TPE, datCredit_train_TFD)
sfResiduals(cox,datCredit_train_TFD,"LN_TPE")
### RESULTS: Clear trend away from 0, therefore reject proportionality.

# [Term]
cox <- coxph(Surv(Start,End, Default_Ind) ~ Term, datCredit_train_TFD)
sfResiduals(cox,datCredit_train_TFD,"Term", c(50,50))
### RESULTS: Clear trend away from 0, therefore reject proportionality.

# [AgeToTerm]
cox <- coxph(Surv(Start,End, Default_Ind) ~ AgeToTerm, datCredit_train_TFD)
sfResiduals(cox,datCredit_train_TFD,"AgeToTerm")
### RESULTS:  NO trend away from 0, although it does have some extreme values
###           Do not reject proportionality.

### CONCLUSION: Add [Principal], [Undrawn_Amt], [AgeToTerm] to the model.

modelVar <- c(modelVar, c("Principal", "Undrawn_Amt", "AgeToTerm"))

# ------ 5. Portfolio variables
varlist <- data.table(vars=c("AgeToTerm_Aggr_Mean",
                             "ArrearsToBalance_Aggr_Prop",
                             "CuringEvents_Aggr_Prop",
                             "InstalmentToBalance_Aggr_Prop",
                             "NewLoans_Aggr_Prop"),
                      vartypes=c("dec", "dec", "dec", "dec", "dec"))

# ------ 5.1 Check correlations of variables to ascertain variables that are similar.

# - Correlation analysis
corGroups <- corrAnalysis(datCredit_train_TFD, varlist$vars[varlist$vartypes != 'cat'], corrThresh = 0.6, method = 'spearman') # Obtain correlation groups

### RESULTS: [AgeToTerm_Aggr_Mean] and [ArrearsToBalance_Aggr_Prop] are highly correlated

vars <- c("AgeToTerm_Aggr_Mean", "ArrearsToBalance_Aggr_Prop")

# Obtain the KS-statistics for goodness of fit tests.
getKS(datCredit_train_TFD, vars)
#                   Variable     KS
# 1        AgeToTerm_Aggr_Mean 0.6509
# 2 ArrearsToBalance_Aggr_Prop 0.6481

### RESULTS:  Non-conclusive results as KS fluctuate too much.

# Obtain the concordance for different variables.
getConcs(datCredit_valid_TFD, vars)
#                       Variable Concordance LR_Statistic
# 1: ArrearsToBalance_Aggr_Prop   0.5367480          118
# 2:        AgeToTerm_Aggr_Mean   0.5091579           24

### RESULTS:  [ArrearsToBalance_Aggr_Prop] has higher concordance than that of
###           [AgeToTerm_Aggr_Mean].

### CONCLUSION: Remove [AgeToTerm_Aggr_Mean] from the data.

varlist <- vecChange(varlist,Remove="AgeToTerm_Aggr_Mean")

# ------ 5.2 Test for proportionality in hazard models

# [ArrearsToBalance_Aggr_Prop]
cox <- coxph(Surv(Start,End, Default_Ind) ~ ArrearsToBalance_Aggr_Prop, datCredit_train_TFD)
sfResiduals(cox,datCredit_train_TFD,"ArrearsToBalance_Aggr_Prop", c(50, 0.0005))
### RESULTS: Minimal trend away from 0, therefore do not reject proportionality.

# [CuringEvents_Aggr_Prop]
cox <- coxph(Surv(Start,End, Default_Ind) ~ CuringEvents_Aggr_Prop, datCredit_train_TFD)
sfResiduals(cox,datCredit_train_TFD,"CuringEvents_Aggr_Prop",c(100, 0.0015))
### RESULTS: Minimal trend away from 0, therefore do not reject proportionality.

# [InstalmentToBalance_Aggr_Prop]
cox <- coxph(Surv(Start,End, Default_Ind) ~ InstalmentToBalance_Aggr_Prop, datCredit_train_TFD)
sfResiduals(cox,datCredit_train_TFD,"InstalmentToBalance_Aggr_Prop",c(100, 0.0015))
### RESULTS: Trends away from 0, therefore reject proportionality.

# [NewLoans_Aggr_Prop]
cox <- coxph(Surv(Start,End, Default_Ind) ~ NewLoans_Aggr_Prop, datCredit_train_TFD)
sfResiduals(cox,datCredit_train_TFD,"NewLoans_Aggr_Prop",c(100, 0.007))
### RESULTS: Minimal trend away from 0, therefore do not reject proportionality.

### CONCLUSION: Add [ArrearsToBalance_Aggr_Prop], [CuringEvents_Aggr_Prop], [NewLoans_Aggr_Prop] to the model.

modelVar <- c(modelVar, c("ArrearsToBalance_Aggr_Prop", "CuringEvents_Aggr_Prop", "NewLoans_Aggr_Prop"))

# ------ 6. Macro Economic
varlist <- data.table(vars=c("M_DTI_Growth","M_Emp_Growth","M_Inflation_Growth",
                             "M_RealGDP_Growth","M_RealIncome_Growth","M_Repo_Rate"),
                      vartypes=c("perc", "perc", "perc", "perc", "perc", "perc"))

# ------ 6.1 Check correlations of variables to ascertain variables that are similar.

# - Correlation analysis
corGroups <- corrAnalysis(datCredit_train_TFD, varlist$vars[varlist$vartypes != 'cat'], corrThresh = 0.6, method = 'spearman') # Obtain correlation groups

### RESULTS: [M_Emp_Growth], [M_RealGDP_Growth] and [M_RealIncome_Growth] are highly correlated
### NOTE: Each variable contains lagging variable, therefore first the lagging variables will
###       tested and compared with the others.

### INVESTIGATE: Lagging [M_Emp_Growth]

vars <- c("M_Emp_Growth_1","M_Emp_Growth_10","M_Emp_Growth_11","M_Emp_Growth_12",
          "M_Emp_Growth_2","M_Emp_Growth_3","M_Emp_Growth_4","M_Emp_Growth_5",
          "M_Emp_Growth_6","M_Emp_Growth_7","M_Emp_Growth_8","M_Emp_Growth_9",
          "M_Emp_Growth")

# Obtain the KS-statistics for goodness of fit tests.
getKS(datCredit_train_TFD, vars)

### RESULTS:  Non-conclusive results as KS fluctuate too much.
#         Variable     KS
# 12  M_Emp_Growth_9 0.6507
# 6   M_Emp_Growth_3 0.6501
# 3  M_Emp_Growth_11 0.6485
# 11  M_Emp_Growth_8 0.6469
# 13    M_Emp_Growth 0.6468
# 2  M_Emp_Growth_10 0.6467
# 4  M_Emp_Growth_12 0.6467
# 1   M_Emp_Growth_1 0.6464
# 10  M_Emp_Growth_7 0.6464
# 7   M_Emp_Growth_4 0.6455
# 9   M_Emp_Growth_6 0.6453
# 5   M_Emp_Growth_2 0.6451
# 8   M_Emp_Growth_5 0.6448

# Obtain the concordance for different variables.
getConcs(datCredit_valid_TFD, vars)
#             Variable Concordance LR_Statistic
# 1:    M_Emp_Growth   0.5267075           45
# 2: M_Emp_Growth_11   0.5265101           21
# 3: M_Emp_Growth_10   0.5225248           14
# 4:  M_Emp_Growth_1   0.5224465           29
# 5: M_Emp_Growth_12   0.5223535           28
# 6:  M_Emp_Growth_2   0.5187988           17
# 7:  M_Emp_Growth_9   0.5171781            8
# 8:  M_Emp_Growth_3   0.5141904            8
# 9:  M_Emp_Growth_4   0.5084908            3
# 10:  M_Emp_Growth_8   0.5076886            4
# 11:  M_Emp_Growth_5   0.5045795            0
# 12:  M_Emp_Growth_7   0.5037239            1
# 13:  M_Emp_Growth_6   0.5001361            0

### RESULTS:  [M_Emp_Growth] has the best concordance.

## INVESTIGATE: Lagging [M_RealGDP_Growth]

vars <- c("M_RealGDP_Growth_1","M_RealGDP_Growth_10","M_RealGDP_Growth_11",
          "M_RealGDP_Growth_12","M_RealGDP_Growth_2","M_RealGDP_Growth_3",
          "M_RealGDP_Growth_4","M_RealGDP_Growth_5","M_RealGDP_Growth_6",
          "M_RealGDP_Growth_7","M_RealGDP_Growth_8","M_RealGDP_Growth_9","M_RealGDP_Growth")

# Obtain the KS-statistics for goodness of fit tests.
getKS(datCredit_train_TFD, vars)
#             Variable     KS
# 13    M_RealGDP_Growth 0.6507
# 3  M_RealGDP_Growth_11 0.6506
# 9   M_RealGDP_Growth_6 0.6500
# 12  M_RealGDP_Growth_9 0.6496
# 11  M_RealGDP_Growth_8 0.6490
# 7   M_RealGDP_Growth_4 0.6489
# 5   M_RealGDP_Growth_2 0.6485
# 4  M_RealGDP_Growth_12 0.6484
# 10  M_RealGDP_Growth_7 0.6470
# 2  M_RealGDP_Growth_10 0.6467
# 1   M_RealGDP_Growth_1 0.6465
# 8   M_RealGDP_Growth_5 0.6464
# 6   M_RealGDP_Growth_3 0.6442

### RESULTS:  Non-conclusive results as KS fluctuate too much.


# Obtain the concordance for different variables.
getConcs(datCredit_valid_TFD, vars)
#                 Variable Concordance LR_Statistic
# 1: M_RealGDP_Growth_12   0.5240902           40
# 2:    M_RealGDP_Growth   0.5211591           18
# 3: M_RealGDP_Growth_11   0.5210253           29
# 4:  M_RealGDP_Growth_1   0.5185701           12
# 5:  M_RealGDP_Growth_2   0.5171132            7
# 6: M_RealGDP_Growth_10   0.5164026           19
# 7:  M_RealGDP_Growth_3   0.5131770            3
# 8:  M_RealGDP_Growth_9   0.5111747           11
# 9:  M_RealGDP_Growth_4   0.5105378            1
# 10:  M_RealGDP_Growth_8   0.5051617            6
# 11:  M_RealGDP_Growth_7   0.5005288            3
# 12:  M_RealGDP_Growth_6   0.4971260            1
# 13:  M_RealGDP_Growth_5   0.4933388            0

### RESULTS:  [M_RealGDP_Growth_12] has the best concordance.

varlist <- vecChange(varlist,Remove="M_RealGDP_Growth",Add=data.table(vars=c("M_RealGDP_Growth_12"), vartypes=c("prc")))

### INVESTIGATE: Lagging [M_RealIncome_Growth]

vars <- c("M_RealIncome_Growth_1","M_RealIncome_Growth_10","M_RealIncome_Growth_11","M_RealIncome_Growth_12",
          "M_RealIncome_Growth_2","M_RealIncome_Growth_3","M_RealIncome_Growth_4","M_RealIncome_Growth_5",
          "M_RealIncome_Growth_6","M_RealIncome_Growth_7","M_RealIncome_Growth_8","M_RealIncome_Growth_9","M_RealIncome_Growth")

# Obtain the KS-statistics for goodness of fit tests.
getKS(datCredit_train_TFD, vars)
#             Variable     KS
# 13    M_RealGDP_Growth 0.6507
# 3  M_RealGDP_Growth_11 0.6506
# 9   M_RealGDP_Growth_6 0.6500
# 12  M_RealGDP_Growth_9 0.6496
# 11  M_RealGDP_Growth_8 0.6490
# 7   M_RealGDP_Growth_4 0.6489
# 5   M_RealGDP_Growth_2 0.6485
# 4  M_RealGDP_Growth_12 0.6484
# 10  M_RealGDP_Growth_7 0.6470
# 2  M_RealGDP_Growth_10 0.6467
# 1   M_RealGDP_Growth_1 0.6465
# 8   M_RealGDP_Growth_5 0.6464
# 6   M_RealGDP_Growth_3 0.6442

### RESULTS:  Non-conclusive results as KS fluctuate too much.


# Obtain the concordance for different variables.
getConcs(datCredit_valid_TFD, vars)
#                   Variable Concordance LR_Statistic
# 1:    M_RealIncome_Growth   0.5341913           68
# 2:  M_RealIncome_Growth_1   0.5318903           58
# 3:  M_RealIncome_Growth_2   0.5282074           46
# 4:  M_RealIncome_Growth_3   0.5245396           33
# 5:  M_RealIncome_Growth_4   0.5202059           23
# 6:  M_RealIncome_Growth_5   0.5173732           15
# 7:  M_RealIncome_Growth_6   0.5155281            9
# 8: M_RealIncome_Growth_12   0.5131676           17
# 9:  M_RealIncome_Growth_7   0.5125671            5
# 10:  M_RealIncome_Growth_8   0.5087736            2
# 11: M_RealIncome_Growth_11   0.5068548            7
# 12:  M_RealIncome_Growth_9   0.5039898            0
# 13: M_RealIncome_Growth_10   0.5013712            1

### RESULTS:  [M_RealIncome_Growth] has the best concordance.

### INVESTIGATE: Which variable should be removed.

vars <- c("M_Emp_Growth","M_RealGDP_Growth_12","M_RealIncome_Growth")

# Obtain the KS-statistics for goodness of fit tests.
getKS(datCredit_train_TFD, vars)
#             Variable     KS
# 1        M_Emp_Growth 0.6472
# 2 M_RealGDP_Growth_12 0.6471
# 3 M_RealIncome_Growth 0.6456

### RESULTS:  Non-conclusive results as KS fluctuate too much.

# Obtain the concordance for different variables.
getConcs(datCredit_valid_TFD, vars)
#                 Variable Concordance LR_Statistic
# 1: M_RealIncome_Growth   0.5341913           68
# 2:        M_Emp_Growth   0.5267075           45
# 3: M_RealGDP_Growth_12   0.5240902           40

### RESULTS:  [M_RealIncome_Growth] has the highest concordance.
### CONCLUSION: Keep M_RealIncome_Growth in the varlist

varlist <- vecChange(varlist,Remove = c("M_Emp_Growth", "M_RealGDP_Growth_12"))

# ------ 6.2 Explore other lagging variables.

### INVESTIGATE: [M_DTI_Growth]

vars <- c("M_DTI_Growth_1","M_DTI_Growth_10","M_DTI_Growth_11","M_DTI_Growth_12",
          "M_DTI_Growth_2","M_DTI_Growth_3","M_DTI_Growth_4","M_DTI_Growth_5",
          "M_DTI_Growth_6","M_DTI_Growth_7","M_DTI_Growth_8","M_DTI_Growth_9","M_DTI_Growth")

# Obtain the KS-statistics for goodness of fit tests.
getKS(datCredit_train_TFD, vars)

### RESULTS:  Non-conclusive results as KS fluctuate too much.
#         Variable     KS
# 7   M_DTI_Growth_4 0.6498
# 8   M_DTI_Growth_5 0.6493
# 12  M_DTI_Growth_9 0.6488
# 6   M_DTI_Growth_3 0.6480
# 10  M_DTI_Growth_7 0.6480
# 1   M_DTI_Growth_1 0.6474
# 3  M_DTI_Growth_11 0.6467
# 2  M_DTI_Growth_10 0.6463
# 13    M_DTI_Growth 0.6463
# 5   M_DTI_Growth_2 0.6460
# 9   M_DTI_Growth_6 0.6453
# 4  M_DTI_Growth_12 0.6443
# 11  M_DTI_Growth_8 0.6439

# Obtain the concordance for different variables.
getConcs(datCredit_valid_TFD, vars)
#           Variable Concordance LR_Statistic
# 1: M_DTI_Growth_11   0.5669068          195
# 2: M_DTI_Growth_12   0.5660056          196
# 3: M_DTI_Growth_10   0.5655813          185
# 4:  M_DTI_Growth_9   0.5563671          166
# 5:  M_DTI_Growth_8   0.5523882          146
# 6:  M_DTI_Growth_7   0.5477875          133
# 7:  M_DTI_Growth_6   0.5447379          121
# 8:  M_DTI_Growth_5   0.5411280          111
# 9:  M_DTI_Growth_4   0.5384984          100
# 10:  M_DTI_Growth_3   0.5361471           96
# 11:  M_DTI_Growth_1   0.5360005          110
# 12:  M_DTI_Growth_2   0.5358584          100
# 13:    M_DTI_Growth   0.5350355          117

### RESULTS:  [M_DTI_Growth_11] has the best concordance.
### CONCLUSION: Replace [M_DTI_Growth] with [M_DTI_Growth_11]

varlist <- vecChange(varlist,Remove="M_DTI_Growth",Add=data.table(vars=c("M_DTI_Growth_11"), vartypes=c("prc")))

### INVESTIGATE: [M_Inflation_Growth]

vars <- c("M_Inflation_Growth_1","M_Inflation_Growth_10","M_Inflation_Growth_11",
          "M_Inflation_Growth_12","M_Inflation_Growth_2","M_Inflation_Growth_3",
          "M_Inflation_Growth_4","M_Inflation_Growth_5","M_Inflation_Growth_6",
          "M_Inflation_Growth_7","M_Inflation_Growth_8","M_Inflation_Growth_9","M_Inflation_Growth")

# Obtain the KS-statistics for goodness of fit tests.
getKS(datCredit_train_TFD, vars)
#               Variable     KS
# 9   M_Inflation_Growth_6 0.6515
# 8   M_Inflation_Growth_5 0.6508
# 12  M_Inflation_Growth_9 0.6499
# 7   M_Inflation_Growth_4 0.6493
# 6   M_Inflation_Growth_3 0.6488
# 5   M_Inflation_Growth_2 0.6486
# 10  M_Inflation_Growth_7 0.6475
# 13    M_Inflation_Growth 0.6466
# 2  M_Inflation_Growth_10 0.6464
# 3  M_Inflation_Growth_11 0.6456
# 11  M_Inflation_Growth_8 0.6444
# 1   M_Inflation_Growth_1 0.6443
# 4  M_Inflation_Growth_12 0.6438
### RESULTS:  Non-conclusive results as KS fluctuate too much.

# Obtain the concordance for different variables.
getConcs(datCredit_valid_TFD, vars)
#                   Variable Concordance LR_Statistic
# 1: M_Inflation_Growth_12   0.5582303          169
# 2: M_Inflation_Growth_11   0.5554290          158
# 3: M_Inflation_Growth_10   0.5541576          166
# 4:  M_Inflation_Growth_8   0.5534695          172
# 5:  M_Inflation_Growth_7   0.5533126          178
# 6:  M_Inflation_Growth_6   0.5529680          182
# 7:  M_Inflation_Growth_9   0.5517474          170
# 8:  M_Inflation_Growth_5   0.5505816          177
# 9:  M_Inflation_Growth_4   0.5466313          164
# 10:  M_Inflation_Growth_3   0.5432397          138
# 11:  M_Inflation_Growth_2   0.5333219          111
# 12:  M_Inflation_Growth_1   0.5275290           83
# 13:    M_Inflation_Growth   0.5202318           55

### RESULTS:  [M_Inflation_Growth_12] has the best concordance, which decrease as
###           the lag increase.

### CONCLUSION: Replace [M_Inflation_Growth] with [M_Inflation_Growth_12]

varlist <- vecChange(varlist,Remove="M_Inflation_Growth",Add=data.table(vars=c("M_Inflation_Growth_12"), vartypes=c("prc")))


### INVESTIGATE: [M_Repo_Rate]

vars <- c("M_Repo_Rate_1","M_Repo_Rate_10","M_Repo_Rate_11","M_Repo_Rate_12",
          "M_Repo_Rate_2","M_Repo_Rate_3","M_Repo_Rate_4","M_Repo_Rate_5",
          "M_Repo_Rate_6","M_Repo_Rate_7","M_Repo_Rate_8","M_Repo_Rate_9", "M_Repo_Rate")

# Obtain the KS-statistics for goodness of fit tests.
getKS(datCredit_train_TFD, vars)
#         Variable     KS
# 13    M_Repo_Rate 0.6519
# 7   M_Repo_Rate_4 0.6505
# 2  M_Repo_Rate_10 0.6499
# 6   M_Repo_Rate_3 0.6499
# 9   M_Repo_Rate_6 0.6485
# 3  M_Repo_Rate_11 0.6475
# 1   M_Repo_Rate_1 0.6474
# 10  M_Repo_Rate_7 0.6466
# 5   M_Repo_Rate_2 0.6461
# 8   M_Repo_Rate_5 0.6449
# 4  M_Repo_Rate_12 0.6447
# 12  M_Repo_Rate_9 0.6444
# 11  M_Repo_Rate_8 0.6420
### RESULTS:  Non-conclusive results as KS fluctuate too much.

# Obtain the concordance for different variables.
getConcs(datCredit_valid_TFD, vars)
#           Variable Concordance LR_Statistic
# 1: M_Repo_Rate_12   0.5419933          154
# 2: M_Repo_Rate_11   0.5408448          144
# 3: M_Repo_Rate_10   0.5387710          131
# 4:  M_Repo_Rate_9   0.5378973          121
# 5:  M_Repo_Rate_8   0.5366254          116
# 6:  M_Repo_Rate_7   0.5363693          109
# 7:  M_Repo_Rate_6   0.5342596           99
# 8:  M_Repo_Rate_5   0.5337216           94
# 9:  M_Repo_Rate_4   0.5313324           84
# 10:  M_Repo_Rate_3   0.5298802           72
# 11:  M_Repo_Rate_2   0.5254836           58
# 12:  M_Repo_Rate_1   0.5233577           46
# 13:    M_Repo_Rate   0.5180882           29

### RESULTS:  [M_Repo_Rate_12] has the best concordance, which decrease as
###           the lag increase.

### CONCLUSION: Replace [M_Repo_Rate] with [M_Repo_Rate_12]

varlist <- vecChange(varlist,Remove="M_Repo_Rate",Add=data.table(vars=c("M_Repo_Rate_12"), vartypes=c("prc")))

# ------ 6.3 Test for proportionality in hazard models

# [M_RealIncome_Growth]
cox <- coxph(Surv(Start,End, Default_Ind) ~ M_RealIncome_Growth, datCredit_train_TFD)
sfResiduals(cox,datCredit_train_TFD,"M_RealIncome_Growth", c(50,0.1))
### RESULTS: Minimal trends away from 0, therefore do not reject proportionality.

# [M_DTI_Growth_11]
cox <- coxph(Surv(Start,End, Default_Ind) ~ M_DTI_Growth_11, datCredit_train_TFD)
sfResiduals(cox,datCredit_train_TFD,"M_DTI_Growth_11", c(100,0.05))
### RESULTS: Clear trend away from 0, therefore reject proportionality.

# [M_Inflation_Growth_12]
cox <- coxph(Surv(Start,End, Default_Ind) ~ M_Inflation_Growth_12, datCredit_train_TFD)
sfResiduals(cox,datCredit_train_TFD,"M_Inflation_Growth_12", c(50,0.1))
### RESULTS: Clear trend away from 0 initially, therefore reject proportionality.

# [M_Repo_Rate_12]
cox <- coxph(Surv(Start,End, Default_Ind) ~ M_Repo_Rate_12, datCredit_train_TFD)
sfResiduals(cox,datCredit_train_TFD,"M_Repo_Rate_12", c(50,0.1))
### RESULTS: Clear trend away from 0 initially, therefore reject proportionality.

### CONCLUSION: Keep [M_RealIncome_Growth]

modelVar <- c(modelVar,"M_RealIncome_Growth")

# ------ 7. Final Model variables
modelVar <- data.table(vars = modelVar)
# [1] "g0_Delinq_Ave"                     "slc_acct_pre_lim_perc_imputed_med"
# [3] "slc_acct_roll_ever_24_imputed_med" "InterestRate_Nom"                 
# [5] "InterestRate_Margin"               "InterestRate_Margin_Aggr_Med_4"   
# [7] "Principal"                         "Undrawn_Amt"                      
# [9] "AgeToTerm"                         "ArrearsToBalance_Aggr_Prop"       
# [11] "CuringEvents_Aggr_Prop"            "NewLoans_Aggr_Prop"               
# [13] "M_RealIncome_Growth"

# ------ 7.1 Check correlations of variables to ascertain variables that are similar.

# - Correlation analysis
corGroups <- corrAnalysis(datCredit_train_TFD, modelVar, corrThresh = 0.6, method = 'spearman') # Obtain correlation groups

### RESULTS: 1) [g0_Delinq_Ave], [InterestRate_Margin_Aggr_Med_4] and [ArrearsToBalance_Aggr_Prop] are highly correlated
###          2) [slc_acct_pre_lim_perc_imputed_med] and [Undrawn_Amt]

### INVESTIGATE: [g0_Delinq_Ave], [InterestRate_Margin_Aggr_Med_4] and [ArrearsToBalance_Aggr_Prop]
vars <- c("g0_Delinq_Ave", "InterestRate_Margin_Aggr_Med_4", "ArrearsToBalance_Aggr_Prop")

# Obtain the KS-statistics for goodness of fit tests.
getKS(datCredit_train_TFD, vars)
#                       Variable     KS
# 3     ArrearsToBalance_Aggr_Prop 0.6502
# 2 InterestRate_Margin_Aggr_Med_4 0.6468
# 1                  g0_Delinq_Ave 0.6431
### RESULTS:  Non-conclusive results as KS fluctuate too much.

# Obtain the concordance for different variables.
getConcs(datCredit_valid_TFD, vars)
#                           Variable Concordance LR_Statistic
# 1: InterestRate_Margin_Aggr_Med_4   0.5611257          104
# 2:                  g0_Delinq_Ave   0.5397435           90
# 3:     ArrearsToBalance_Aggr_Prop   0.5367480          118

### RESULTS:  [InterestRate_Margin_Aggr_Med_4] has the best concordance, significantly
###           better than the other two options.

modelVar <- vecChange(modelVar,Remove=c("g0_Delinq_Ave", "ArrearsToBalance_Aggr_Prop"))

### INVESTIGATE: [slc_acct_pre_lim_perc_imputed_med] and [Undrawn_Amt]
vars <- c("slc_acct_pre_lim_perc_imputed_med", "Undrawn_Amt")

# Obtain the KS-statistics for goodness of fit tests.
getKS(datCredit_train_TFD, vars)
#                           Variable     KS
# 1 slc_acct_pre_lim_perc_imputed_med 0.6497
# 2                       Undrawn_Amt 0.6486
### RESULTS:  Non-conclusive results as KS fluctuate too much.

# Obtain the concordance for different variables.
getConcs(datCredit_valid_TFD, vars)
#                             Variable Concordance LR_Statistic
# 1:                       Undrawn_Amt   0.6883511         2126
# 2: slc_acct_pre_lim_perc_imputed_med   0.6360895         2865

### RESULTS:  [Undrawn_Amt] has higher concordance than [slc_acct_pre_lim_perc_imputed_med].

modelVar <- vecChange(modelVar,Remove=c("slc_acct_pre_lim_perc_imputed_med"))

### INVESTIGATE: Explore variables

# Obtain the KS-statistics for goodness of fit tests.
getKS(datCredit_train_TFD, modelVar$vars)
#                           Variable     KS
# 3                   InterestRate_Nom 0.6500
# 4                InterestRate_Margin 0.6488
# 5     InterestRate_Margin_Aggr_Med_4 0.6481
# 7                        Undrawn_Amt 0.6480
# 8                          AgeToTerm 0.6479
# 10                NewLoans_Aggr_Prop 0.6473
# 11               M_RealIncome_Growth 0.6468
# 9             CuringEvents_Aggr_Prop 0.6466
# 1  slc_acct_pre_lim_perc_imputed_med 0.6465
# 6                          Principal 0.6457
# 2  slc_acct_roll_ever_24_imputed_med 0.6306


# Obtain the concordance for different variables.
getConcs(datCredit_valid_TFD, modelVar$vars)
#                             Variable Concordance LR_Statistic
# 1: slc_acct_roll_ever_24_imputed_med   0.8588785        18196
# 2:                       Undrawn_Amt   0.6883511         2126
# 3: slc_acct_pre_lim_perc_imputed_med   0.6360895         2865
# 4:                         AgeToTerm   0.6179618           96
# 5:                         Principal   0.5827355          335
# 6:    InterestRate_Margin_Aggr_Med_4   0.5611257          104
# 7:                NewLoans_Aggr_Prop   0.5571220          194
# 8:                  InterestRate_Nom   0.5497724          165
# 9:               InterestRate_Margin   0.5452677          174
# 10:            CuringEvents_Aggr_Prop   0.5359976           82
# 11:               M_RealIncome_Growth   0.5341913           68

# ------ Final test of TFD Model
# Create the formula dynamically
cox_formula_TFD <- as.formula(paste("Surv(Start, End, Default_Ind) ~", paste(modelVar$vars, collapse = " + ")))

# Fit the Cox model
cox_TFD <- coxph(cox_formula_TFD, data = datCredit_train_TFD)
c <- coefficients(cox_TFD)
c <- data.table(var = names(c), coef = c, expCoef = round(exp(c),2))
#                             var               coef                                     expCoef
# 1: slc_acct_roll_ever_24_imputed_med    1.8476538050625                                        6.34
# 2:                       Undrawn_Amt   -0.0000007237391                                        1.00
# 3: slc_acct_pre_lim_perc_imputed_med  -54.3095078384995                                        0.00
# 4:                         AgeToTerm    0.1132746140864                                        1.12
# 5:                         Principal    0.0000002068101                                        1.00
# 6:    InterestRate_Margin_Aggr_Med_4 -102.1416180753126                                        0.00
# 7:                NewLoans_Aggr_Prop   91.5560480545481 5784774777758415903262608220248024662406.00
# 8:                  InterestRate_Nom    3.4494865119585                                       31.48
# 9:               InterestRate_Margin   -0.0548564444060                                        0.95
# 10:            CuringEvents_Aggr_Prop -318.9388154248800                                        0.00
# 11:               M_RealIncome_Growth   -0.7269168824253                                        0.48

concordance(cox_TFD)
# 0.9278

tdROC(datCredit_valid_TFD, cox_TFD, month_Start = 0, month=3, numDigits=1);gc()
tdROC(datCredit_valid_TFD, cox_TFD, month_Start = 4, month=12, numDigits=1);gc()




cox.zph(cox_TFD)
