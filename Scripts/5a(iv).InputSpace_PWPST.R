# ============================================ INPUT SPACE =========================================
# Divide data into thematic groups and perform data analysis on them to compile an input space for the TPWPST model
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
#   - datCredit_train_PWPST | Prepared from script 3b
#   - datCredit_valid_PWPST | Prepared from script 3b

#
# -- Outputs:
#   - Input_Space
# ------------------------------------------------------------------------------------------------------

# ------ 1. Preliminaries
# - Confirm prepared datasets are loaded into memory
if (!exists('datCredit_train_PWPST')) unpack.ffdf(paste0(genPath,"creditdata_train_PWPST"), tempPath);gc()
if (!exists('datCredit_valid_PWPST')) unpack.ffdf(paste0(genPath,"creditdata_valid_PWPST"), tempPath);gc()

# BS: will work these changes into script c once master's is done.
datCredit_train_PWPST[, slc_acct_arr_dir_3_Change_Ind := ifelse(slc_acct_arr_dir_3 != "SAME", 1,0)]
datCredit_valid_PWPST[, slc_acct_arr_dir_3_Change_Ind := ifelse(slc_acct_arr_dir_3 != "SAME", 1,0)]

### HOMEWORK: Create [Removed], [slc_acct_roll_ever_24_imputed_med_f]
datCredit_train_PWPST[,Removed := ifelse(Date==PerfSpell_Max_Date,T,F)]
datCredit_valid_PWPST[,Removed := ifelse(Date==PerfSpell_Max_Date,T,F)]


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
concTable <- function(data_train, data_valid, variables) {
  # Use lapply to efficiently compute concordances for all univariate models
  results <- lapply(variables, function(var) {
    formula <- as.formula(paste0("Surv(Start, End, Default_Ind) ~ PerfSpell_Num +", var))
    tryCatch({
      model <- coxph(formula,id=PerfSpell_Key, data = data_train)# Fit Cox model
      c <- concordance(model, newdata=data_valid)
      conc <- as.numeric(c[1])# Extract concordance
      sd <- sqrt(c$var)# Extract concordance variability
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
### NOTE: Bivariate models are created with one always being [PerfSpell_Num]

# Table B statistics of single variable cox ph models
# Used to compare the goodness of fit of variables
csTable <- function(data,variables,seedVal=1,numIt=3){
  
  # Simulate null distribution if seedVal is not NA
  # Initialize results
  results <- data.frame(Variable = variables, B_Statistic = NA_real_)
  
  # Simulate null distribution if seedVal is not NA
  if (!is.na(seedVal)) {
    set.seed(seedVal, kind = "Mersenne-Twister")
    #null_distribution <- rexp(nrow(data))
    
    # Vectorized calculation for B statistics
    results$B_Statistic <- sapply(variables, function(var) {
      formula <- as.formula(paste0("Surv(Start, End, Default_Ind) ~ PerfSpell_Num +", var))
      tryCatch({
        model <- coxph(formula, data = data, id=PerfSpell_Key)  # Fit a univariate Cox model
        GoF_CoxSnell_KS(model, data, GraphInd=F)$Stat  # Calculate B statistic
      }, warning = function(w) {
        cat("Warning: ", w$message, " for variable: ", var, "\n")
        NA
      }, error = function(e) {
        cat("Error: ", e$message, " for variable: ", var, "\n")
        NA
      })
    })
    
    # Sort results by B statistic in descending order
    results <- results[order(-results$B_Statistic, na.last = TRUE), ]
    
    # Return results and range of B statistics
    return(list(Table = results, Range = diff(range(results$B_Statistic, na.rm = TRUE))))
    
  } else {
    # Perform iterative B calculation when seedVal is NA
    matResults <- matrix(NA, nrow = length(variables), ncol = numIt,
                         dimnames = list(variables,
                                         paste0("Iteration_", 1:numIt))) %>%
      as.data.table()
    
    for (it in seq_len(numIt)) {
      #null_distribution <- rexp(nrow(data))
      
      # Vectorized iteration for B statistics
      matResults[, it] <- sapply(variables, function(var) {
        formula <- as.formula(paste0("Surv(Start, End, Default_Ind) ~ PerfSpell_Num +", var))
        tryCatch({
          model <- coxph(formula, id=PerfSpell_Key, data = data)
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
    
    # Return matrix of B statistics
    return(matResults)
  }
}
### NOTE: Bivariate models are created with one always being [PerfSpell_Num]
#============================================================================================
# ------ 1. Delinquency measures
varlist <- data.table(vars=c("g0_Delinq","g0_Delinq_fac","PerfSpell_g0_Delinq_Num",
                             "Arrears" ,"TimeInDelinqState","g0_Delinq_Any_Aggr_Prop",
                             "g0_Delinq_Ave","slc_acct_arr_dir_3",
                             "slc_acct_roll_ever_24_imputed_med"),
                      vartypes=c("acc", "cat", "acc", "acc", "acc", "dte", "dte",
                                 "cat", "acc"))

#=========================================================================================



# ------ 1.1 Which time window length is the best in calculating Delinquency volatility?

# Initialize variables to be tested
vars <- c("g0_Delinq_SD_4", "g0_Delinq_SD_5", "g0_Delinq_SD_6", "g0_Delinq_SD_9", "g0_Delinq_SD_12")

# Goodness of fit test
csTable(datCredit_train_PWPST,vars)
#           Variable B_Statistic
# 5 g0_Delinq_SD_12      0.6463
# 1  g0_Delinq_SD_4      0.6458
# 2  g0_Delinq_SD_5      0.6324
# 4  g0_Delinq_SD_9      0.6236
# 3  g0_Delinq_SD_6      0.6126

### RESULTS:  [g0_Delinq_SD_12] fits the data the best, slightly better than [g0_Delinq_SD_4]

# Accuracy test
concTable(datCredit_train_PWPST, datCredit_valid_PWPST, vars)
#           Variable Concordance          SD LR_Statistic
# 1:  g0_Delinq_SD_4   0.9806753 0.001401249       107861
# 2:  g0_Delinq_SD_5   0.9733834 0.001754087       117219
# 3:  g0_Delinq_SD_6   0.9538233 0.002316051       116953
# 4:  g0_Delinq_SD_9   0.9212898 0.002923247       110487
# 5: g0_Delinq_SD_12   0.8878487 0.003321776        91619

### RESULTS: As the SD period increase, there is a slight decrease in concordance.
### NOTE: Concordance is extremely high with low variability

### Conclusion: Larger window are less influenced by large changes therefore significant changes
###             are less pronounced. Include [g0_Delinq_SD_4] in the varlist.

varlist <- vecChange(varlist,Add=data.table(vars=c("g0_Delinq_SD_4"), vartypes=c("acc")))

# ------ 1.2 Which variables are highly correlated?

# Correlation analysis
corrAnalysis(datCredit_train_PWPST, varlist[vartypes!="cat"]$vars, corrThresh = 0.6, method = 'spearman') # Obtain correlation groups

### RESULTS:  1) [g0_Delinq_Any_Aggr_Prop] and [g0_Delinq_Ave] with a correlation of 1
###           2) [PerfSpell_g0_Delinq_Num] and [slc_acct_roll_ever_24_imputed_med]
### NOTE: Group 1) are also highly correlated with [g0_Delinq_Any_Aggr_Prop_Lag_3],
###       which is to be expected.
###           3) [g0_Delinq] and [Arrears]

### CONCLUSION: A single variable from each group must be retained while the rest are removed.



# ------ 1.2.1 Which variable should be kept from group 1) [g0_Delinq_Any_Aggr_Prop] and [g0_Delinq_Ave]

vars <- c("g0_Delinq_Any_Aggr_Prop", "g0_Delinq_Ave")

# Goodness of fit
csTable(datCredit_train_PWPST,vars,seedVal = NA)
#                   Variables Iteration_1 Iteration_2 Iteration_3 Iteration_4 Iteration_5 Average
# 1:           g0_Delinq_Ave      0.6796      0.6827      0.6842      0.6808      0.6785 0.68116
# 2: g0_Delinq_Any_Aggr_Prop      0.6824      0.6803      0.6812      0.6793      0.6792 0.68048
# 3:                   Range      0.0028      0.0024      0.0030      0.0015      0.0007 0.00208


### RESULTS: [g0-Delinq_Ave] seems to have the better goodness of fit over 5 iterations.

# Accuracy
concTable(datCredit_train_PWPST, datCredit_valid_PWPST,vars)
#                   Variable Concordance          SD LR_Statistic
# 1:           g0_Delinq_Ave   0.6863403 0.003767392         5347
# 2: g0_Delinq_Any_Aggr_Prop   0.6859647 0.003783753         5305

### RESULTS: [g0-Delinq_Ave] has a slightly better concordance.

### CONCLUSION: [g0-Delinq_Ave] seems to outperform [g0_Delinq_Any_Aggr_Prop] and therefore is kept in the model
###             and [g0_Delinq_Any_Aggr_Prop_Lag_3] is removed along with [g0_Delinq_Any_Aggr_Prop] due to the high correlation
###             I expect similar results.

varlist <- vecChange(varlist,Remove=c("g0_Delinq_Any_Aggr_Prop"))



# ------ 1.2.2 Which variable should be kept from group 2) [PerfSpell_g0_Delinq_Num] and [slc_acct_roll_ever_24_imputed_med]

vars <- c("PerfSpell_g0_Delinq_Num", "slc_acct_roll_ever_24_imputed_med")

# Goodness of fit
csTable(datCredit_train_PWPST,vars,seedVal = NA)
#                           Variables Iteration_1 Iteration_2 Iteration_3 Iteration_4 Iteration_5 Average
# 1:           PerfSpell_g0_Delinq_Num      0.6769      0.6760      0.6762      0.6726      0.6756 0.67546
# 2: slc_acct_roll_ever_24_imputed_med      0.6568      0.6587      0.6573      0.6562      0.6556 0.65692
# 3:                             Range      0.0201      0.0173      0.0189      0.0164      0.0200 0.01854
# 
### RESULTS: [PerfSpell_g0_Delinq_Num] seems to have the better goodness of fit over 5 iterations.

# Accuracy
concTable(datCredit_train_PWPST, datCredit_valid_PWPST,vars)
#                             Variable Concordance           SD LR_Statistic
# 1:           PerfSpell_g0_Delinq_Num   0.9268651 0.0007233635        13438
# 2: slc_acct_roll_ever_24_imputed_med   0.8515251 0.0035445945        43093

### RESULTS: [PerfSpell_g0_Delinq_Num] seems to have a significant better concordance with the concordances having low SD's.

### CONCLUSION: Keep [PerfSpell_g0_Delinq_Num] in the model and remove [slc_acct_roll_ever_24_imputed_med]

varlist <- vecChange(varlist,Remove="slc_acct_roll_ever_24_imputed_med")



# ------ 1.3 Which version of [g0_Delinq] should be kept in the model?

# [g0_Delinq]
cox <- coxph(Surv(Start,End,Default_Ind) ~ g0_Delinq + PerfSpell_Num, id=PerfSpell_Key, datCredit_train_PWPST)
summary(cox);rm(cox)
### RESULTS: Beta is unstable with high variability.

# [g0_Delinq_fac]
cox <- coxph(Surv(Start,End,Default_Ind) ~ g0_Delinq_fac + PerfSpell_Num, id=PerfSpell_Key, datCredit_train_PWPST)
summary(cox);rm(cox)
### RESULTS: coef is unstable with high variability.

### INVESTIGATE: Is quasi-complete seperation present?
datCredit_train_PWPST[g0_Delinq==3 & Default_Ind==0, .N]
### RESULTS: 0
datCredit_train_PWPST[g0_Delinq!=3 & Default_Ind==1, .N]
### RESULTS: 0

### CONCLUSION: Quasi-complete separation is present for [g0_Delinq]=3

varlist <- vecChange(varlist,Remove=c("g0_Delinq","g0_Delinq_fac"))

### INVESTIGATE: Should an indicator version of [g0_Delinq] be included in the model.

# An indicator for when the value is greater than 0
datCredit_train_PWPST[,g0_Delinq_Ind := ifelse(g0_Delinq > 0, 1, 0)]
cox <- coxph(Surv(Start,End,Default_Ind) ~ g0_Delinq_Ind + PerfSpell_Num, id=PerfSpell_Key, datCredit_train_PWPST)
summary(cox); rm(cox)
### RESULTS: coef is unstable with high variability.
datCredit_train_PWPST[,g0_Delinq_Ind := NULL]

### INVESTIGATE: Should an lagged version of [g0_Delinq] be included in the model.

# A lag version for [g0_Delinq]
datCredit_train_PWPST[,g0_Delinq_Lag_1 := shift(g0_Delinq,fill=0),by=LoanID]
cox <- coxph(Surv(Start,End,Default_Ind) ~ g0_Delinq_Lag_1 + PerfSpell_Num, id=PerfSpell_Key, datCredit_train_PWPST)
### RESULTS: exp(coef) tends to Inf
datCredit_train_PWPST[,g0_Delinq_Lag_1 := NULL]

# Arrears
cox <- coxph(Surv(Start,End,Default_Ind) ~ Arrears + PerfSpell_Num, id=PerfSpell_Key, datCredit_train_PWPST)
summary(cox); rm(cox)
### RESULTS: Obtain a stable model

### CONCLUSION: Unable to add [Delinq_0] to the model, since the various forms'
###             coef is unstable, however, [Arrears] compiles a seemingly stable
###             model, therefore it can serve as a proxy for g0_Delinq.

# ------ 1.4 How does [PerfSpell_g0_Delinq_Num] compare to [g0_Delinq_SD_4]?

# Initialize variables
vars <- c("PerfSpell_g0_Delinq_Num", "g0_Delinq_SD_4")

# Goodness of fit
csTable(datCredit_train_PWPST,vars,seedVal = NA)
#                 Variables Iteration_1 Iteration_2 Iteration_3    Average
# 1: PerfSpell_g0_Delinq_Num      0.6791      0.6777      0.6762 0.67766667
# 2:          g0_Delinq_SD_4      0.6466      0.6428      0.6459 0.64510000
# 3:                   Range      0.0325      0.0349      0.0303 0.03256667

### RESULTS: [PerfSpell_g0_Delinq_Num] seems to have a much better goodness of fit.

# Accuracy
concTable(datCredit_train_PWPST, datCredit_valid_PWPST,vars)
#                   Variable Concordance           SD LR_Statistic
# 1:          g0_Delinq_SD_4   0.9806753 0.0014012488       107861
# 2: PerfSpell_g0_Delinq_Num   0.9268651 0.0007233635        13438

### RESULTS: [g0_Delinq_SD_4] seems to have better concordance.

### CONCLUSION: Keep [g0_Delinq_SD_4] in the model and remove [PerfSpell_g0_Delinq_Num].

varlist <- vecChange(varlist,Remove="PerfSpell_g0_Delinq_Num")



# ------ 1.4 What is the performance of current thematic variables in univariate models?

# Build thematic model based on remaining delinquency variables.
vars <- c("g0_Delinq_SD_4","TimeInDelinqState","g0_Delinq_Ave",
          "slc_acct_arr_dir_3", "Arrears")

# Test Goodness of fit
csTable(datCredit_train_PWPST,vars)
#             Variable B_Statistic
# 3      g0_Delinq_Ave      0.6791
# 5            Arrears      0.6775
# 1     g0_Delinq_SD_4      0.6458
# 2  TimeInDelinqState          NA
# 4 slc_acct_arr_dir_3          NA

### RESULTS: Fits are close to on another.
### NOTE: [TimeInDelinqState] ran out of iterations and did not converge
### NOTE: [slc_acct_arr_dir_3] exp overflow due to covariates

# ------ 1.4.1 Why does [TimeInDelinqState] not converge?

cox <- coxph(Surv(Start,End,Default_Ind) ~ TimeInDelinqState + PerfSpell_Num, id=PerfSpell_Key, datCredit_train_PWPST)
### RESULTS:  Beta tends to Inf. After some inspection on the data it relates to all
###           Defaulting events starting in a new delinquency state,
###           i.e. [TimeInDelinqState] = 1. Therefore quasi-complete separation
###           seems to be present.

### INVESTIGATE: WHETHER QUASI-COMPLETE SEPERATION IS PRESENT
datCredit_train_PWPST[TimeInDelinqState!=1 & Default_Ind==1, .N]
### RESULTS: 0
datCredit_train_PWPST[TimeInDelinqState==1 & Default_Ind==0, .N]
### RESULTS: 196899

### CONCLUSION: Quasi-complete separation seems to be present for [g0_Delinq]=3
###             and should therefore be removed.

# Test a lagged version of TimeInDelinqState
cox <- coxph(Surv(Start,End,Default_Ind) ~ TimeInDelinqState_Lag_1 + PerfSpell_Num
             , id=PerfSpell_Key, datCredit_train_PWPST)
summary(cox); rm(cox)

### CONCLUSION: Replace old variable with new variable

### RESULTS: Variable is significant with a high concordance of 0.942, thus replace
###           [TimeInDelinqState] with [TimeInDelinqState_Lag_1]

varlist <- vecChange(varlist,Remove="TimeInDelinqState",
                     Add=data.table(vars="TimeInDelinqState_Lag_1",vartypes="acc"))


# ------ 1.4.2 Why does [slc_acct_arr_dir_3] exp overflow?

describe(datCredit_train_PWPST[,slc_acct_arr_dir_3])
# Value            CURING MISSING_DATA      ROLLING         SAME
# Proportion        0.015        0.126        0.022        0.837

describe(datCredit_train_PWPST[Default_Ind==1,slc_acct_arr_dir_3])
# Value            CURING MISSING_DATA      ROLLING         SAME
# Proportion        0.007        0.160        0.783        0.049

### RESULTS:  Although the ROLLING level is not dominant in datCredit_train_PWPST,
###           we can clearly see that it is highly predictive of a default event
###           occurring, given its high proportion if Default_Ind == 1.

# Create an indicator function
datCredit_train_PWPST[, slc_acct_arr_dir_3_ROLLING_Ind := 
                      ifelse(slc_acct_arr_dir_3 == "ROLLING", 1,0)]
cox <- coxph(Surv(Start,End,Default_Ind) ~ slc_acct_arr_dir_3_ROLLING_Ind + 
               PerfSpell_Num, datCredit_train_PWPST, id=PerfSpell_Key)
### RESULTS: exp overflows

# Make the indicator categorical
cox <- coxph(Surv(Start,End,Default_Ind) ~ factor(slc_acct_arr_dir_3_ROLLING_Ind) + 
               PerfSpell_Num, datCredit_train_PWPST, id=PerfSpell_Key)
### RESULTS: exp overflows
datCredit_train_PWPST[,slc_acct_arr_dir_3_ROLLING_Ind := NULL]

# Use an indicator variable for a change in account
cox <- coxph(Surv(Start,End,Default_Ind) ~ slc_acct_arr_dir_3_Change_Ind +
            PerfSpell_Num, datCredit_valid_PWPST, id=PerfSpell_Key)
summary(cox); rm(cox)

### CONCLUSION: Replace old variable with new variable

varlist <- vecChange(varlist,Remove=c("slc_acct_arr_dir_3") ,
                     Add=data.table(vars=c("slc_acct_arr_dir_3_Change_Ind"),
                                    vartypes=c("bin")))



# ------ 1.4.3 What is the performance of current thematic variables in univariate models?

vars <- c("g0_Delinq_SD_4","g0_Delinq_Ave", "Arrears",
          "TimeInDelinqState_Lag_1","slc_acct_arr_dir_3_Change_Ind")

# Goodness of fit
csTable(datCredit_train_PWPST,vars)
#                       Variable B_Statistic
# 2                 g0_Delinq_Ave      0.6791
# 3                       Arrears      0.6775
# 5 slc_acct_arr_dir_3_Change_Ind      0.6741
# 4       TimeInDelinqState_Lag_1      0.6643
# 1                g0_Delinq_SD_4      0.6458

### RESULTS: All variables seems to be a good fit.

# Accuracy
concTable(datCredit_train_PWPST, datCredit_valid_PWPST,vars)
#                         Variable Concordance          SD LR_Statistic
# 1:                g0_Delinq_SD_4   0.9806753 0.001401249       107861
# 2:                       Arrears   0.9450365 0.001432871        20316
# 3:       TimeInDelinqState_Lag_1   0.9219226 0.002705651        43522
# 4: slc_acct_arr_dir_3_Change_Ind   0.8891861 0.001995888        43030
# 5:                 g0_Delinq_Ave   0.6863403 0.003767392         5347

### RESULTS: [g0_Delinq_Ave] has significant less concordance than the other variables.

### CONCLUSION: Leave all variables in the model (including [g0_Delinq_Ave], since it seems to have the best fit for the data).



# ------ 1.5 What is the performance of current thematic cox ph model?

# Build cox model based on all thematic variables
coxDelinq <- coxph(Surv(Start,End,Default_Ind) ~ g0_Delinq_Ave +
                    g0_Delinq_SD_4 + slc_acct_arr_dir_3_Change_Ind +
                    TimeInDelinqState_Lag_1 + Arrears + PerfSpell_Num,
                    id=PerfSpell_Key, data=datCredit_train_PWPST)

summary(coxDelinq)

### RESULTS: All variables are significant (p-value < 0.001)

# Test goodness of fit
GoF_CoxSnell_KS(coxDelinq,datCredit_train_PWPST,GraphInd = FALSE) # 0.6436

### RESULTS: Goodness of fit is a bit low

# Test accuracy
concordance(coxDelinq, newdata=datCredit_valid_PWPST) # Concordance= 0.9952 se= 0.0003395

### RESUTLS: extremely good prediction accuracy.

### CONCLUSION: Keep all variables in the model

# House keeping
rm(coxDelinq)

# (0,3) (4,12) (13,24) (0,12) (0,36)
#timedROC(datCredit_valid_PWPST, coxDelinq_valid, month_Start=0, month_End=36,
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
#                    g0_Delinq_Any_Aggr_Prop_Lag_5, datCredit_train_PWPST)
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

modelVar <- c("g0_Delinq_SD_4", "Arrears", "g0_Delinq_Ave",
              "TimeInDelinqState_Lag_1", "slc_acct_arr_dir_3_Change_Ind")

#============================================================================================
# ------ 2. Engineered measures
varlist <- data.table(vars=c("slc_acct_pre_lim_perc_imputed_med",
                             "slc_acct_prepaid_perc_dir_12_imputed_med",
                             "pmnt_method_grp") ,
                      vartypes=c("prc", "dec", "cat"))

#=========================================================================================

# ------ 2.1 Which variables should be removed due to high correlation?

# - Correlation analysis
corGroups <- corrAnalysis(datCredit_train_PWPST, varlist[vartypes!="cat"]$vars, corrThresh = 0.6, method = 'spearman') # Obtain correlation groups

### RESULTS:  1) [slc_acct_pre_lim_perc_imputed_med] and [slc_acct_prepaid_perc_dir_12_imputed_med] with a correlation of 0.82

# Initialize variables to be tested
vars <- c("slc_acct_pre_lim_perc_imputed_med", "slc_acct_prepaid_perc_dir_12_imputed_med")

# Compare goodness of fit of different variables
csTable(datCredit_train_PWPST,vars)
#                                   Variable B_Statistic
# 1        slc_acct_pre_lim_perc_imputed_med      0.6801
# 2 slc_acct_prepaid_perc_dir_12_imputed_med          NA

### RESULTS:  [slc_acct_prepaid_perc_dir_12_imputed_med] ran out of iterations and did not converge

### INVESTIGATE: Why did slc_acct_prepaid_perc_dir_12_imputed_med] ran out of iterations and not converge?

hist(datCredit_train_PWPST[,slc_acct_prepaid_perc_dir_12_imputed_med])
### RESULTS: A significant portion of the values are 0 with possible extreme values.

# Proportion of values being 0
datCredit_train_PWPST[slc_acct_prepaid_perc_dir_12_imputed_med==0,.N]/datCredit_train_PWPST[,.N]
# 0.7463936

hist(datCredit_train_PWPST[Default_Ind==1, slc_acct_prepaid_perc_dir_12_imputed_med])
### Similar distribution with most of the values being 0, but all values are below 3.5

# Proportion of values being 0 | Default_Ind == 1
datCredit_train_PWPST[slc_acct_prepaid_perc_dir_12_imputed_med==0 & Default_Ind == 1,.N]/datCredit_train_PWPST[Default_Ind == 1,.N]
# 0.9984765

# Proportion of values defaulted | slc_acct_prepaid_perc_dir_12_imputed_med == 0
datCredit_train_PWPST[slc_acct_prepaid_perc_dir_12_imputed_med==0 & Default_Ind == 1,.N]/datCredit_train_PWPST[slc_acct_prepaid_perc_dir_12_imputed_med==0,.N]
# 0.003805984

# RESULTS: (Default_Ind == 1) => (slc_acct_prepaid_perc_dir_12_imputed_med == 0)

### CONCLUSION: Quasi-complete separation seems to be present for [slc_acct_prepaid_perc_dir_12_imputed_med],
###             therefore remove it from the varlist

varlist <- vecChange(varlist,Remove="slc_acct_prepaid_perc_dir_12_imputed_med")

# ------ 2.2 What is the predictive power of the variables left in varlist?

# Initialize variables to be tested
vars <- c("pmnt_method_grp", "slc_acct_pre_lim_perc_imputed_med")

csTable(datCredit_train_PWPST,vars)
#                             Variable B_Statistic
# 2 slc_acct_pre_lim_perc_imputed_med      0.6781
# 1                   pmnt_method_grp      0.6745

### RESULTS:  [slc_acct_pre_lim_perc_imputed_med] fits the better and there is an improvement for [pmnt_method_grp]

concTable(datCredit_train_PWPST, datCredit_valid_PWPST,vars)
#                             Variable Concordance          SD LR_Statistic
# 1:                   pmnt_method_grp   0.7546841 0.003472790        14162
# 2: slc_acct_pre_lim_perc_imputed_med   0.7241381 0.002584756        10364

### RESULTS: [pmnt_method_grp]  has a much better concordance than [slc_acct_pre_lim_perc_imputed_med],
###           but both concordances are reasonably high.

### Conclusion: All values should be kept in the model



# ------ 2.3 How predictive is a single model based on the thematic variables?

# Initialize variables
vars <- c("pmnt_method_grp", "slc_acct_pre_lim_perc_imputed_med")

# Build Cox model based on each variable
coxEngineered <- coxph(Surv(Start, End, Default_Ind) ~ pmnt_method_grp +
                               slc_acct_pre_lim_perc_imputed_med, id=LoanID,
                             data=datCredit_train_PWPST)
summary(coxEngineered)

### RESULTS: Coefficients for pmnt_method_grp are noticably high.

# Goodness of fit
GoF_CoxSnell_KS(coxEngineered,datCredit_train_PWPST, GraphInd=FALSE) # 0.675

### RESULTS: slight decrease in goodness of fit, although it could be due to random sampling

# Accuracy
concordance(coxEngineered,newdata=datCredit_valid_PWPST) # Concordance= 0.784 se= 0.002925

# # Time dependent AUC
# timedROC(datCredit_valid_PWPST, coxEngineered_valid, month_Start=0, month_End=36,
#          fld_ID="LoanID", fld_Event="Default_Ind",fld_StartTime="Start",
#          fld_EndTime="End", numDigits=0, Graph=FALSE) # 0.7040435

### CONCLUSION: Due to increase in Concordance keep all current variables in the model.

# House keeping
rm(coxEngineered); gc()

#============================================================================================

modelVar <- c(modelVar,"slc_acct_pre_lim_perc_imputed_med","pmnt_method_grp")

#============================================================================================
# ------ 3. Interest Rate
varlist <- data.table(vars=c("InterestRate_Nom", "InterestRate_Margin"),
                      vartypes=c("prc","prc"))

#=========================================================================================

# ------ 3.1 How correlated are the two variables?

# - Correlation analysis
corrAnalysis(datCredit_train_PWPST, varlist[vartypes!="cat"]$vars, corrThresh = 0.6, method = 'spearman') # Obtain correlation groups

### RESULTS: Correlation of 0.41, therefore no significant correlations.

# ------ 3.1.1 How does the variables compare in terms of their predictive power?

# Obtain the B-statistics for goodness of fit tests.
csTable(datCredit_train_PWPST, varlist$vars, seedVal=NA)
#               Variables Iteration_1 Iteration_2 Iteration_3     Average
# 1:    InterestRate_Nom      0.6828      0.6791      0.6823 0.681400000
# 2: InterestRate_Margin      0.6819      0.6797      0.6807 0.680766667
# 3:               Range      0.0009      0.0006      0.0016 0.001033333

### RESULTS:  [InterestRate_Nom] is a better fit than [InterestRate_Margin]

# Obtain the concordance for different variables.
concTable(datCredit_train_PWPST, datCredit_valid_PWPST, varlist$vars)
#               Variable Concordance          SD LR_Statistic
# 1:    InterestRate_Nom   0.7009491 0.003957952         5460
# 2: InterestRate_Margin   0.6515720 0.004432220         4725

### RESULTS:  [InterestRate_Nom] has a much higher concordance than [InterestRate_Margin].

### CONCLUSION: Variables are quite similar in predictive power, however are dissimilar
###              in correlation, therefore both can be kept in the model.



# ------ 3.2 How does [InterestRate_Margin] compare with [InterestRate_Margin_imputed_bin]?

# Initialize variables
vars <- c("InterestRate_Margin", "InterestRate_Margin_imputed_bin")

# Obtain the B-statistics for goodness of fit tests.
csTable(datCredit_train_PWPST, vars, seedVal=NA)
#                           Variables Iteration_1 Iteration_2 Iteration_3     Average
# 1:             InterestRate_Margin      0.6810      0.6829      0.6815 0.681800000
# 2: InterestRate_Margin_imputed_bin      0.6788      0.6790      0.6754 0.677733333
# 3:                           Range      0.0022      0.0039      0.0061 0.004066667

### RESULTS: [InterestRate_Margin] has a better fit than [InterestRate_Margin_imputed_bin]

# Obtain the concordance for different variables.
concTable(datCredit_train_PWPST, datCredit_valid_PWPST, vars)
#                           Variable Concordance          SD LR_Statistic
# 1: InterestRate_Margin_imputed_bin   0.6520048 0.004361896         4905
# 2:             InterestRate_Margin   0.6515720 0.004432220         4725

### RESULTS:  [InterestRate_Margin_imputed_bin] has the higher concordance.

### CONCLUSIONS: Do not replace [InterestRate_Margin] with [InterestRate_Margin_imputed_bin],
###               they have similar concordances (differ with 0.0004328), but the B-statistic
###               differ more and [InterestRate_Margin] is consistently better.


# ------ 3.3 How does [InterestRate_Margin], InterestRate_Margin_Aggr_Med and [M_Repo_Rate] compare with one another?

vars <- c("InterestRate_Margin", "M_Repo_Rate", "InterestRate_Margin_Aggr_Med")

# Evaluate their correlation
corrAnalysis(datCredit_train_PWPST,vars)
### RESULTS: no significant correlations were detected between the variables.

# Obtain the B-statistics for goodness of fit tests.
csTable(datCredit_train_PWPST, vars, seedVal=NA)
#                       Variables Iteration_1 Iteration_2 Iteration_3   Average
# 1:                  M_Repo_Rate      0.6787      0.6854      0.6817 0.6819333
# 2:          InterestRate_Margin      0.6804      0.6815      0.6803 0.6807333
# 3: InterestRate_Margin_Aggr_Med      0.6797      0.6779      0.6777 0.6784333
# 4:                        Range      0.0017      0.0075      0.0040 0.0044000
### RESULTS:  [M_Repo_Rate] seems to be have the best fit.

# Obtain the concordance for different variables.
concTable(datCredit_train_PWPST, datCredit_valid_PWPST, vars)
#                       Variable Concordance          SD LR_Statistic
# 1: InterestRate_Margin_Aggr_Med   0.6909229 0.003841607         5340
# 2:                  M_Repo_Rate   0.6785476 0.004009584         4892
# 3:          InterestRate_Margin   0.6515720 0.004432220         4725

### RESULTS:  [InterestRate_Margin_Aggr_Med] has the highest concordance.

### CONCLUSIONS: Use [M_Repo_Rate] in the model, since it seems to have the best fit and second best concordance.

varlist <- vecChange(varlist,Remove=c("InterestRate_Margin_Aggr_Med", "InterestRate_Margin"), Add=data.table(vars="M_Repo_Rate", vartypes="prc"))


# ------ 3.2 Should we add lagging variables to the model?

vars <- c("M_Repo_Rate","M_Repo_Rate_1","M_Repo_Rate_12","M_Repo_Rate_3",
          "M_Repo_Rate_6","M_Repo_Rate_9")

# Obtain the B-statistics for goodness of fit tests.
csTable(datCredit_train_PWPST, vars, seedVal=NA)
#         Variables Iteration_1 Iteration_2 Iteration_3     Average
# 1: M_Repo_Rate_12      0.6827      0.6825      0.6837 0.682966667
# 2:    M_Repo_Rate      0.6812      0.6821      0.6819 0.681733333
# 3:  M_Repo_Rate_3      0.6805      0.6812      0.6814 0.681033333
# 4:  M_Repo_Rate_9      0.6770      0.6814      0.6816 0.680000000
# 5:  M_Repo_Rate_1      0.6800      0.6784      0.6762 0.678200000
# 6:  M_Repo_Rate_6      0.6767      0.6790      0.6768 0.677500000
# 7:          Range      0.0060      0.0041      0.0075 0.005866667

### RESULTS:  [InterestRate_Margin_Aggr_Med_12] and [InterestRate_Margin_Aggr_Med_3] seems to be better fits than [InterestRate_Margin_Aggr_Med]

# Obtain the concordance for different variables.
concTable(datCredit_train_PWPST, datCredit_valid_PWPST, vars)
#         Variable Concordance          SD LR_Statistic
# 1:  M_Repo_Rate_6   0.6848303 0.003886111         5294
# 2:  M_Repo_Rate_9   0.6845921 0.003878557         5325
# 3:  M_Repo_Rate_3   0.6840043 0.003913224         5170
# 4:  M_Repo_Rate_1   0.6815728 0.003965619         4999
# 5: M_Repo_Rate_12   0.6811567 0.003867908         5377
# 6:    M_Repo_Rate   0.6785476 0.004009584         4892

### RESULTS:  Difficult to make a conclusion based on the contrasting goodness of
###           fit and accuracy evidence. [M_Repo_Rate_9] is the most likely on to
###           include given its high concordance and 4th highest B-statistic.

### CONCLUSION: Add [M_Repo_Rate_9] to the model concidering the results.

varlist <- vecChange(varlist,Add=data.table(vars=c("M_Repo_Rate_9"),
                                            vartypes=c("prc")))



# ------ 3.3 What is the performance of current thematic variables in univariate models?

vars <- c("InterestRate_Nom", "M_Repo_Rate", "M_Repo_Rate_9")

# Goodness of fit
csTable(datCredit_train_PWPST,vars)
#           Variable B_Statistic
# 1 InterestRate_Nom      0.6807
# 2      M_Repo_Rate      0.6788
# 3    M_Repo_Rate_9      0.6780

### RESULTS: Variables seem to be a good fit.

# Accuracy
concTable(datCredit_train_PWPST, datCredit_valid_PWPST,vars)
#           Variable Concordance          SD LR_Statistic
# 1: InterestRate_Nom   0.7009491 0.003957952         5460
# 2:    M_Repo_Rate_9   0.6845921 0.003878557         5325
# 3:      M_Repo_Rate   0.6785476 0.004009584         4892

### RESULTS: Concordances are quite high.
 

# ------ 3.4 What is the performance of current thematic cox ph model?

# Build cox model based on all thematic variables
coxInterest <- coxph(Surv(Start,End,Default_Ind) ~ InterestRate_Nom +
                      M_Repo_Rate + M_Repo_Rate_9 + PerfSpell_Num,
                     id=PerfSpell_Key, data=datCredit_train_PWPST)
summary(coxInterest)

### RESULTS: All variables are significant.

# Kolmogorov-Smirnof of coxDelinq
GoF_CoxSnell_KS(coxInterest,datCredit_train_PWPST,GraphInd = FALSE) # 0.6808
### RESULTS: Good fit

# Accuracy
concordance(coxInterest, newdata=datCredit_valid_PWPST) # Concordance= 0.7092 se= 0.003822
### RESULTS: Good accuracy

timedROC(datCredit_valid_PWPST, coxInterest_valid, month_Start=0, month_End=36,
         fld_ID="LoanID", fld_Event="Default_Ind",fld_StartTime="Start",
         fld_EndTime="End", numDigits=0, Graph=FALSE)
# AUC: 0.5214215

# House keeping
rm(coxInterest)
#===========================================================================================

modelVar <- c(modelVar,"InterestRate_Nom", "M_Repo_Rate", "M_Repo_Rate_9")

#============================================================================================
# ------ 4. General
varlist <- data.table(vars=c("Balance","Instalment",
                             "Principal","Undrawn_Amt",
                             "LN_TPE","Term"),
                      vartypes=c("int", "fin", "fin","int","cat", "int"))

#=========================================================================================
# ------ 4.1 Which variables are highly correlated in the group?

# - Correlation analysis
corrAnalysis(datCredit_train_PWPST, varlist[vartypes!="cat"]$vars, corrThresh = 0.6, method = 'spearman') # Obtain correlation groups

### RESULTS: 1) [Balance], [Installment] and [Principal] are highly correlated.

# ------ 4.2 Which variable in the correlated group should be kept in the model?

# Initialize variable
vars <- c("Balance", "Instalment", "Principal")

# Goodness of fit test of variables
csTable(datCredit_train_PWPST, vars, seedVal=NA)
#     Variables Iteration_1 Iteration_2 Iteration_3     Average
# 1:    Balance      0.6828      0.6798      0.6837 0.682100000
# 2: Instalment      0.6818      0.6821      0.6787 0.680866667
# 3:  Principal      0.6789      0.6806      0.6801 0.679866667
# 4:      Range      0.0039      0.0023      0.0050 0.003733333

### RESULTS:  Difficult to make a conclusion with varying B-statistics.

# Obtain the concordance for different variables.
concTable(datCredit_train_PWPST, datCredit_valid_PWPST, vars)
#       Variable Concordance          SD LR_Statistic
# 1:  Principal   0.6673059 0.004317894         4955
# 2: Instalment   0.6528624 0.004529993         4621
# 3:    Balance   0.6416711 0.004499720         4475

## RESULTS: [Principal] has the highest concordance

### CONCOLUSION:  Considering Principal has the best concordance, but the worse fit
###               complicates the decision.

# ------ 4.2.1 Can [Balance] and [Instalment] be replaced with [InstalmentToBalance_Aggr_Prop]?

# Initialize variable
vars <- c("Balance", "Instalment", "InstalmentToBalance_Aggr_Prop")

# Goodness of fit test of variables
csTable(datCredit_train_PWPST, vars, seedVal=NA)
#                       Variables Iteration_1 Iteration_2 Iteration_3     Average
# 1:                       Balance      0.6787      0.6811      0.6791 0.679633333
# 2: InstalmentToBalance_Aggr_Prop      0.6799      0.6830      0.6756 0.679500000
# 3:                    Instalment      0.6779      0.6788      0.6815 0.679400000
# 4:                         Range      0.0020      0.0042      0.0059 0.004033333

### RESULTS:  Inconclusive results from iterations.

# Accuracy of variables
concTable(datCredit_train_PWPST, datCredit_valid_PWPST, vars)
#                         Variable Concordance          SD LR_Statistic
# 1: InstalmentToBalance_Aggr_Prop   0.6708871 0.004148979         4637
# 2:                    Instalment   0.6528624 0.004529993         4621
# 3:                       Balance   0.6416711 0.004499720         4475

### RESULTS:  [InstalmentToBalance_Aggr_Prop] has a better concordance than [Instalment] and [Balance].

### CONCLUSION: [InstalmentToBalance_Aggr_Prop] can replace [Instalment] and [Balance].


# ------ 4.2.1 How does [InstalmentToBalance_Aggr_Prop] compare with [Principal]?
# Initialize variable
vars <- c("Principal", "InstalmentToBalance_Aggr_Prop")

# Goodness of fit test of variables
csTable(datCredit_train_PWPST, vars, seedVal=NA)
#                       Variables Iteration_1 Iteration_2 Iteration_3     Average
# 1:                     Principal      0.6803      0.6797      0.6850 0.681666667
# 2: InstalmentToBalance_Aggr_Prop      0.6786      0.6817      0.6781 0.679466667
# 3:                         Range      0.0017      0.0020      0.0069 0.003533333

### RESULTS:  Inconclusive results from iterations.

# Accuracy of variables
concTable(datCredit_train_PWPST, datCredit_valid_PWPST, vars)
#                         Variable Concordance          SD LR_Statistic
# 1: InstalmentToBalance_Aggr_Prop   0.6708871 0.004148979         4637
# 2:                     Principal   0.6673059 0.004317894         4955

### RESULTS:  [InstalmentToBalance_Aggr_Prop] has a better concordance than [Principal].

### CONCLUSION: [InstalmentToBalance_Aggr_Prop] can replace [Principal].

varlist <- vecChange(varlist, Remove=c("Balance", "Instalment", "Principal"), Add = data.table(vars="InstalmentToBalance_Aggr_Prop", vartypes="prc"))



# ------ 4.2.2 How does [ArrearsToBalance_Aggr_Prop] compare with [Arrears] and [Balance]?

# Initialize variables
vars <- c("ArrearsToBalance_Aggr_Prop", "Arrears", "Balance")


# Goodness of fit test of variables
csTable(datCredit_train_PWPST, vars, seedVal=NA)
#                     Variables Iteration_1 Iteration_2 Iteration_3   Average
# 1: ArrearsToBalance_Aggr_Prop      0.6824      0.6800      0.6813 0.6812333
# 2:                    Arrears      0.6820      0.6804      0.6786 0.6803333
# 3:                    Balance      0.6822      0.6797      0.6785 0.6801333
# 4:                      Range      0.0004      0.0007      0.0028 0.0013000

### RESULTS:  [ArrearsToBalance_Aggr_Prop] seem to be the best fit for the model.

# Accuracy of variables
concTable(datCredit_train_PWPST, datCredit_valid_PWPST, vars)
#                       Variable Concordance          SD LR_Statistic
# 1:                    Arrears   0.9450365 0.001432871        20316
# 2: ArrearsToBalance_Aggr_Prop   0.6814453 0.003817135         5413
# 3:                    Balance   0.6416711 0.004499720         4475

### RESULTS:  [Arrears] has much better concordance than [ArrearsToBalance_Aggr_Prop] and [Balance].

### CONCLUSION: [ArrearsToBalance_Aggr_Prop] should not be included in the model.



# ------ 4.3 How does [Term] compare to [AgeToTerm_Aggr_Mean]?

# Initialize variables
vars <- c("Term", "AgeToTerm_Aggr_Mean")

# Goodness of fit
csTable(datCredit_train_PWPST,vars,seedVal = NA)
#             Variables Iteration_1 Iteration_2 Iteration_3     Average
# 1:                Term      0.6817      0.6843      0.6801 0.682033333
# 2: AgeToTerm_Aggr_Mean      0.6786      0.6816      0.6782 0.679466667
# 3:               Range      0.0031      0.0027      0.0019 0.002566667

### RESULTS: [Term] seems to consistently outfit [AgeToTerm_Aggr_Mean].

# Accuracy
concTable(datCredit_train_PWPST, datCredit_valid_PWPST,vars)
#               Variable Concordance          SD LR_Statistic
# 1: AgeToTerm_Aggr_Mean   0.6766470 0.003819651         4962
# 2:                Term   0.6254012 0.003779138         4362

### RESULTS: [AgeToTerm_Aggr_Mean] has a better concordance than [Term].

### CONCLUSION: Include [AgeToTerm] since it contains more information than [Term].

varlist <- vecChange(varlist,Remove=c("Term"), Add=data.table(vars="AgeToTerm", vartypes="dec"))

# ------ 4.4 What is the predictive power of the variables in univariate models?

vars <- c("InstalmentToBalance_Aggr_Prop", "Undrawn_Amt", "LN_TPE", "AgeToTerm")

# Goodness of fit
csTable(datCredit_train_PWPST,vars, seedVal=NA)
#                         Variables Iteration_1 Iteration_2 Iteration_3   Average
# 1:                        LN_TPE      0.6823      0.6819      0.6786 0.6809333
# 2: InstalmentToBalance_Aggr_Prop      0.6811      0.6766      0.6823 0.6800000
# 3:                     AgeToTerm      0.6803      0.6783      0.6812 0.6799333
# 4:                   Undrawn_Amt      0.6796      0.6766      0.6808 0.6790000
# 5:                         Range      0.0027      0.0053      0.0037 0.0039000

### RESULTS: Each bivariate model has suitable B-statistics.

# Accuracy
concTable(datCredit_train_PWPST, datCredit_valid_PWPST,vars)
#                         Variable Concordance          SD LR_Statistic
# 1:                   Undrawn_Amt   0.7539662 0.002225233         6580
# 2: InstalmentToBalance_Aggr_Prop   0.6708871 0.004148979         4637
# 3:                     AgeToTerm   0.6437042 0.004368801         4728
# 4:                        LN_TPE   0.6325628 0.003571573         4523

### RESULTS:  Each bivariate model has suitable concordances.

# ------ 4.5 What is the predictive power of the variables in a single model?

# Goodness of fit
coxGen <- coxph(as.formula(paste0("Surv(Start,End,Default_Ind) ~ PerfSpell_Num + ",
                                  paste(vars,collapse=" + "))), id=PerfSpell_Num,
                      datCredit_train_PWPST)

summary(coxGen)# [InstalmentToBalance_Aggr_Prop] and [Undrawn_Amt] are not significant (p-value < 0.001)

GoF_CoxSnell_KS(coxGen,datCredit_train_PWPST,GraphInd = FALSE) # 0.6762
### RESULTS: Reduce in goodness of fit compared to bivariate models.

# Accuracy
concordance(coxGen, newdata=datCredit_valid_PWPST) # Concordance= 0.7338 se= 0.003416
### RESULTS: Good Concordance

# House keeping
rm(coxGen)

### CONCLUSION: Leave all variables in the model.

#============================================================================================

modelVar <- c(modelVar, "InstalmentToBalance_Aggr_Prop", "Undrawn_Amt", "LN_TPE",
              "AgeToTerm")

#============================================================================================
# ------ 5. Portfolio Level
varlist <- data.table(vars=c("CuringEvents_Aggr_Prop","NewLoans_Aggr_Prop"),
                      vartypes=c("dec", "dec"))

#============================================================================================
# ------ 5.1 Are the values correlated?

# - Correlation analysis
corrAnalysis(datCredit_train_PWPST, varlist[vartypes!="cat"]$vars, corrThresh = 0.6, method = 'spearman') # Obtain correlation groups

### RESULTS: The variables are not correlated significantly.

# ------ 5.2 What is the performance of the variables in bivariate models?

# Initialize variable
vars <- c("CuringEvents_Aggr_Prop","NewLoans_Aggr_Prop")

# Goodness of fit test of variables
csTable(datCredit_train_PWPST, vars, seedVal=NA)
#                 Variables Iteration_1 Iteration_2 Iteration_3      Average
# 1:     NewLoans_Aggr_Prop      0.6796      0.6811      0.6796 0.6801000000
# 2: CuringEvents_Aggr_Prop      0.6802      0.6802      0.6788 0.6797333333
# 3:                  Range      0.0006      0.0009      0.0008 0.0007666667

### RESULTS:  Both seem to have a good goodness of fit.

# Obtain the concordance for different variables.
concTable(datCredit_train_PWPST, datCredit_valid_PWPST, vars)
#                 Variable Concordance          SD LR_Statistic
# 1:     NewLoans_Aggr_Prop   0.6445252 0.004497458         4742
# 2: CuringEvents_Aggr_Prop   0.6274140 0.004543555         4550

## RESULTS: The concordances are good

### CONCOLUSION:  Keep both variables in the model

# ------ 5.2 What is the performance of the variables in a single model?

# Goodness of fit
coxPort <- coxph(as.formula(paste0("Surv(Start,End,Default_Ind) ~ PerfSpell_Num + ",
                                        paste(vars,collapse=" + "))), id=PerfSpell_Key,
                      datCredit_train_PWPST)

summary(coxPort)
### RESULS: [CuringEvents_Aggr_Prop] has a high se(coef) ((18)realative to coef (-51).

### CONCLUSION: Remove [CuringEvents_Aggr_Prop]

#============================================================================================

modelVar <- c(modelVar,"NewLoans_Aggr_Prop")

#============================================================================================
# ------ 6. Macro Economic
varlist <- data.table(vars=c("M_DTI_Growth","M_Emp_Growth","M_Inflation_Growth",
                             "M_RealGDP_Growth","M_RealIncome_Growth"),
                      vartypes=c("prc", "prc", "prc", "prc", "prc"))
#============================================================================================

# ------ 6.1 Which economic variables are highly correlated with one another?

# - Correlation analysis
corrAnalysis(datCredit_train_PWPST, varlist$vars[varlist$vartypes != 'cat'],
             corrThresh = 0.6, method = 'spearman') # Obtain correlation groups

### RESULTS: [M_Emp_Growth], [M_RealGDP_Growth] and [M_RealIncome_Growth] are highly correlated



# ------ 6.2 Which correlated variable should remain in the model?

# Initialize variables
vars <- c("M_Emp_Growth", "M_RealGDP_Growth", "M_RealIncome_Growth")

# Obtain the B-statistics for goodness of fit tests.
csTable(datCredit_train_PWPST, vars, seedVal=NA)
#               Variables Iteration_1 Iteration_2 Iteration_3     Average
# 1: M_RealIncome_Growth      0.6800      0.6812      0.6801 0.680433333
# 2:    M_RealGDP_Growth      0.6780      0.6794      0.6816 0.679666667
# 3:        M_Emp_Growth      0.6807      0.6784      0.6794 0.679500000
# 4:               Range      0.0027      0.0028      0.0022 0.002566667

### RESULTS:  Non-conclusive results as rankings fluctuated too much over 3 iterations.

# Obtain the concordance for different variables.
concTable(datCredit_train_PWPST, datCredit_valid_PWPST, vars)
#               Variable Concordance          SD LR_Statistic
# 1:    M_RealGDP_Growth   0.6467414 0.004210124         4372
# 2: M_RealIncome_Growth   0.6149946 0.004579253         4372
# 3:        M_Emp_Growth   0.6071295 0.004643992         4358

### RESULTS:  [M_RealIncome_Growth] has the highest concordance by a significant
###           margin.

### CONCLUSION: Keep M_RealIncome_Growth in the model.

varlist = vecChange(varlist, Remove=c("M_Emp_Growth", "M_RealGDP_Growth"))



# ------ 6.3 What are the predictive powers of the current variables in the models?

# Initialize variables
vars <- c("M_DTI_Growth", "M_Inflation_Growth", "M_RealIncome_Growth")

# Compare goodness of fit of different variables
csTable(datCredit_train_PWPST,vars,seedVal = NA)
#             Variables Iteration_1 Iteration_2 Iteration_3   Average
# 1:  M_Inflation_Growth      0.6794      0.6800      0.6829 0.6807667
# 2: M_RealIncome_Growth      0.6790      0.6831      0.6799 0.6806667
# 3:        M_DTI_Growth      0.6802      0.6794      0.6794 0.6796667
# 4:               Range      0.0012      0.0037      0.0035 0.0028000
# 
### RESULTS:  [M_Inflation_Growth] seems to have the best fit.

# Compare concordance of different variables
concTable(datCredit_train_PWPST, datCredit_valid_PWPST,vars)
#               Variable Concordance          SD LR_Statistic
# 1:        M_DTI_Growth   0.6797836 0.003779943         5451
# 2:  M_Inflation_Growth   0.6719732 0.004056093         4933
# 3: M_RealIncome_Growth   0.6149946 0.004579253         4372

### RESULTS: [M_DTI_Growth] has the highest concordance.
### CONCLUSION: The predictive powers may be enhanced with lags.


# ------ 6.4 Which lags for current variables be included in the models?

# ------ 6.4.1 Which lags for [M_DTI_Growth] be included in the models?
# Initialize variables
vars <- c("M_DTI_Growth_1","M_DTI_Growth_12","M_DTI_Growth_3", "M_DTI_Growth_6",
          "M_DTI_Growth_9","M_DTI_Growth")

# Obtain the B-statistics for goodness of fit tests.
csTable(datCredit_train_PWPST, vars, seedVal=NA)
#         Variables Iteration_1 Iteration_2 Iteration_3     Average
# 1:  M_DTI_Growth_1      0.6798      0.6808      0.6804 0.680333333
# 2:  M_DTI_Growth_3      0.6838      0.6774      0.6794 0.680200000
# 3: M_DTI_Growth_12      0.6813      0.6811      0.6776 0.680000000
# 4:  M_DTI_Growth_6      0.6834      0.6815      0.6750 0.679966667
# 5:    M_DTI_Growth      0.6792      0.6800      0.6788 0.679333333
# 6:  M_DTI_Growth_9      0.6791      0.6805      0.6774 0.679000000
# 7:           Range      0.0047      0.0041      0.0054 0.004733333

### RESULTS: Only [M_DTI_Growth_1], [M_DTI_Growth_3] and [M_DTI_Growth_12] appears to have the best fit.

# Obtain the concordance for different variables.
concTable(datCredit_train_PWPST, datCredit_valid_PWPST, vars)
#           Variable Concordance          SD LR_Statistic
# 1:  M_DTI_Growth_6   0.6818990 0.003900483         5406
# 2:  M_DTI_Growth_9   0.6818670 0.003950138         5485
# 3:  M_DTI_Growth_3   0.6815223 0.003852030         5344
# 4: M_DTI_Growth_12   0.6812717 0.003983199         5449
# 5:  M_DTI_Growth_1   0.6806493 0.003803199         5410
# 6:    M_DTI_Growth   0.6797836 0.003779943         5451

### RESULTS: Only [M_DTI_Growth_6], [M_DTI_Growth_9] and [M_DTI_Growth_3] appears to have the best fit.

### CONCLUSION: Include the following lags in the model: [M_DTI_Growth_3]

varlist <- vecChange(varlist,Add=data.table(vars=c("M_DTI_Growth_3"),
                                            vartypes=c("prc")))

# ------ 6.4.2 Should lags for [M_Inflation_Growth] be included in the models?
# Initialize variables
vars <- c("M_Inflation_Growth_1","M_Inflation_Growth_12","M_Inflation_Growth_3",
          "M_Inflation_Growth_6","M_Inflation_Growth_9","M_Inflation_Growth")

# Obtain the B-statistics for goodness of fit tests.
csTable(datCredit_train_PWPST, vars, seedVal=NA)
#                 Variables Iteration_1 Iteration_2 Iteration_3     Average
# 1:  M_Inflation_Growth_3      0.6794      0.6820      0.6832 0.681533333
# 2: M_Inflation_Growth_12      0.6836      0.6790      0.6772 0.679933333
# 3:  M_Inflation_Growth_6      0.6800      0.6817      0.6781 0.679933333
# 4:  M_Inflation_Growth_9      0.6794      0.6786      0.6802 0.679400000
# 5:  M_Inflation_Growth_1      0.6754      0.6815      0.6805 0.679133333
# 6:    M_Inflation_Growth      0.6797      0.6791      0.6774 0.678733333
# 7:                 Range      0.0082      0.0034      0.0060 0.005866667

### RESULTS: [M_Inflation_Growth_3], [M_Inflation_Growth_12] and [M_Inflation_Growth_6] appear to have the best fits.

# Obtain the concordance for different variables.
concTable(datCredit_train_PWPST, datCredit_valid_PWPST, vars)
#                 Variable Concordance          SD LR_Statistic
# 1:  M_Inflation_Growth_6   0.6847471 0.003911142         5363
# 2:  M_Inflation_Growth_3   0.6828872 0.003918695         5256
# 3:  M_Inflation_Growth_9   0.6752656 0.004030759         5345
# 4:  M_Inflation_Growth_1   0.6745120 0.004014520         5030
# 5:    M_Inflation_Growth   0.6719732 0.004056093         4933
# 6: M_Inflation_Growth_12   0.6679151 0.004118109         5144

### RESULTS: [M_Inflation_Growth_3], [M_Inflation_Growth_6] and [M_Inflation_Growth_9] appear to be the most accurate.

### CONCLUSION: Add [M_Inflation_Growth_3] and [M_Inflation_Growth_6] to the model since it seems
###             to have the first best goodness of fit and second best accuracy

varlist <- vecChange(varlist,Add=data.table(vars=c("M_Inflation_Growth_3","M_Inflation_Growth_6"),
                                            vartypes=c("prc","prc")))

 # ------ 6.4.3 Should lags for [M_RealIncome_Growth] be included in the models?
vars <- c("M_RealIncome_Growth_1","M_RealIncome_Growth_12","M_RealIncome_Growth_3",
          "M_RealIncome_Growth_6","M_RealIncome_Growth_9","M_RealIncome_Growth")

# Obtain the B-statistics for goodness of fit tests.
csTable(datCredit_train_PWPST, vars, seedVal=NA)
# Variables Iteration_1 Iteration_2 Iteration_3   Average
# 1:    M_RealIncome_Growth      0.6791      0.6820      0.6834 0.6815000
# 2:  M_RealIncome_Growth_6      0.6818      0.6805      0.6816 0.6813000
# 3:  M_RealIncome_Growth_1      0.6822      0.6790      0.6806 0.6806000
# 4:  M_RealIncome_Growth_9      0.6829      0.6776      0.6811 0.6805333
# 5: M_RealIncome_Growth_12      0.6804      0.6783      0.6811 0.6799333
# 6:  M_RealIncome_Growth_3      0.6793      0.6777      0.6787 0.6785667
# 7:                  Range      0.0038      0.0044      0.0047 0.0043000

### RESULTS: [M_RealIncome_Growth_6], [M_RealIncome_Growth] and [M_RealIncome_Growth_1] appear to have the best fits.

# Obtain the concordance for different variables.
concTable(datCredit_train_PWPST, datCredit_valid_PWPST, vars)
#                 Variable Concordance          SD LR_Statistic
# 1: M_RealIncome_Growth_12   0.6758466 0.003950397         4669
# 2:  M_RealIncome_Growth_9   0.6628651 0.004058959         4481
# 3:  M_RealIncome_Growth_6   0.6513517 0.004159279         4393
# 4:  M_RealIncome_Growth_3   0.6413519 0.004227750         4355
# 5:    M_RealIncome_Growth   0.6149946 0.004579253         4372
# 6:  M_RealIncome_Growth_1   0.6118485 0.004585836         4359

### RESULTS: [M_RealIncome_Growth_12], [M_RealIncome_Growth_9] and [M_RealIncome_Growth_6] appear to be the most accurate,
###           although all of the have extremely weak concordances.

### CONCLUSION: Include lag [M_RealIncome_Growth_6].

varlist <- vecChange(varlist,Add=data.table(vars=c("M_Inflation_Growth_9","M_Inflation_Growth_6"),
                                            vartypes=c("prc","prc")))

# ------ 6.5 What is the predictive performance of the current univariate thematic models?

# Initialize macro-economic thematic variables
vars <- c("M_DTI_Growth", "M_Inflation_Growth","M_DTI_Growth_3", "M_Inflation_Growth_9",
          "M_Inflation_Growth_6", "M_Inflation_Growth_3", "M_RealIncome_Growth")

# Test Goodness of fit
csTable(datCredit_train_PWPST,vars)
#               Variable B_Statistic
# 5 M_Inflation_Growth_6      0.6831
# 6 M_Inflation_Growth_3      0.6822
# 4 M_Inflation_Growth_9      0.6812
# 1         M_DTI_Growth      0.6809
# 7  M_RealIncome_Growth      0.6792
# 2   M_Inflation_Growth      0.6789
# 3       M_DTI_Growth_3      0.6779

### RESULTS: All variables seem to be relatively good fits with the lowest being 0.6449.

# Test Accuracy
concTable(datCredit_train_PWPST, datCredit_valid_PWPST, vars)
#                 Variable Concordance          SD LR_Statistic
# 1: M_Inflation_Growth_6   0.6847471 0.003911142         5363
# 2: M_Inflation_Growth_3   0.6828872 0.003918695         5256
# 3:       M_DTI_Growth_3   0.6815223 0.003852030         5344
# 4:         M_DTI_Growth   0.6797836 0.003779943         5451
# 5: M_Inflation_Growth_9   0.6752656 0.004030759         5345
# 6:   M_Inflation_Growth   0.6719732 0.004056093         4933
# 7:  M_RealIncome_Growth   0.6149946 0.004579253         4372

### RESULTS: All variables should be kept in the model. Although [M_RealIncome_Growth]
###           has a low concordance, it is kept in the model as a proxy for M_Emp_Growth.

# ------ 6.6 What is the predictive performance of the current thematic model?
# Initialize macro-economic thematic variables
vars <- c("M_DTI_Growth", "M_Inflation_Growth","M_DTI_Growth_3", "M_Inflation_Growth_9",
          "M_Inflation_Growth_6", "M_Inflation_Growth_3", "M_RealIncome_Growth")

# Build model based on macro economic variables
coxMacro <- coxph(as.formula(paste0("Surv(Start,End,Default_Ind) ~ PerfSpell_Num + ",
                                    paste(vars,collapse=" + "))), id=PerfSpell_Num,
                  datCredit_train_PWPST)
summary(coxMacro)

 ### RESULTS: [M_Inflation_Growth_6] and [M_DTI_Growth_3] have high p-values, with [M_DTI_Growth_3] 
###           the lowest Goodness of Fit and therefore should be removed.

vars <- c("M_DTI_Growth", "M_Inflation_Growth","M_DTI_Growth_3", "M_Inflation_Growth_9",
          "M_RealIncome_Growth")

# Build model based on macro economic variables
coxMacro <- coxph(as.formula(paste0("Surv(Start,End,Default_Ind) ~ ",
                                    paste(vars,collapse=" + "))), id=LoanID,
                  datCredit_train_PWPST)
summary(coxMacro)

### RESULTS:  All values have significant p-value.

# Test Goodness of fit of the single model based on all variables.
GoF_CoxSnell_KS(coxMacro,datCredit_train_PWPST,GraphInd = F) # 0.6826
# Goodness of fit is equal to greatest Goodness of fit bivariate model.

# Test accuracy of the single model based on all variables.
concordance(coxMacro, newdata=datCredit_valid_PWPST)# Concordance= 0.5686 se= 0.004245
# Concordance is much worse than that of bivariate models.

#============================================================================================

modelVar <- c(modelVar,"M_DTI_Growth", "M_Inflation_Growth", "M_Inflation_Growth_6",
              "M_RealIncome_Growth")

#============================================================================================
#============================================================================================
# ------ 7. Multi-spell variables
#varlist <- data.table(vars=as.character(),vartypes=as.character())

#============================================================================================

# # ------ 7.1 Should we include a variable containing the counter length of a previous performance spell?
# 
# # Create a variable for the duration of the previous spell
# # Step 1: Calculate the maximum value of PerfSpell_Counter for each PerfSpell_Key group
# datCredit_train_PWPST[, PerfSpell_Max_Counter := max(PerfSpell_Counter), by = PerfSpell_Key]
# datCredit_valid_PWPST[, PerfSpell_Max_Counter := max(PerfSpell_Counter), by = PerfSpell_Key]
# 
# # Step 2: For rows where PerfSpell_Counter is 1, calculate a lagged (shifted) version of max_Counter
# # Grouping is done by LoanID to ensure the shift occurs within each LoanID group
# datCredit_train_PWPST[PerfSpell_Counter == 1, PerfSpell_Prev_Length := shift(PerfSpell_Max_Counter), by = LoanID]
# datCredit_valid_PWPST[PerfSpell_Counter == 1, PerfSpell_Prev_Length := shift(PerfSpell_Max_Counter), by = LoanID]
# 
# # Step 3: Fill down (forward fill) missing values in shift_Counter within each LoanID group
# # This ensures that NA values are replaced by the last non-NA value in the same LoanID group
# datCredit_train_PWPST[, PerfSpell_Prev_Length := nafill(PerfSpell_Prev_Length, type = "locf"), by = LoanID]
# datCredit_valid_PWPST[PerfSpell_Counter == 1, PerfSpell_Prev_Length := shift(PerfSpell_Max_Counter), by = LoanID]
# 
# # Step 4: Set shift_Counter to 0 for rows where PerfSpell_Num equals 1
# # This is likely to handle a specific condition for the beginning of a PerfSpell_Num sequence
# datCredit_train_PWPST[PerfSpell_Num == 1, PerfSpell_Prev_Length := 0]
# datCredit_valid_PWPST[PerfSpell_Counter == 1, PerfSpell_Prev_Length := shift(PerfSpell_Max_Counter), by = LoanID]
# 
# datCredit_train_PWPST[is.na(datCredit_train_PWPST$PerfSpell_Prev_Length) &
#                         order(PerfSpell_Key), list(Date,PerfSpell_Key, Removed,
#                                                    Default_Ind, 
#                                                    PerfSpell_Prev_Length)]
# 
# ### NOTE: Some rows still contain NA values.

#============================================================================================

# ------ 7. Final Model
# ------ 7.1 How does the univariate models compare with one another?
# Initialize variables
vars <- c("g0_Delinq_SD_4", "Arrears", "g0_Delinq_Ave", "TimeInDelinqState_Lag_1",      
          "slc_acct_arr_dir_3_Change_Ind", "slc_acct_pre_lim_perc_imputed_med", 
          "pmnt_method_grp", "InterestRate_Nom", "M_Repo_Rate", "M_Repo_Rate_9", 
          "InstalmentToBalance_Aggr_Prop", "Undrawn_Amt", "LN_TPE", "AgeToTerm",
          "NewLoans_Aggr_Prop", "M_DTI_Growth","M_Inflation_Growth",
          "M_Inflation_Growth_6", "M_RealIncome_Growth")

# Test Goodness of fit
csTable(datCredit_train_PWPST,vars)
#                             Variable B_Statistic
# 11     InstalmentToBalance_Aggr_Prop      0.6838
# 19               M_RealIncome_Growth      0.6831
# 9                        M_Repo_Rate      0.6824
# 6  slc_acct_pre_lim_perc_imputed_med      0.6814
# 17                M_Inflation_Growth      0.6814
# 10                     M_Repo_Rate_9      0.6809
# 13                            LN_TPE      0.6802
# 16                      M_DTI_Growth      0.6800
# 8                   InterestRate_Nom      0.6798
# 18              M_Inflation_Growth_6      0.6792
# 14                         AgeToTerm      0.6788
# 12                       Undrawn_Amt      0.6786
# 3                      g0_Delinq_Ave      0.6783
# 2                            Arrears      0.6782
# 15                NewLoans_Aggr_Prop      0.6779
# 5      slc_acct_arr_dir_3_Change_Ind      0.6741
# 7                    pmnt_method_grp      0.6730
# 4            TimeInDelinqState_Lag_1      0.6643
# 1                     g0_Delinq_SD_4      0.6458

### RESULTS: All variables seem to have comparible Goodness of fits, although
###           [TimeInDelinqState_Lag_1] and [g0_Delinq_SD_4] have noticable bad ones.

# Test accuracy
# concTable(datCredit_train_PWPST, datCredit_valid_PWPST, vars)
#                 Variable Concordance          SD LR_Statistic
# 1:                    g0_Delinq_SD_4   0.9806753 0.001401249       107861
# 2:                           Arrears   0.9450365 0.001432871        20316
# 3:           TimeInDelinqState_Lag_1   0.9219226 0.002705651        43522
# 4:     slc_acct_arr_dir_3_Change_Ind   0.8891861 0.001995888        43030
# 5:                   pmnt_method_grp   0.7546841 0.003472790        14162
# 6:                       Undrawn_Amt   0.7539662 0.002225233         6580
# 7: slc_acct_pre_lim_perc_imputed_med   0.7241381 0.002584756        10364
# 8:                  InterestRate_Nom   0.7009491 0.003957952         5460
# 9:                     g0_Delinq_Ave   0.6863403 0.003767392         5347
# 10:              M_Inflation_Growth_6   0.6847471 0.003911142         5363
# 11:                     M_Repo_Rate_9   0.6845921 0.003878557         5325
# 12:                      M_DTI_Growth   0.6797836 0.003779943         5451
# 13:                       M_Repo_Rate   0.6785476 0.004009584         4892
# 14:                M_Inflation_Growth   0.6719732 0.004056093         4933
# 15:     InstalmentToBalance_Aggr_Prop   0.6708871 0.004148979         4637
# 16:                NewLoans_Aggr_Prop   0.6445252 0.004497458         4742
# 17:                         AgeToTerm   0.6437042 0.004368801         4728
# 18:                            LN_TPE   0.6325628 0.003571573         4523
# 19:               M_RealIncome_Growth   0.6149946 0.004579253         4372

### RESULTS: All variables seem to have comparible concordances.

# ------ 7.2 What is the predictive performance of a single model?
# Initialize variables
vars <- c("g0_Delinq_SD_4", "Arrears", "g0_Delinq_Ave", "TimeInDelinqState_Lag_1",      
          "slc_acct_arr_dir_3_Change_Ind", "slc_acct_pre_lim_perc_imputed_med", 
          "pmnt_method_grp", "InterestRate_Nom", "M_Repo_Rate", "M_Repo_Rate_9", 
          "InstalmentToBalance_Aggr_Prop", "Undrawn_Amt", "LN_TPE", "AgeToTerm",
          "NewLoans_Aggr_Prop", "M_DTI_Growth","M_Inflation_Growth",
          "M_Inflation_Growth_6", "M_RealIncome_Growth")

vars2 <- c("g0_Delinq_SD_4", "Arrears", "g0_Delinq_Ave", "TimeInDelinqState_Lag_1",      
           "slc_acct_arr_dir_3_Change_Ind", "slc_acct_pre_lim_perc_imputed_med", 
           "InterestRate_Nom", "M_Repo_Rate", "M_Repo_Rate_9", 
           "InstalmentToBalance_Aggr_Prop", "LN_TPE", "AgeToTerm",
           "NewLoans_Aggr_Prop", "M_DTI_Growth","M_Inflation_Growth",
           "M_Inflation_Growth_6", "M_RealIncome_Growth", "PerfSpell_Grp")

# Test Goodness of fit
csTable_PWPST <- csTable(datCredit_train_PWPST,vars2,TimeDef="PWP_ST")
#                             Variable B_Statistic
# 11                            LN_TPE      0.6854
# 6  slc_acct_pre_lim_perc_imputed_med      0.6835
# 9                      M_Repo_Rate_9      0.6834
# 10     InstalmentToBalance_Aggr_Prop      0.6833
# 17               M_RealIncome_Growth      0.6833
# 8                        M_Repo_Rate      0.6832
# 16              M_Inflation_Growth_6      0.6825
# 12                         AgeToTerm      0.6822
# 7                   InterestRate_Nom      0.6817
# 13                NewLoans_Aggr_Prop      0.6817
# 14                      M_DTI_Growth      0.6816
# 2                            Arrears      0.6815
# 15                M_Inflation_Growth      0.6807
# 3                      g0_Delinq_Ave      0.6800
# 5      slc_acct_arr_dir_3_Change_Ind      0.6797
# 4            TimeInDelinqState_Lag_1      0.6713
# 1                     g0_Delinq_SD_4      0.6468

# Test accuracy
concTable_PWPST <- concTable(datCredit_train_PWPST, datCredit_valid_PWPST, vars2, TimeDef="PWP_ST")
#                             Variable Concordance           SD LR_Statistic
# 1:                    g0_Delinq_SD_4   0.9929890 0.0004440688       109216
# 2:                           Arrears   0.9301293 0.0015083583         9288
# 3:           TimeInDelinqState_Lag_1   0.9163394 0.0026692102        39902
# 4:     slc_acct_arr_dir_3_Change_Ind   0.9045385 0.0015469709        39628
# 5: slc_acct_pre_lim_perc_imputed_med   0.7369588 0.0025027221        10239
# 6:                  InterestRate_Nom   0.7213484 0.0039525630         5506
# 7:                     g0_Delinq_Ave   0.7151798 0.0038156470         5652
# 8:                      M_DTI_Growth   0.7128062 0.0038242592         5623
# 9:                     M_Repo_Rate_9   0.7114630 0.0039792989         5390
# 10:              M_Inflation_Growth_6   0.7065827 0.0040608535         5300
# 11:                       M_Repo_Rate   0.7036419 0.0041361893         5025
# 12:                M_Inflation_Growth   0.7003644 0.0041659212         4906
# 13:     InstalmentToBalance_Aggr_Prop   0.6923587 0.0043362701         4612
# 14:                NewLoans_Aggr_Prop   0.6600105 0.0047112105         4360
# 15:               M_RealIncome_Growth   0.6458617 0.0045948006         4078
# 16:                            LN_TPE   0.6450998 0.0036884899         4212
# 17:                         AgeToTerm   0.6307607 0.0038752408         4163

# Create table object
Table_PWPST <- left_join(csTable_PWPST$Table,concTable_PWPST,by="Variable")

# 
# Build model based on variables
cox_PWPST <- coxph(as.formula(paste0("Surv(Start,End,Default_Ind) ~ ",paste(vars2,collapse=" + "))), 
                   id=PerfSpell_Key, datCredit_train_PWPST, ties="efron")
# NOTE: Default option for handling ties is recently "Efron's method", but we specify this for backwards compatability
c <- coefficients(cox_PWPST)
c <- data.table(Variable=names(c),Coefficient=c)
#                         Variable         Coefficient
# 1:                     PerfSpell_Num    0.15096051311539
# 2:                    g0_Delinq_SD_4    6.75209715853658
# 3:                           Arrears    0.00002264476443
# 4:                     g0_Delinq_Ave   -1.19661791411541
# 5:           TimeInDelinqState_Lag_1   -0.03252747945563
# 6:     slc_acct_arr_dir_3_Change_Ind    1.21831026120410
# 7: slc_acct_pre_lim_perc_imputed_med  -17.22076821317840
# 8:       pmnt_method_grpMISSING_DATA    0.03660829390555
# 9:    pmnt_method_grpSalary/Suspense    0.77678877192294
# 10:          pmnt_method_grpStatement    0.39000075830291
# 11:                  InterestRate_Nom    7.33906786444637
# 12:                       M_Repo_Rate   -1.14215556492784
# 13:                     M_Repo_Rate_9   -4.75512297497422
# 14:     InstalmentToBalance_Aggr_Prop -129.73549435563544
# 15:                       Undrawn_Amt   -0.00000003076138
# 16:                         LN_TPEWHL   -0.20291518212577
# 17:                         AgeToTerm   -0.06736246885203
# 18:                NewLoans_Aggr_Prop    2.79644547316863
# 19:                      M_DTI_Growth    2.10615496221582
# 20:                M_Inflation_Growth    1.86190351858316
# 21:              M_Inflation_Growth_6    3.12806371982533
# 22:               M_RealIncome_Growth    0.27999127862750
# Variable         Coefficient

# Goodnes of fit
GoF_CoxSnell_KS(cox_PWPST,datCredit_train_PWPST, GraphInd=TRUE, legPos=c(0.6,0.4)) # 0.6408

### RESULTS: Goodness of fit for the model seems to be a bit low.

# Accuracy

# Build model based on variables
concordance(cox_PWPST, newdata=datCredit_valid_PWPST)
# Concordance= 0.9942 se= 0.0005169

tROC.multi(datCredit_valid_PWPST, cox_PWPST, month_Start=0, month_End=12, sLambda=0.05,
           estMethod="NN-0/1", numDigits=2,fld_ID="PerfSpell_Key", fld_Event="Default_Ind",
           eventVal=1, fld_StartTime="Start", fld_EndTime="End",Graph=TRUE,
           graphName="timedROC-Graph_PWPST",
           genFigPath="C:/Users/R8873885/OneDrive - FRG/Documents/LifetimePD-TermStructure-RecurrentEvents/Figures/PWPST/tdROC/")


### RESULTS: Accuracy for the model seems to be a bit high.

pack.ffdf(paste0(genObjPath,"PWPST_Univariate_Models"), Table_PWPST)
pack.ffdf(paste0(genPath,"PWPST_Cox_Model"), cox_PWPST)

# - Cleanup
rm(datCredit_train_PWPST, datCredit_valid_PWPST, cox_PWPST); gc()
