# ============================================ INPUT SPACE: AG =========================================
# Divide data into thematic groups and perform data analysis on them to compile an input space for the AG model
# ------------------------------------------------------------------------------------------------------
# PROJECT TITLE: Default survival modelling
# SCRIPT AUTHOR(S): Bernard Scheepers (BS)
# ------------------------------------------------------------------------------------------------------
# -- Script dependencies:
#   - 0.Setup.R
#   - 1.Data_Import.R
#   - 2a.Data_Prepare_Credit_Basic.R
#   - 2b.Data_Prepare_Credit_Advanced.R
#   - 2c.Data_Prepare_Credit_Advanced2.R
#   - 2d.Data_Enrich.R
#   - 2e.Data_Prepare_Macro.R
#   - 2f.Data_Fusion1.R
#   - 3a(i).Data_Transform.R
#   - 3c(ii).Data_Fusion2_AG.R

# -- Inputs:
#   - datCredit_train_AG | Prepared from script 3c
#   - datCredit_valid_AG | Prepared from script 3c

#
# -- Outputs:
#   - Input_Space
# ------------------------------------------------------------------------------------------------------

# ------ 1. Preliminaries
# - Confirm prepared datasets are loaded into memory
if (!exists('datCredit_train_AG')) unpack.ffdf(paste0(genPath,"creditdata_train_AG"), tempPath);gc()
if (!exists('datCredit_valid_AG')) unpack.ffdf(paste0(genPath,"creditdata_valid_AG"), tempPath);gc()


#============================================================================================
# ------ 1. Delinquency measures
varlist <- data.table(vars=c("PerfSpell_g0_Delinq_Num","Arrears" ,"TimeInDelinqState_Lag_1",
                             "g0_Delinq_Any_Aggr_Prop","g0_Delinq_Ave","slc_acct_arr_dir_3_Change_Ind"),
                      vartypes=c("acc", "acc", "acc", "dte", "dte", "bin"))

#==========================================================================================



# ------ 1.1 Which time window length is the best in calculating Delinquency volatility?

# Initialize variables to be tested
vars <- c("g0_Delinq_SD_4", "g0_Delinq_SD_5", "g0_Delinq_SD_6", "g0_Delinq_SD_9", "g0_Delinq_SD_12")

# Goodness of fit test
csTable(datCredit_train_AG, vars, TimeDef="AG") # [g0_Delinq_SD_12] has the lowest KS-statistic

aicTable(datCredit_train_AG, vars, TimeDef="AG", genPath=genObjPath) # [g0_Delinq_SD_5] has the lowest AIC

# Accuracy test
concTable(datCredit_train_AG, datCredit_valid_AG, vars, TimeDef="AG", genPath=genObjPath) # [g0_Delinq_SD_4] has the highest concordance

### Conclusion: Different results for each test, although [go_Delinq_SD_4] holistically ranked the best.

varlist <- vecChange(varlist,Remove="PerfSpell_g0_Delinq_SD",Add=data.table(vars=c("g0_Delinq_SD_4"), vartypes=c("acc")))

# ------ 1.2 Which variables are highly correlated?

# Correlation analysis
corrAnalysis(datCredit_train_AG, varlist[vartypes!="cat"]$vars, corrThresh = 0.6, method = 'spearman') # Obtain correlation groups

### RESULTS:  [g0_Delinq_Any_Aggr_Prop] and [g0_Delinq_Ave] with a correlation of 100%

### CONCLUSION: A single variable from each group must be retained while the rest are removed.



# ------ 1.2.1 Which variable should be kept from group 1) [g0_Delinq_Any_Aggr_Prop] and [g0_Delinq_Ave]

vars <- c("g0_Delinq_Any_Aggr_Prop", "g0_Delinq_Ave")

# Goodness of fit
csTable(datCredit_train_AG, vars, TimeDef="AG") # No noticeable difference

aicTable(datCredit_train_AG, vars, TimeDef="AG", genPath=genObjPath) # [g0_Delinq_Ave] has the lowest AIC

# Accuracy
concTable(datCredit_train_AG, datCredit_valid_AG, vars, TimeDef="AG", genPath=NA) # [g0_Delinq_Ave] has the highest concordance

### CONCLUSION: [g0_Delinq_Ave] seems to outperform [g0_Delinq_Any_Aggr_Prop] and therefore is kept in the model.

varlist <- vecChange(varlist,Remove=c("g0_Delinq_Any_Aggr_Prop"))

# ------ 1.3 How does [PerfSpell_g0_Delinq_Num] compare to [g0_Delinq_SD_4]?

# Initialize variables
vars <- c("PerfSpell_g0_Delinq_Num", "g0_Delinq_SD_4")

# Goodness of fit
csTable(datCredit_train_AG, vars, TimeDef="AG") # [PerfSpell_g0_Delinq_Num] has the lowest KS-statistic

aicTable(datCredit_train_AG, vars, TimeDef="AG", genPath=genObjPath) # [g0_Delinq_SD_4] has the lowest AIC

# Accuracy
concTable(datCredit_train_AG, datCredit_valid_AG, vars, TimeDef="AG", genPath=genObjPath) # [g0_Delinq_SD_4] has the highest concordance

### CONCLUSION: Keep [g0_Delinq_SD_4] in the model and remove [PerfSpell_g0_Delinq_Num]

varlist <- vecChange(varlist,Remove="PerfSpell_g0_Delinq_Num")



# ------ 1.4 What is the performance of current thematic variables in univariate models?

vars <- c("g0_Delinq_SD_4","g0_Delinq_Ave", "Arrears",
          "TimeInDelinqState_Lag_1","slc_acct_arr_dir_3_Change_Ind")

# Goodness of fit test
csTable(datCredit_train_AG, vars, TimeDef="AG")
#                   Variable      D
# 2                 g0_Delinq_Ave 0.6866
# 3                       Arrears 0.6862
# 5 slc_acct_arr_dir_3_Change_Ind 0.6805
# 4       TimeInDelinqState_Lag_1 0.6717
# 1                g0_Delinq_SD_4 0.6488

aicTable(datCredit_train_AG, vars, TimeDef="AG", genPath=genObjPath)
#                     Variable      AIC      pValue
# 1:                g0_Delinq_SD_4 140797.9 0.011818182
# 2:       TimeInDelinqState_Lag_1 208701.9 0.011080982
# 3: slc_acct_arr_dir_3_Change_Ind 211618.7 0.009816939
# 4:                       Arrears 240024.5 0.008952998
# 5:                 g0_Delinq_Ave 247568.4 0.008735097

# Accuracy test
concTable(datCredit_train_AG, datCredit_valid_AG, vars, TimeDef="AG", genPath=genObjPath)
#                         Variable Concordance           SD LR_Statistic
# 1:                g0_Delinq_SD_4   0.9940415 0.0004646086       113312
# 2:                       Arrears   0.9537163 0.0017187368        14086
# 3:       TimeInDelinqState_Lag_1   0.9432592 0.0016388105        45408
# 4: slc_acct_arr_dir_3_Change_Ind   0.8994939 0.0015165245        42492
# 5:                 g0_Delinq_Ave   0.6896797 0.0042802653         6542

### CONCLUSION: Leave all variables in the model.



# ------ 1.5 What is the performance of current thematic cox ph model?

# Build cox model based on all thematic variables
coxDelinq <- coxph(Surv(Start,End,Default_Ind) ~ g0_Delinq_Ave +
                     g0_Delinq_SD_4 + slc_acct_arr_dir_3_Change_Ind +
                     TimeInDelinqState_Lag_1 + Arrears, id=LoanID,
                   data=datCredit_train_AG)

summary(coxDelinq)

### RESULTS: All variables are significant (p-value < 0.001)

# Test goodness of fit
GoF_CoxSnell_KS(coxDelinq,datCredit_train_AG,GraphInd = FALSE) # 0.6456

AIC(coxDelinq) # 134698.2

# Test accuracy
concordance(coxDelinq, newdata=datCredit_valid_AG) # Concordance= 0.9965 se= 0.0001942

### CONCLUSION: Keep all variables in the model

# House keeping
rm(coxDelinq)

#============================================================================================

modelVar <- c("g0_Delinq_SD_4", "Arrears", "g0_Delinq_Ave",
              "TimeInDelinqState_Lag_1", "slc_acct_arr_dir_3_Change_Ind",
              "slc_acct_roll_ever_24_imputed_mean")

#============================================================================================





#============================================================================================
# ------ 2. Engineered measures
varlist <- data.table(vars=c("slc_acct_pre_lim_perc_imputed_med",
                             "slc_acct_prepaid_perc_dir_12_imputed_med",
                             "pmnt_method_grp") ,
                      vartypes=c("prc", "dec", "cat"))
#=========================================================================================

# ------ 2.1 Which variables should be removed due to high correlation?

# - Correlation analysis
corGroups <- corrAnalysis(datCredit_train_AG, varlist[vartypes!="cat"]$vars, corrThresh = 0.6, method = 'spearman') # Obtain correlation groups

### RESULTS:  1) [slc_acct_pre_lim_perc_imputed_med] and [slc_acct_prepaid_perc_dir_12_imputed_med] with a correlation of 0.83

# Initialize variables to be tested
vars <- c("slc_acct_pre_lim_perc_imputed_med", "slc_acct_prepaid_perc_dir_12_imputed_med")

# Goodness of fit test
csTable(datCredit_train_AG, vars, TimeDef="AG") # [slc_acct_prepaid_perc_dir_12_imputed_med] has a lower KS-statistic

aicTable(datCredit_train_AG, vars, TimeDef="AG", genPath=genObjPath) # [slc_acct_pre_lim_perc_imputed_med] has the lowest AIC

# Accuracy test
concTable(datCredit_train_AG, datCredit_valid_AG, vars, TimeDef="AG", genPath=genObjPath) # [slc_acct_pre_lim_perc_imputed_med] has the highest concordance

### RESULTS: [slc_acct_prepaid_perc_dir_12_imputed_med] should be removed from the model.

varlist <- vecChange(varlist,Remove="slc_acct_prepaid_perc_dir_12_imputed_med")



# ------ 2.2 What is the predictive power of the variables left in varlist?

# Initialize variables to be tested
vars <- c("pmnt_method_grp", "slc_acct_pre_lim_perc_imputed_med")

csTable(datCredit_train_AG, vars, TimeDef="AG") # [slc_acct_pre_lim_perc_imputed_med ] has the lowest KS-statistic

aicTable(datCredit_train_AG, vars, TimeDef="AG", genPath=genObjPath) # [pmnt_method_grp] has the lowest AIC

# Accuracy test
concTable(datCredit_train_AG, datCredit_valid_AG, vars, TimeDef="AG", genPath=genObjPath) # [pmnt_method_grp] has the highest concordance

### Conclusion: All values should be kept in the model



# ------ 2.3 How predictive is a single model based on the thematic variables?

# Initialize variables
vars <- c("pmnt_method_grp", "slc_acct_pre_lim_perc_imputed_med")

# Build Cox model based on each variable
coxEngineered <- coxph(Surv(Start, End, Default_Ind) ~ pmnt_method_grp +
                         slc_acct_pre_lim_perc_imputed_med, id=LoanID,
                       data=datCredit_train_AG)
summary(coxEngineered)

### RESULTS: All variables are significant.

# Goodness of fit
GoF_CoxSnell_KS(coxEngineered,datCredit_train_AG, GraphInd=FALSE) # 0.6791

AIC(coxEngineered) # 237024.8

# Accuracy
concordance(coxEngineered,newdata=datCredit_valid_AG) # Concordance= 0.8042 se= 0.002868

# House keeping
rm(coxEngineered); gc()

#============================================================================================

modelVar <- c(modelVar,"slc_acct_pre_lim_perc_imputed_med","pmnt_method_grp")

#============================================================================================



#===========================================================================================
# ------ 3. Interest Rate
varlist <- data.table(vars=c("InterestRate_Nom", "InterestRate_Margin"),
                      vartypes=c("prc","prc"))
#==========================================================================================



# ------ 3.1 How correlated are the two variables?

# - Correlation analysis
corrAnalysis(datCredit_train_AG, varlist[vartypes!="cat"]$vars, corrThresh = 0.6, method = 'spearman') # Obtain correlation groups

### RESULTS: Correlation of 0.42, therefore no significant correlations.

# Goodness of fit test
csTable(datCredit_train_AG, varlist$vars, TimeDef="AG") # Insignificant difference.

aicTable(datCredit_train_AG, varlist$vars, TimeDef="AG", genPath=genObjPath) # [InterestRate_Nom] has the lowest AIC

# Accuracy test
concTable(datCredit_train_AG, datCredit_valid_AG,  varlist$vars, TimeDef="AG", genPath=genObjPath) # [InterestRate_Nom] has the highest concordance

### CONCLUSION: Keep both variables in the model.



# ------ 3.2 How does [InterestRate_Margin] compare with [InterestRate_Margin_imputed_bin]?

# Initialize variables
vars <- c("InterestRate_Margin", "InterestRate_Margin_imputed_bin")

# Goodness of fit test
csTable(datCredit_train_AG, vars, TimeDef="AG") # No noticeable differences.

aicTable(datCredit_train_AG, vars, TimeDef="AG", genPath=genObjPath) # [InterestRate_Margin_imputed_bin] has the lowest AIC

# Accuracy test
concTable(datCredit_train_AG, datCredit_valid_AG, vars, TimeDef="AG", genPath=genObjPath) # [InterestRate_Margin] has the highest concordance

### CONCLUSIONS: Keep [InterestRate_Margin] in the model.


# ------ 3.3 How does [InterestRate_Margin], [InterestRate_Margin_Aggr_Med] and [M_Repo_Rate] compare with one another?

vars <- c("InterestRate_Margin", "M_Repo_Rate", "InterestRate_Margin_Aggr_Med")

# Goodness of fit test
csTable(datCredit_train_AG, vars, TimeDef="AG") # No noticeable differences

aicTable(datCredit_train_AG, vars, TimeDef="AG", genPath=genObjPath) # [InterestRate_Margin_Aggr_Med] has the lowest AIC

# Accuracy test
concTable(datCredit_train_AG, datCredit_valid_AG, vars, TimeDef="AG", genPath=genObjPath) # [InterestRate_Margin_Aggr_Med] has the highest concordance

### CONCLUSIONS: Keep [InterestRate_Margin_Aggr_Med] in the model.

varlist <- vecChange(varlist,Remove="InterestRate_Margin", Add=data.table(vars="InterestRate_Margin_Aggr_Med", vartypes="prc"))



# ------ 3.4 Should we add lagging variables to the model?

vars <- c("InterestRate_Margin_Aggr_Med","InterestRate_Margin_Aggr_Med_1",
          "InterestRate_Margin_Aggr_Med_3","InterestRate_Margin_Aggr_Med_2",
          "InterestRate_Margin_Aggr_Med_9")

# Goodness of fit test
csTable(datCredit_train_AG, vars, TimeDef="AG") # No noticeable difference

aicTable(datCredit_train_AG, vars, TimeDef="AG", genPath=genObjPath) # [InterestRate_Margin_Aggr_Med] has the lowest AIC

# Accuracy test
concTable(datCredit_train_AG, datCredit_valid_AG, vars, TimeDef="AG", genPath=genObjPath) # [InterestRate_Margin_Aggr_Med] has the highest concordance

### CONCLUSION: Do not add a lag to the model, since [InterestRate_Margin_Aggr_Med] outperforms all other lags.




# ------ 3.5 What is the performance of current thematic variables in univariate models?

vars <- c("InterestRate_Nom", "InterestRate_Margin_Aggr_Med")

# Goodness of fit test
csTable(datCredit_train_AG, vars, TimeDef="AG") # Not a noticeable difference


aicTable(datCredit_train_AG, vars, TimeDef="AG", genPath=genObjPath) # [InterestRate_Nom] has the lowest AIC
#                   Variable      AIC      pValue
# 1:             InterestRate_Nom 247811.6 0.008640697
# 2: InterestRate_Margin_Aggr_Med 248042.3 0.008899713

# Accuracy test
concTable(datCredit_train_AG, datCredit_valid_AG, vars, TimeDef="AG", genPath=genObjPath) # [InterestRate_Margin_Aggr_Med_3] has the highest concordance
#                         Variable Concordance          SD LR_Statistic
# 1:             InterestRate_Nom   0.6907904 0.004433932         6299
# 2: InterestRate_Margin_Aggr_Med   0.6907121 0.004409733         6068

### CONCLUSION: Keep all variables in the model.



# ------ 3.6 What is the performance of current thematic cox ph model?

# Build cox model based on all thematic variables
coxInterest <- coxph(Surv(Start,End,Default_Ind) ~ InterestRate_Nom +
                       InterestRate_Margin_Aggr_Med,
                     id=LoanID, data=datCredit_train_AG)
summary(coxInterest)

### RESULTS: All variables are significant.

# Test goodness of fit
GoF_CoxSnell_KS(coxInterest,datCredit_train_AG,GraphInd = FALSE) # 0.6867

AIC(coxInterest) # 252576.8

# Accuracy
concordance(coxInterest, newdata=datCredit_valid_AG)# Concordance= 0.6455 se= 0.004273

# House keeping
rm(coxInterest)

#===========================================================================================

modelVar <- c(modelVar,"InterestRate_Nom", "InterestRate_Margin_Aggr_Med")

#============================================================================================





#============================================================================================
# ------ 4. General
varlist <- data.table(vars=c("Balance","Instalment",
                             "Principal","LN_TPE","Term"),
                      vartypes=c("int", "fin", "fin","cat", "int"))
#=========================================================================================



# ------ 4.1 Which variables are highly correlated in the group?

# - Correlation analysis
corrAnalysis(datCredit_train_AG, varlist[vartypes!="cat"]$vars, corrThresh = 0.6, method = 'spearman') # Obtain correlation groups

### RESULTS: 1) [Balance], [Installment] and [Principal] are highly correlated.



# ------ 4.2 Which variable in the correlated group should be kept in the model?

# Initialize variable
vars <- c("Balance", "Instalment", "Principal")

# Goodness of fit test
csTable(datCredit_train_AG, vars, TimeDef="AG") # No noticeable differences.

aicTable(datCredit_train_AG, vars, TimeDef="AG", genPath=genObjPath) # [Principal] has the lowest AIC

# Accuracy test
concTable(datCredit_train_AG, datCredit_valid_AG, vars, TimeDef="AG", genPath=genObjPath) # [Principal] has the highest concordance

### CONCLUSION:  Keep Principal in the model.


# ------ 4.2.1 Can [Balance] and [Installment] be replaced with [InstalmentToBalance_Aggr_Prop]?

# Initialize variable
vars <- c("Balance", "Instalment", "InstalmentToBalance_Aggr_Prop")

# Goodness of fit test of variables
aicTable(datCredit_train_AG, vars, TimeDef="AG", genPath=genObjPath) # [InstalmentToBalance_Aggr_Prop] has the lowest AIC

# Accuracy test
concTable(datCredit_train_AG, datCredit_valid_AG, vars, TimeDef="AG", genPath=genObjPath) # [InstalmentToBalance_Aggr_Prop] has the highest concordance

### CONCLUSION: [InstalmentToBalance_Aggr_Prop] can replace [Instalment] and [Balance].

varlist <- vecChange(varlist,Remove=c("Instalment","Balance"), Add=data.table(vars="InstalmentToBalance_Aggr_Prop", vartypes="prc"))


# ------ 4.2.2 How does [InstalmentToBalance_Aggr_Prop] compare to [BalanceToPrincipal]?

# Initialize variable
vars <- c("InstalmentToBalance_Aggr_Prop", "BalanceToPrincipal")

# Goodness of fit test of variables
csTable(datCredit_train_AG, vars, TimeDef="AG") # No noticeable differences.

aicTable(datCredit_train_AG, vars, TimeDef="AG", genPath=genObjPath) # [InstalmentToBalance_Aggr_Prop ] has the lowest AIC

# Accuracy test
concTable(datCredit_train_AG, datCredit_valid_AG, vars, TimeDef="AG", genPath=genObjPath) # [BalanceToPrincipal] has the highest concordance

### CONCLUSION: [BalanceToPrincipal] has a significantly greater concordance, therefore it should replace 

varlist <- vecChange(varlist,Remove=c("InstalmentToBalance_Aggr_Prop"), Add=data.table(vars="BalanceToPrincipal", vartypes="prc"))


# ------ 4.2.2 How does [BalanceToPrincipal] compare to [Principal]?

# Initialize variable
vars <- c("BalanceToPrincipal", "Principal")

# Goodness of fit test of variables
csTable(datCredit_train_AG, vars, TimeDef="AG") # No noticeable differences.

aicTable(datCredit_train_AG, vars, TimeDef="AG", genPath=genObjPath) # [BalanceToPrincipal] has the lowest AIC

# Accuracy test
concTable(datCredit_train_AG, datCredit_valid_AG, vars, TimeDef="AG", genPath=genObjPath) # [BalanceToPrincipal] has the highest concordance

### CONCLUSION: [BalanceToPrincipal] should replace [Principal]

varlist <- vecChange(varlist,Remove=c("Principal"), Add=data.table(vars="BalanceToPrincipal", vartypes="prc"))



# ------ 4.3 How does [Term] compare to [AgeToTerm_Aggr_Mean]?

# Initialize variables
vars <- c("Term", "AgeToTerm_Aggr_Mean")

# Goodness of fit
aicTable(datCredit_train_AG, vars, TimeDef="AG", genPath=genObjPath) # [AgeToTerm_Aggr_Mean] has the lowest AIC

# Accuracy test
concTable(datCredit_train_AG, datCredit_valid_AG, vars, TimeDef="AG", genPath=genObjPath) # [AgeToTerm_Aggr_Mean] has the highest concordance

### CONCLUSION: Remove [Term].

varlist <- vecChange(varlist,Remove=c("Term"), Add=data.table(vars=c("AgeToTerm_Aggr_Mean"), vartypes=c("prc")))



# ------ 4.4 What is the predictive power of the variables in univariate models?

vars <- c("LN_TPE", "BalanceToPrincipal", "AgeToTerm_Aggr_Mean")

# Goodness of fit test
csTable(datCredit_train_AG, vars, TimeDef="AG") # [AgeToTerm_Aggr_Mean] has the lowest KS-statisitc
#         Variable      D
# 4 AgeToTerm_Aggr_Mean 0.6867
# 3  BalanceToPrincipal 0.6864
# 1           Principal 0.6863
# 2              LN_TPE 0.6863

aicTable(datCredit_train_AG, vars, TimeDef="AG", genPath=genObjPath) # [AgeToTerm_Aggr_Mean] has the lowest AIC
#               Variable      AIC      pValue
# 1: AgeToTerm_Aggr_Mean 247934.1 0.008753135
# 2:  BalanceToPrincipal 248913.4 0.008831758
# 3:              LN_TPE 249221.7 0.008852371

# Accuracy test
concTable(datCredit_train_AG, datCredit_valid_AG, vars, TimeDef="AG", genPath=genObjPath) # [AgeToTerm_Aggr_Mean] has the highest concordance
#               Variable Concordance          SD LR_Statistic
# 1:  BalanceToPrincipal   0.8001090 0.004305891         5197
# 2: AgeToTerm_Aggr_Mean   0.6704183 0.004044660         6176
# 3:              LN_TPE   0.5945567 0.003070909         4889

### RESULTS:  All variables should be kept in the model.



# ------ 4.5 What is the predictive power of the variables in a single model?

# Goodness of fit
coxGen <- coxph(Surv(Start,End,Default_Ind) ~ BalanceToPrincipal + AgeToTerm_Aggr_Mean 
                + LN_TPE, id=LoanID,datCredit_train_AG)

summary(coxGen)# All variables are significant (p-value < 0.001)

GoF_CoxSnell_KS(coxGen,datCredit_train_AG,GraphInd = FALSE) # 0.6889

AIC(coxGen) # 252509.4

# Accuracy
concordance(coxGen, newdata=datCredit_valid_AG) # Concordance= 0.6449 se= 0.004189

# House keeping
rm(coxGen)

### CONCLUSION: Leave all variables in the model.


#============================================================================================

modelVar <- c(modelVar, "LN_TPE", "BalanceToPrincipal",
              "AgeToTerm_Aggr_Mean")

#============================================================================================





#============================================================================================
# ------ 5. Portfolio Level
varlist <- data.table(vars=c("CuringEvents_Aggr_Prop","NewLoans_Aggr_Prop"),
                      vartypes=c("dec", "dec"))
#============================================================================================



# ------ 5.1 Are the values correlated?

# - Correlation analysis
corrAnalysis(datCredit_train_AG, varlist[vartypes!="cat"]$vars, corrThresh = 0.6, method = 'spearman') # Obtain correlation groups

### RESULTS: The variables are not correlated significantly.



# ------ 5.2 What is the performance of the variables in univariate models?

# Initialize variable
vars <- c("CuringEvents_Aggr_Prop","NewLoans_Aggr_Prop")

# Goodness of fit test
csTable(datCredit_train_AG, vars, TimeDef="AG") # No noticeable difference

aicTable(datCredit_train_AG, vars, TimeDef="AG", genPath=genObjPath) # [NewLoans_Aggr_Prop] has the lowest AIC

# Accuracy test
concTable(datCredit_train_AG, datCredit_valid_AG, vars, TimeDef="AG", genPath=genObjPath) # [NewLoans_Aggr_Prop] has the highest concordance

### CONCOLUSION:  Keep both variables in the model



# ------ 5.3 What is the performance of the variables in a single model?

# Goodness of fit
coxPort <- coxph(Surv(Start,End,Default_Ind) ~ CuringEvents_Aggr_Prop 
                 + NewLoans_Aggr_Prop, id=LoanID,datCredit_train_AG)

summary(coxPort)
### RESULS: [CuringEvents_Aggr_Prop] has a high se(coef) (16 realative to coef -39)

### CONCLUSION: Remove [CuringEvents_Aggr_Prop]

varlist <- vecChange(varlist, Remove="CuringEvents_Aggr_Prop")



#============================================================================================

modelVar <- c(modelVar,"NewLoans_Aggr_Prop")

#============================================================================================





#============================================================================================
# ------ 6. Macro Economic
varlist <- data.table(vars=c("M_DTI_Growth","M_Emp_Growth","M_Inflation_Growth",
                             "M_RealGDP_Growth","M_RealIncome_Growth"),
                      vartypes=c("prc", "prc", "prc", "prc", "prc"))
#============================================================================================



# ------ 6.1 Which economic variables are highly correlated with one another?

# - Correlation analysis
corrAnalysis(datCredit_train_AG, varlist$vars[varlist$vartypes != 'cat'],
             corrThresh = 0.6, method = 'spearman') # Obtain correlation groups

### RESULTS: [M_Emp_Growth], [M_RealGDP_Growth] and [M_RealIncome_Growth] are highly correlated



# ------ 6.2 Which correlated variable should remain in the model?

# Initialize variables
vars <- c("M_Emp_Growth", "M_RealGDP_Growth", "M_RealIncome_Growth")

# Goodness of fit test
csTable(datCredit_train_AG, vars, TimeDef="AG") # No noticeable difference

aicTable(datCredit_train_AG, vars, TimeDef="AG", genPath=genObjPath) # [M_RealGDP_Growth] has the lowest AIC

# Accuracy test
concTable(datCredit_train_AG, datCredit_valid_AG, vars, TimeDef="AG", genPath=genObjPath) # [M_RealGDP_Growth] has the highest concordance

### CONCLUSION: Keep [M_RealGDP_Growth] in the model.

varlist = vecChange(varlist, Remove=c("M_Emp_Growth", "M_RealIncome_Growth"))



# ------ 6.3 What are the predictive powers of the current variables in the models?

# Initialize variables
vars <- c("M_DTI_Growth", "M_Inflation_Growth", "M_RealGDP_Growth")

# Goodness of fit test
csTable(datCredit_train_AG, vars, TimeDef="AG") # No noticeable difference.

aicTable(datCredit_train_AG, vars, TimeDef="AG", genPath=genObjPath) # [M_DTI_Growth] has the lowest AIC

# Accuracy test
concTable(datCredit_train_AG, datCredit_valid_AG, vars, TimeDef="AG", genPath=genObjPath) # [M_DTI_Growth] has the highest concordance

### RESULTS: [M_DTI_Growth] seems to be the best suited.



# ------ 6.4 Which lags for current variables be included in the models?


# ------ 6.4.1 Which lags for [M_DTI_Growth] be included in the models?
# Initialize variables
vars <- c("M_DTI_Growth_1","M_DTI_Growth_12","M_DTI_Growth_3", "M_DTI_Growth_6",
          "M_DTI_Growth_9","M_DTI_Growth")

# Goodness of fit test
csTable(datCredit_train_AG, vars, TimeDef="AG") # No noticeable difference

aicTable(datCredit_train_AG, vars, TimeDef="AG", genPath=genObjPath)
### RESULTS: [M_DTI_Growth], [M_DTI_Growth_1], [M_DTI_Growth_9] has the lowest AIC

# Accuracy test
concTable(datCredit_train_AG, datCredit_valid_AG, vars, TimeDef="AG", genPath=genObjPath)

### RESULTS: Only [M_DTI_Growth_1], [M_DTI_Growth_9] and [M_DTI_Growth_6] appears to have the best fit.

### CONCLUSION: Include the following lags in the model: M_DTI_Growth_9]

varlist <- vecChange(varlist,Add=data.table(vars=c("M_DTI_Growth_9"),
                                            vartypes=c("prc")))


# ------ 6.4.2 Should lags for [M_Inflation_Growth] be included in the models?
# Initialize variables
vars <- c("M_Inflation_Growth_1","M_Inflation_Growth_12","M_Inflation_Growth_3",
          "M_Inflation_Growth_6","M_Inflation_Growth_9","M_Inflation_Growth")

# Goodness of fit test
csTable(datCredit_train_AG, vars, TimeDef="AG") # No noticeable difference

aicTable(datCredit_train_AG, vars, TimeDef="AG", genPath=genObjPath)

### RESULTS: [M_Inflation_Growth_6], [M_Inflation_Growth_9] and [M_Inflation_Growth_3] appear to have the best fits.

# Accuracy test
concTable(datCredit_train_AG, datCredit_valid_AG, vars, TimeDef="AG", genPath=genObjPath)

### RESULTS: [M_Inflation_Growth_3], [M_Inflation_Growth_6] and [M_Inflation_Growth_9] appear to be the most accurate.

### CONCLUSION: Add [M_Inflation_Growth_6] to the model since it seems
###             to have the best goodness of fit and accuracy.

varlist <- vecChange(varlist,Add=data.table(vars=c("M_Inflation_Growth_6"),
                                            vartypes=c("prc")))


# ------ 6.4.3 Should lags for [M_RealIncome_Growth] be included in the models?
vars <- c("M_RealIncome_Growth_1","M_RealIncome_Growth_12","M_RealIncome_Growth_3",
          "M_RealIncome_Growth_6","M_RealIncome_Growth_9","M_RealIncome_Growth")

# Goodness of fit test
csTable(datCredit_train_AG, vars, TimeDef="AG") # No noticeable difference

aicTable(datCredit_train_AG, vars, TimeDef="AG", genPath=genObjPath)

### RESULTS: [M_RealIncome_Growth_12], [M_RealIncome_Growth_9] and [M_RealIncome_Growth_6] appear to have the best fits.

# Accuracy test
concTable(datCredit_train_AG, datCredit_valid_AG, vars, TimeDef="AG", genPath=genObjPath)

### RESULTS: [M_RealIncome_Growth_12], [M_RealIncome_Growth_9] and [M_RealIncome_Growth_6] appear to be the most accurate.

### CONCLUSION: Include [M_RealIncome_Growth_12] and [M_RealIncome_Growth_6] into the model.

varlist <- vecChange(varlist,Add=data.table(vars=c("M_RealIncome_Growth_12", "M_RealIncome_Growth_6"),
                                            vartypes=c("prc", "prc")))



# ------ 6.5 What is the predictive performance of the current univariate thematic models?

# Initialize macro-economic thematic variables
vars <- c("M_DTI_Growth", "M_Inflation_Growth", "M_RealGDP_Growth","M_DTI_Growth_1", "M_DTI_Growth_9",
          "M_Inflation_Growth_6", "M_Inflation_Growth_12", "M_RealIncome_Growth",
          "M_RealIncome_Growth_12", "M_RealIncome_Growth_6")

# Goodness of fit test
csTable(datCredit_train_AG, vars, TimeDef="AG") # No noticeable differences

aicTable(datCredit_train_AG, vars, TimeDef="AG", genPath=genObjPath)
#                 Variable      AIC      pValue
# 1:           M_DTI_Growth 247526.8 0.008854540
# 2:         M_DTI_Growth_1 247589.1 0.008858048
# 3:         M_DTI_Growth_9 247700.0 0.008867799
# 4:   M_Inflation_Growth_6 248073.6 0.008747257
# 5:     M_Inflation_Growth 248382.9 0.008744750
# 6: M_RealIncome_Growth_12 248576.1 0.008759174
# 7:  M_Inflation_Growth_12 248590.6 0.008812558
# 8:  M_RealIncome_Growth_6 249133.0 0.008778271
# 9:       M_RealGDP_Growth 249247.7 0.008797323
# 10:    M_RealIncome_Growth 249339.0 0.008820367

# Accuracy test
concTable(datCredit_train_AG, datCredit_valid_AG, vars, TimeDef="AG", genPath=genObjPath)
#                 Variable Concordance          SD LR_Statistic
# 1:         M_DTI_Growth_9   0.6879258 0.004409342         6410
# 2:         M_DTI_Growth_1   0.6858043 0.004251746         6521
# 3:           M_DTI_Growth   0.6857039 0.004231540         6584
# 4:   M_Inflation_Growth_6   0.6835126 0.004487438         6037
# 5: M_RealIncome_Growth_12   0.6653178 0.004485600         5534
# 6:     M_Inflation_Growth   0.6651781 0.004613387         5727
# 7:  M_Inflation_Growth_12   0.6628177 0.004564789         5520
# 8:  M_RealIncome_Growth_6   0.6258651 0.004616276         4977
# 9:       M_RealGDP_Growth   0.6149635 0.004656062         4863
# 10:    M_RealIncome_Growth   0.5873021 0.004778452         4771

### RESULTS: All variables should be kept in the model. Although [M_RealIncome_Growth]
###           has a low concordance, it is kept in the model as a proxy for [M_Emp_Growth].



# ------ 6.6 What is the predictive performance of the current thematic model?

# Initialize macro-economic thematic variables
vars <- c("M_DTI_Growth", "M_Inflation_Growth", "M_RealGDP_Growth","M_DTI_Growth_1", "M_DTI_Growth_9",
          "M_Inflation_Growth_6", "M_Inflation_Growth_12", "M_RealIncome_Growth",
          "M_RealIncome_Growth_12", "M_RealIncome_Growth_6")

# Build model based on macro economic variables
coxMacro <- coxph(as.formula(paste0("Surv(Start,End,Default_Ind) ~ ",
                                    paste(vars,collapse=" + "))), id=LoanID,
                  datCredit_train_AG)
summary(coxMacro)

### RESULTS: [M_RealIncome_Growth_12] and [M_RealIncome_Growth_6] have p-values greater than 0.05,
###           therefore they should be removed. Furthermore, [M_Inflation_Growth_12] and [M_DTI_Growth_1]
###           have coefficient signs opposite to that of their non-lagging counterparts, indicating a possible
###           illogical relationship.

# Initialize macro-economic thematic variables
vars <- c("M_DTI_Growth", "M_Inflation_Growth", "M_RealGDP_Growth","M_DTI_Growth_9",
          "M_Inflation_Growth_6","M_RealIncome_Growth")

# Build model based on macro economic variables
coxMacro <- coxph(as.formula(paste0("Surv(Start,End,Default_Ind) ~ ",
                                    paste(vars,collapse=" + "))), id=LoanID,
                  datCredit_train_AG)
summary(coxMacro)

vars <- c("M_DTI_Growth", "M_Inflation_Growth", "M_Inflation_Growth_6", "M_RealIncome_Growth")

# Build model based on macro economic variables
coxMacro <- coxph(as.formula(paste0("Surv(Start,End,Default_Ind) ~ ",
                                    paste(vars,collapse=" + "))), id=LoanID,
                  datCredit_train_AG)
summary(coxMacro)

### RESULTS:  All values have significant p-value.

# Test goodness of fit
GoF_CoxSnell_KS(coxMacro,datCredit_train_AG,GraphInd = FALSE) # 0.6523

AIC(coxMacro) # 187854.8

# Test accuracy
concordance(coxMacro, newdata=datCredit_valid_AG) # Concordance= 0.5876 se= 0.004408



#============================================================================================

modelVar <- c(modelVar,"M_DTI_Growth", "M_DTI_Growth_9", "M_Inflation_Growth",
              "M_Inflation_Growth_6","M_RealIncome_Growth")

#============================================================================================



# ------ 7. Check semi-final model for statistical significance and parsimony

# - Initialize variables
vars2 <- c("g0_Delinq_SD_4", "Arrears", "g0_Delinq_Ave", "TimeInDelinqState_Lag_1",      
           "slc_acct_arr_dir_3_Change_Ind", "slc_acct_roll_ever_24_imputed_mean",
           "pmnt_method_grp","slc_acct_pre_lim_perc_imputed_med","M_DTI_Growth_9",
           "InterestRate_Margin_Aggr_Med","InterestRate_Nom", "LN_TPE","NewLoans_Aggr_Prop",
           "BalanceToPrincipal","AgeToTerm_Aggr_Mean","M_DTI_Growth","M_Inflation_Growth",
           "M_Inflation_Growth_6", "M_RealIncome_Growth")

# - Build model based on variables
cox_AG <- coxph(as.formula(paste0("Surv(Start,End,Default_Ind) ~ PerfSpell_Num + ", paste(vars2,collapse=" + "))),
                 id=LoanID, datCredit_train_AG, ties="efron")
# NOTE: Default option for handling ties is recently "Efron's method", but we specify this for backwards compatability
summary(cox_AG)
### RESULTS: [InterestRate_Margin_Aggr_Med] and [NewLoans_Aggr_Prop] have significant
###           coefficient standard errors, leading to the variables being removed.

vars2 <- vars2[!(vars2 %in% c("InterestRate_Margin_Aggr_Med", "NewLoans_Aggr_Prop"))]

# - Build model based on remaining variables
cox_AG <- coxph(as.formula(paste0("Surv(Start,End,Default_Ind) ~ PerfSpell_Num + ", paste(vars2,collapse=" + "))),
                 id=LoanID, datCredit_train_AG, ties="efron")
summary(cox_AG)
### RESULTS: [M_DTI_Growth] and [M_RealIncome_Growth] have insignificant
###           p-values, therefore the variables should be removed.

vars2 <- vars2[!(vars2 %in% c("M_DTI_Growth", "M_RealIncome_Growth"))]

# - Build model based on remaining variables
cox_AG <- coxph(as.formula(paste0("Surv(Start,End,Default_Ind) ~ PerfSpell_Num + ", paste(vars2,collapse=" + "))),
                 id=LoanID, datCredit_train_AG, ties="efron")
summary(cox_AG)
### RESULTS: [M_Inflation_Growth_6] has an insignificant p-value, therefore the variables should be removed.

vars2 <- vars2[!(vars2 %in% c("M_Inflation_Growth_6"))]

# - Build model based on remaining variables
cox_AG <- coxph(as.formula(paste0("Surv(Start,End,Default_Ind) ~ PerfSpell_Num + ", paste(vars2,collapse=" + "))),
                 id=LoanID, datCredit_train_AG, ties="efron")
summary(cox_AG)
### RESULTS: All variables are significant with low standard errors for their coefficients.

#============================================================================================





# ------ 8. Final model

# - Confirm prepared datasets are loaded into memory
if (!exists('datCredit_train_AG')) unpack.ffdf(paste0(genPath,"creditdata_train_AG"), tempPath);gc()
if (!exists('datCredit_valid_AG')) unpack.ffdf(paste0(genPath,"creditdata_valid_AG"), tempPath);gc()

# - Initialize variables
vars2 <- c("PerfSpell_Num_binned","g0_Delinq_SD_4", "Arrears", "g0_Delinq_Ave", "TimeInDelinqState_Lag_1",
           "slc_acct_arr_dir_3_Change_Ind", "slc_acct_roll_ever_24_imputed_mean","LN_TPE",
           "slc_acct_pre_lim_perc_imputed_med","pmnt_method_grp","M_Inflation_Growth",
           "InterestRate_Nom", "BalanceToPrincipal","AgeToTerm_Aggr_Mean","M_DTI_Growth_9")

# - Build model based on variables
cox_AG <- coxph(as.formula(paste0("Surv(Start,End,Default_Ind) ~ ", paste(vars2,collapse=" + "))),
                 id=LoanID, datCredit_train_AG, ties="efron")
summary(cox_AG); AIC(cox_AG); concordance(cox_AG)
### RESULTS: AIC: 130465 Harrell's c: 0.9967

c <- coefficients(cox_AG)
(c <- data.table(Variable=names(c),Coefficient=c))
#                             Variable    Coefficient
#1:               PerfSpell_Num_binned -0.24508682626
#2:                     g0_Delinq_SD_4  6.81174643986
#3:                            Arrears  0.00001029733
#4:                      g0_Delinq_Ave -6.00455861524
#5:            TimeInDelinqState_Lag_1 -0.00719281717
#6:      slc_acct_arr_dir_3_Change_Ind  0.96410503046
#7: slc_acct_roll_ever_24_imputed_mean  0.63416390620
#8:                          LN_TPEWHL -0.11479813795
#9:  slc_acct_pre_lim_perc_imputed_med -2.65764500753
#10:        pmnt_method_grpMISSING_DATA  1.08217256453
#11:     pmnt_method_grpSalary/Suspense  0.68522656615
#12:           pmnt_method_grpStatement  0.36291864635
#13:                 M_Inflation_Growth  4.46077469559
#14:                   InterestRate_Nom  2.56822098859
#15:                 BalanceToPrincipal  0.28129085495
#16:                AgeToTerm_Aggr_Mean -4.39417001418
#17:                     M_DTI_Growth_9 -2.27078230596


# -Test Goodness of fit using bootstrapped B-statistics (1-KS statistic) over single-factor models
csTable_AG <- csTable(datCredit_train_AG,vars2, TimeDef="AG", seedVal=1, numIt=10,
                       fldSpellID="PerfSpell_Key", fldLstRowInd="PerfSpell_Exit_Ind", fldEventInd="Default_Ind", genPath=genPath)
### RESULTS: Top 3 single-factor models: AgeToTerm_Aggr_Mean  + g0_Delinq_Ave  + M_Inflation_Growth   

aicTable_AG <- aicTable(datCredit_train_AG, variables=vars2,fldSpellID="PerfSpell_Key",
                         TimeDef="AG", genPath=genPath)
### RESULTS: Top 3 single-factor models: g0_Delinq_SD_4 + slc_acct_roll_ever_24_imputed_mean + TimeInDelinqState_Lag_1

# Test accuracy using Harrell's c-statistic over single-factor models
concTable_AG <- concTable(datCredit_train_AG, datCredit_valid_AG, variables=vars2, 
                           fldSpellID="PerfSpell_Key", TimeDef="AG", genPath=genPath)
### RESULTS: Top x single-factor models (>90%):
# g0_Delinq_SD_4  + Arrears + TimeInDelinqState_Lag_1 + slc_acct_roll_ever_24_imputed_mean + slc_acct_arr_dir_3_Change_Ind  

# - Combine results into a single object
Table_AG <- concTable_AG[,1:2] %>% left_join(aicTable_AG, by ="Variable") %>% left_join(data.table(csTable_AG$Results), by="Variable")

# - Test Goodnes-of-fit using Cox-Snell, having measured distance between residual distribution and unit exponential using KS-statistic
GoF_CoxSnell_KS(cox_AG,datCredit_train_AG, GraphInd=TRUE, legPos=c(0.6,0.4), panelTitle="Andersen-Gill (AG) model",
                fileName = paste0(genFigPath, "AG/KS_Test_CoxSnellResiduals_Exp_AG", ".png"), dpi=280) # 0.6167
### RESULTS: Goodness of fit for the model seems to be a bit low.


concordance(cox_AG) # Concordance= 0.9971 se= 0.0001901

# Save objects
pack.ffdf(paste0(genObjPath,"AG_Univariate_Models"), Table_AG)
pack.ffdf(paste0(genPath,"AG_Cox_Model"), cox_AG)

# - Cleanup
rm(datCredit_train_TFD, datCredit_valid_TFD, cox_TFD, c); gc()
