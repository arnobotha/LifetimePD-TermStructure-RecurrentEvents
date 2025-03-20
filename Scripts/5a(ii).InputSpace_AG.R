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
if (!exists('datCredit_valid_AG')) unpack.ffdf(paste0(genPath,"creditdata_final_AG"), tempPath);gc()


#============================================================================================
# ------ 1. Delinquency measures
varlist <- data.table(vars=c("g0_Delinq","g0_Delinq_fac","PerfSpell_g0_Delinq_Num",
                             "Arrears" ,"TimeInDelinqState","g0_Delinq_Any_Aggr_Prop",
                             "g0_Delinq_Ave","slc_acct_arr_dir_3",
                             "slc_acct_roll_ever_24_imputed_mean"),
                      vartypes=c("acc", "cat", "acc", "acc", "acc", "dte", "dte",
                                 "cat", "acc"))

#=========================================================================================


### AB: I suspect very strongly that PWP_ST can be copied, pasted here, and rerun accordingly





#============================================================================================
# ------ 7. Final Model

# - Confirm prepared datasets are loaded into memory
if (!exists('datCredit_train_AG')) unpack.ffdf(paste0(genPath,"creditdata_train_AG"), tempPath);gc()
if (!exists('datCredit_valid_AG')) unpack.ffdf(paste0(genPath,"creditdata_valid_AG"), tempPath);gc()

# - Initialize variables | AB-variant
vars2 <- c("PerfSpell_g0_Delinq_Num", "Arrears", "g0_Delinq_Ave", "TimeInDelinqState_Lag_1",      
           "slc_acct_arr_dir_3_Change_Ind", "slc_acct_pre_lim_perc_imputed_med", 
           "InterestRate_Margin_Aggr_Med", "InterestRate_Margin_Aggr_Med_3", 
           "M_DTI_Growth","M_Inflation_Growth", "M_Inflation_Growth_6", "M_RealIncome_Growth")

# - Build model based on variables
cox_AG <- coxph(as.formula(paste0("Surv(Start,End,Default_Ind) ~ ", paste(vars2,collapse=" + "))),
                 id=LoanID, datCredit_train_AG, ties="efron")
# NOTE: Default option for handling ties is recently "Efron's method", but we specify this for backwards compatability
summary(cox_AG)
### RESULTS: All variables are statistically significant

c <- coefficients(cox_AG)
(c <- data.table(Variable=names(c),Coefficient=c))
#                           Variable      Coefficient
#1:           PerfSpell_g0_Delinq_Num    0.01088066805
#2:                           Arrears    0.00002245914
#3:                     g0_Delinq_Ave   -6.29712945291
#4:           TimeInDelinqState_Lag_1   -0.16289352028
#5:     slc_acct_arr_dir_3_Change_Ind    3.16314704845
#6: slc_acct_pre_lim_perc_imputed_med   -8.77841646561
#7:      InterestRate_Margin_Aggr_Med  250.40823968411
#8:    InterestRate_Margin_Aggr_Med_3 -333.82470866537
#9:                      M_DTI_Growth   -5.45866171416
#10:                M_Inflation_Growth   10.26239576449
#11:              M_Inflation_Growth_6    6.86623111352
#12:               M_RealIncome_Growth  -17.83127157512

# -T Test Goodness of fit using bootstrapped B-statistics (1-KS statistic) over single-factor models
csTable_AG <- csTable(datCredit_train_AG,vars2, TimeDef="AG", seedVal=1, numIt=10,
                       fldSpellID="PerfSpell_Key", fldLstRowInd="PerfSpell_Exit_Ind", fldEventInd="Default_Ind", genPath=genPath)
### RESULTS: Top 3 single-factor models: M_RealIncome_Growth  + M_Inflation_Growth  + Arrears  
# Results are, however, very close one another such that meaningful analysis is not possible

# Test accuracy using Harrell's c-statistic over single-factor models
concTable_AG <- concTable(datCredit_train_AG, datCredit_valid_AG, variables=vars2, 
                           fldSpellID="PerfSpell_Key", TimeDef="AG", genPath=genPath)
### RESULTS: Top x single-factor models (>80%):
# Arrears + PerfSpell_g0_Delinq_Num + TimeInDelinqState_Lag_1 + slc_acct_arr_dir_3_Change_Ind  

# - Combine results into a single object
Table_AG <- left_join(csTable_AG$Results, concTable_AG, by="Variable")

# - Test Goodnes-of-fit using Cox-Snell, having measured distance between residual distribution and unit exponential using KS-statistic
GoF_CoxSnell_KS(cox_AG,datCredit_train_AG, GraphInd=TRUE, legPos=c(0.6,0.4), panelTitle="Andersen-Gill (AG) model",
                fileName = paste0(genFigPath, "AG/KS_Test_CoxSnellResiduals_Exp_AG", ".png"), dpi=280) # 0.6372
### RESULTS: Goodness of fit for the model seems to be a bit low.

# Save objects
pack.ffdf(paste0(genObjPath,"AG_Univariate_Models"), Table_AG)
pack.ffdf(paste0(genPath,"AG_Cox_Model"), cox_AG)

# - Cleanup
rm(datCredit_train_AG, datCredit_valid_AG, cox_AG, c); gc()
