# ========================= INVESTIGATING SURVIVALROC() ==========================
# As a poof of concept, we shall compare various functions from various packages
# in conducting time-dependent ROC-analyses on the same fitted Cox regression model
# --------------------------------------------------------------------------------
# PROJECT TITLE: Default Survival Modelling
# SCRIPT AUTHOR(S): Dr Arno Botha (AB)
# VERSION: 1.0 (November-2024)
# --------------------------------------------------------------------------------
# -- Script dependencies:
#   - 0.Setup.R
# ================================================================================


### AB: Compare the differences amongst TFD1-3 and establish whether or not these 
# are significant by building different Cox regression models from them.
# Method of comparison: time-dependent ROC-analysis.


head(datCredit_TFD[,list(LoanID, PerfSpell_Num, Start, End, Default_Ind)])
testIDs <- datCredit_TFD[which(datCredit_TFD$Default_ind == 1), LoanID][1:3, LoanID]