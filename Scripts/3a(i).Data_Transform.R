# =========================================== Data enrichment ===========================================
# Prepare extractions from the base credit dataset towards creating various time definitions, themselves
#  used in subsequent Cox regression-based modelling of the default term-structure
# ------------------------------------------------------------------------------------------------------
# PROJECT TITLE: Default survival modelling
# SCRIPT AUTHOR(S): Bernard Scheepers, Dr Arno Botha (AB)

# DESCRIPTION:
# This script performs the following high-level tasks:
#   1) Extract a smaller subset of the main credit dataset
#   2) Engineer additional features in assisting with analyses and modelling later
#   3) Extract, prepare, and save subsets of the dataset according to the following time definitions: 
#       a) Time to first default; b) Anderson-Gill (AG) method; c) Prentice-Williams-Peterson (PWP) method 
#         based on Total-time approach; d) PWP-method based on Spell-time approach.
#
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

# -- Inputs:
#   - datCredit_real | Prepared from script 2b.
#
# -- Outputs:
#   - datCredit_TFD
#   - datCredit_AG
#   - datCredit_PWP_TT
#   - datCredit_PWP_ST
# ------------------------------------------------------------------------------------------------------




# --------- 1. Preliminaries
# -- Confirm prepared database is loaded into memory
if (!exists('datCredit_real')) unpack.ffdf(paste0(genPath,"creditdata_final4a"), tempPath)

Defcol <- c("DefSpell_Num", "TimeInDefSpell") # Default variables to remove
ncol <- ncol(datCredit_real)
datCredit_real[, (Defcol) := NULL]

# Sanity Check - TRUE
ncol(datCredit_real) + length(Defcol) == ncol




# --------- 2. Additional (light) feature engineering towards implementing various time definitions
# - Max date of each performance spell (used as a stratifier)
datCredit_real[!is.na(PerfSpell_Key), PerfSpell_Max_Date := max(Date, na.rm=T), by=list(PerfSpell_Key)]
datCredit_real[is.na(PerfSpell_Key), PerfSpell_Max_Date:= NA]

# - Min date of each performance spell (used as a stratifier)
datCredit_real[!is.na(PerfSpell_Key), PerfSpell_Min_Date := min(Date, na.rm=T), by=list(PerfSpell_Key)]
datCredit_real[is.na(PerfSpell_Key), PerfSpell_Min_Date:= NA]

# - Creating an indicator variable variable for when a loan exists a performance spell
datCredit_real[,PerfSpell_Exit_Ind := with(datCredit_train_TFD, ave(seq_along(PerfSpell_Key), PerfSpell_Key, FUN = function(x) x == max(x)))]

# - Creating new spell resolution types
# Performance spells
datCredit_real <- datCredit_real %>% mutate(PerfSpellResol_Type_Hist2 = case_when(PerfSpellResol_Type_Hist=="Defaulted" ~ "Defaulted",
                                                                                  PerfSpellResol_Type_Hist=="Censored" ~ "Censored",
                                                                                  TRUE ~ "Settled & Other"))
# Checking the proportions of the newly created variable
describe(datCredit_real$PerfSpellResol_Type_Hist2)




# --------- 3. Prepare data according to time definitions


# --- 3.1 Time to First Default time definition
# NOTE: The first spell condition is only enforced during the resampling of data into the training set, since we purposefully
# would like the validation set to include multiple spells in order to validate certain assumptions.
datCredit_TFD <- subset(datCredit_real, !is.na(PerfSpell_Num)) %>%
  mutate(Start = TimeInPerfSpell-1, End = TimeInPerfSpell,
         Default_Ind = DefaultStatus1)
datCredit_TFD[is.na(PerfSpell_Num),.N] == 0  # TRUE, field created successfully
# - Save snapshots to disk (zip) for quick disk-based retrieval later
pack.ffdf(paste0(genPath,"creditdata_final_TFD"), datCredit_TFD)
# - Remove from memory as an expedient
rm(datCredit_TFD); gc()


# BS: I think for the article we have to ascertain whether the Age_Adj is the best method to capute the timing component.
### AB: If we want to know the "true age" of a loan at any t, then [Age_Adj] is the best we have. [Age] is the base but flawed variable that we corrected into [Age_Adj]
# --- 3.2 Anderson-Gill (AG) time definition
datCredit_AG <- subset(datCredit_real, !is.na(PerfSpell_Num)) %>%
  mutate(Start = Age_Adj-1, End = Age_Adj,
         Default_Ind = DefaultStatus1)
datCredit_AG[is.na(PerfSpell_Num),.N] == 0  # TRUE, field created successfully
# - Save snapshots to disk (zip) for quick disk-based retrieval later
pack.ffdf(paste0(genPath,"creditdata_final_AG"), datCredit_AG)
# - Remove from memory as an expedient
rm(datCredit_AG); gc()


# BS: I think for the article we have to ascertain whether the Age_Adj is the best method to capute the timing component.
# - 3.3 Prentice-Williams-Peterson (PWP) Total-time  definition
datCredit_PWPTT <- subset(datCredit_real, !is.na(PerfSpell_Num)) %>%
  mutate(Start = Age_Adj-1, End = Age_Adj,
         Default_Ind = DefaultStatus1)
datCredit_PWPTT[is.na(PerfSpell_Num),.N] == 0  # TRUE, field created successfully
# - Save snapshots to disk (zip) for quick disk-based retrieval later
pack.ffdf(paste0(genPath,"creditdata_final_PWPTT"), datCredit_PWPTT)
# - Remove from memory as an expedient
rm(datCredit_PWPTT); gc()


# - 3.4 Prentice-Williams-Peterson (PWP) Spell-time definition
datCredit_PWPST <- subset(datCredit_real, !is.na(PerfSpell_Num)) %>%
  mutate(Start = TimeInPerfSpell-1, End = TimeInPerfSpell,
         Default_Ind = DefaultStatus1)
datCredit_PWPST[is.na(PerfSpell_Num),.N] == 0  # TRUE, field created successfully
# - Save snapshots to disk (zip) for quick disk-based retrieval later
pack.ffdf(paste0(genPath,"creditdata_final_PWPST"), datCredit_PWPST)
# - Remove from memory as an expedient
rm(datCredit_PWPST); gc()




# --------- Housekeeping
suppressWarnings(rm(datCredit_real, exclusions_all, exclusions_credit)); gc()
