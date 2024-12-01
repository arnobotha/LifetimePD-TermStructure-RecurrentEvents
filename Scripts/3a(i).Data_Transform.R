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

Defcol <- c("DefSpell_Num","TimeInDefSpell") # Default variables to remove
ncol <- ncol(datCredit_real)
datCredit_real[, (Defcol) := NULL]

# Sanity Check - TRUE
ncol(datCredit_real) + length(Defcol) == ncol




# --------- 2. Additional (light) feature engineering
# - Max date of each performance spell (used as a stratifier)
datCredit_real[!is.na(PerfSpell_Key), PerfSpell_Max_Date := max(Date, na.rm=T), by=list(PerfSpell_Key)]
datCredit_real[is.na(PerfSpell_Key), PerfSpell_Max_Date:= NA]

# - Min date of each performance spell (used as a stratifier)
datCredit_real[!is.na(PerfSpell_Key), PerfSpell_Min_Date := min(Date, na.rm=T), by=list(PerfSpell_Key)]
datCredit_real[is.na(PerfSpell_Key), PerfSpell_Min_Date:= NA]

# - Creating a variable for the first observation of a loan (used as stratification variable)
datCredit_real[, Date_First := Date[1], by=LoanID]
### AB: This variable would likely be the same as Date_Origination[1]. We can check later, verify, and possibly streamline by
# using Date_Origination[1] towards more efficient programming, especially when we are considering to publish!

# - Creating an indicator variable variable for when a loan exists a performance spell
datCredit_real[,PerfSpell_Exit_Ind := ifelse(Date==PerfSpell_Max_Date,1,0)]
### AB: I can appreciate the novelty of this variable, well done. However, its naming sucks, so I have changed it
# across the entire codebase (up to script 4b, at least, since I know you were working on the later ones)
# to something more descriptive.

# - Creating new spell resolution types
# Performance spells
datCredit_real <- datCredit_real %>% mutate(PerfSpellResol_Type_Hist2 = case_when(PerfSpellResol_Type_Hist=="Defaulted" ~ "Defaulted",
                                                                                  PerfSpellResol_Type_Hist=="Censored" ~ "Censored",
                                                                                  TRUE ~ "Settled & Other"))
# Checking the proportions of the newly created variable
datCredit_real$PerfSpellResol_Type_Hist2 %>% table() %>% prop.table()




# --------- 3. Prepare data according to time definitions


# --- 3.1a Time to First Default time definition
# NOTE: The first spell condition is only enforced during the resampling of data into the training set, since we purposefully 
# would like the validation set to include multiple spells in order to validate certain assumptions.

# Dplyr-actions include in order the following tasks:
# 1) Rename and adjust start and ending times for survival modelling
# 2) Subset for performance spells only
datCredit_TFD <- datCredit_real %>% mutate( Start = ifelse(is.na(PerfSpell_Counter),NA,PerfSpell_Counter-1), # Records the start of time interval
                                            End = ifelse(is.na(PerfSpell_Counter),NA,PerfSpell_Counter), # Records the end of time interval
                                            Default_Ind = ifelse(!is.na(PerfSpell_Num) & DefaultStatus1==1,1,
                                                                  ifelse(!is.na(PerfSpell_Num),0,NA))) %>%
                                            filter(!is.na(PerfSpell_Num)) # Filter for only performance spells
# Sanity check - Should be TRUE
datCredit_TFD[is.na(PerfSpell_Num),.N] == 0  # TRUE, field created successfully
# - Save snapshots to disk (zip) for quick disk-based retrieval later
pack.ffdf(paste0(genPath,"creditdata_final_TFD"), datCredit_TFD)
# - Remove from memory as an expedient
rm(datCredit_TFD); gc()



# --- 3.1b Time to First Default time definition | Alternate version from AB
# NOTE: The first spell condition is only enforced during the resampling of data into the training set, since we purposefully 
# would like the validation set to include multiple spells in order to validate certain assumptions.

# - Dplyr-actions include in order the following tasks:
# 1) Subset for performance spells only
# 2) Rename and adjust start and ending times for survival modelling
datCredit_TFD2 <- subset(datCredit_real, !is.na(PerfSpell_Num)) %>% 
  mutate(Start = PerfSpell_Counter-1, End = PerfSpell_Counter,
         Default_Ind = DefaultStatus1)
### AB: Using the "_Counter" variety is incorrect since it refers simply to the row index within a spell, i.e.,
# 1, ..., <Spell_Length>. Instead, one should rather use the "timeInPerfSpell" as the spell period since it 
# accounts for left-truncated spells, which are quite prevalent in our dataset.
# We should fix this after the close-out and rerun all subsequent analyses/models downstream from this point.
# It is an important update and the results should shift on some axis, though will perhaps not change in 
# distributional shape.
# Sanity check - Should be TRUE
datCredit_TFD2[is.na(PerfSpell_Num),.N] == 0  # TRUE, field created successfully
# - Save snapshots to disk (zip) for quick disk-based retrieval later
pack.ffdf(paste0(genPath,"creditdata_final_TFD2"), datCredit_TFD2)
# - Remove from memory as an expedient
rm(datCredit_TFD2); gc()



# --- 3.1c Time to First Default time definition | Alternate and corrected version from AB
# NOTE: The first spell condition is only enforced during the resampling of data into the training set, since we purposefully 
# would like the validation set to include multiple spells in order to validate certain assumptions.

# - Dplyr-actions include in order the following tasks:
# 1) Subset for performance spells only
# 2) Rename and adjust start and ending times for survival modelling
datCredit_TFD3 <- subset(datCredit_real, !is.na(PerfSpell_Num)) %>% 
  mutate(Start = TimeInPerfSpell-1, End = TimeInPerfSpell,
         Default_Ind = DefaultStatus1)
datCredit_TFD3[is.na(PerfSpell_Num),.N] == 0  # TRUE, field created successfully
# - Save snapshots to disk (zip) for quick disk-based retrieval later
pack.ffdf(paste0(genPath,"creditdata_final_TFD3"), datCredit_TFD3)
# - Remove from memory as an expedient
rm(datCredit_TFD3); gc()



# --- 3.2 Anderson-Gill (AG) time definition
datCredit_AG <- datCredit_real %>% mutate(Start = ifelse(is.na(PerfSpell_Counter),NA,Counter-1), # Records the start of time interval
                                          End = ifelse(is.na(PerfSpell_Counter),NA,Counter), # Records the end of time interval
                                          Default_Ind = ifelse(!is.na(PerfSpell_Num) & DefaultStatus1==1,1,
                                                                ifelse(!is.na(PerfSpell_Num),0,NA))) %>%
                                          filter(!is.na(PerfSpell_Num)) # Filter for only the performance spells
# Sanity check - Should be TRUE
datCredit_AG[is.na(PerfSpell_Num),.N] == 0 # TRUE, field created successfully
# - Save snapshots to disk (zip) for quick disk-based retrieval later
pack.ffdf(paste0(genPath,"creditdata_final_AG"), datCredit_AG)
# - Remove from memory as an expedient
rm(datCredit_AG); gc()



# - 3.3 Prentice-Williams-Peterson (PWP) Total-time  definition
datCredit_PWPTT <- datCredit_real %>% mutate( Start = ifelse(is.na(PerfSpell_Counter),NA,Counter-1), # Records the start of time interval
                                              End = ifelse(is.na(PerfSpell_Counter),NA,Counter), # Records the end of time interval
                                              Default_Ind = ifelse(!is.na(PerfSpell_Num) & DefaultStatus1==1,1,
                                                                     ifelse(!is.na(PerfSpell_Num),0,NA))) %>%
                                              filter(!is.na(PerfSpell_Num)) # Filter for only the performance spells
# Sanity check - Should be TRUE
datCredit_PWPTT[is.na(PerfSpell_Num),.N] == 0 # TRUE, field created successfully
# - Save snapshots to disk (zip) for quick disk-based retrieval later
pack.ffdf(paste0(genPath,"creditdata_final_PWP_TT"), datCredit_PWPTT)
# - Remove from memory as an expedient
rm(datCredit_PWPTT); gc()


# - 3.4 Prentice-Williams-Peterson (PWP) Spell-time definition
datCredit_PWPST <- datCredit_real %>% mutate( Start = ifelse(is.na(PerfSpell_Counter),NA,PerfSpell_Counter-1), # Records the start of time interval
                                              End = ifelse(is.na(PerfSpell_Counter),NA,PerfSpell_Counter), # Records the end of time interval
                                              Default_Ind = ifelse(!is.na(PerfSpell_Num) & DefaultStatus1==1,1,
                                                                   ifelse(!is.na(PerfSpell_Num),0,NA))) %>%
                                              filter(!is.na(PerfSpell_Num)) # Filter for only the performance spells
# Sanity check - Should be TRUE
datCredit_PWPST[is.na(PerfSpell_Num),.N] == 0 # TRUE, field created successfully
# - Save snapshots to disk (zip) for quick disk-based retrieval later
pack.ffdf(paste0(genPath,"creditdata_final_PWP_ST"), datCredit_PWPST)
# - Remove from memory as an expedient
rm(datCredit_PWPST); gc()




# --------- Housekeeping
suppressWarnings(rm(datCredit_real, exclusions_all, exclusions_credit)); gc()
