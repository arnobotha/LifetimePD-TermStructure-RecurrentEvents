# =========================================== Data enrichment ===========================================
# Build Cox PH models based on the Anderson Gill and Prentice-Williams-Peterson (Total Time and Spell 
# Duration Time) time definitions
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
#   - 2d(i).Data_Fusion.R

# -- Inputs:
#   - datCredit_real | Prepared from script 2b.
#
# -- Outputs:
#   - datCredit_TFD
#   - datCredit_AG
#   - datCredit_PWP_TT
#   - datCredit_PWP_ST
# ------------------------------------------------------------------------------------------------------

# -- Confirm prepared database is loaded into memory
if (!exists('datCredit_real')) unpack.ffdf(paste0(genPath,"creditdata_final4"), tempPath)

Defcol <- c("DefSpell_Key","DefaultStatus1_lead_12","DefaultStatus1_lead_12_max","DefaultStatus2",
            "DefSpell_Num","DefSpell_Counter","TimeInDefSpell","DefSpell_LeftTrunc","DefSpell_Censored", "DefSpell_Age",
            "DefSpellResol_Type_Hist","HasLeftTruncDefSpell","DefSpell_LastStart","DefSpellResol_TimeEnd", "ReceiptPV") # Default variables to remove
ncol <- ncol(datCredit_real)
datCredit_real[, (Defcol) := NULL]

# Sanity Check - TRUE
ncol(datCredit_real) + length(Defcol) == ncol

# -- 1. Some feature engineering
# - Max date of each performance spell (used as a stratifier)
datCredit_real[!is.na(PerfSpell_Key), PerfSpell_Max_Date := max(Date, na.rm=T), by=list(PerfSpell_Key)]
datCredit_real[is.na(PerfSpell_Key), PerfSpell_Max_Date:= NA]

# - Min date of each performance spell (used as a stratifier)
datCredit_real[!is.na(PerfSpell_Key), PerfSpell_Min_Date := min(Date, na.rm=T), by=list(PerfSpell_Key)]
datCredit_real[is.na(PerfSpell_Key), PerfSpell_Min_Date:= NA]

# - Creating a variable for the first observation of a loan (used as stratification variable)
datCredit_real[, Date_First := Date[1], by=LoanID]

# - Creating new spell resolution types
# Performance spells
datCredit_real <- datCredit_real %>% mutate(PerfSpellResol_Type_Hist2 = case_when(PerfSpellResol_Type_Hist=="Defaulted" ~ "Defaulted",
                                                                                  PerfSpellResol_Type_Hist=="Censored" ~ "Censored",
                                                                                  TRUE ~ "Settled & Other"))
# Checking the proportions of the newly created variable
datCredit_real$PerfSpellResol_Type_Hist2 %>% table() %>% prop.table()

# -- 2. Prepare data according to time definitions
# - 2.1 Time to First Default time definition
datCredit_TFD <- datCredit_real %>% mutate( Start = ifelse(is.na(PerfSpell_Counter),NA,PerfSpell_Counter-1), # Records the start of time interval
                                            End = ifelse(is.na(PerfSpell_Counter),NA,PerfSpell_Counter), # Records the end of time interval
                                            Default_Ind = ifelse(!is.na(PerfSpell_Num) & DefaultStatus1==1,1,
                                                                  ifelse(!is.na(PerfSpell_Num),0,NA))) %>%
                                            filter(!is.na(PerfSpell_Num)) # Filter for only the first performance spells
# Sanity check - Should be TRUE
datCredit_TFD[is.na(PerfSpell_Num),.N] == 0  # TRUE, field created successfully

# - 2.2 Anderson Gill time definition
datCredit_AG <- datCredit_real %>% mutate(Start = ifelse(is.na(PerfSpell_Counter),NA,Counter-1), # Records the start of time interval
                                          End = ifelse(is.na(PerfSpell_Counter),NA,Counter), # Records the end of time interval
                                          Default_Ind = ifelse(!is.na(PerfSpell_Num) & DefaultStatus1==1,1,
                                                                ifelse(!is.na(PerfSpell_Num),0,NA))) %>%
                                          filter(!is.na(PerfSpell_Num)) # Filter for only the performance spells
# Sanity check - Should be TRUE
datCredit_AG[is.na(PerfSpell_Num),.N] == 0 # TRUE, field created successfully

# - 2.3 Prentice-Williams-Peterson Total Time time definition
datCredit_PWPTT <- datCredit_real %>% mutate( Start = ifelse(is.na(PerfSpell_Counter),NA,Counter-1), # Records the start of time interval
                                              End = ifelse(is.na(PerfSpell_Counter),NA,Counter), # Records the end of time interval
                                              Default_Ind = ifelse(!is.na(PerfSpell_Num) & DefaultStatus1==1,1,
                                                                     ifelse(!is.na(PerfSpell_Num),0,NA))) %>%
                                              filter(!is.na(PerfSpell_Num)) # Filter for only the performance spells
# Sanity check - Should be TRUE
datCredit_PWPTT[is.na(PerfSpell_Num),.N] == 0 # TRUE, field created successfully

# - 2.4 Prentice-Williams-Peterson Spell Time time definition
datCredit_PWPST <- datCredit_real %>% mutate( Start = ifelse(is.na(PerfSpell_Counter),NA,PerfSpell_Counter-1), # Records the start of time interval
                                              End = ifelse(is.na(PerfSpell_Counter),NA,PerfSpell_Counter), # Records the end of time interval
                                              Default_Ind = ifelse(!is.na(PerfSpell_Num) & DefaultStatus1==1,1,
                                                                   ifelse(!is.na(PerfSpell_Num),0,NA))) %>%
                                              filter(!is.na(PerfSpell_Num)) # Filter for only the performance spells
# Sanity check - Should be TRUE
datCredit_PWPST[is.na(PerfSpell_Num),.N] == 0 # TRUE, field created successfully

# -- 3. Save snapshot to disk (zip) for quick disk-based retrieval later
pack.ffdf(paste0(genPath,"creditdata_final_TFD"), datCredit_TFD)
pack.ffdf(paste0(genPath,"creditdata_final_AG"), datCredit_AG)
pack.ffdf(paste0(genPath,"creditdata_final_PWP_TT"), datCredit_PWPTT)
pack.ffdf(paste0(genPath,"creditdata_final_PWP_ST"), datCredit_PWPST)

# -- 4. Housekeeping
rm(datCredit_real)
suppressWarnings(rm(exclusions_all, exclusions_credit)); gc()