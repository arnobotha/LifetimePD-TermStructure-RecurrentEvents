# ========================= INVESTIGATING SURVIVALROC() ==========================
# As a proof of concept, we shall compare various functions from various packages
# in conducting time-dependent ROC-analyses on the same fitted Cox regression model,
# having used the prepared credit data
# --------------------------------------------------------------------------------
# PROJECT TITLE: Default Survival Modelling
# SCRIPT AUTHOR(S): Dr Arno Botha (AB)
# --------------------------------------------------------------------------------
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
#   - 3c.Data_Fusion2.R
#   - 5a(i).InputSpace_TFD.R
#   - 5a(iv).InputSpace_PWPST.R

# -- Inputs:
#   - datCredit_train_TFD | Prepared from script 3b
#   - datCredit_valid_TFD | Prepared from script 3b

#
# -- Outputs:
#   - Input_Space
# ================================================================================





# ----------------- 1. Load data

# - General parameters
sMaxSpellAge <- 300 # max for [PerfSpell_Age], as determined in earlier analyses (script 4a(i))

# - Confirm prepared datasets are loaded into memory
if (!exists('datCredit_train_TFD')) unpack.ffdf(paste0(genPath,"creditdata_train_TFD"), tempPath);gc()
if (!exists('datCredit_valid_TFD')) unpack.ffdf(paste0(genPath,"creditdata_valid_TFD"), tempPath);gc()
if (!exists('datCredit_train_PWPST')) unpack.ffdf(paste0(genPath,"creditdata_train_PWPST"), tempPath);gc()
if (!exists('datCredit_valid_PWPST')) unpack.ffdf(paste0(genPath,"creditdata_valid_PWPST"), tempPath);gc()
setDT(datCredit_train_TFD, key="PerfSpell_Key")
setDT(datCredit_valid_TFD, key="PerfSpell_Key")
setDT(datCredit_train_PWPST, key="PerfSpell_Key")
setDT(datCredit_valid_PWPST, key="PerfSpell_Key")





# ----------------- 2. Fit a Cox regression model on the resampled prepared data

# ------ Time to first Default (TFD) definition
# - Initialize variables | AB-variant
vecVars_TFD <- c("PerfSpell_g0_Delinq_Num", "Arrears", "g0_Delinq_Ave", "TimeInDelinqState_Lag_1",      
                 "slc_acct_arr_dir_3_Change_Ind", "slc_acct_pre_lim_perc_imputed_med", 
                 "InterestRate_Margin_Aggr_Med", "InterestRate_Margin_Aggr_Med_3", "M_DTI_Growth","M_Inflation_Growth",
                 "M_Inflation_Growth_6", "M_RealIncome_Growth")

# - Fit a Cox Proportional Hazards model with time-varying covariates, and clustered observations
# NOTE: Assume dependence (by specifying ID-field) amongst certain observations clustered around ID-values
cox_TFD <- coxph(as.formula(paste0("Surv(Start,End,Default_Ind) ~ ", paste(vecVars_TFD,collapse=" + "))), 
                 ties="efron", id=LoanID, datCredit_train_TFD)
summary(cox_TFD)



# ------ Prentice-Williams-Peterson (PWP) Total-time definition
# - Initialize variables | AB-variant
vecVars_PWPST <- c("g0_Delinq_SD_4", "Arrears", "g0_Delinq_Ave", "TimeInDelinqState_Lag_1",      
                 "slc_acct_arr_dir_3_Change_Ind", "slc_acct_pre_lim_perc_imputed_med", 
                 "InterestRate_Nom", "M_Repo_Rate", "M_Repo_Rate_9", 
                 "InstalmentToBalance_Aggr_Prop", "LN_TPE", "AgeToTerm",
                 "NewLoans_Aggr_Prop", "M_DTI_Growth","M_Inflation_Growth",
                 "M_Inflation_Growth_6", "M_RealIncome_Growth", "PerfSpell_Grp")

# # Fit a Cox Proportional Hazards model with time-varying covariates, and clustered observations
# # NOTE: Assume dependence (by specifying ID-field) amongst certain observations clustered around ID-values
cox_PWP <- coxph(as.formula(paste0("Surv(Start,End,Default_Ind) ~ ",paste(vecVars_PWPST,collapse=" + "))), 
                 ties="efron", id=PerfSpell_Key, datCredit_train_PWPST)
summary(cox_PWP)





# ----------------- 3. Actual & Expected Term structure for a chosen cohort

# ------ Time to first Default (TFD) definition

# --- Preliminaries
numSpellKeys <- length(unique(datCredit_valid_TFD$PerfSpell_Key))
vSpellKeys <- unique(datCredit_valid_TFD$PerfSpell_Key)



# --- Iterate across spell keys and calculate survival-related quantities using survQuants()
ptm <- proc.time() #IGNORE: for computation time calculation
#cl.port <- makeCluster(round(6)); registerDoParallel(cl.port) # multi-threading setup
cat("New Job: Estimating various survival quantities at the loan-period level for a given dataset ..",
    file=paste0(genPath,"survQuants_log.txt"), append=F)

datSurv_TFD <- foreach(j=1:numSpellKeys, .combine='rbind', .verbose=F, .inorder=T,
                   .packages=c('data.table', 'survival'), .export=c('survQuants')) %do%
  { # ----------------- Start of Inner Loop -----------------
    # - Testing conditions
    # j <- 1
    prepDat <- survQuants(datGiven=subset(datCredit_valid_TFD, PerfSpell_Key == vSpellKeys[j]), coxGiven = cox_TFD,
                          it=j, numKeys=numSpellKeys, genPath=genPath)
  } # ----------------- End of Inner Loop -----------------
#stopCluster(cl.port); 
proc.time() - ptm; # 54h

# - Save snapshots to disk (zip) for quick disk-based retrieval later
pack.ffdf(paste0(genPath,"datSurv_TFD"), datSurv_TFD)







### SCRATCH -----------------------------------
test <- subset(datCredit_valid_TFD, LoanID %in% unique(datCredit_valid_TFD[PerfSpell_Age > 5 & PerfSpell_Num > 2,LoanID])[1:3],
               select=c("LoanID", "Date", vecVars_TFD, "PerfSpell_Key", "PerfSpell_Num","PerfSpell_Counter","Start", "End", "Default_Ind"))
numSpellKeys <- length(unique(test$PerfSpell_Key))
vSpellKeys <- unique(test$PerfSpell_Key)
j <- 1
# Compare outputs
(prepDat <- survQuants(datGiven=subset(datCredit_valid_TFD, PerfSpell_Key == vSpellKeys[j]), coxGiven = cox_TFD,
                      it=j, numKeys=numSpellKeys, genPath=genPath))

prepDat2 <- survQuants.data(datGiven=subset(datCredit_valid_TFD, PerfSpell_Key == vSpellKeys[j]),
                            
                      it=j, numKeys=numSpellKeys, genPath=genPath)


### SCRATCH-END -----------------------------------


# --- Graphing the event density / probability mass function f(t)

# - Confirm prepared datasets are loaded into memory
if (!exists('datSurv_TFD')) unpack.ffdf(paste0(genPath,"datSurv_TFD"), tempPath);gc()
if (!exists('datSurv')) unpack.ffdf(paste0(genPath,"datSurv_KM_FirstSpell"), tempPath);gc()
setDT(datSurv_TFD, key="End")

# - Determine population average survival and event rate across loans per time period
datAggr <- datSurv_TFD[, list(EventRate = mean(EventProb,na.rm=T), Freq=.N),by=list(End)]
plot(datAggr[End <= sMaxSpellAge, End], datAggr[End <= sMaxSpellAge, EventRate], type="b")

# - Fitting natural cubic regression splines
sDf_Act <- 10; sDf_Exp <- 10
smthEventRate_Act <- lm(EventRate ~ ns(Time, df=sDf_Act), data=datSurv[Time <= sMaxSpellAge,])
smthEventRate_Exp <- lm(EventRate ~ ns(End, df=sDf_Exp), data=datAggr[End <= sMaxSpellAge])
summary(smthEventRate_Act); summary(smthEventRate_Exp)

# - Render predictions based on fitted smoother, with standard errors for confidence intervals
vPredSmth_Act <- predict(smthEventRate_Act, newdata=datSurv, se=T)
vPredSmth_Exp <- predict(smthEventRate_Exp, newdata=datAggr, se=T)

# - Add smoothed estimate to graphing object
datSurv[, EventRate_spline := vPredSmth_Act$fit]
datAggr[, EventRate_spline := vPredSmth_Exp$fit]

# - Create graphing data object
datGraph <- rbind(datSurv[,list(Time, EventRate, Type="a_Actual")],
                  datSurv[,list(Time, EventRate=EventRate_spline, Type="b_Actual_spline")],
                  datAggr[, list(Time=End, EventRate, Type="c_Expected")],
                  datAggr[, list(Time=End, EventRate=EventRate_spline, Type="d_Expected_spline")]
                  )

# - Create different groupings for graphing purposes
datGraph[Type %in% c("a_Actual","c_Expected"), EventRatePoint := EventRate ]
datGraph[Type %in% c("b_Actual_spline","d_Expected_spline"), EventRateLine := EventRate ]

# - Aesthetic engineering
vCol_Point <- brewer.pal(8, "Pastel1")[c(1,2)]
vCol_Line <- brewer.pal(8, "Set1")[c(1,2)]
mainEventName <- "Default"

# - Graphing parameters
vCol <- brewer.pal(10, "Paired")[c(3,4,1,2)]
vLabel2 <- c("b_Actual_spline"=paste0("Actual natural spline (df=",sDf_Act,")"), 
             "d_Expected_spline"=paste0("Scored natural spline (df=", sDf_Exp,")"),
             "a_Actual"="Actual", "c_Expected"="Scored")
vSize <- c(0.5,0.5,0.5,0.5)
vLineType <- c("dashed", "solid", "dashed", "solid")

# - Create main graph 
(gsurv_ft <- ggplot(datGraph[Time <= sMaxSpellAge,], aes(x=Time, y=EventRate, group=Type)) + theme_minimal() +
    labs(y=bquote(plain(Event~probability~~italic(f(t))*" ["*.(mainEventName)*"]"*"")), 
         x=bquote(Discrete~time~italic(t)*" (months) in spell: First-spell")) + 
    theme(text=element_text(family=chosenFont),legend.position="bottom") + 
    # Main graph
    geom_point(aes(y=EventRatePoint, colour=Type, shape=Type), size=1) + 
    geom_line(aes(y=EventRateLine, colour=Type, linetype=Type, linewidth=Type)) + 
    # Scales and options
    scale_colour_manual(name="", values=vCol, labels=vLabel2) + 
    scale_linewidth_manual(name="", values=vSize, labels=vLabel2) + 
    scale_linetype_manual(name="", values=vLineType, labels=vLabel2) + 
    scale_shape_discrete(name="", labels=vLabel2) + 
    scale_y_continuous(breaks=breaks_pretty(), label=percent) + 
    scale_x_continuous(breaks=breaks_pretty(n=8), label=comma)
  )

# - Save plot
dpi <- 180 # reset
ggsave(gsurv_ft, file=paste0(genFigPath, "TFD/EventProb-", mainEventName,"_SpellLevel_FirstSpell_ActVsExp.png"),
       width=1200/dpi, height=1000/dpi,dpi=dpi, bg="white")





### AB: First test concept before merging back

# - Merge results back to main credit dataset
datCredit_valid_TFD <- merge(datCredit_valid_TFD, datSurv, by=c("PerfSpell_Key", "End"), all.x=T) %>%
  relocate(End, CHaz, Survival, Hazard, EventProb, .after=Start)
setDT(datCredit_valid_TFD, key=c("LoanID","Date"))





# ---- SCRATCH

ptm <- proc.time() #IGNORE: for computation time calculation
survFit_pred <- survfit(cox_TFD, centered=F, newdata=datCredit_valid_TFD, id=PerfSpell_Key)
proc.time() - ptm # 86m



vOnes <- which(survFit_pred$time == 1)
matSurv <- matrix(NA, nrow=max(survFit_pred$time), ncol=length(vOnes))
k <-2
t <- 1
for (i in 1:length(survFit_pred$time)) {
  if (i < vOnes[k] ) {
    matSurv[t,k-1] <- survFit_pred$surv[i]
    t <- t+1
  } else {
    k <- k+1
    t <- 1
    matSurv[t,k-1] <- survFit_pred$surv[i]
  }
}


testCase1 <- matSurv[,1]; testCase1 <- testCase1[!is.na(testCase1)]
plot(testCase1)
testCase2 <- datSurv_TFD[PerfSpell_Key==unique(datCredit_valid_TFD$PerfSpell_Key)[1],]
plot(testCase2$Survival)
### RESULTS: Not exactly the same

# Test event density (average) and compare

