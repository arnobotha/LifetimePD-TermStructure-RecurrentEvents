# ========================= ANALYTICS: TERM-STRUCTURE ==========================
# Derive the term-structure of default risk across time from the various Cox
# regression models.
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
vecVars_TFD <- c("g0_Delinq_SD_4", "Arrears", "g0_Delinq_Ave",      
                 "slc_acct_arr_dir_3_Change_Ind", "slc_acct_roll_ever_24_imputed_mean",
                 "slc_acct_pre_lim_perc_imputed_med", "M_DTI_Growth",
                 "M_RealIncome_Growth","pmnt_method_grp",
                 "InterestRate_Nom", "Principal")

# - Fit a Cox Proportional Hazards model with time-varying covariates, and clustered observations
# NOTE: Assume dependence (by specifying ID-field) amongst certain observations clustered around ID-values
cox_TFD <- coxph(as.formula(paste0("Surv(Start,End,Default_Ind) ~ ", paste(vecVars_TFD,collapse=" + "))), 
                 ties="efron", id=LoanID, datCredit_train_TFD, model=T) # Keep model frame (model=T)
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
                 ties="efron", id=PerfSpell_Key, datCredit_train_PWPST, model=T) # Keep model frame (model=T)
summary(cox_PWP)





# ----------------- 3. Term-structure of default risk
# Implement the preferred (and tested) approach towards deriving the term structure


# ------ Time to first Default (TFD) definition

# --- Preliminaries
numSpellKeys <- length(unique(datCredit_valid_TFD$PerfSpell_Key))
vSpellKeys <- unique(datCredit_valid_TFD$PerfSpell_Key)



# --- Iterate across spell keys and calculate survival-related quantities using survQuants()
ptm <- proc.time() #IGNORE: for computation time calculation
cl.port <- makeCluster(8); registerDoParallel(cl.port) # multi-threading setup
cat("New Job: Estimating various survival quantities at the loan-period level for a given dataset ..",
    file=paste0(genPath,"survQuants_log.txt"), append=F)

datSurv_TFD <- foreach(j=1:numSpellKeys, .combine='rbind', .verbose=F, .inorder=T,
                   .packages=c('data.table', 'survival'), .export=c('survQuants')) %dopar%
  { # ----------------- Start of Inner Loop -----------------
    # - Testing conditions
    # j <- 1
    prepDat <- survQuants(datGiven=subset(datCredit_valid_TFD, PerfSpell_Key == vSpellKeys[j]), coxGiven = cox_TFD,
                          it=j, numKeys=numSpellKeys, genPath=genPath)
  } # ----------------- End of Inner Loop -----------------
stopCluster(cl.port); 
proc.time() - ptm; # 42.9h (single-threaded); 8.1h (multi-threaded)

# - Save snapshots to disk (zip) for quick disk-based retrieval later
pack.ffdf(paste0(genPath,"datSurv_TFD"), datSurv_TFD)


# --- Graphing the event density / probability mass function f(t)

# - Confirm prepared datasets are loaded into memory
if (!exists('datSurv_TFD')) unpack.ffdf(paste0(genPath,"datSurv_TFD"), tempPath);gc()
if (!exists('datSurv')) unpack.ffdf(paste0(genPath,"datSurv_KM_FirstSpell"), tempPath);gc()
setDT(datSurv_TFD, key="End")

# - Determine population average survival and event rate across loans per time period
datAggr <- datSurv_TFD[, list(EventRate = mean(EventProb,na.rm=T), Freq=.N),by=list(End)]
plot(datAggr[End <= 300, End], datAggr[End <= 300, EventRate], type="b")
plot(datSurv[Time <= 300, Time], datSurv[Time <= 300, EventRate], type="b")

# - General parameters
sMaxSpellAge <- 240 # max for [PerfSpell_Age], as determined in earlier analyses (script 4a(i))

# - Fitting natural cubic regression splines
sDf_Act <- 12; sDf_Exp <- 18
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
datGraph[, FacetLabel := "Time to First Default (TFD) model"]

# - Aesthetic engineering
chosenFont <- "Cambria"
vCol_Point <- brewer.pal(8, "Pastel1")[c(1,2)]
vCol_Line <- brewer.pal(8, "Set1")[c(1,2)]
mainEventName <- "Default"

# - Calculate MAE between event rates
datFusion <- merge(datSurv[Time <= sMaxSpellAge], 
                   datAggr[End <= sMaxSpellAge,list(Time=End, EventRate_Exp=EventRate)], by="Time")
MAE_eventProb <- mean(abs(datFusion$EventRate - datFusion$EventRate_Exp), na.rm=T)

# - Graphing parameters
vCol <- brewer.pal(10, "Paired")[c(3,4,1,2)]
vLabel2 <- c("b_Actual_spline"=paste0("Actual spline (df=",sDf_Act,")"), 
             "d_Expected_spline"=paste0("Scored spline (df=", sDf_Exp,")"),
             "a_Actual"="Actual", "c_Expected"="Scored")
vSize <- c(0.5,0.6,0.5,0.6)
vLineType <- c("dashed", "solid", "dashed", "solid")

# - Create main graph 
(gsurv_ft <- ggplot(datGraph[Time <= sMaxSpellAge,], aes(x=Time, y=EventRate, group=Type)) + theme_minimal() +
    labs(y=bquote(plain(Event~probability~~italic(f(t))*" ["*.(mainEventName)*"]"*"")), 
         x=bquote(Discrete~time~italic(t)*" (months) in spell: First-spell")) + 
    theme(text=element_text(family=chosenFont),legend.position = "bottom",
          axis.text.x=element_text(angle=90), #legend.text=element_text(family=chosenFont), 
          strip.background=element_rect(fill="snow2", colour="snow2"),
          strip.text=element_text(size=8, colour="gray50"), strip.text.y.right=element_text(angle=90)) + 
    # Main graph
    geom_point(aes(y=EventRatePoint, colour=Type, shape=Type), size=1.25) + 
    geom_line(aes(y=EventRate, colour=Type, linetype=Type, linewidth=Type)) + 
    # Annotations
    annotate("text", y=0.01,x=100, label=paste0("MAE: ", percent(MAE_eventProb, accuracy=0.0001)), family=chosenFont,
             size = 3) + 
    # Scales and options
    facet_grid(FacetLabel ~ .) + 
    scale_colour_manual(name="", values=vCol, labels=vLabel2) + 
    scale_linewidth_manual(name="", values=vSize, labels=vLabel2) + 
    scale_linetype_manual(name="", values=vLineType, labels=vLabel2) + 
    scale_shape_discrete(name="", labels=vLabel2) + 
    scale_y_continuous(breaks=breaks_pretty(), label=percent) + 
    scale_x_continuous(breaks=breaks_pretty(n=8), label=comma) + 
    guides(color=guide_legend(nrow=2))
)

# - Save plot
dpi <- 260 # reset
ggsave(gsurv_ft, file=paste0(genFigPath, "TFD/EventProb-", mainEventName,"_SpellLevel_FirstSpell_ActVsExp_TFD.png"),
       width=1200/dpi, height=1000/dpi,dpi=dpi, bg="white")


