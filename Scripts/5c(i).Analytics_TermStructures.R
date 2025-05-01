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
#   - 4a(ii).Analysis_PerformanceSpells_HazardRates
#   - 5a(i).InputSpace_TFD.R
#   - 5a(ii).InputSpace_AG.R
#   - 5a(iii).InputSpace_PWPST.R

# -- Inputs:
#   - datCredit_train_TFD | Prepared from script 3b
#   - datCredit_valid_TFD | Prepared from script 3b\
#
# -- Outputs:
#   - <Analytics> | Term-structure graphs
# ================================================================================





# ----------------- 1. Load data

# - Confirm prepared datasets are loaded into memory
if (!exists('datCredit_train_TFD')) unpack.ffdf(paste0(genPath,"creditdata_train_TFD"), tempPath);gc()
if (!exists('datCredit_valid_TFD')) unpack.ffdf(paste0(genPath,"creditdata_valid_TFD"), tempPath);gc()
if (!exists('datCredit_train_AG')) unpack.ffdf(paste0(genPath,"creditdata_train_AG"), tempPath);gc()
if (!exists('datCredit_valid_AG')) unpack.ffdf(paste0(genPath,"creditdata_valid_AG"), tempPath);gc()
if (!exists('datCredit_train_PWPST')) unpack.ffdf(paste0(genPath,"creditdata_train_PWPST"), tempPath);gc()
if (!exists('datCredit_valid_PWPST')) unpack.ffdf(paste0(genPath,"creditdata_valid_PWPST"), tempPath);gc()
setDT(datCredit_train_TFD, key="PerfSpell_Key")
setDT(datCredit_valid_TFD, key="PerfSpell_Key")
setDT(datCredit_train_AG, key="PerfSpell_Key")
setDT(datCredit_valid_AG, key="PerfSpell_Key")
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
summary(cox_TFD); AIC(cox_TFD); concordance(cox_TFD)



# ------ Time to first Default (TFD) definition
# - Initialize variables | AB-variant
vecVars_AG <- c("PerfSpell_Num","g0_Delinq_SD_4", "Arrears", "g0_Delinq_Ave", "TimeInDelinqState_Lag_1",
                "slc_acct_arr_dir_3_Change_Ind", "slc_acct_roll_ever_24_imputed_mean","LN_TPE",
                "slc_acct_pre_lim_perc_imputed_med","pmnt_method_grp","M_Inflation_Growth",
                "InterestRate_Nom", "BalanceToPrincipal","AgeToTerm_Aggr_Mean","M_DTI_Growth_9")

# - Fit a Cox Proportional Hazards model with time-varying covariates, and clustered observations
# NOTE: Assume dependence (by specifying ID-field) amongst certain observations clustered around ID-values
cox_AG <- coxph(as.formula(paste0("Surv(Start,End,Default_Ind) ~ ", paste(vecVars_AG,collapse=" + "))), 
                 ties="efron", id=LoanID, datCredit_train_AG, model=T) # Keep model frame (model=T)
summary(cox_AG); AIC(cox_AG); concordance(cox_AG)



# ------ Prentice-Williams-Peterson (PWP) Total-time definition
# - Initialize variables
vecVars_PWPST <- c( # Delinquency-themed
  "g0_Delinq_SD_4", "slc_acct_roll_ever_24_imputed_mean", "g0_Delinq_Ave", "Arrears", "PerfSpell_Num",
  # Portfolio-level variables
  "AgeToTerm_Aggr_Mean",
  # Loan-level variables
  "BalanceToPrincipal", "pmnt_method_grp", "InterestRate_Nom", "slc_acct_arr_dir_3_Change_Ind",
  # Macroeconomic variables
  "M_DTI_Growth_9", "M_Inflation_Growth_6", "M_Repo_Rate_6")

# # Fit a Cox Proportional Hazards model with time-varying covariates, and clustered observations
# # NOTE: Assume dependence (by specifying ID-field) amongst certain observations clustered around ID-values
cox_PWPST <- coxph(as.formula(paste0("Surv(Start,End,Default_Ind) ~ ", paste(vecVars_PWPST,collapse=" + "), 
                                     " + strata(PerfSpell_Num_binned)")),
                   id=PerfSpell_Key, datCredit_train_PWPST, ties="efron", model=T) # Keep model frame (model=T)
summary(cox_PWPST); AIC(cox_PWPST); concordance(cox_PWPST)





# ----------------- 3. Term-structure of default risk
# Implement the preferred (and tested) approach towards deriving the term structure


# ------ Time to first Default (TFD) definition

# --- Preliminaries
numSpellKeys <- length(unique(datCredit_valid_TFD$PerfSpell_Key))
vSpellKeys <- unique(datCredit_valid_TFD$PerfSpell_Key) 

# - Testing conditions
# Get a much smaller sample on which we can test the multithreading setup within reasonable timeframes
#vSpellKeys <- data.table(PerfSpell_Key=unique(datCredit_valid_TFD$PerfSpell_Key)) %>% slice_sample(prop=0.05)
#numSpellKeys <- vSpellKeys[,.N]
# vSpellKeys <- vSpellKeys$PerfSpell_Key


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
suppressWarnings(rm(datSurv))

# - Confirm prepared datasets are loaded into memory
if (!exists('datSurv_TFD')) unpack.ffdf(paste0(genPath,"datSurv_TFD"), tempPath);gc()
if (!exists('datSurv')) unpack.ffdf(paste0(genPath,"datSurv_KM_FirstSpell"), tempPath);gc()
setDT(datSurv_TFD, key="End")

# - Determine population average survival and event rate across loans per time period
describe(datSurv_TFD$EventProb); hist(datSurv_TFD[End<=300, EventProb])
#test <- unique(subset(datSurv_TFD, EventProb > 1)$PerfSpell_Key)[1]
#j <- which(vSpellKeys==test)
datAggr <- datSurv_TFD[, list(EventRate = mean(EventProb,na.rm=T), Freq=.N),by=list(End)]
plot(datAggr[End <= 300, End], datAggr[End <= 300, EventRate], type="b") # correct shape
plot(datSurv[Time <= 300, Time], datSurv[Time <= 300, EventRate], type="b") # correct shape
cumsum(datAggr$EventRate) # can exceed 1, incorrect
cumsum(datSurv$EventRate) # approaches 1, correct

# - General parameters
sMaxSpellAge <- 240 # max for [PerfSpell_Age], as determined in earlier analyses (script 4a(i))
sMaxSpellAge_graph <- 240 # max for [PerfSpell_Age] for graphing purposes

# - Fitting natural cubic regression splines
sDf_Act <- 12; sDf_Exp <- 12
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
datGraph[, FacetLabel := "Time to First Default (TFD)"]

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
             "d_Expected_spline"=paste0("Expected spline (df=", sDf_Exp,")"),
             "a_Actual"="Actual", "c_Expected"="Expected")
vSize <- c(0.5,0.6,0.5,0.6)
vLineType <- c("dashed", "solid", "dashed", "solid")

# - Create main graph 
(gsurv_ft <- ggplot(datGraph[Time <= sMaxSpellAge_graph,], aes(x=Time, y=EventRate, group=Type)) + theme_minimal() +
    labs(y=bquote(plain(Event~probability~~italic(f(t))*" ["*.(mainEventName)*"]"*"")), 
         x=bquote(Discrete~time~italic(t)*" (months) in spell: First-spell")) + 
    theme(text=element_text(family=chosenFont),legend.position = "bottom",
          strip.background=element_rect(fill="snow2", colour="snow2"),
          strip.text=element_text(size=8, colour="gray50"), strip.text.y.right=element_text(angle=90)) + 
    # Main graph
    geom_point(aes(y=EventRatePoint, colour=Type, shape=Type), size=1.25) + 
    geom_line(aes(y=EventRate, colour=Type, linetype=Type, linewidth=Type)) + 
    # Annotations
    annotate("text", y=0.0025,x=100, label=paste0("MAE: ", percent(MAE_eventProb, accuracy=0.0001)), family=chosenFont,
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

# - Create dataset containing only the actual term-structure towards creating an inset graph
datGraph2 <- datGraph %>% subset(Type %in% c("a_Actual", "b_Actual_spline") & Time <= sMaxSpellAge_graph)

# - Create inset graph
(gsurv_ft_act <- ggplot(datGraph2, aes(x=Time, y=EventRate, group=Type)) + theme_bw() +
    labs(y="", x="", title="Actual Term-structure \n(Kaplan-Meier)") + 
    theme(legend.position.inside=c(0.75,0.40), text=element_text(size=12, family="Cambria"),
          #specific for plot-in-plot
          axis.text.y=element_text(margin=unit(c(0,0,0,0), "mm"), size=9),
          axis.text.x=element_text(margin=unit(c(0,0,0,0), "mm"), size=9),
          axis.ticks=element_blank(), axis.title.x=element_blank(), #axis.title.y=element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_rect(color="black", fill="white"),
          plot.background = element_rect(color="white"), plot.margin = unit(c(0,0,0,0),"mm"),
          plot.title = element_text(hjust=0.55,vjust=-10,margin=margin(t=-12), size=10)) + 
    # Main graph
    geom_point(aes(y=EventRatePoint, colour=Type, shape=Type), size=0.5) + 
    geom_line(aes(y=EventRate, colour=Type, linetype=Type, linewidth=Type)) + 
    # Scales and options
    scale_colour_manual(name="", values=vCol, labels=vLabel2) + 
    scale_linewidth_manual(name="", values=vSize, labels=vLabel2) + 
    scale_linetype_manual(name="", values=vLineType, labels=vLabel2) + 
    scale_shape_discrete(name="", labels=vLabel2) + 
    scale_y_continuous(breaks=breaks_pretty(), label=percent) + 
    scale_x_continuous(breaks=breaks_pretty(n=8), label=comma) + 
    guides(color="none", linewidth="none", linetype="none", shape="none"))


# - Combining the two above plots onto a single graph
ymin <- diff(ggplot_build(gsurv_ft)$layout$panel_params[[1]]$y.range) * 0.3
ymax <- max(ggplot_build(gsurv_ft)$layout$panel_params[[1]]$y.range) * 0.99
(plot.full <- gsurv_ft + annotation_custom(grob = ggplotGrob(gsurv_ft_act), xmin=15, xmax=205, 
                                           ymin=ymin, ymax=ymax))

# - Save plot
dpi <- 260 # reset
ggsave(plot.full, file=paste0(genFigPath, "TFD/EventProb-", mainEventName,"_SpellLevel_FirstSpell_ActVsExp_TFD.png"),
       width=1200/dpi, height=1000/dpi,dpi=dpi, bg="white")




# ------ Andersen-Gill (AG) definition

# --- Preliminaries
numSpellKeys <- length(unique(datCredit_valid_AG$PerfSpell_Key))
vSpellKeys <- unique(datCredit_valid_AG$PerfSpell_Key) 

# --- Iterate across spell keys and calculate survival-related quantities using survQuants()
ptm <- proc.time() #IGNORE: for computation time calculation
cl.port <- makeCluster(8); registerDoParallel(cl.port) # multi-threading setup
cat("New Job: Estimating various survival quantities at the loan-period level for a given dataset ..",
    file=paste0(genPath,"survQuants_log.txt"), append=F)

datSurv_AG <- foreach(j=1:numSpellKeys, .combine='rbind', .verbose=F, .inorder=T,
                       .packages=c('data.table', 'survival'), .export=c('survQuants')) %dopar%
  { # ----------------- Start of Inner Loop -----------------
    # - Testing conditions
    # j <- 1
    prepDat <- survQuants(datGiven=subset(datCredit_valid_AG, PerfSpell_Key == vSpellKeys[j]), coxGiven = cox_AG,
                          it=j, numKeys=numSpellKeys, genPath=genPath)
  } # ----------------- End of Inner Loop -----------------
stopCluster(cl.port); 
proc.time() - ptm; # 9.7h (multi-threaded)

# - Save snapshots to disk (zip) for quick disk-based retrieval later
pack.ffdf(paste0(genPath,"datSurv_AG"), datSurv_AG)



# --- Graphing the event density / probability mass function f(t)
suppressWarnings(rm(datSurv))

# - Confirm prepared datasets are loaded into memory
if (!exists('datSurv_AG')) unpack.ffdf(paste0(genPath,"datSurv_AG"), tempPath);gc()
if (!exists('datSurv')) unpack.ffdf(paste0(genPath,"datSurv_KM_MultiSpell"), tempPath);gc()
setDT(datSurv_AG, key="End")

# - Determine population average survival and event rate across loans per time period
describe(datSurv_AG$EventProb); hist(datSurv_AG[End<=300, EventProb])
#test <- unique(subset(datSurv_AG, EventProb > 1)$PerfSpell_Key)[1]
#j <- which(vSpellKeys==test)
datAggr <- datSurv_AG[, list(EventRate = mean(EventProb,na.rm=T), Freq=.N),by=list(End)]
plot(datAggr[End <= 300, End], datAggr[End <= 300, EventRate], type="b") # correct shape
plot(datSurv[Time <= 300, Time], datSurv[Time <= 300, EventRate], type="b") # correct shape
cumsum(datAggr$EventRate) # can exceed 1, incorrect
cumsum(datSurv$EventRate) # approaches 1, correct

# - General parameters
sMaxSpellAge <- 240 # max for [PerfSpell_Age], as determined in earlier analyses (script 4a(i))
sMaxSpellAge_graph <- 240 # max for [PerfSpell_Age] for graphing purposes

# - Fitting natural cubic regression splines
sDf_Act <- 12; sDf_Exp <- 12
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
datGraph[, FacetLabel := "Andersen-Gill (AG)"]

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
             "d_Expected_spline"=paste0("Expected spline (df=", sDf_Exp,")"),
             "a_Actual"="Actual", "c_Expected"="Expected")
vSize <- c(0.5,0.6,0.5,0.6)
vLineType <- c("dashed", "solid", "dashed", "solid")

# - Create main graph 
(gsurv_ft <- ggplot(datGraph[Time <= sMaxSpellAge_graph,], aes(x=Time, y=EventRate, group=Type)) + theme_minimal() +
    labs(y=bquote(plain(Event~probability~~italic(f(t))*" ["*.(mainEventName)*"]"*"")), 
         x=bquote(Discrete~time~italic(t)*" (months) in spell: First-spell")) + 
    theme(text=element_text(family=chosenFont),legend.position = "bottom",
          strip.background=element_rect(fill="snow2", colour="snow2"),
          strip.text=element_text(size=8, colour="gray50"), strip.text.y.right=element_text(angle=90)) + 
    # Main graph
    geom_point(aes(y=EventRatePoint, colour=Type, shape=Type), size=1.25) + 
    geom_line(aes(y=EventRate, colour=Type, linetype=Type, linewidth=Type)) + 
    # Annotations
    annotate("text", y=0.0025,x=100, label=paste0("MAE: ", percent(MAE_eventProb, accuracy=0.0001)), family=chosenFont,
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

# - Create dataset containing only the actual term-structure towards creating an inset graph
datGraph2 <- datGraph %>% subset(Type %in% c("a_Actual", "b_Actual_spline") & Time <= sMaxSpellAge_graph)

# - Create inset graph
(gsurv_ft_act <- ggplot(datGraph2, aes(x=Time, y=EventRate, group=Type)) + theme_bw() +
    labs(y="", x="", title="Actual Term-structure \n(Kaplan-Meier)") + 
    theme(legend.position.inside=c(0.75,0.40), text=element_text(size=12, family="Cambria"),
          #specific for plot-in-plot
          axis.text.y=element_text(margin=unit(c(0,0,0,0), "mm"), size=9),
          axis.text.x=element_text(margin=unit(c(0,0,0,0), "mm"), size=9),
          axis.ticks=element_blank(), axis.title.x=element_blank(), #axis.title.y=element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_rect(color="black", fill="white"),
          plot.background = element_rect(color="white"), plot.margin = unit(c(0,0,0,0),"mm"),
          plot.title = element_text(hjust=0.55,vjust=-10,margin=margin(t=-12), size=10)) + 
    # Main graph
    geom_point(aes(y=EventRatePoint, colour=Type, shape=Type), size=0.5) + 
    geom_line(aes(y=EventRate, colour=Type, linetype=Type, linewidth=Type)) + 
    # Scales and options
    scale_colour_manual(name="", values=vCol, labels=vLabel2) + 
    scale_linewidth_manual(name="", values=vSize, labels=vLabel2) + 
    scale_linetype_manual(name="", values=vLineType, labels=vLabel2) + 
    scale_shape_discrete(name="", labels=vLabel2) + 
    scale_y_continuous(breaks=breaks_pretty(), label=percent) + 
    scale_x_continuous(breaks=breaks_pretty(n=8), label=comma) + 
    guides(color="none", linewidth="none", linetype="none", shape="none"))


# - Combining the two above plots onto a single graph
ymin <- diff(ggplot_build(gsurv_ft)$layout$panel_params[[1]]$y.range) * 0.3
ymax <- max(ggplot_build(gsurv_ft)$layout$panel_params[[1]]$y.range) * 0.99
(plot.full <- gsurv_ft + annotation_custom(grob = ggplotGrob(gsurv_ft_act), xmin=25, xmax=205, 
                                           ymin=ymin, ymax=ymax))

# - Save plot
dpi <- 260 # reset
ggsave(plot.full, file=paste0(genFigPath, "AG/EventProb-", mainEventName,"_SpellLevel_FirstSpell_ActVsExp_AG.png"),
       width=1200/dpi, height=1000/dpi,dpi=dpi, bg="white")





# ------ Prentice-Williams-Peterson (PWP) definition

# --- Preliminaries
numSpellKeys <- length(unique(datCredit_valid_PWPST$PerfSpell_Key))
vSpellKeys <- unique(datCredit_valid_PWPST$PerfSpell_Key)


# --- Iterate across spell keys and calculate survival-related quantities using survQuants()
ptm <- proc.time() #IGNORE: for computation time calculation
cl.port <- makeCluster(8); registerDoParallel(cl.port) # multi-threading setup
cat("New Job: Estimating various survival quantities at the loan-period level for a given dataset ..",
    file=paste0(genPath,"survQuants_log.txt"), append=F)

datSurv_PWPST <- foreach(j=1:numSpellKeys, .combine='rbind', .verbose=F, .inorder=T,
                       .packages=c('data.table', 'survival'), .export=c('survQuants')) %dopar%
  { # ----------------- Start of Inner Loop -----------------
    # - Testing conditions
    # j <- 1
    prepDat <- survQuants(datGiven=subset(datCredit_valid_PWPST, PerfSpell_Key == vSpellKeys[j]), coxGiven = cox_PWPST,
                          it=j, numKeys=numSpellKeys, genPath=genPath)
  } # ----------------- End of Inner Loop -----------------
stopCluster(cl.port); 
proc.time() - ptm; # 10.6h (multi-threaded)

# - Save snapshots to disk (zip) for quick disk-based retrieval later
pack.ffdf(paste0(genPath,"datSurv_PWPST"), datSurv_PWPST)



# --- Graphing the event density / probability mass function f(t)
suppressWarnings(rm(datSurv))

# - Confirm prepared datasets are loaded into memory
if (!exists('datSurv_PWPST')) unpack.ffdf(paste0(genPath,"datSurv_PWPST"), tempPath);gc()
if (!exists('datSurv')) unpack.ffdf(paste0(genPath,"datSurv_KM_MultiSpell"), tempPath);gc()
setDT(datSurv_PWPST, key="End")

# - Determine population average survival and event rate across loans per time period
describe(datSurv_PWPST$EventProb); hist(datSurv_PWPST[End<=300, EventProb])
datAggr <- datSurv_PWPST[, list(EventRate = mean(EventProb,na.rm=T), Freq=.N),by=list(End)]
plot(datAggr[End <= 300, End], datAggr[End <= 300, EventRate], type="b")
plot(datSurv[Time <= 300, Time], datSurv[Time <= 300, EventRate], type="b")


# - General parameters
sMaxSpellAge <- 240 # max for [PerfSpell_Age], as determined in earlier analyses (script 4a(i))
sMaxSpellAge_graph <- 240 # max for [PerfSpell_Age] for graphing purposes

# - Fitting natural cubic regression splines
sDf_Act <- 12; sDf_Exp <- 12
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
datGraph[, FacetLabel := "Prentice-Williams-Peterson (PWP) spell-time"]

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
             "d_Expected_spline"=paste0("Expected spline (df=", sDf_Exp,")"),
             "a_Actual"="Actual", "c_Expected"="Expected")
vSize <- c(0.5,0.6,0.5,0.6)
vLineType <- c("dashed", "solid", "dashed", "solid")

# - Create main graph 
(gsurv_ft <- ggplot(datGraph[Time <= sMaxSpellAge_graph,], aes(x=Time, y=EventRate, group=Type)) + theme_minimal() +
    labs(y=bquote(plain(Event~probability~~italic(f(t))*" ["*.(mainEventName)*"]"*"")), 
         x=bquote(Discrete~time~italic(t)*" (months) in spell: Multi-spell")) + 
    theme(text=element_text(family=chosenFont),legend.position = "bottom",
          strip.background=element_rect(fill="snow2", colour="snow2"),
          strip.text=element_text(size=8, colour="gray50"), strip.text.y.right=element_text(angle=90)) + 
    # Main graph
    geom_point(aes(y=EventRatePoint, colour=Type, shape=Type), size=1.25) + 
    geom_line(aes(y=EventRate, colour=Type, linetype=Type, linewidth=Type)) + 
    # Annotations
    annotate("text", y=0.0025,x=100, label=paste0("MAE: ", percent(MAE_eventProb, accuracy=0.0001)), family=chosenFont,
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

# - Create dataset containing only the actual term-structure towards creating an inset graph
datGraph2 <- datGraph %>% subset(Type %in% c("a_Actual", "b_Actual_spline") & Time <= sMaxSpellAge_graph)

# - Create inset graph
(gsurv_ft_act <- ggplot(datGraph2, aes(x=Time, y=EventRate, group=Type)) + theme_bw() +
  labs(y="", x="", title="Actual Term-structure \n(Kaplan-Meier)") + 
  theme(legend.position.inside=c(0.75,0.40), text=element_text(size=12, family="Cambria"),
        #specific for plot-in-plot
        axis.text.y=element_text(margin=unit(c(0,0,0,0), "mm"), size=9),
        axis.text.x=element_text(margin=unit(c(0,0,0,0), "mm"), size=9),
        axis.ticks=element_blank(), axis.title.x=element_blank(), #axis.title.y=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(color="black", fill="white"),
        plot.background = element_rect(color="white"), plot.margin = unit(c(0,0,0,0),"mm"),
        plot.title = element_text(hjust=0.55,vjust=-10,margin=margin(t=-12), size=10)) + 
  # Main graph
  geom_point(aes(y=EventRatePoint, colour=Type, shape=Type), size=0.5) + 
  geom_line(aes(y=EventRate, colour=Type, linetype=Type, linewidth=Type)) + 
  # Scales and options
  scale_colour_manual(name="", values=vCol, labels=vLabel2) + 
  scale_linewidth_manual(name="", values=vSize, labels=vLabel2) + 
  scale_linetype_manual(name="", values=vLineType, labels=vLabel2) + 
  scale_shape_discrete(name="", labels=vLabel2) + 
  scale_y_continuous(breaks=breaks_pretty(), label=percent) + 
  scale_x_continuous(breaks=breaks_pretty(n=8), label=comma) + 
  guides(color="none", linewidth="none", linetype="none", shape="none"))


# - Combining the two above plots onto a single graph
ymin <- diff(ggplot_build(gsurv_ft)$layout$panel_params[[1]]$y.range) * 0.3
ymax <- max(ggplot_build(gsurv_ft)$layout$panel_params[[1]]$y.range) * 0.99
(plot.full <- gsurv_ft + annotation_custom(grob = ggplotGrob(gsurv_ft_act), xmin=26, xmax=205, 
                                           ymin=ymin, ymax=ymax))

# - Save plot
dpi <- 260 # reset
ggsave(plot.full, file=paste0(genFigPath, "PWP ST/EventProb-", mainEventName,"_SpellLevel_FirstSpell_ActVsExp_PWPST.png"),
       width=1200/dpi, height=1000/dpi,dpi=dpi, bg="white")
