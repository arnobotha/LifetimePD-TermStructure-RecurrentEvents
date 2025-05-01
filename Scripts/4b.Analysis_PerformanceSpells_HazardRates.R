# ================================= HAZARD RATE ANALYSIS =========================================
# Analysis performed on the hazard rate of different performance spells. The 
# analysis aids in the decision making about a suitable cut-off point beyond 
# which all performance spells are grouped together.
# ------------------------------------------------------------------------------------------------------
# PROJECT TITLE: Classifier Diagnostics
# SCRIPT AUTHOR(S): Dr Arno Botha, Bernard Scheepers
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

# -- Inputs:
#   - datCredit_real | Prepared from script 2f.
#
# -- Outputs:
#   - Cumulative baseline hazard rates by performance spells (2 different groupings)
#   - Event probabilities
#   - Hazard rates
#   - Kaplan-Meyer analysis by performance spells
#   - datSurv objects | Respective to each setting, containiing survival, cumulative hezard, & event probabilities
# ------------------------------------------------------------------------------------------------------




# -------- 0 Preliminarnies

# - Generic parameters
mainEventName <- "Default"





# -------- 1 Kaplan-Meier analysis across all performance spells | Left-truncation included

# - Confirm prepared datasets are loaded into memory
if (!exists('datCredit_PWPST')) unpack.ffdf(paste0(genPath,"creditdata_final_PWPST"), tempPath);gc()

# --- Fit Kaplan-Meier (KM) nonparametric model
# Compute Kaplan-Meier survival estimates (product-limit) for main-event | Spell-level with right-censoring & left-truncation
km_All <- survfit(Surv(time=Start, time2=End, event=Default_Ind==1,type="counting") ~ 1, 
                  id=PerfSpell_Key, data=datCredit_PWPST)
summary(km_All)$table # overall summary statistics
### RESULTS: 16306 events with a median survival time of 322 and 95% confidence interval of [255,262]
plot(km_All, conf.int = T) # survival curve

# -- Calculate various summary statistics
# - Chosen percentile of survival times, used later as a prediction time, e.g., 75%
# NOTE: Median would also work, though sample sizes typically dwindle at that point, which can complicate analyses
survPercentile <- 0.75
cat(paste0("Kaplan-Meier: Survival time at the xth percentile of the distribution of unique survival times: \n\t",
           percTimepoint <- min(km_All$time[km_All$surv <= survPercentile],na.rm=T)))
# - Calculate the Truncated Mean (or "Restricted Mean Survival Time" [RMST]) of survival times
# NOTE: The classical mean survival time is only defined when S(t) can reach 0
# NOTE2: Calculating RMST involves integrating the KM-curve up to the last (hence truncated) observed unique event time
#         I.e., sum the survival probabilities, as weighted by the time difference between two successive survival times
cat(paste0("Kaplan-Meier: Truncated/Restricted Mean Survival Time: \n\t", 
           sum(km_All$surv*diff(c(0,km_All$time)), na.rm=T)))
# - Retrieve the mean survival probability at the previously found percentile-based prediction point
cat(paste0("Kaplan-Meier: Retrieved mean survival probability at a given time point (t=", percTimepoint, "): \n\t",
           percent(km_All$surv[max(which(km_All$time <= percTimepoint))], accuracy=0.01)))
# - Retrieve the mean survival probability at the last unique event time as prediction time, i.e., 
# the "pooled mean survival probability"
cat(paste0("Kaplan-Meier: Retrieved mean survival probability at last unique event point as time point (t=", max(km_All$time), "): \n\t",
           percent(km_All$surv[max(which(km_All$time <= max(km_All$time)))], accuracy=0.01)))



# --- Cumulative Lifetime Distribution F(t)=1-S(t)

# - Graphing parameters
vCol <- brewer.pal(8, "Set1")[c(1)] # for F(t)
medLife <- summary(km_All)$table["median"]
median_survival <- data.frame(x = medLife, y = 0.5)
maxTime <- max(km_All$time)
chosenFont <- "Cambria"
dpi <- 190

# - Graphing logic: Entire domain of time to event
gsurv_Ft <- ggsurvplot(km_All, fun="event", conf.int=T, surv.scale = "percent", legend="none", 
                               break.time.by=round(maxTime)/8, palette=vCol,
                               xlab = bquote(Discrete~time~italic(t)*" (months) in spell: Multi-spell"),
                               ylab = bquote(CLD~"["*.(mainEventName)*"]"*~italic(F(t))*": Kaplan-Meier (incl. left-truncation)"),
                               xlim=c(0, maxTime+1), censor=F, 
                               ggtheme = theme_bw(base_family=chosenFont), tables.theme = theme_cleantable(),
                               tables.height=0.10, tables.y.text=F, tables.y.text.col=T, risk.table = "abs_pct", risk.table.pos = "out",
                               cumevents=T, cumevents.title="Cumulative number of events", 
                               cumcensor=T, cumcensor.title="Cumulative number of censored observations (incl. competing risks)",
                               risk.table.title = "Number in (% of) sample at risk of main event",
                               font.family=chosenFont, fontsize=2.5,
                               surv.median.line = "hv", size=0.5)
# Add median annotation for illustrative purposes
gsurv_Ft$plot <- gsurv_Ft$plot  +    
  annotate("text", x = medLife, y = 0.5,
           label = paste0("50th Percentile: ", medLife, " months"),
           vjust=-0.5, hjust=-0.1, color = "black", size = 2.5, family = chosenFont)
# Add Custom percentile line segment for illustrative purposes   
gsurv_Ft$plot <- gsurv_Ft$plot  +    
    geom_segment(x = 0, xend=percTimepoint, y=1-survPercentile, yend = 1-survPercentile, linetype = "dashed", color = "black") +
    geom_segment(x = percTimepoint, xend=percTimepoint, y=0, yend = 1-survPercentile, linetype = "dashed", color = "black") +
    annotate("text", x = percTimepoint, y = 1 - survPercentile, label = paste0(comma((1-survPercentile)*100), "th Percentile: ", round(percTimepoint, 2), " months"),
             vjust=-0.5, hjust=-0.1, color = "black", size = 2.5, family = chosenFont)
### RESULTS: Based on survival analysis, about 100% of the dataset was consumed at t=300; We shall henceforth remove outliers purely
# for graphing purposes
maxPeriod <- 300

# - Save graph
ggsave(gsurv_Ft$plot, file=paste0(genFigPath, "FULL SET/CumulLifetimeProb_", mainEventName,"_SpellLevel_MultiSpell_LatentComp_InclLeftTrunc.png"), width=1200/dpi, height=1000/dpi, dpi=dpi, bg="white")
dpi <- 185 # need to decrease size for risk tables' text
ggsave(print(gsurv_Ft,newpage=F), file=paste0(genFigPath,"FULL SET/CumulLifetimeProb_", mainEventName,"_SpellLevel_MultiSpell_LatentComp_InclLeftTrunc_RiskTable.png"),width=1200/dpi, height=1200/dpi,dpi=dpi, bg="white")


# - Graphing logic: Limited domain of time to event [maxPeriod]
gsurv_Ft2 <- ggsurvplot(km_All, fun="event", conf.int=T, surv.scale = "percent", legend="none", 
                       break.time.by=round(maxPeriod)/8, palette=vCol,
                       xlab = bquote(Discrete~time~italic(t)*" (months) in spell: Multi-spell"),
                       ylab = bquote(CLD~"["*.(mainEventName)*"]"*~italic(F(t))*": Kaplan-Meier (incl. left-truncation)"),
                       xlim=c(0, maxPeriod), censor=F, 
                       ggtheme = theme_bw(base_family=chosenFont), tables.theme = theme_cleantable(),
                       tables.height=0.10, tables.y.text=F, tables.y.text.col=T, risk.table = "abs_pct", risk.table.pos = "out",
                       cumevents=T, cumevents.title="Cumulative number of events", 
                       cumcensor=T, cumcensor.title="Cumulative number of censored observations (incl. competing risks)",
                       risk.table.title = "Number in (% of) sample at risk of main event",
                       font.family=chosenFont, fontsize=2.5,
                       surv.median.line = "hv", size=0.5)
# Add median annotation for illustrative purposes
gsurv_Ft2$plot <- gsurv_Ft2$plot  +    
  annotate("text", x = medLife, y = 0.51,
           label = paste0("50th Percentile: ", medLife, " months"),
           vjust=-0.5, hjust=1, color = "black", size = 2.5, family = chosenFont)
# Add Custom percentile line segment for illustrative purposes   
gsurv_Ft2$plot <- gsurv_Ft2$plot  +    
  geom_segment(x = 0, xend=percTimepoint, y=1-survPercentile, yend = 1-survPercentile, linetype = "dashed", color = "black") +
  geom_segment(x = percTimepoint, xend=percTimepoint, y=0, yend = 1-survPercentile, linetype = "dashed", color = "black") +
  annotate("text", x = percTimepoint, y = 1 - survPercentile+0.01, label = paste0(comma((1-survPercentile)*100), "th Percentile: ", round(percTimepoint, 2), " months"),
           vjust=-0.5, hjust=1, color = "black", size = 2.5, family = chosenFont)

# - Save graph
ggsave(print(gsurv_Ft2,newpage=F), file=paste0(genFigPath,"FULL SET/CumulLifetimeProb_", mainEventName,"_SpellLevel_MultiSpell_LatentComp_InclLeftTrunc_RiskTable_LimPeriod.png"),width=1200/dpi, height=1200/dpi,dpi=dpi, bg="white")



# --- Create data and graphing objects
# - Create survival table using surv_summary()
(datSurv <- surv_summary(km_All)) # Survival table
datSurv <- datSurv %>% rename(Time=time, AtRisk_n=`n.risk`, Event_n=`n.event`, Censored_n=`n.censor`, SurvivalProb_KM=`surv`) %>%
  mutate(Hazard_Actual = Event_n/AtRisk_n, Hazard_Actual2 = 1 - SurvivalProb_KM/shift(SurvivalProb_KM,n=1,fill=1)) %>% 
  mutate(Hazard_Actual2 = ifelse(is.na(Hazard_Actual2), 0, Hazard_Actual2)) %>% # Handle NaN-values
  mutate(CHaz = cumsum(Hazard_Actual),
         CHaz2 = -log(SurvivalProb_KM),# Created as a sanity check
         SurvivalProb_KM_Disc = cumprod(1-Hazard_Actual)) %>% # Created as a sanity check
  mutate(EventRate = Hazard_Actual*shift(SurvivalProb_KM, n=1, fill=1)) %>%  # probability mass function f(t)=h(t).S(t-1)
  filter(Event_n > 0 | Censored_n >0) %>% as.data.table()

# - conduct sanity checks
all.equal(datSurv$Hazard_Actual, datSurv$Hazard_Actual2) # Should be TRUE
all.equal(datSurv$SurvivalProb_KM, datSurv$SurvivalProb_KM_Disc) # Should be TRUE
all.equal(km_All$cumhaz, datSurv$CHaz) # Should be TRUE
all.equal(datSurv$CHaz, datSurv$CHaz2)
plot(datSurv$CHaz - datSurv$CHaz2, type="b")
### RESULTS: CHaz2 is very similar to CHaz, derived fundamentally from KM-estimate of S(t), though kept for comparative purposes

# - Remove sanity checks
datSurv[, c("SurvivalProb_KM_Disc","CHaz2", "Hazard_Actual2") := NULL]

# - Distributional analyses
# Hazard rate
describe(datSurv$Hazard_Actual); hist(datSurv$Hazard_Actual, breaks="FD")
plot(datSurv$Hazard_Actual, type="b")
### RESULTS: No apparent shape though increasingly affected by outliers as time increases,
# so we should put less stock in those right-most results on the x-axis.

# Event rate
describe(datSurv$EventRate); hist(datSurv$EventRate, breaks="FD")
plot(datSurv$EventRate, type="b")
### RESULTS: Shows a U-shaped distribution, as expected



# --- Graphing the hazard rate h(t)

# Fitting Locally Estimated Scatterplot Smoothing (LOESS)
sSpan <- 0.2
smthHazard_Act <- loess(datSurv$Hazard_Actual ~ datSurv[,list(x=1:.N)]$x, span=sSpan)
summary(smthHazard_Act)

# - Render predictions based on fitted smoother, with standard errors for confidence intervals
vPredSmth <- predict(smthHazard_Act, newdata=datSurv, se=T)

# - Add smoothed hazard to graphing object
datSurv[, Hazard_spline := vPredSmth$fit]
alpha <- 0.05 # significance level for confidence intervals
critVal <- qt(1-alpha/2, df=vPredSmth$df) # use t-distribution for small sample sizes (<30)
datSurv[, Hazard_spline_upper := Hazard_spline + critVal*vPredSmth$se.fit]
datSurv[, Hazard_spline_lower := Hazard_spline - critVal*vPredSmth$se.fit]

# - Graphing parameters
vCol2 <- brewer.pal(10, "Paired")[c(10,9)]
vlabel <- c(bquote('LOESS-smoothed [span: '*.(sSpan)*'] hazard'*~italic(h(t))*' with 95% CI'))
datSurv[, Group := "1"] # only necessary for aesthetics when using geom_smooth() during plotting

# - Create main graph 
(gsurv_ht <- ggplot(datSurv[Time<= maxPeriod, ], aes(x=Time, y=Hazard_Actual)) + theme_minimal() +
    # Main graph
    geom_point(aes(y=Hazard_Actual), colour=vCol2[2]) + 
    geom_line(aes(y=Hazard_Actual), colour=vCol2[2], linetype="solid") + 
    #geom_smooth(aes(colour=Group, fill=Group), se=T, method="loess", span=sSpan, alpha=0.25, linetype="dotted") +
    # Smoothed quantity
    geom_line(aes(y=Hazard_spline), colour="black", linetype="dotted") +
    geom_ribbon(aes(x=Time, ymin=Hazard_spline_lower, ymax=Hazard_spline_upper), fill=vCol2[1], alpha=0.25)+
    # Scales and options
    labs(y=bquote(plain(Discrete~hazard~~italic(h(t))*" ["*.(mainEventName)*"]"*":  Kaplan-Meier (incl. left-truncation)")), 
         x=bquote(Discrete~time~italic(t)*" (months) in spell: Multi-spell")) + 
    theme(text=element_text(family=chosenFont),legend.position="bottom") + 
    scale_colour_manual(name="", values=vCol2, labels=vlabel) + 
    scale_fill_manual(name="", values=vCol2, labels=vlabel) + 
    scale_y_continuous(breaks=breaks_pretty(), label=comma) + 
    scale_x_continuous(breaks=breaks_pretty(n=8), label=comma))
### RESULTS: The hazard appears to be U-shaped, as expected

# - Save plot
dpi <- 180 # reset
ggsave(gsurv_ht, file=paste0(genFigPath, "FULL SET/DiscreteHazard_", mainEventName,"_SpellLevel_MultiSpell-LatentComp_InclLeftTrunc-AG.png"),
       width=1200/dpi, height=1000/dpi,dpi=dpi, bg="white")



# --- Graphing the event density / probability mass function f(t)

# Fitting Locally Estimated Scatterplot Smoothing (LOESS)
sSpan <- 0.3
smthEventRate_Act <- loess(datSurv$EventRate ~ datSurv[,list(x=1:.N)]$x, span=sSpan)
summary(smthEventRate_Act)

# - Render predictions based on fitted smoother, with standard errors for confidence intervals
vPredSmth <- predict(smthEventRate_Act, newdata=datSurv, se=T)

# - Add smoothed estimate to graphing object
datSurv[, EventRate_spline := vPredSmth$fit]
alpha <- 0.05 # significance level for confidence intervals
critVal <- qt(1-alpha/2, df=vPredSmth$df) # use t-distribution for small sample sizes (<30)
datSurv[, EventRate_spline_upper := EventRate_spline + critVal*vPredSmth$se.fit]
datSurv[, EventRate_spline_lower := EventRate_spline - critVal*vPredSmth$se.fit]

# - Graphing parameters
vCol3 <- brewer.pal(10, "Paired")[c(4,3)]
vlabel <- c(bquote('LOESS-smoothed [span: '*.(sSpan)*'] estimate'*~italic(f(t))*' with 95% CI'))

# - Create main graph 
(gsurv_ft <- ggplot(datSurv[Time <= maxPeriod,], aes(x=Time, y=EventRate)) + theme_minimal() +
    # Main graph
    geom_point(aes(y=EventRate), colour=vCol3[2]) + 
    geom_line(aes(y=EventRate), colour=vCol3[2], linetype="solid") + 
    # Smoothed quantity
    geom_line(aes(y=EventRate_spline), colour="black", linetype="dotted") +
    geom_ribbon(aes(x=Time, ymin=EventRate_spline_lower, ymax=EventRate_spline_upper), fill=vCol3[1], alpha=0.15)+
    # Scales and options
    labs(y=bquote(plain(Event~probability~~italic(f(t))*" ["*.(mainEventName)*"]"*":  Kaplan-Meier (incl. left-truncation)")), 
         x=bquote(Discrete~time~italic(t)*" (months) in spell: Multi-spell")) + 
    theme(text=element_text(family=chosenFont),legend.position="bottom") + 
    scale_colour_manual(name="", values=vCol3, labels=vlabel) + 
    scale_fill_manual(name="", values=vCol3, labels=vlabel) + 
    scale_y_continuous(breaks=breaks_pretty(), label=percent) + 
    scale_x_continuous(breaks=breaks_pretty(n=8), label=comma) + 
    guides(colour = "none"))

# - Save plot
dpi <- 180 # reset
ggsave(gsurv_ft, file=paste0(genFigPath, "FULL SET/EventProb-", mainEventName,"_SpellLevel_MultiSpell_LatentComp_InclLeftTrunc.png"),
       width=1200/dpi, height=1000/dpi,dpi=dpi, bg="white")

# - Save snapshots to disk (zip) for quick disk-based retrieval later
pack.ffdf(paste0(genPath,"datSurv_KM_MultiSpell"), datSurv)

# - Housekeeping
rm(gsurv_Ft, gsurv_Ft2, gsurv_ht, gsurv_ft, km_All, median_survival, datSurv, vlabel, datCredit_PWPST, 
   smthEventRate_Act, smthHazard_Act, vPredSmth)





# -------- 2 Kaplan-Meier analysis on first performance spell | Left-truncation included

# - Confirm prepared datasets are loaded into memory
if (!exists('datCredit_TFD')) unpack.ffdf(paste0(genPath,"creditdata_final_TFD"), tempPath);gc()

# - Initialize a dataset with only the first performance spell
dat <- datCredit_TFD[PerfSpell_Num==1,]

# --- Fit Kaplan-Meier (KM) nonparametric model
# Compute Kaplan-Meier survival estimates (product-limit) for main-event | Spell-level with right-censoring & left-truncation
km_first <- survfit(Surv(time=Start, time2=End, event=Default_Ind==1,type="counting") ~ 1, 
                  id=PerfSpell_Key, data=dat)
summary(km_first)$table # overall summary statistics
### RESULTS: 16306 events with a median survival time of 322 and 95% confidence interval of [255,262]
plot(km_first, conf.int = T) # survival curve

# -- Calculate various summary statistics
# - Chosen percentile of survival times, used later as a prediction time, e.g., 75%
# NOTE: Median would also work, though sample sizes typically dwindle at that point, which can complicate analyses
survPercentile <- 0.75
cat(paste0("Kaplan-Meier: Survival time at the xth percentile of the distribution of unique survival times: \n\t",
           percTimepoint <- min(km_first$time[km_first$surv <= survPercentile],na.rm=T)))
# - Calculate the Truncated Mean (or "Restricted Mean Survival Time" [RMST]) of survival times
# NOTE: The classical mean survival time is only defined when S(t) can reach 0
# NOTE2: Calculating RMST involves integrating the KM-curve up to the last (hence truncated) observed unique event time
#         I.e., sum the survival probabilities, as weighted by the time difference between two successive survival times
cat(paste0("Kaplan-Meier: Truncated/Restricted Mean Survival Time: \n\t", 
           sum(km_first$surv*diff(c(0,km_first$time)), na.rm=T)))
# - Retrieve the mean survival probability at the previously found percentile-based prediction point
cat(paste0("Kaplan-Meier: Retrieved mean survival probability at a given time point (t=", percTimepoint, "): \n\t",
           percent(km_first$surv[max(which(km_first$time <= percTimepoint))], accuracy=0.01)))
# - Retrieve the mean survival probability at the last unique event time as prediction time, i.e., 
# the "pooled mean survival probability"
cat(paste0("Kaplan-Meier: Retrieved mean survival probability at last unique event point as time point (t=", max(km_first$time), "): \n\t",
           percent(km_first$surv[max(which(km_first$time <= max(km_first$time)))], accuracy=0.01)))



# --- Cumulative Lifetime Distribution F(t)=1-S(t)

# - Graphing parameters
vCol <- brewer.pal(8, "Set1")[c(1)] # for F(t)
medLife <- summary(km_first)$table["median"]
median_survival <- data.frame(x = medLife, y = 0.5)
maxTime <- max(km_first$time)
chosenFont <- "Cambria"
dpi <- 190

# - Graphing logic: Entire domain of time to event
gsurv_Ft <- ggsurvplot(km_first, fun="event", conf.int=T, surv.scale = "percent", legend="none", 
                       break.time.by=round(maxTime)/8, palette=vCol,
                       xlab = bquote(Discrete~time~italic(t)*" (months) in spell: First-spell"),
                       ylab = bquote(CLD~"["*.(mainEventName)*"]"*~italic(F(t))*": Kaplan-Meier (incl. left-truncation)"),
                       xlim=c(0, maxTime+1), censor=F, 
                       ggtheme = theme_bw(base_family=chosenFont), tables.theme = theme_cleantable(),
                       tables.height=0.10, tables.y.text=F, tables.y.text.col=T, risk.table = "abs_pct", risk.table.pos = "out",
                       cumevents=T, cumevents.title="Cumulative number of events", 
                       cumcensor=T, cumcensor.title="Cumulative number of censored observations (incl. competing risks)",
                       risk.table.title = "Number in (% of) sample at risk of main event",
                       font.family=chosenFont, fontsize=2.5,
                       surv.median.line = "hv", size=0.5)
# Add median annotation for illustrative purposes
gsurv_Ft$plot <- gsurv_Ft$plot  +    
  annotate("text", x = medLife, y = 0.5,
           label = paste0("50th Percentile: ", medLife, " months"),
           vjust=-0.5, hjust=-0.1, color = "black", size = 2.5, family = chosenFont)
# Add Custom percentile line segment for illustrative purposes   
gsurv_Ft$plot <- gsurv_Ft$plot  +    
  geom_segment(x = 0, xend=percTimepoint, y=1-survPercentile, yend = 1-survPercentile, linetype = "dashed", color = "black") +
  geom_segment(x = percTimepoint, xend=percTimepoint, y=0, yend = 1-survPercentile, linetype = "dashed", color = "black") +
  annotate("text", x = percTimepoint, y = 1 - survPercentile, label = paste0(comma((1-survPercentile)*100), "th Percentile: ", round(percTimepoint, 2), " months"),
           vjust=-0.5, hjust=-0.1, color = "black", size = 2.5, family = chosenFont)
### RESULTS: Based on survival analysis, about 100% of the dataset was consumed at t=300; We shall henceforth remove outliers purely
# for graphing purposes
maxPeriod <- 300

# - Save graph
ggsave(gsurv_Ft$plot, file=paste0(genFigPath, "FULL SET/CumulLifetimeProb_", mainEventName,"_SpellLevel_FirstSpell_LatentComp_InclLeftTrunc.png"), width=1200/dpi, height=1000/dpi, dpi=dpi, bg="white")
dpi <- 185 # need to decrease size for risk tables' text
ggsave(print(gsurv_Ft,newpage=F), file=paste0(genFigPath,"FULL SET/CumulLifetimeProb_", mainEventName,"_SpellLevel_FirstSpell_LatentComp_InclLeftTrunc_RiskTable.png"),width=1200/dpi, height=1200/dpi,dpi=dpi, bg="white")


# --- Create data and graphing objects
# - Create survival table using surv_summary()
(datSurv <- surv_summary(km_first)) # Survival table
datSurv <- datSurv %>% rename(Time=time, AtRisk_n=`n.risk`, Event_n=`n.event`, Censored_n=`n.censor`, SurvivalProb_KM=`surv`) %>%
  mutate(Hazard_Actual = Event_n/AtRisk_n, Hazard_Actual2 = 1 - SurvivalProb_KM/shift(SurvivalProb_KM,n=1,fill=1)) %>% 
  mutate(Hazard_Actual2 = ifelse(is.na(Hazard_Actual2), 0, Hazard_Actual2)) %>% # Handle NaN-values
  mutate(CHaz = cumsum(Hazard_Actual),
         CHaz2 = -log(SurvivalProb_KM),# Created as a sanity check
         SurvivalProb_KM_Disc = cumprod(1-Hazard_Actual)) %>% # Created as a sanity check
  mutate(EventRate = Hazard_Actual*shift(SurvivalProb_KM, n=1, fill=1)) %>%  # probability mass function f(t)=h(t).S(t-1)
  filter(Event_n > 0 | Censored_n >0) %>% as.data.table()

# - conduct sanity checks
all.equal(datSurv$Hazard_Actual, datSurv$Hazard_Actual2) # Should be TRUE
all.equal(datSurv$SurvivalProb_KM, datSurv$SurvivalProb_KM_Disc) # Should be TRUE
all.equal(km_first$cumhaz, datSurv$CHaz) # Should be TRUE
all.equal(datSurv$CHaz, datSurv$CHaz2)
plot(datSurv$CHaz - datSurv$CHaz2, type="b")
### RESULTS: CHaz2 is very similar to CHaz, derived fundamentally from KM-estimate of S(t), though kept for comparative purposes

# - Remove sanity checks
datSurv[, c("SurvivalProb_KM_Disc","CHaz2", "Hazard_Actual2") := NULL]

# - Distributional analyses
# Hazard rate
describe(datSurv$Hazard_Actual); hist(datSurv$Hazard_Actual, breaks="FD")
plot(datSurv$Hazard_Actual, type="b")
### RESULTS: No apparent shape though increasingly affected by outliers as time increases,
# so we should put less stock in those right-most results on the x-axis.

# Event rate
describe(datSurv$EventRate); hist(datSurv$EventRate, breaks="FD")
plot(datSurv$EventRate, type="b")
### RESULTS: Shows a U-shaped distribution, as expected



# --- Graphing the hazard rate h(t)

# Fitting Locally Estimated Scatterplot Smoothing (LOESS)
sSpan <- 0.2
smthHazard_Act <- loess(datSurv$Hazard_Actual ~ datSurv[,list(x=1:.N)]$x, span=sSpan)
summary(smthHazard_Act)

# - Render predictions based on fitted smoother, with standard errors for confidence intervals
vPredSmth <- predict(smthHazard_Act, newdata=datSurv, se=T)

# - Add smoothed hazard to graphing object
datSurv[, Hazard_spline := vPredSmth$fit]
alpha <- 0.05 # significance level for confidence intervals
critVal <- qt(1-alpha/2, df=vPredSmth$df) # use t-distribution for small sample sizes (<30)
datSurv[, Hazard_spline_upper := Hazard_spline + critVal*vPredSmth$se.fit]
datSurv[, Hazard_spline_lower := Hazard_spline - critVal*vPredSmth$se.fit]

# - Graphing parameters
vCol2 <- brewer.pal(10, "Paired")[c(10,9)]
vlabel <- c(bquote('LOESS-smoothed [span: '*.(sSpan)*'] hazard'*~italic(h(t))*' with 95% CI'))
datSurv[, Group := "1"] # only necessary for aesthetics when using geom_smooth() during plotting

# - Create main graph 
(gsurv_ht <- ggplot(datSurv[Time<= maxPeriod, ], aes(x=Time, y=Hazard_Actual)) + theme_minimal() +
    # Main graph
    geom_point(aes(y=Hazard_Actual), colour=vCol2[2]) + 
    geom_line(aes(y=Hazard_Actual), colour=vCol2[2], linetype="solid") + 
    #geom_smooth(aes(colour=Group, fill=Group), se=T, method="loess", span=sSpan, alpha=0.25, linetype="dotted") +
    # Smoothed quantity
    geom_line(aes(y=Hazard_spline), colour="black", linetype="dotted") +
    geom_ribbon(aes(x=Time, ymin=Hazard_spline_lower, ymax=Hazard_spline_upper), fill=vCol2[1], alpha=0.25)+
    # Scales and options
    labs(y=bquote(plain(Discrete~hazard~~italic(h(t))*" ["*.(mainEventName)*"]"*":  Kaplan-Meier (incl. left-truncation)")), 
         x=bquote(Discrete~time~italic(t)*" (months) in spell: First-spell")) + 
    theme(text=element_text(family=chosenFont),legend.position="bottom") + 
    scale_colour_manual(name="", values=vCol2, labels=vlabel) + 
    scale_fill_manual(name="", values=vCol2, labels=vlabel) + 
    scale_y_continuous(breaks=breaks_pretty(), label=comma) + 
    scale_x_continuous(breaks=breaks_pretty(n=8), label=comma))
### RESULTS: The hazard appears to be U-shaped, as expected

# - Save plot
dpi <- 180 # reset
ggsave(gsurv_ht, file=paste0(genFigPath, "FULL SET/DiscreteHazard_", mainEventName,"_SpellLevel_FirstSpell-LatentComp_InclLeftTrunc-AG.png"),
       width=1200/dpi, height=1000/dpi,dpi=dpi, bg="white")



# --- Graphing the event density / probability mass function f(t)

# Fitting Locally Estimated Scatterplot Smoothing (LOESS)
sSpan <- 0.3
smthEventRate_Act <- loess(datSurv$EventRate ~ datSurv[,list(x=1:.N)]$x, span=sSpan)
summary(smthEventRate_Act)

# - Render predictions based on fitted smoother, with standard errors for confidence intervals
vPredSmth <- predict(smthEventRate_Act, newdata=datSurv, se=T)

# - Add smoothed estimate to graphing object
datSurv[, EventRate_spline := vPredSmth$fit]
alpha <- 0.05 # significance level for confidence intervals
critVal <- qt(1-alpha/2, df=vPredSmth$df) # use t-distribution for small sample sizes (<30)
datSurv[, EventRate_spline_upper := EventRate_spline + critVal*vPredSmth$se.fit]
datSurv[, EventRate_spline_lower := EventRate_spline - critVal*vPredSmth$se.fit]

# - Graphing parameters
vCol3 <- brewer.pal(10, "Paired")[c(4,3)]
vlabel <- c(bquote('LOESS-smoothed [span: '*.(sSpan)*'] estimate'*~italic(f(t))*' with 95% CI'))

# - Create main graph 
(gsurv_ft <- ggplot(datSurv[Time <= maxPeriod,], aes(x=Time, y=EventRate)) + theme_minimal() +
    # Main graph
    geom_point(aes(y=EventRate), colour=vCol3[2]) + 
    geom_line(aes(y=EventRate), colour=vCol3[2], linetype="solid") + 
    # Smoothed quantity
    geom_line(aes(y=EventRate_spline), colour="black", linetype="dotted") +
    geom_ribbon(aes(x=Time, ymin=EventRate_spline_lower, ymax=EventRate_spline_upper), fill=vCol3[1], alpha=0.15)+
    # Scales and options
    labs(y=bquote(plain(Event~probability~~italic(f(t))*" ["*.(mainEventName)*"]"*":  Kaplan-Meier (incl. left-truncation)")), 
         x=bquote(Discrete~time~italic(t)*" (months) in spell: First-spell")) + 
    theme(text=element_text(family=chosenFont),legend.position="bottom") + 
    scale_colour_manual(name="", values=vCol3, labels=vlabel) + 
    scale_fill_manual(name="", values=vCol3, labels=vlabel) + 
    scale_y_continuous(breaks=breaks_pretty(), label=percent) + 
    scale_x_continuous(breaks=breaks_pretty(n=8), label=comma) + 
    guides(colour = "none"))

# - Save plot
dpi <- 180 # reset
ggsave(gsurv_ft, file=paste0(genFigPath, "FULL SET/EventProb-", mainEventName,"_SpellLevel_FirstSpell_LatentComp_InclLeftTrunc.png"),
       width=1200/dpi, height=1000/dpi,dpi=dpi, bg="white")

# - Save snapshots to disk (zip) for quick disk-based retrieval later
pack.ffdf(paste0(genPath,"datSurv_KM_FirstSpell"), datSurv)

# - Housekeeping
rm(gsurv_Ft, gsurv_ht, gsurv_ft, km_first, median_survival, datSurv, vlabel, datCredit_TFD, 
   smthEventRate_Act, smthHazard_Act, vPredSmth, dat)






# -------- 3 Kaplan-Meier analysis on first performance spell | Left-truncation excluded

# - Confirm prepared datasets are loaded into memory
if (!exists('datCredit_TFD')) unpack.ffdf(paste0(genPath,"creditdata_final_TFD"), tempPath);gc()

# - Initialize a dataset with only the first performance spell and excluding left-truncated spells
dat <- datCredit_TFD[PerfSpell_Num==1 & PerfSpell_LeftTrunc==0,]

# --- Fit Kaplan-Meier (KM) nonparametric model
# Compute Kaplan-Meier survival estimates (product-limit) for main-event | Spell-level with right-censoring & left-truncation
km_first <- survfit(Surv(time=Start, time2=End, event=Default_Ind==1,type="counting") ~ 1, 
                    id=PerfSpell_Key, data=dat)
summary(km_first)$table # overall summary statistics
### RESULTS: 16306 events with a median survival time of 322 and 95% confidence interval of [255,262]
plot(km_first, conf.int = T) # survival curve

# -- Calculate various summary statistics
# - Chosen percentile of survival times, used later as a prediction time, e.g., 75%
# NOTE: Median would also work, though sample sizes typically dwindle at that point, which can complicate analyses
survPercentile <- 0.75
cat(paste0("Kaplan-Meier: Survival time at the xth percentile of the distribution of unique survival times: \n\t",
           percTimepoint <- min(km_first$time[km_first$surv <= survPercentile],na.rm=T)))
# - Calculate the Truncated Mean (or "Restricted Mean Survival Time" [RMST]) of survival times
# NOTE: The classical mean survival time is only defined when S(t) can reach 0
# NOTE2: Calculating RMST involves integrating the KM-curve up to the last (hence truncated) observed unique event time
#         I.e., sum the survival probabilities, as weighted by the time difference between two successive survival times
cat(paste0("Kaplan-Meier: Truncated/Restricted Mean Survival Time: \n\t", 
           sum(km_first$surv*diff(c(0,km_first$time)), na.rm=T)))
# - Retrieve the mean survival probability at the previously found percentile-based prediction point
cat(paste0("Kaplan-Meier: Retrieved mean survival probability at a given time point (t=", percTimepoint, "): \n\t",
           percent(km_first$surv[max(which(km_first$time <= percTimepoint))], accuracy=0.01)))
# - Retrieve the mean survival probability at the last unique event time as prediction time, i.e., 
# the "pooled mean survival probability"
cat(paste0("Kaplan-Meier: Retrieved mean survival probability at last unique event point as time point (t=", max(km_first$time), "): \n\t",
           percent(km_first$surv[max(which(km_first$time <= max(km_first$time)))], accuracy=0.01)))



# --- Cumulative Lifetime Distribution F(t)=1-S(t)

# - Graphing parameters
vCol <- brewer.pal(8, "Set1")[c(1)] # for F(t)
medLife <- summary(km_first)$table["median"]
median_survival <- data.frame(x = medLife, y = 0.5)
maxTime <- max(km_first$time)
chosenFont <- "Cambria"
dpi <- 190

# - Graphing logic: Entire domain of time to event
gsurv_Ft <- ggsurvplot(km_first, fun="event", conf.int=T, surv.scale = "percent", legend="none", 
                       break.time.by=round(maxTime)/8, palette=vCol,
                       xlab = bquote(Discrete~time~italic(t)*" (months) in spell: First-spell"),
                       ylab = bquote(CLD~"["*.(mainEventName)*"]"*~italic(F(t))*": Kaplan-Meier (excl. left-truncation)"),
                       xlim=c(0, maxTime+1), censor=F, 
                       ggtheme = theme_bw(base_family=chosenFont), tables.theme = theme_cleantable(),
                       tables.height=0.10, tables.y.text=F, tables.y.text.col=T, risk.table = "abs_pct", risk.table.pos = "out",
                       cumevents=T, cumevents.title="Cumulative number of events", 
                       cumcensor=T, cumcensor.title="Cumulative number of censored observations (incl. competing risks)",
                       risk.table.title = "Number in (% of) sample at risk of main event",
                       font.family=chosenFont, fontsize=2.5,
                       surv.median.line = "hv", size=0.5)
# Add median annotation for illustrative purposes
gsurv_Ft$plot <- gsurv_Ft$plot  +    
  annotate("text", x = medLife, y = 0.5,
           label = paste0("50th Percentile: ", medLife, " months"),
           vjust=-0.5, hjust=0.5, color = "black", size = 2.5, family = chosenFont)
# Add Custom percentile line segment for illustrative purposes   
gsurv_Ft$plot <- gsurv_Ft$plot  +    
  geom_segment(x = 0, xend=percTimepoint, y=1-survPercentile, yend = 1-survPercentile, linetype = "dashed", color = "black") +
  geom_segment(x = percTimepoint, xend=percTimepoint, y=0, yend = 1-survPercentile, linetype = "dashed", color = "black") +
  annotate("text", x = percTimepoint, y = 1 - survPercentile, label = paste0(comma((1-survPercentile)*100), "th Percentile: ", round(percTimepoint, 2), " months"),
           vjust=-0.6, hjust=1, color = "black", size = 2.5, family = chosenFont)
### RESULTS: Based on survival analysis, about 100% of the dataset was consumed at t=300; We shall henceforth remove outliers purely
# for graphing purposes
maxPeriod <- 300

# - Save graph
ggsave(gsurv_Ft$plot, file=paste0(genFigPath, "FULL SET/CumulLifetimeProb_", mainEventName,"_SpellLevel_FirstSpell_LatentComp_ExclLeftTrunc.png"), width=1200/dpi, height=1000/dpi, dpi=dpi, bg="white")
dpi <- 185 # need to decrease size for risk tables' text
ggsave(print(gsurv_Ft,newpage=F), file=paste0(genFigPath,"FULL SET/CumulLifetimeProb_", mainEventName,"_SpellLevel_FirstSpell_LatentComp_ExclLeftTrunc_RiskTable.png"),width=1200/dpi, height=1200/dpi,dpi=dpi, bg="white")


# --- Create data and graphing objects
# - Create survival table using surv_summary()
(datSurv <- surv_summary(km_first)) # Survival table
datSurv <- datSurv %>% rename(Time=time, AtRisk_n=`n.risk`, Event_n=`n.event`, Censored_n=`n.censor`, SurvivalProb_KM=`surv`) %>%
  mutate(Hazard_Actual = Event_n/AtRisk_n, Hazard_Actual2 = 1 - SurvivalProb_KM/shift(SurvivalProb_KM,n=1,fill=1)) %>% 
  mutate(Hazard_Actual2 = ifelse(is.na(Hazard_Actual2), 0, Hazard_Actual2)) %>% # Handle NaN-values
  mutate(CHaz = cumsum(Hazard_Actual),
         CHaz2 = -log(SurvivalProb_KM),# Created as a sanity check
         SurvivalProb_KM_Disc = cumprod(1-Hazard_Actual)) %>% # Created as a sanity check
  mutate(EventRate = Hazard_Actual*shift(SurvivalProb_KM, n=1, fill=1)) %>%  # probability mass function f(t)=h(t).S(t-1)
  filter(Event_n > 0 | Censored_n >0) %>% as.data.table()

# - conduct sanity checks
all.equal(datSurv$Hazard_Actual, datSurv$Hazard_Actual2) # Should be TRUE
all.equal(datSurv$SurvivalProb_KM, datSurv$SurvivalProb_KM_Disc) # Should be TRUE
all.equal(km_first$cumhaz, datSurv$CHaz) # Should be TRUE
all.equal(datSurv$CHaz, datSurv$CHaz2)
plot(datSurv$CHaz - datSurv$CHaz2, type="b")
### RESULTS: CHaz2 is very similar to CHaz, derived fundamentally from KM-estimate of S(t), though kept for comparative purposes

# - Remove sanity checks
datSurv[, c("SurvivalProb_KM_Disc","CHaz2", "Hazard_Actual2") := NULL]

# - Distributional analyses
# Hazard rate
describe(datSurv$Hazard_Actual); hist(datSurv$Hazard_Actual, breaks="FD")
plot(datSurv$Hazard_Actual, type="b")
### RESULTS: No apparent shape though increasingly affected by outliers as time increases,
# so we should put less stock in those right-most results on the x-axis.

# Event rate
describe(datSurv$EventRate); hist(datSurv$EventRate, breaks="FD")
plot(datSurv$EventRate, type="b")
### RESULTS: Shows a U-shaped distribution, as expected



# --- Graphing the hazard rate h(t)

# Fitting Locally Estimated Scatterplot Smoothing (LOESS)
sSpan <- 0.2
smthHazard_Act <- loess(datSurv$Hazard_Actual ~ datSurv[,list(x=1:.N)]$x, span=sSpan)
summary(smthHazard_Act)

# - Render predictions based on fitted smoother, with standard errors for confidence intervals
vPredSmth <- predict(smthHazard_Act, newdata=datSurv, se=T)

# - Add smoothed hazard to graphing object
datSurv[, Hazard_spline := vPredSmth$fit]
alpha <- 0.05 # significance level for confidence intervals
critVal <- qt(1-alpha/2, df=vPredSmth$df) # use t-distribution for small sample sizes (<30)
datSurv[, Hazard_spline_upper := Hazard_spline + critVal*vPredSmth$se.fit]
datSurv[, Hazard_spline_lower := Hazard_spline - critVal*vPredSmth$se.fit]

# - Graphing parameters
vCol2 <- brewer.pal(10, "Paired")[c(10,9)]
vlabel <- c(bquote('LOESS-smoothed [span: '*.(sSpan)*'] hazard'*~italic(h(t))*' with 95% CI'))
datSurv[, Group := "1"] # only necessary for aesthetics when using geom_smooth() during plotting

# - Create main graph 
(gsurv_ht <- ggplot(datSurv[Time<= maxPeriod, ], aes(x=Time, y=Hazard_Actual)) + theme_minimal() +
    # Main graph
    geom_point(aes(y=Hazard_Actual), colour=vCol2[2]) + 
    geom_line(aes(y=Hazard_Actual), colour=vCol2[2], linetype="solid") + 
    #geom_smooth(aes(colour=Group, fill=Group), se=T, method="loess", span=sSpan, alpha=0.25, linetype="dotted") +
    # Smoothed quantity
    geom_line(aes(y=Hazard_spline), colour="black", linetype="dotted") +
    geom_ribbon(aes(x=Time, ymin=Hazard_spline_lower, ymax=Hazard_spline_upper), fill=vCol2[1], alpha=0.25)+
    # Scales and options
    labs(y=bquote(plain(Discrete~hazard~~italic(h(t))*" ["*.(mainEventName)*"]"*":  Kaplan-Meier (excl. left-truncation)")), 
         x=bquote(Discrete~time~italic(t)*" (months) in spell: First-spell")) + 
    theme(text=element_text(family=chosenFont),legend.position="bottom") + 
    scale_colour_manual(name="", values=vCol2, labels=vlabel) + 
    scale_fill_manual(name="", values=vCol2, labels=vlabel) + 
    scale_y_continuous(breaks=breaks_pretty(), label=comma) + 
    scale_x_continuous(breaks=breaks_pretty(n=8), label=comma))
### RESULTS: The hazard appears to be slightly U-shaped, as expected

# - Save plot
dpi <- 180 # reset
ggsave(gsurv_ht, file=paste0(genFigPath, "FULL SET/DiscreteHazard_", mainEventName,"_SpellLevel_FirstSpell-LatentComp_ExclLeftTrunc-AG.png"),
       width=1200/dpi, height=1000/dpi,dpi=dpi, bg="white")



# --- Graphing the event density / probability mass function f(t)

# Fitting Locally Estimated Scatterplot Smoothing (LOESS)
sSpan <- 0.3
smthEventRate_Act <- loess(datSurv$EventRate ~ datSurv[,list(x=1:.N)]$x, span=sSpan)
summary(smthEventRate_Act)

# - Render predictions based on fitted smoother, with standard errors for confidence intervals
vPredSmth <- predict(smthEventRate_Act, newdata=datSurv, se=T)

# - Add smoothed estimate to graphing object
datSurv[, EventRate_spline := vPredSmth$fit]
alpha <- 0.05 # significance level for confidence intervals
critVal <- qt(1-alpha/2, df=vPredSmth$df) # use t-distribution for small sample sizes (<30)
datSurv[, EventRate_spline_upper := EventRate_spline + critVal*vPredSmth$se.fit]
datSurv[, EventRate_spline_lower := EventRate_spline - critVal*vPredSmth$se.fit]

# - Graphing parameters
vCol3 <- brewer.pal(10, "Paired")[c(4,3)]
vlabel <- c(bquote('LOESS-smoothed [span: '*.(sSpan)*'] estimate'*~italic(f(t))*' with 95% CI'))

# - Create main graph 
(gsurv_ft <- ggplot(datSurv[Time <= maxPeriod,], aes(x=Time, y=EventRate)) + theme_minimal() +
    # Main graph
    geom_point(aes(y=EventRate), colour=vCol3[2]) + 
    geom_line(aes(y=EventRate), colour=vCol3[2], linetype="solid") + 
    # Smoothed quantity
    geom_line(aes(y=EventRate_spline), colour="black", linetype="dotted") +
    geom_ribbon(aes(x=Time, ymin=EventRate_spline_lower, ymax=EventRate_spline_upper), fill=vCol3[1], alpha=0.15)+
    # Scales and options
    labs(y=bquote(plain(Event~probability~~italic(f(t))*" ["*.(mainEventName)*"]"*":  Kaplan-Meier (excl. left-truncation)")), 
         x=bquote(Discrete~time~italic(t)*" (months) in spell: First-spell")) + 
    theme(text=element_text(family=chosenFont),legend.position="bottom") + 
    scale_colour_manual(name="", values=vCol3, labels=vlabel) + 
    scale_fill_manual(name="", values=vCol3, labels=vlabel) + 
    scale_y_continuous(breaks=breaks_pretty(), label=percent) + 
    scale_x_continuous(breaks=breaks_pretty(n=8), label=comma) + 
    guides(colour = "none"))

# - Save plot
dpi <- 180 # reset
ggsave(gsurv_ft, file=paste0(genFigPath, "FULL SET/EventProb-", mainEventName,"_SpellLevel_FirstSpell_LatentComp_ExclLeftTrunc.png"),
       width=1200/dpi, height=1000/dpi,dpi=dpi, bg="white")

# - Save snapshots to disk (zip) for quick disk-based retrieval later
pack.ffdf(paste0(genPath,"datSurv_KM_FirstSpell_ExclLeftTrunc"), datSurv)

# - Housekeeping
rm(gsurv_Ft, gsurv_ht, gsurv_ft, km_first, median_survival, datSurv, vlabel, datCredit_TFD, 
   smthEventRate_Act, smthHazard_Act, vPredSmth, dat)





# -------- 3 Kaplan-Meier analysis across multiple spells, having capped certain spell numbers | Left-truncation included

# - Confirm prepared datasets are loaded into memory
if (!exists('datCredit_PWPST')) unpack.ffdf(paste0(genPath,"creditdata_final_TFD_smp2"), tempPath);gc()

# - Initialize data by excluding the 9th performance spell
dat <- datCredit_PWPST[PerfSpell_Num !=9,]

# - Bin [PerfSpell_Num] based on previous analysis (script 4a(i)) towards grouping later spells together
### NOTE: This variable is created anyways in script 3c(iv), though within the subsampld set
dat[,PerfSpell_Grp := fifelse(PerfSpell_Num <= 3, PerfSpell_Num, 4)]


# --- Fit Kaplan-Meier (KM) nonparametric model
km_PerfSpells <- survfit(Surv(time=Start, time2=End, event=Default_Ind==1,type="counting") ~ PerfSpell_Grp, 
                              id=PerfSpell_Key, data=dat)
summary(km_PerfSpells)$table # overall summary statistics
### RESULTS: The 7th and 8th performance spell have no upper limit for their respective median 95th confidence interval
(km_PerfSpells_tableSummary <- surv_summary(km_PerfSpells)) # Survival table


# --- Cumulative Lifetime Distribution F(t) = 1 - S(t)

# - Graphing parameters
vCol <- brewer.pal(8, "Dark2")[1:4] # for F(t)

# Create facet groupings for performance spells
dat[, PerfSpell_Label := fifelse(PerfSpell_Num <= 3, "Perf. Spells 1-3",
                                       fifelse(PerfSpell_Num <= 5,
                                               "Perf. Spells 4-5", "Perf. Spells 6-8"))]
mainEventName <- "Default"
chosenFont <- "Cambria"
dpi <- 150

# - Graphing logic
gsurv_Ft_PerfSpells <- ggsurvplot(km_PerfSpells, data=dat, fun="event", conf.int=T, surv.scale="percent",
                                palette=vCol, xlab = bquote(Discrete~time~italic(t)*" (months) in spells"),
                                ylab = bquote(CLD~"["*.(mainEventName)*"]"*~italic(F(t))*": Kaplan-Meier (By Performance Spells)"),
                                legend.title="Performance Spells", legend="bottom",
                                legend.labs=c("Spell 1", "Spell 2","Spell 3", "Spell 4+"),
                                censor=F, ggtheme = theme_bw(base_family=chosenFont))

# - Save graph and object
ggsave(print(gsurv_Ft_PerfSpells,newpage=F), file=paste0(genFigPath, "FULL SET/gsurv_Ft_PerfSpells-", mainEventName,"_Surv-KaplanMeier-SpellLevel-MultiSpell-LatentComp_Correct.png"), width=1200/dpi, height=1000/dpi, dpi=dpi, bg="white")


# --- Discrete baseline hazard function h(t)
# - Create plotting data object
haz_dat <- data.table(Time=km_PerfSpells$time, AtRisk_n=km_PerfSpells$n.risk, 
                      Event_n = km_PerfSpells$n.event, Censored_n=km_PerfSpells$n.censor,
                      hazard=km_PerfSpells$n.event/km_PerfSpells$n.risk, 
                      CumulHazard = km_PerfSpells$cumhaz, #Nelson-Aalen estimator
                      Group=km_PerfSpells_tableSummary$PerfSpell_Grp,Surv_KM = km_PerfSpells$surv) %>% 
  filter(Event_n > 0 | Censored_n >0) %>% group_by(Group) %>%
  # Discrete-time variants
  mutate(CumulHazard_Disc = -cumsum(log(1-hazard)), Surv_KM_Disc = cumprod(1-hazard)) %>% 
  mutate(Event_KM_Disc = 1-Surv_KM_Disc) %>% ungroup() %>% as.data.table()
haz_dat[, Surv_KM_Disc_prev:= shift(Surv_KM_Disc, n=1, type="lag"), by=list(Group)]
# - create alternative versions for sanity checks
haz_dat[Time==Time[1], hazard2 := 1 - Surv_KM_Disc, by=list(Group)]
haz_dat[Time>Time[1], hazard2 := 1 - Surv_KM_Disc/Surv_KM_Disc_prev, by=list(Group)]
# - conduct sanity checks
all.equal(haz_dat$hazard, haz_dat$hazard2) # Should be TRUE
# Again the last 8 values in the first performance spell, but the first observation of [PerfSpell_Num]==8
all.equal(haz_dat$Surv_KM, haz_dat$Surv_KM_Disc) # Should be TRUE
all.equal(haz_dat$CumulHazard, haz_dat$CumulHazard_Disc) # usually FALSE
plot(haz_dat$Time, haz_dat$CumulHazard - haz_dat$CumulHazard_Disc, type="b")
### RESULTS: The discrepancy is very small difference due to estimator method differences

# --- Graphing survival and related quantities from fitted KM-model | h(t)

# - Graphing parameters
vCol2_Line <- brewer.pal(8, "Dark2")[c(1:4)] # for h(t) lines
vCol2_Point <- brewer.pal(8, "Pastel2")[c(1:4)] # for h(t) points
# Create facet groupings for performance spells
haz_dat[, GroupLabel := factor(Group, levels = 1:4,
                               labels = c(
                                 bquote(~italic(h[0](t))*": Performance Spell 1"),
                                 bquote(~italic(h[0](t))*": Performance Spell 2"),
                                 bquote(~italic(h[0](t))*": Performance Spell 3"),
                                 bquote(~italic(h[0](t))*": Performance Spell 4+")))]

# - Fitting natural cubic splines
sDf1 <- 5; sDf2 <- 4; sDf3 <- 5; sDf4 <- 5;
smthHazard_Act1 <- lm(haz_dat[Group==1,hazard] ~ ns(haz_dat[Group==1,Time],df=sDf1))
smthHazard_Act2 <- lm(haz_dat[Group==2,hazard] ~ ns(haz_dat[Group==2,Time],df=sDf2))
smthHazard_Act3 <- lm(haz_dat[Group==3,hazard] ~ ns(haz_dat[Group==3,Time],df=sDf3))
smthHazard_Act4 <- lm(haz_dat[Group==4,hazard] ~ ns(haz_dat[Group==4,Time],df=sDf4))
summary(smthHazard_Act1); summary(smthHazard_Act2)
summary(smthHazard_Act3); summary(smthHazard_Act4)
plot(predict(smthHazard_Act1, newdata=data.table(stop=unique(haz_dat[Group==1, Time]))))
plot(predict(smthHazard_Act2, newdata=data.table(stop=unique(haz_dat[Group==2, Time]))))
plot(predict(smthHazard_Act3, newdata=data.table(stop=unique(haz_dat[Group==3, Time]))))
plot(predict(smthHazard_Act4, newdata=data.table(stop=unique(haz_dat[Group==4, Time]))))

# - Render predictions
haz_dat[Group==1, Spline := predict(smthHazard_Act1, newdata=data.table(stop=unique(haz_dat[Group==1, Time])))]
haz_dat[Group==2, Spline := predict(smthHazard_Act2, newdata=data.table(stop=unique(haz_dat[Group==2, Time])))]
haz_dat[Group==3, Spline := predict(smthHazard_Act3, newdata=data.table(stop=unique(haz_dat[Group==3, Time])))]
haz_dat[Group==4, Spline := predict(smthHazard_Act4, newdata=data.table(stop=unique(haz_dat[Group==4, Time])))]

# - Aesthetic engineering
vLabel <- c("1"=paste0("Group 1 natural spline (df=", sDf1, ")"), "2"=paste0("Group 2 natural spline (df=", sDf2, ")"),
            "3"=paste0("Group 3 natural spline (df=", sDf3, ")"), "4"=paste0("Group 1 natural spline (df=", sDf4, ")"))

# - Graph object for shorter time, informed by previous graphs
(gsurv_ht_PerfSpell <- ggplot(haz_dat, aes(x=Time)) + theme_minimal() +
    labs(y=bquote(plain(Estimated~baseline~hazard*" ["*.(mainEventName)*"]"*~italic(h[0](t))*":  Kaplan-Meier (spell-level)")), 
         x=bquote(Discrete~time~italic(t)*" (months) in each spell: Single-spell")) + 
    theme(text=element_text(family=chosenFont),legend.position="bottom",
          strip.background=element_rect(fill="snow2", colour="snow2"),
          strip.text=element_text(size=8, colour="gray50")) +
    geom_point(aes(y=hazard, colour=Group), alpha=0.25) + 
    geom_line(aes(y=Spline, colour=Group), linetype="solid") +
    facet_wrap(~GroupLabel, scales="free", labeller=label_parsed) + 
    scale_colour_manual(name="", values=vCol2_Line, labels=vLabel) +
    scale_y_continuous(breaks=breaks_pretty(), label=comma) + 
    scale_x_continuous(breaks=breaks_pretty(n=8), label=comma))
### RESULTS:  As the performance spell increase, the volatility in the hazard rate increases accordingly. Note that
###           Spell 5 on wards seem to peak near the end of their life time, while the spells 2-4 peak at the start of the 
###           life time (Spell 1 would also have a similar display as the latter group if not for the left-truncated spells).

# -- Save plots
dpi <- 233 # reset
ggsave(gsurv_ht_PerfSpell, file=paste0(genFigPath,"FULL SET/gsurv_ht_PerfSpell-", mainEventName,"_Hazard-KaplanMeier-SpellLevel-MultiSpell-LatentComp-incLeftTrunc_Correct.png"),
       width=1200/dpi, height=1000/dpi,dpi=dpi, bg="white")

# Housekeeping
rm(dat, vlabel, datCredit_PWPST, datCredit_TFD, Facet_lbl, gsurv_Ft_PerfSpells, gsurv_ht_PerfSpell, haz_dat, km_PerfSpells, 
   km_PerfSpells_tableSummary, smthHazard_Act1, smthHazard_Act2, smthHazard_Act3, smthHazard_Act4)

