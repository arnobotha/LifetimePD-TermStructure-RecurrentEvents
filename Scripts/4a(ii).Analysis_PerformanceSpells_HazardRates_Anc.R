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
#   - Kaplan-Meyer analysis by performance spells
# ------------------------------------------------------------------------------------------------------

# - Confirm prepared datasets are loaded into memory
if (!exists('datCredit_TFD')) unpack.ffdf(paste0(genPath,"creditdata_final_TFD"), tempPath);gc()
if (!exists('datCredit_PWPST')) unpack.ffdf(paste0(genPath,"creditdata_final_PWPST"), tempPath);gc()

# Apply performance spell grouping from script 4a(i)
# Load datasets
datCredit_TFD[,PerfSpell_Grp := fifelse(PerfSpell_Num <= 3, PerfSpell_Num, 4)]
datCredit_PWPST[,PerfSpell_Grp := fifelse(PerfSpell_Num <= 3, PerfSpell_Num, 4)]

# Sanity check
all.equal(datCredit_TFD[PerfSpell_Num < 4,PerfSpell_Num], datCredit_TFD[PerfSpell_Grp < 4,PerfSpell_Grp])# Should be true
all.equal(datCredit_PWPST[PerfSpell_Num < 4,PerfSpell_Num], datCredit_PWPST[PerfSpell_Grp < 4,PerfSpell_Grp])# Should be true

# -------- 1 Kaplan-Meier analysis on first performance spell when left-truncation is incorporated

# Initialize a dataset with only the first performance spell
dat <- datCredit_TFD[PerfSpell_Num==1,]

# -- 1.1 Cumulative Lifetime Distribution on all months

# --- Fit Kaplan-Meier (KM) nonparametric model
# Compute Kaplan-Meier survival estimates (product-limit) for main-event | Spell-level with right-censoring & left-truncation
km_TFD_All <- survfit(Surv(time=Start, time2=End, event=Default_Ind==1,type="counting") ~ 1, 
                  id=PerfSpell_Key, data=dat)
summary(km_TFD_All)$table # overall summary statistics
### RESULTS: 16306 events with a median survival time of 322 and 95% confidence interval of [255,262]
(km_TFD_All_tableSummary <- surv_summary(km_TFD_All)) # Survival table

# -- Calculate various summary statistics
# - Chosen percentile of survival times, used later as a prediction time, e.g., 75%
# NOTE: Median would also work, though sample sizes typically dwindle at that point, which can complicate analyses
survPercentile <- 0.75
cat(paste0("Kaplan-Meier: Survival time at the xth percentile of the distribution of unique survival times: \n\t",
           percTimepoint <- min(km_TFD_All$time[km_TFD_All$surv <= survPercentile],na.rm=T)))
# - Calculate the Truncated Mean (or "Restricted Mean Survival Time" [RMST]) of survival times
# NOTE: The classical mean survival time is only defined when S(t) can reach 0
# NOTE2: Calculating RMST involves integrating the KM-curve up to the last (hence truncated) observed unique event time
#         I.e., sum the survival probabilities, as weighted by the time difference between two successive survival times
cat(paste0("Kaplan-Meier: Truncated/Restricted Mean Survival Time: \n\t", 
           sum(km_TFD_All$surv*diff(c(0,km_TFD_All$time)), na.rm=T)))
# - Retrieve the mean survival probability at the previously found percentile-based prediction point
cat(paste0("Kaplan-Meier: Retrieved mean survival probability at a given time point (t=", percTimepoint, "): \n\t",
           percent(km_TFD_All$surv[max(which(km_TFD_All$time <= percTimepoint))], accuracy=0.01)))
# - Retrieve the mean survival probability at the last unique event time as prediction time, i.e., 
# the "pooled mean survival probability"
cat(paste0("Kaplan-Meier: Retrieved mean survival probability at last unique event point as time point (t=", max(km_TFD_All$time), "): \n\t",
           percent(km_TFD_All$surv[max(which(km_TFD_All$time <= max(km_TFD_All$time)))], accuracy=0.01)))


# --- Graphing survival and related quantities from fitted KM-model | F(t)

# -- Graphing parameters
vCol <- brewer.pal(8, "Set1")[c(1)] # for F(t)
medLife <- summary(km_TFD_All)$table["median"]
median_survival <- data.frame(x = medLife, y = 0.5)
mainEventName <- "Default"
chosenFont <- "Cambria"
dpi <- 190

# -- Cumulative event/lifetime probability: F(t)=1-S(t)
gsurv_Ft_TFD_All <- ggsurvplot(km_TFD_All, fun="event", conf.int=T, surv.scale = "percent", legend="none", 
                               break.time.by=round(max(km_TFD_All$time)/8), palette=vCol,
                               xlab = bquote(Discrete~time~italic(t)*" (months) in spell: First-spell"),
                               ylab = bquote(CLD~"["*.(mainEventName)*"]"*~italic(F(t))*": Kaplan-Meier (first spell incl. left-truncation)"),
                               xlim=c(0, max(km_TFD_All$time)+1), censor=F, 
                               ggtheme = theme_bw(base_family=chosenFont), tables.theme = theme_cleantable(),
                               tables.height=0.10, tables.y.text=F, tables.y.text.col=T, risk.table = "abs_pct", risk.table.pos = "out",
                               cumevents=T, cumevents.title="Cumulative number of events", 
                               cumcensor=T, cumcensor.title="Cumulative number of censored observations (incl. competing risks)",
                               risk.table.title = "Number in (% of) sample at risk of main event",
                               font.family=chosenFont, fontsize=2.5,
                               surv.median.line = "hv")
# Add median annotation for illustrative purposes
gsurv_Ft_TFD_All$plot <- gsurv_Ft_TFD_All$plot  +    
  annotate("text", x = medLife, y = 0.5,
           label = paste0("50th Percentile: ", medLife, " months"),
           vjust=-0.5, hjust=-0.1, color = "black", size = 3, family = chosenFont)
# Add Custom percentile line segment for illustrative purposes   
(gsurv_Ft_TFD_All$plot <- gsurv_Ft_TFD_All$plot  +    
    geom_segment(x = 0, xend=percTimepoint, y=1-survPercentile, yend = 1-survPercentile, linetype = "dashed", color = "black") +
    geom_segment(x = percTimepoint, xend=percTimepoint, y=0, yend = 1-survPercentile, linetype = "dashed", color = "black") +
    annotate("text", x = percTimepoint, y = 1 - survPercentile, label = paste0(comma((1-survPercentile)*100), "th Percentile: ", round(percTimepoint, 2), " months"),
             vjust=-0.5, hjust=-0.1, color = "black", size = 3, family = chosenFont))

# - Save graph and object
ggsave(gsurv_Ft_TFD_All$plot, file=paste0(genFigPath, "FULL SET/gsurv_Ft_TFD_All-", mainEventName,"_Surv-KaplanMeier-SpellLevel-MultiSpell-LatentComp-InclLeftTrunc_Correct.png"), width=1200/dpi, height=1000/dpi, dpi=dpi, bg="white")
pack.ffdf(paste0(genObjPath,"gsurv_Ft_TFD_All_RiskTable-", mainEventName,"_Surv-KaplanMeier-SpellLevel-MultiSpell-LatentComp-InclLeftTrunc_Correct.png"), km_TFD_All_tableSummary)
dpi <- 150 # need to decrease size for risk tables' text
ggsave(print(gsurv_Ft_TFD_All,newpage=F), file=paste0(genFigPath,"FULL SET/gsurv_Ft_TFD_All_RiskTable-", mainEventName,"_Surv-KaplanMeier-SpellLevel-MultiSpell-LatentComp-InclLeftTrunc_Correct.png"),width=1200/dpi, height=1000/dpi,dpi=dpi, bg="white")



# -- 1.2 Discrete baseline hazard function on first performance spell when left-truncation is incorporated
# - h(t) | Empirical estimation method
#  create plotting data object
haz_dat <- data.table(Time=km_TFD_All$time, AtRisk_n=km_TFD_All$n.risk, 
                      Event_n = km_TFD_All$n.event, Censored_n=km_TFD_All$n.censor,
                      hazard=km_TFD_All$n.event/km_TFD_All$n.risk, 
                      CumulHazard = km_TFD_All$cumhaz, #Nelson-Aalen estimator
                      Group="1",Surv_KM = km_TFD_All$surv) %>% 
  filter(Event_n > 0 | Censored_n >0) %>%
  # Discrete-time variants
  mutate(CumulHazard_Disc = -cumsum(log(1-hazard)), Surv_KM_Disc = cumprod(1-hazard)) %>% 
  mutate(Event_KM_Disc = 1-Surv_KM_Disc) %>% as.data.table()
haz_dat[, Surv_KM_Disc_prev:= shift(Surv_KM_Disc, n=1, type="lag"), by=list(Group)]
# - create alternative versions for sanity checks
haz_dat[Time==Time[1], hazard2 := 1 - Surv_KM_Disc]
haz_dat[Time>Time[1], hazard2 := 1 - Surv_KM_Disc/Surv_KM_Disc_prev]
# - conduct sanity checks
all.equal(haz_dat$hazard, haz_dat$hazard2) # Should be TRUE
### RESULTS: 8 last values are NAN, since the 8th last value has a [hazard] = 1 which results in [Surv_KM_Disc]=0
all.equal(haz_dat$Surv_KM, haz_dat$Surv_KM_Disc) # Should be TRUE
all.equal(haz_dat$CumulHazard, haz_dat$CumulHazard_Disc) # usually FALSE
plot(km_TFD_All$time, haz_dat$CumulHazard - haz_dat$CumulHazard_Disc, type="b")
### RESULTS: The discrepancy is very small difference due to estimator method differences

# --- Graphing survival and related quantities from fitted KM-model | h(t)

# -- Cubic spline for hazard rate
haz_dat$Spline <- spline_estimation(haz_dat$Time, haz_dat$hazard, 5, 3)

# -- Graphing parameters
vCol2 <- brewer.pal(10, "Paired")[c(10,9)] # for h(t)
vlabel <- c("Cubic Spline", bquote(~italic(h[0](t))))

# - Graph object for shorter time, informed by previous graphs
(gsurv_ht_TFD_All <- ggplot(haz_dat, aes(x=Time)) + theme_minimal() +
                    geom_point(aes(y=hazard, colour=vCol2[2])) + 
                    geom_line(aes(y=Spline, colour=vCol2[1]), linetype="solid") +
                    labs(y=bquote(plain(Estimated~~italic(h[0](t))*" ["*.(mainEventName)*"]"*":  Kaplan-Meier (All Months)")), 
                         x=bquote(Discrete~time~italic(t)*" (months) in spell: Single-spell")) + 
                    theme(text=element_text(family=chosenFont),legend.position="bottom") + 
                    scale_colour_manual(name="", values=vCol2, labels=vlabel) + 
                    scale_fill_manual(name="", values=vCol2, labels=vlabel) + 
                    scale_y_continuous(breaks=breaks_pretty(), label=percent) + 
                    scale_x_continuous(breaks=breaks_pretty(n=8), label=comma))
### RESULTS: The hazard appears to be near-constant over time, except after month 225 when a peak forms

# -- Save plots
dpi <- 215 # reset
ggsave(gsurv_ht_TFD_All, file=paste0(genFigPath,"FULL SET/gsurv_ht_TFD_All-", mainEventName,"_Hazard-KaplanMeier-SpellLevel-MultiSpell-LatentComp-InclLeftTrunc_Correct.png"),
       width=1200/dpi, height=1000/dpi,dpi=dpi, bg="white")

# Housekeeping
rm(dat, gsurv_Ft_TFD_All, gsurv_ht_TFD_All, km_TFD_All, km_TFD_All_tableSummary, median_survival, haz_dat, vlabel)

# Reduce months to 300
# Load datasets
datCredit_TFD <- subset(datCredit_TFD, End <= 300)
datCredit_PWPST <- subset(datCredit_PWPST, End <= 300)

# Sanity check
max(datCredit_TFD$Start) == 299
max(datCredit_PWPST$Start) == 299





# -------- 2 Kaplan-Meier analysis on first performance spell when left-truncation is incorporated

# Initialize a dataset with only the first performance spell
dat <- datCredit_TFD[PerfSpell_Num==1,]

# -- 2.1 Cumulative Lifetime Distribution on first performance spell when left-truncation is incorporated

# --- Fit Kaplan-Meier (KM) nonparametric model
# Compute Kaplan-Meier survival estimates (product-limit) for main-event | Spell-level with right-censoring & left-truncation
km_TFD <- survfit(Surv(time=Start, time2=End, event=Default_Ind==1,type="counting") ~ 1, 
                     id=PerfSpell_Key, data=dat)
summary(km_TFD)$table # overall summary statistics
### RESULTS: 16306 events with a median survival time of 258 and 95% confidence interval of [255,262]
(km_TFD_tableSummary <- surv_summary(km_TFD)) # Survival table
### RESULTS: Median survival time of 258 has standard error of 1%, which is low.

# -- Calculate various summary statistics
# - Chosen percentile of survival times, used later as a prediction time, e.g., 75%
# NOTE: Median would also work, though sample sizes typically dwindle at that point, which can complicate analyses
survPercentile <- 0.75
cat(paste0("Kaplan-Meier: Survival time at the xth percentile of the distribution of unique survival times: \n\t",
           percTimepoint <- min(km_TFD$time[km_TFD$surv <= survPercentile],na.rm=T)))
# - Calculate the Truncated Mean (or "Restricted Mean Survival Time" [RMST]) of survival times
# NOTE: The classical mean survival time is only defined when S(t) can reach 0
# NOTE2: Calculating RMST involves integrating the KM-curve up to the last (hence truncated) observed unique event time
#         I.e., sum the survival probabilities, as weighted by the time difference between two successive survival times
cat(paste0("Kaplan-Meier: Truncated/Restricted Mean Survival Time: \n\t", 
           sum(km_TFD$surv*diff(c(0,km_TFD$time)), na.rm=T)))
# - Retrieve the mean survival probability at the previously found percentile-based prediction point
cat(paste0("Kaplan-Meier: Retrieved mean survival probability at a given time point (t=", percTimepoint, "): \n\t",
           percent(km_TFD$surv[max(which(km_TFD$time <= percTimepoint))], accuracy=0.01)))
# - Retrieve the mean survival probability at the last unique event time as prediction time, i.e., 
# the "pooled mean survival probability"
cat(paste0("Kaplan-Meier: Retrieved mean survival probability at last unique event point as time point (t=", max(km_TFD$time), "): \n\t",
           percent(km_TFD$surv[max(which(km_TFD$time <= max(km_TFD$time)))], accuracy=0.01)))


# --- Graphing survival and related quantities from fitted KM-model | F(t)

# -- Graphing parameters
vCol <- brewer.pal(8, "Set1")[c(1)] # for F(t)
medLife <- summary(km_TFD)$table["median"]
median_survival <- data.frame(x = medLife, y = 0.5)
mainEventName <- "Default"
chosenFont <- "Cambria"
dpi <- 190

# -- Cumulative event/lifetime probability: F(t)=1-S(t)
gsurv_Ft_TFD <- ggsurvplot(km_TFD, fun="event", conf.int=T, surv.scale = "percent", legend="none", 
                          break.time.by=round(max(km_TFD$time)/8), palette=vCol,
                          xlab = bquote(Discrete~time~italic(t)*" (months) in spell: First-spell"),
                          ylab = bquote(CLD~"["*.(mainEventName)*"]"*~italic(F(t))*": Kaplan-Meier (incl. left-truncation)"),
                          xlim=c(0, max(km_TFD$time)+1), censor=F, 
                          ggtheme = theme_bw(base_family=chosenFont), tables.theme = theme_cleantable(),
                          tables.height=0.10, tables.y.text=F, tables.y.text.col=T, risk.table = "abs_pct", risk.table.pos = "out",
                          cumevents=T, cumevents.title="Cumulative number of events", 
                          cumcensor=T, cumcensor.title="Cumulative number of censored observations (incl. competing risks)",
                          risk.table.title = "Number in (% of) sample at risk of main event",
                          font.family=chosenFont, fontsize=2.5,
                          surv.median.line = "hv")
# Add median annotation for illustrative purposes
gsurv_Ft_TFD$plot <- gsurv_Ft_TFD$plot  +    
                      annotate("text", x = medLife, y = 0.5,
                      label = paste0("50th Percentile: ", medLife, " months"),
                      vjust=-0.5, hjust=1, color = "black", size = 3, family = chosenFont)
# Add Custom percentile line segment for illustrative purposes   
(gsurv_Ft_TFD$plot <- gsurv_Ft_TFD$plot  +    
                      geom_segment(x = 0, xend=percTimepoint, y=1-survPercentile, yend = 1-survPercentile, linetype = "dashed", color = "black") +
                      geom_segment(x = percTimepoint, xend=percTimepoint, y=0, yend = 1-survPercentile, linetype = "dashed", color = "black") +
                      annotate("text", x = percTimepoint, y = 1 - survPercentile, label = paste0(comma((1-survPercentile)*100), "th Percentile: ", round(percTimepoint, 2), " months"),
                               vjust=-0.5, hjust=1, color = "black", size = 3, family = chosenFont))

# - Save graph and object
ggsave(gsurv_Ft_TFD$plot, file=paste0(genFigPath, "FULL SET/gsurv_Ft_TFD-", mainEventName,"_Surv-KaplanMeier-SpellLevel-MultiSpell-LatentComp-InclLeftTrunc_Correct.png"), width=1200/dpi, height=1000/dpi, dpi=dpi, bg="white")
pack.ffdf(paste0(genObjPath,"gsurv_Ft_TFD_RiskTable-", mainEventName,"_Surv-KaplanMeier-SpellLevel-MultiSpell-LatentComp-InclLeftTrunc_Correct.png"), km_TFD_tableSummary)
dpi <- 150 # need to decrease size for risk tables' text (152 is the cut-off value for the table annotation)
ggsave(print(gsurv_Ft_TFD,newpage=F), file=paste0(genFigPath,"FULL SET/gsurv_Ft_TFD_RiskTable-", mainEventName,"_Surv-KaplanMeier-SpellLevel-MultiSpell-LatentComp-InclLeftTrunc_Correct.png"),width=1200/dpi, height=1000/dpi,dpi=dpi, bg="white")



# -- 2.2 Discrete baseline hazard function on first performance spell when left-truncation is incorporated
# - h(t) | Empirical estimation method
#  create plotting data object
haz_dat <- data.table(Time=km_TFD$time, AtRisk_n=km_TFD$n.risk, 
                      Event_n = km_TFD$n.event, Censored_n=km_TFD$n.censor,
                      hazard=km_TFD$n.event/km_TFD$n.risk, 
                      CumulHazard = km_TFD$cumhaz, #Nelson-Aalen estimator
                      Group="1",Surv_KM = km_TFD$surv) %>% 
  filter(Event_n > 0 | Censored_n >0) %>%
  # Discrete-time variants
  mutate(CumulHazard_Disc = -cumsum(log(1-hazard)), Surv_KM_Disc = cumprod(1-hazard)) %>% 
  mutate(Event_KM_Disc = 1-Surv_KM_Disc) %>% as.data.table()
haz_dat[, Surv_KM_Disc_prev:= shift(Surv_KM_Disc, n=1, type="lag"), by=list(Group)]
# - create alternative versions for sanity checks
haz_dat[Time==Time[1], hazard2 := 1 - Surv_KM_Disc]
haz_dat[Time>Time[1], hazard2 := 1 - Surv_KM_Disc/Surv_KM_Disc_prev]
# - conduct sanity checks
all.equal(haz_dat$hazard, haz_dat$hazard2) # Should be TRUE
### RESULTS: 8 last values are NAN, since the 8th last value has a [hazard] = 1 which results in [Surv_KM_Disc]=0
all.equal(haz_dat$Surv_KM, haz_dat$Surv_KM_Disc) # Should be TRUE
all.equal(haz_dat$CumulHazard, haz_dat$CumulHazard_Disc) # usually FALSE
plot(km_TFD$time, haz_dat$CumulHazard - haz_dat$CumulHazard_Disc, type="b")
### RESULTS: The discrepancy is very small difference due to estimator method differences

# --- Graphing survival and related quantities from fitted KM-model | h(t)

# -- Cubic spline for hazard rate
# HW: Pack Spline haz_date
haz_dat$Spline <- spline_estimation(haz_dat$Time, haz_dat$hazard, 10, 3)

# -- Graphing parameters
vCol2 <- brewer.pal(10, "Paired")[c(10,9)] # for h(t)
vlabel <- c("Cubic Spline", bquote(~italic(h[0](t))))

# - Graph object for shorter time, informed by previous graphs
(gsurv_ht_TFD <- ggplot(haz_dat, aes(x=Time)) + theme_minimal() +
                  geom_point(aes(y=hazard, colour=vCol2[2])) + 
                  geom_line(aes(y=Spline, colour=vCol2[1]), linetype="solid") +
                  labs(y=bquote(plain(Estimated~~italic(h[0](t))*" ["*.(mainEventName)*"]"*": Kaplan-Meier (incl. left-trunc.)")), 
                       x=bquote(Discrete~time~italic(t)*" (months) in spell: Single-spell")) + 
                  theme(text=element_text(family=chosenFont),legend.position="bottom") + 
                  scale_colour_manual(name="", values=vCol2, labels=vlabel) + 
                  scale_fill_manual(name="", values=vCol2, labels=vlabel) + 
                  scale_y_continuous(breaks=breaks_pretty(), label=percent) + 
                  scale_x_continuous(breaks=breaks_pretty(n=8), label=comma))
### RESULTS: The hazard appears to be near-constant over time, except after month 225 when a peak forms

# -- Save plots
dpi <- 215 # reset
ggsave(gsurv_ht_TFD, file=paste0(genFigPath,"FULL SET/gsurv_ht_TFD-", mainEventName,"_Hazard-KaplanMeier-SpellLevel-MultiSpell-LatentComp-InclLeftTrunc_Correct.png"),
       width=1200/dpi, height=1000/dpi,dpi=dpi, bg="white")

# Housekeeping
rm(dat, gsurv_Ft_TFD, gsurv_ht_TFD, km_TFD, km_TFD_tableSummary, median_survival, haz_dat, vlabel)


# -- 1.3 Cumulative Lifetime Distribution function on first performance spell when left-truncation is not incorporated

# Initialize a dataset with only the first performance spell
dat <- datCredit_TFD[PerfSpell_Num==1 & PerfSpell_LeftTrunc==0,]

# --- Fit Kaplan-Meier (KM) nonparametric model
# Compute Kaplan-Meier survival estimates (product-limit) for main-event | Spell-level with right-censoring & left-truncation
km_TFD_noLeftTrunc <- survfit(Surv(time=Start, time2=End, event=Default_Ind==1,type="counting") ~ 1, 
                  id=PerfSpell_Key, data=dat)
summary(km_TFD_noLeftTrunc)$table # overall summary statistics
### RESULTS: 7489 events with no median survival time
(km_TFD_noLeftTrunc_tableSummary <- surv_summary(km_TFD_noLeftTrunc)) # Survival table

# -- Calculate various summary statistics
# - Chosen percentile of survival times, used later as a prediction time, e.g., 75%
# NOTE: Median would also work, though sample sizes typically dwindle at that point, which can complicate analyses
survPercentile <- 0.75
cat(paste0("Kaplan-Meier: Survival time at the xth percentile of the distribution of unique survival times: \n\t",
           percTimepoint <- min(km_TFD_noLeftTrunc$time[km_TFD_noLeftTrunc$surv <= survPercentile],na.rm=T)))
# - Calculate the Truncated Mean (or "Restricted Mean Survival Time" [RMST]) of survival times
# NOTE: The classical mean survival time is only defined when S(t) can reach 0
# NOTE2: Calculating RMST involves integrating the KM-curve up to the last (hence truncated) observed unique event time
#         I.e., sum the survival probabilities, as weighted by the time difference between two successive survival times
cat(paste0("Kaplan-Meier: Truncated/Restricted Mean Survival Time: \n\t", 
           sum(km_TFD_noLeftTrunc$surv*diff(c(0,km_TFD_noLeftTrunc$time)), na.rm=T)))
# - Retrieve the mean survival probability at the previously found percentile-based prediction point
cat(paste0("Kaplan-Meier: Retrieved mean survival probability at a given time point (t=", percTimepoint, "): \n\t",
           percent(km_TFD_noLeftTrunc$surv[max(which(km_TFD_noLeftTrunc$time <= percTimepoint))], accuracy=0.01)))
# - Retrieve the mean survival probability at the last unique event time as prediction time, i.e., 
# the "pooled mean survival probability"
cat(paste0("Kaplan-Meier: Retrieved mean survival probability at last unique event point as time point (t=", max(km_TFD_noLeftTrunc$time), "): \n\t",
           percent(km_TFD_noLeftTrunc$surv[max(which(km_TFD_noLeftTrunc$time <= max(km_TFD_noLeftTrunc$time)))], accuracy=0.01)))


# --- Graphing survival and related quantities from fitted KM-model | F(t)

# -- Graphing parameters
vCol <- brewer.pal(8, "Set1")[c(1)] # for F(t)
maxProb <- min(km_TFD_noLeftTrunc_tableSummary$surv)
mainEventName <- "Default"
chosenFont <- "Cambria"
dpi <- 170

# -- Cumulative event/lifetime probability: F(t)=1-S(t)
gsurv_Ft_TFD_noLeftTrunc <- ggsurvplot(km_TFD_noLeftTrunc, fun="event", conf.int=T, surv.scale = "percent", legend="none", 
                             break.time.by=round(max(km_TFD_noLeftTrunc$time)/8), palette=vCol,
                             xlab = bquote(Discrete~time~italic(t)*" (months) in spell: First-spell"),
                             ylab = bquote(CLD~"["*.(mainEventName)*"]"*~italic(F(t))*": Kaplan-Meier (first spell exc. left-truncation)"),
                             xlim=c(0, max(km_TFD_noLeftTrunc$time)+1), censor=F, 
                             ggtheme = theme_bw(base_family=chosenFont), tables.theme = theme_cleantable(),
                             tables.height=0.10, tables.y.text=F, tables.y.text.col=T, risk.table = "abs_pct", risk.table.pos = "out",
                             cumevents=T, cumevents.title="Cumulative number of events", 
                             cumcensor=T, cumcensor.title="Cumulative number of censored observations (incl. competing risks)",
                             risk.table.title = "Number in (% of) sample at risk of main event",
                             font.family=chosenFont, fontsize=2.5)
# Add Custom percentile line segment for illustrative purposes   
(gsurv_Ft_TFD_noLeftTrunc$plot <- gsurv_Ft_TFD_noLeftTrunc$plot  +    
                                  geom_segment(x = 0, xend=192, y=1-maxProb, yend = 1-maxProb, linetype = "dashed", color = "black") +
                                  geom_segment(x = 192, xend=192, y=0, yend = 1-maxProb, linetype = "dashed", color = "black") +
                                  annotate("text", x = 192, y = 1 - maxProb, label = paste0(comma((1-maxProb)*100), "th Percentile: 192 months"),
                                           vjust=-0.5, hjust=1, color = "black", size = 3, family = chosenFont))

# - Save graph and object
ggsave(gsurv_Ft_TFD_noLeftTrunc$plot, file=paste0(genFigPath, "FULL SET/gsurv_Ft_TFD_noLeftTrunc-", mainEventName,"_Surv-KaplanMeier-SpellLevel-MultiSpell-LatentComp-ExclLeftTrunc_Correct.png"), width=1200/dpi, height=1000/dpi, dpi=dpi, bg="white")
pack.ffdf(paste0(genObjPath,"gsurv_Ft_TFD_noLeftTrunc_RiskTable-", mainEventName,"_Surv-KaplanMeier-SpellLevel-MultiSpell-LatentComp-InclLeftTrunc_Correct"), km_TFD_noLeftTrunc_tableSummary)
dpi <- 150 # need to decrease size for risk tables' text
ggsave(print(gsurv_Ft_TFD_noLeftTrunc,newpage=F), file=paste0(genFigPath,"FULL SET/gsurv_Ft_TFD_noLeftTrunc_RiskTable-", mainEventName,"_Surv-KaplanMeier-SpellLevel-MultiSpell-LatentComp-InclLeftTrunc_Correct.png"),width=1200/dpi, height=1000/dpi,dpi=dpi, bg="white")



# -- 1.4 Discrete baseline hazard function on first performance spell when left-truncation is not incorporated
# - h(t) | Empirical estimation method
# create plotting data object
haz_dat <- data.table(Time=km_TFD_noLeftTrunc$time, AtRisk_n=km_TFD_noLeftTrunc$n.risk, 
                      Event_n = km_TFD_noLeftTrunc$n.event, Censored_n=km_TFD_noLeftTrunc$n.censor,
                      hazard=km_TFD_noLeftTrunc$n.event/km_TFD_noLeftTrunc$n.risk, 
                      CumulHazard = km_TFD_noLeftTrunc$cumhaz, #Nelson-Aalen estimator
                      Group="1",Surv_KM = km_TFD_noLeftTrunc$surv) %>% 
  filter(Event_n > 0 | Censored_n >0) %>%
  # Discrete-time variants
  mutate(CumulHazard_Disc = -cumsum(log(1-hazard)), Surv_KM_Disc = cumprod(1-hazard)) %>% 
  mutate(Event_KM_Disc = 1-Surv_KM_Disc) %>% as.data.table()
haz_dat[, Surv_KM_Disc_prev:= shift(Surv_KM_Disc, n=1, type="lag"), by=list(Group)]
# - create alternative versions for sanity checks
haz_dat[Time==Time[1], hazard2 := 1 - Surv_KM_Disc]
haz_dat[Time>Time[1], hazard2 := 1 - Surv_KM_Disc/Surv_KM_Disc_prev]
# - conduct sanity checks
all.equal(haz_dat$hazard, haz_dat$hazard2) # Should be TRUE
all.equal(haz_dat$Surv_KM, haz_dat$Surv_KM_Disc) # Should be TRUE
all.equal(haz_dat$CumulHazard, haz_dat$CumulHazard_Disc) # usually FALSE
plot(km_TFD_noLeftTrunc$time, haz_dat$CumulHazard - haz_dat$CumulHazard_Disc, type="b")
### RESULTS: The discrepancy is very small difference due to estimator method differences

# --- Graphing survival and related quantities from fitted KM-model | h(t)

# -- Cubic spline for hazard rate
haz_dat$Spline <- spline_estimation(haz_dat$Time, haz_dat$hazard, 10, 3)

# -- Graphing parameters
vCol2 <- brewer.pal(10, "Paired")[c(10,9)] # for h(t)
vlabel <- c("Cubic Spline", bquote(~italic(h[0](t))))

# - Graph object for shorter time, informed by previous graphs
(gsurv_ht_TFD_noLeftTrunc <- ggplot(haz_dat, aes(x=Time)) + theme_minimal() +
                            geom_point(aes(y=hazard, colour=vCol2[2])) + 
                            geom_line(aes(y=Spline, colour=vCol2[1]), linetype="solid") +
                            labs(y=bquote(plain(Estimated~italic(h[0](t))*"["*.(mainEventName)*"]: Kaplan-Meier (excl. left-trunc.)")), 
                                 x=bquote(Discrete~time~italic(t)*" (months) in spell: First spell")) + 
                            theme(text=element_text(family=chosenFont),legend.position="bottom") + 
                            scale_colour_manual(name="", values=vCol2, labels=vlabel) + 
                            scale_fill_manual(name="", values=vCol2, labels=vlabel) + 
                            scale_y_continuous(breaks=breaks_pretty(), label=percent) + 
                            scale_x_continuous(breaks=breaks_pretty(n=8), label=comma))
### RESULTS: The hazard appears to be near-constant over time, with some notable oscillation over some prediction periods.
# However, when viewed in tandem with S(t), itself almost a straight downward-sloping line, it makes sense for hazard
# to be near-flat. The oscillation also seems more pronounced towards later prediction periods than earlier ones.

# -- Save plots
dpi <- 220 # reset
ggsave(gsurv_ht_TFD_noLeftTrunc, file=paste0(genFigPath,"FULL SET/gsurv_ht_TFD_noLeftTrunc-", mainEventName,"_Hazard-KaplanMeier-SpellLevel-MultiSpell-LatentComp-excLeftTrunc_Correct.png"),
       width=1200/dpi, height=1000/dpi,dpi=dpi, bg="white")

# Housekeeping
rm(dat,gsurv_ht_TFD_noLeftTrunc, km_TFD_noLeftTrunc, km_TFD_noLeftTrunc_tableSummary, haz_dat)




# -- 1.5 Cumulative Lifetime Distribution function on different performance spells
# Initialize data by one loan in the 9th performance spell that was censored
dat <- datCredit_PWPST[PerfSpell_Num !=9,]

# --- Fit Cox-model based on [PerfSpell_Num]
km_PerfSpells <- survfit(Surv(time=Start, time2=End, event=Default_Ind==1,type="counting") ~ PerfSpell_Grp, 
                              id=PerfSpell_Key, data=dat)
summary(km_PerfSpells)$table # overall summary statistics
### RESULTS: The 7th and 8th performance spell have no upper limit for their respective median 95th confidence interval
(km_PerfSpells_tableSummary <- surv_summary(km_PerfSpells)) # Survival table

# --- Graphing survival and related quantities from fitted KM-model | F(t)

# -- Graphing parameters
vCol <- brewer.pal(8, "Dark2")[1:4] # for F(t)

# Create facet groupings for performance spells
dat[, PerfSpell_Label := fifelse(PerfSpell_Num <= 3, "Perf. Spells 1-3",
                                       fifelse(PerfSpell_Num <= 5,
                                               "Perf. Spells 4-5", "Perf. Spells 6-8"))]
mainEventName <- "Default"
chosenFont <- "Cambria"
dpi <- 150
# -- Cumulative event/lifetime probability: F(t)=1-S(t)
(gsurv_Ft_PerfSpells <- ggsurvplot(km_PerfSpells, data=dat, fun="event", conf.int=T, surv.scale="percent",
                                palette=vCol, xlab = bquote(Discrete~time~italic(t)*" (months) in spells"),
                                ylab = bquote(CLD~"["*.(mainEventName)*"]"*~italic(F(t))*": Kaplan-Meier (By Performance Spells)"),
                                legend.title="Performance Spells", legend="bottom",
                                legend.labs=c("Spell 1", "Spell 2","Spell 3", "Spell 4+"),
                                censor=F, ggtheme = theme_bw(base_family=chosenFont)))

# - Save graph and object
ggsave(print(gsurv_Ft_PerfSpells,newpage=F), file=paste0(genFigPath, "FULL SET/gsurv_Ft_PerfSpells-", mainEventName,"_Surv-KaplanMeier-SpellLevel-MultiSpell-LatentComp_Correct.png"), width=1200/dpi, height=1000/dpi, dpi=dpi, bg="white")
pack.ffdf(paste0(genObjPath,"KM_PerfSpells_RiskTable"), km_PerfSpells_tableSummary)


# -- 1.4 Discrete baseline hazard function on different performance spells
# - h(t) | Empirical estimation method
# create plotting data object
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

# -- Graphing parameters
vCol2_Line <- brewer.pal(8, "Dark2")[c(1:4)] # for h(t) lines
vCol2_Point <- brewer.pal(8, "Pastel2")[c(1:4)] # for h(t) points
# Create facet groupings for performance spells
haz_dat[, GroupLabel := factor(Group, levels = 1:4,
                               labels = c(
                                 bquote(~italic(h[0](t))*": Performance Spell 1"),
                                 bquote(~italic(h[0](t))*": Performance Spell 2"),
                                 bquote(~italic(h[0](t))*": Performance Spell 3"),
                                 bquote(~italic(h[0](t))*": Performance Spell 4+")))]

# -- Cubic spline for hazard rate
Spline_1 <- spline_estimation(haz_dat[Group==1,Time], haz_dat[Group==1,hazard], 10, 3)
Spline_2 <- spline_estimation(haz_dat[Group==2,Time], haz_dat[Group==2,hazard], 10, 3)
Spline_3 <- spline_estimation(haz_dat[Group==3,Time], haz_dat[Group==3,hazard], 10, 3)
Spline_4 <- spline_estimation(haz_dat[Group==4,Time], haz_dat[Group==4,hazard], 10, 3)
haz_dat[Group==1, `:=` (Spline=Spline_1, colLine=vCol2_Line[1], colPoint=vCol2_Point[1])]
haz_dat[Group==2, `:=` (Spline=Spline_2, colLine=vCol2_Line[2], colPoint=vCol2_Point[2])]
haz_dat[Group==3, `:=` (Spline=Spline_3, colLine=vCol2_Line[3], colPoint=vCol2_Point[3])]
haz_dat[Group==4, `:=` (Spline=Spline_4, colLine=vCol2_Line[4], colPoint=vCol2_Point[4])]

# - Graph object for shorter time, informed by previous graphs
(gsurv_ht_PerfSpell <- ggplot(haz_dat, aes(x=Time)) + theme_minimal() +
                              geom_point(aes(y=hazard, colour=colPoint)) + 
                              geom_line(aes(y=Spline, colour=colLine), linetype="solid") +
                              labs(y=bquote(plain(Estimated~hazard*" ["*.(mainEventName)*"]"*~italic(h[0](t))*":  Kaplan-Meier (spell-level)")), 
                                   x=bquote(Discrete~time~italic(t)*" (months) in each spell: Single-spell")) + 
                              theme(text=element_text(family=chosenFont),legend.position="bottom",
                                    strip.background=element_rect(fill="snow2", colour="snow2"),
                                    strip.text=element_text(size=8, colour="gray50")) +
                              facet_wrap(~GroupLabel, scales="free", labeller=label_parsed) + 
                              scale_colour_identity() +
                              scale_y_continuous(breaks=breaks_pretty(), label=percent) + 
                              scale_x_continuous(breaks=breaks_pretty(n=8), label=comma))
### RESULTS:  As the performance spell increase, the volatility in the hazard rate increases accordingly. Note that
###           Spell 5 on wards seem to peak near the end of their life time, while the spells 2-4 peak at the start of the 
###           life time (Spell 1 would also have a similar display as the latter group if not for the left-truncated spells).

# -- Save plots
dpi <- 233 # reset
ggsave(gsurv_ht_PerfSpell, file=paste0(genFigPath,"FULL SET/gsurv_ht_PerfSpell-", mainEventName,"_Hazard-KaplanMeier-SpellLevel-MultiSpell-LatentComp-incLeftTrunc_Correct.png"),
       width=1200/dpi, height=1000/dpi,dpi=dpi, bg="white")

# Housekeeping
rm(dat, vlabel, datCredit_PWPST, datCredit_TFD, Facet_lbl, gsurv_Ft_PerfSpells, gsurv_ht_PerfSpell, haz_dat, km_PerfSpells, km_PerfSpells_tableSummary)

# Save Datasets with new changes
# Load datasets
if (!exists('datCredit_train_TFD')) unpack.ffdf(paste0(genPath,"creditdata_train_TFD"), tempPath);gc()
if (!exists('datCredit_valid_TFD')) unpack.ffdf(paste0(genPath,"creditdata_valid_TFD"), tempPath);gc()
if (!exists('datCredit_valid_PWPST')) unpack.ffdf(paste0(genPath,"creditdata_valid_PWPST"), tempPath);gc()
if (!exists('datCredit_train_PWPST')) unpack.ffdf(paste0(genPath,"creditdata_train_PWPST"), tempPath);gc()

# Reduce months to 300
# Load datasets
datCredit_train_TFD <- subset(datCredit_train_TFD, End <= 300)
datCredit_valid_TFD <- subset(datCredit_valid_TFD, End <= 300)
datCredit_train_PWPST <- subset(datCredit_train_PWPST, End <= 300)
datCredit_valid_PWPST <- subset(datCredit_valid_PWPST, End <= 300)

# Load datasets
datCredit_train_TFD[,PerfSpell_Grp := fifelse(PerfSpell_Num <= 3, PerfSpell_Num, 4)]
datCredit_valid_TFD[,PerfSpell_Grp := fifelse(PerfSpell_Num <= 3, PerfSpell_Num, 4)]
datCredit_train_PWPST[,PerfSpell_Grp := fifelse(PerfSpell_Num <= 3, PerfSpell_Num, 4)]
datCredit_valid_PWPST[,PerfSpell_Grp := fifelse(PerfSpell_Num <= 3, PerfSpell_Num, 4)]

# --- 5.2 Saving the dataset amendments scheme
# - Training dataset
pack.ffdf(paste0(genPath,"creditdata_train_TFD"), datCredit_train_TFD)
pack.ffdf(paste0(genPath,"creditdata_train_PWPST"), datCredit_train_PWPST)

# - Validation dataset
pack.ffdf(paste0(genPath,"creditdata_valid_TFD"), datCredit_valid_TFD)
pack.ffdf(paste0(genPath,"creditdata_valid_PWPST"), datCredit_valid_PWPST)

# Housekeeping
rm(datCredit_train_PWPST, datCredit_train_TFD, datCredit_valid_PWPST, datCredit_valid_TFD)
