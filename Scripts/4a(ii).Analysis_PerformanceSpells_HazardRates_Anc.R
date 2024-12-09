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

# -------- 1 Kaplan-Meier analysis on first performance spell when left-truncation is incorporated

# Initialize a dataset with only the first performance spell
dat <- datCredit_TFD[PerfSpell_Num==1,]

# -- 1.1 Cumulative Lifetime Distribution on first performance spell when left-truncation is incorporated

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
dpi <- 180

# -- Cumulative event/lifetime probability: F(t)=1-S(t)
gsurv_Ft_TFD <- ggsurvplot(km_TFD, fun="event", conf.int=T, surv.scale = "percent", legend="none", 
                          break.time.by=round(max(km_TFD$time)/8), palette=vCol,
                          xlab = bquote(Discrete~time~italic(t)*" (months) in spell: first-spell"),
                          ylab = bquote(CLD~"["*.(mainEventName)*"]"*~italic(F(t))*": Kaplan-Meier (first spell incl. left-truncation)"),
                          xlim=c(0, max(km_TFD$time)+1), censor=F, 
                          ggtheme = theme_bw(base_family=chosenFont), tables.theme = theme_cleantable(),
                          tables.height=0.10, tables.y.text=F, tables.y.text.col=T, risk.table = "abs_pct", risk.table.pos = "out",
                          cumevents=T, cumevents.title="Cumulative number of events", 
                          cumcensor=T, cumcensor.title="Cumulative number of censored observations (incl. competing risks)",
                          risk.table.title = "Number in (% of) sample at risk of main event",
                          font.family=chosenFont, fontsize=2.5,
                          surv.median.line = "hv")
# Add median annotation for illustrative purposes
(gsurv_Ft_TFD$plot <- gsurv_Ft_TFD$plot  +    
                      annotate("text", x = medLife, y = 0.5,
                      label = paste0("50th Percentile: ", medLife, " months"),
                      hjust = -0.1, color = "black", size = 3, family = chosenFont))
# Add Custom percentile line segment for illustrative purposes   
(gsurv_Ft_TFD$plot <- gsurv_Ft_TFD$plot  +    
                      geom_segment(x = 0, xend=percTimepoint, y=1-survPercentile, yend = 1-survPercentile, linetype = "dashed", color = "black") +
                      geom_segment(x = percTimepoint, xend=percTimepoint, y=0, yend = 1-survPercentile, linetype = "dashed", color = "black") +
                      annotate("text", x = percTimepoint, y = 1 - survPercentile, label = paste0(comma((1-survPercentile)*100), "th Percentile: ", round(percTimepoint, 2), " months"),
                                hjust = -0.1, color = "black", size = 3, family = chosenFont))

# - Save graph and object
ggsave(gsurv_Ft_TFD$plot, file=paste0(genFigPath, "FULL SET/gsurv_Ft_TFD-", mainEventName,"_Surv-KaplanMeier-SpellLevel-MultiSpell-LatentComp-InclLeftTrunc_Correct.png"), width=1200/dpi, height=1000/dpi, dpi=dpi, bg="white")
pack.ffdf(paste0(genObjPath,"gsurv_Ft_TFD_RiskTable-", mainEventName,"_Surv-KaplanMeier-SpellLevel-MultiSpell-LatentComp-InclLeftTrunc_Correct.png"), km_TFD_tableSummary)
dpi <- 150 # need to decrease size for risk tables' text
ggsave(print(gsurv_Ft_TFD,newpage=F), file=paste0(genFigPath,"FULL SET/gsurv_Ft_TFD_RiskTable-", mainEventName,"_Surv-KaplanMeier-SpellLevel-MultiSpell-LatentComp-InclLeftTrunc_Correct.png"),width=1200/dpi, height=1000/dpi,dpi=dpi, bg="white")

# -- 1.2 Discrete baseline hazard function on first performance spell when left-truncation is incorporated
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

# -- Graphing parameters
vCol2 <- brewer.pal(10, "Paired")[c(10,9)] # for h(t)
sSpan <- 0.1; # span for LOESS-smoother in h(t)
vlabel <- paste0("Loess-smoothed hazard [span: ", sSpan, "]") # for h(t)

# - Graph object for shorter time, informed by previous graphs
(gsurv_ht_TFD <- ggplot(haz_dat[Time<=300,], aes(x=Time,y=hazard)) + theme_minimal() +
                  geom_line(linetype="solid", colour=vCol2[1]) + geom_point(colour=vCol2[1]) + 
                  geom_smooth(aes(colour=Group, fill=Group), se=T, method="loess",
                              span=sSpan, alpha=0.25, linetype="dotted") +
                  labs(y=bquote(plain(Estimated~hazard*" ["*.(mainEventName)*"]"*~italic(h(t))*":  Kaplan-Meier (first spell incl. left-truncation)")), 
                       x=bquote(Discrete~time~italic(t)*" (months) in spell: Multi-spell")) + 
                  theme(text=element_text(family=chosenFont),legend.position="bottom") + 
                  scale_colour_manual(name="", values=vCol2[2], labels=vlabel) + 
                  scale_fill_manual(name="", values=vCol2[2], labels=vlabel) + 
                  scale_y_continuous(breaks=breaks_pretty(), label=percent) + 
                  scale_x_continuous(breaks=breaks_pretty(n=8), label=comma))
### RESULTS: The hazard appears to be near-constant over time, except after month 225 when a peak forms

# -- Save plots
dpi <- 180 # reset
ggsave(gsurv_ht_TFD, file=paste0(genFigPath,"FULL SET/gsurv_ht_TFD-", mainEventName,"_Hazard-KaplanMeier-SpellLevel-MultiSpell-LatentComp-InclLeftTrunc_Correct.png"),
       width=1200/dpi, height=1000/dpi,dpi=dpi, bg="white")

# Housekeeping
rm(dat, gsurv_Ft_TFD, gsurv_ht_TFD, km_TFD, km_TFD_tableSummary, median_survival, haz_dat)


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
dpi <- 180

# -- Cumulative event/lifetime probability: F(t)=1-S(t)
gsurv_Ft_TFD_noLeftTrunc <- ggsurvplot(km_TFD_noLeftTrunc, fun="event", conf.int=T, surv.scale = "percent", legend="none", 
                             break.time.by=round(max(km_TFD_noLeftTrunc$time)/8), palette=vCol,
                             xlab = bquote(Discrete~time~italic(t)*" (months) in spell: first-spell"),
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
pack.ffdf(paste0(genObjPath,"gsurv_Ft_TFD_noLeftTrunc_RiskTable-", mainEventName,"_Surv-KaplanMeier-SpellLevel-MultiSpell-LatentComp-InclLeftTrunc_Correct.png"), km_TFD_noLeftTrunc_tableSummary)
dpi <- 150 # need to decrease size for risk tables' text
ggsave(print(gsurv_Ft_TFD_noLeftTrunc,newpage=F), file=paste0(genFigPath,"FULL SET/gsurv_Ft_TFD_RiskTable-", mainEventName,"_Surv-KaplanMeier-SpellLevel-MultiSpell-LatentComp-InclLeftTrunc_Correct.png"),width=1200/dpi, height=1000/dpi,dpi=dpi, bg="white")

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

# -- Graphing parameters
vCol2 <- brewer.pal(10, "Paired")[c(10,9)] # for h(t)
sSpan <- 0.1; # span for LOESS-smoother in h(t)
vlabel <- paste0("Loess-smoothed hazard [span: ", sSpan, "]") # for h(t)

# - Graph object for shorter time, informed by previous graphs
(gsurv_ht_TFD_noLeftTrunc <- ggplot(haz_dat[Time<=300,], aes(x=Time,y=hazard)) + theme_minimal() +
                            geom_line(linetype="solid", colour=vCol2[1]) + geom_point(colour=vCol2[1]) + 
                            geom_smooth(aes(colour=Group, fill=Group), se=T, method="loess",
                                        span=sSpan, alpha=0.25, linetype="dotted") +
                            labs(y=bquote(plain(Estimated~hazard*" ["*.(mainEventName)*"]"*~italic(h(t))*":  Kaplan-Meier (first spell incl. left-truncation)")), 
                                 x=bquote(Discrete~time~italic(t)*" (months) in spell: Multi-spell")) + 
                            theme(text=element_text(family=chosenFont),legend.position="bottom") + 
                            scale_colour_manual(name="", values=vCol2[2], labels=vlabel) + 
                            scale_fill_manual(name="", values=vCol2[2], labels=vlabel) + 
                            scale_y_continuous(breaks=breaks_pretty(), label=percent) + 
                            scale_x_continuous(breaks=breaks_pretty(n=8), label=comma))
### RESULTS: The hazard appears to be near-constant over time, with some notable oscillation over some prediction periods.
# However, when viewed in tandem with S(t), itself almost a straight downward-sloping line, it makes sense for hazard
# to be near-flat. The oscillation also seems more pronounced towards later prediction periods than earlier ones.

# -- Save plots
dpi <- 180 # reset
ggsave(gsurv_ht_TFD_noLeftTrunc, file=paste0(genFigPath,"FULL SET/gsurv_ht_TFD_noLeftTrunc-", mainEventName,"_Hazard-KaplanMeier-SpellLevel-MultiSpell-LatentComp-excLeftTrunc_Correct.png"),
       width=1200/dpi, height=1000/dpi,dpi=dpi, bg="white")

# Housekeeping
rm(dat,gsurv_ht_TFD_noLeftTrunc, km_TFD_noLeftTrunc, km_TFD_noLeftTrunc_tableSummary, haz_dat)




# -- 1.4 Cumulative Lifetime Distribution function on different performance spells
# Initialize data by one loan in the 9th performance spell that was censored
dat <- datCredit_PWPST[PerfSpell_Num !=9,]

# --- Fit Cox-model based on [PerfSpell_Num]
km_PerfSpells <- survfit(Surv(time=Start, time2=End, event=Default_Ind==1,type="counting") ~ PerfSpell_Num, 
                              id=PerfSpell_Key, data=dat)
summary(km_PerfSpells)$table # overall summary statistics
### RESULTS: The 7th and 8th performance spell have no upper limit for their respective median 95th confidence interval
(km_PerfSpells_tableSummary <- surv_summary(km_PerfSpells)) # Survival table

# --- Graphing survival and related quantities from fitted KM-model | F(t)

# -- Graphing parameters
vCol <- brewer.pal(8, "Dark2") # for F(t)

# Create facet groupings for performance spells
dat[, PerfSpell_grp := fifelse(PerfSpell_Num==1,"Spells 1",
                               fifelse(PerfSpell_Num <= 3, "Spells 2-3",
                                       fifelse(PerfSpell_Num <= 5,
                                               "Spells 4-5", "Spells 6-8")))]
mainEventName <- "Default"
chosenFont <- "Cambria"
dpi <- 150
# -- Cumulative event/lifetime probability: F(t)=1-S(t)
(gsurv_Ft_PerfSpells <- ggsurvplot(km_PerfSpells, data=dat, fun="event", conf.int=T, surv.scale="percent",
                                palette=vCol, xlab = bquote(Discrete~time~italic(t)*" (months) in spells"),
                                ylab = bquote(CLD~"["*.(mainEventName)*"]"*~italic(F(t))*": Kaplan-Meier (By Performance Spells)"),
                                legend.title="Performance Spells", legend="bottom",
                                legend.labs=c("Spell 1", "Spell 2","Spell 3", "Spell 4",
                                              "Spell 5","Spell 6","Spell 7", "Spell 8"),
                                censor=F, ggtheme = theme_bw(base_family=chosenFont),
                                facet.by="PerfSpell_grp", short.panel.labs=T,scales="free"))

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
                      Group=km_PerfSpells_tableSummary$PerfSpell_Num,Surv_KM = km_PerfSpells$surv) %>% 
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
vCol2_Line <- brewer.pal(8, "Dark2") # for h(t) lines
vCol2_CI <- brewer.pal(8, "Pastel2") # for h(t) confidence intervals
sSpan <- 0.1; # span for LOESS-smoother in h(t)
# Create facet groupings for performance spells
Facet_grp <- paste0("Perf. Spell ", 1:8)
Facet_lbl <- c("1" = Facet_grp[1], "2" = Facet_grp[2], "3" = Facet_grp[3],
                    "4" = Facet_grp[4], "5" = Facet_grp[5], "6" = Facet_grp[6],
                    "7" = Facet_grp[7], "8" = Facet_grp[8])

# - Graph object for shorter time, informed by previous graphs
(gsurv_ht_PerfSpell <- ggplot(haz_dat[Time<=300,], aes(x=Time,y=hazard, group=Group, color=Group)) + theme_minimal() +
                              geom_line(linetype="solid") + geom_point() + 
                              geom_smooth(se=T, method="loess", span=sSpan, alpha=0.25, linetype="dotted") +
                              labs(y=bquote(plain(Estimated~hazard*" ["*.(mainEventName)*"]"*~italic(h(t))*":  Kaplan-Meier (spell-level)")), 
                                   x=bquote(Discrete~time~italic(t)*" (months) in each spell: Single-spell")) + 
                              theme(text=element_text(family=chosenFont),legend.position="inside",
                                    legend.position.inside=c(0.82,0.14), legend.background=element_rect(fill="snow2"),
                                    strip.background=element_rect(fill="snow2", colour="snow2"),
                                    strip.text=element_text(size=8, colour="gray50")) +
                              facet_wrap(~Group, scales="free", labeller=as_labeller(Facet_lbl)) + 
                              scale_colour_manual(name=paste0("Loess-smoothed hazards[", sSpan,"]"),
                                                  values=vCol2_Line, labels=Facet_grp) + 
                              guides(color=guide_legend(ncol=2)) +
                              scale_fill_manual(name="", values=vCol2_CI) + 
                              scale_y_continuous(breaks=breaks_pretty(), label=percent) + 
                              scale_x_continuous(breaks=breaks_pretty(n=8), label=comma))
### RESULTS:  As the performance spell increase, the volatility in the hazard rate increases accordingly. Note that
###           Spell 5 on wards seem to peak near the end of their life time, while the spells 2-4 peak at the start of the 
###           life time (Spell 1 would also have a similar display as the latter group if not for the left-truncated spells).

# -- Save plots
dpi <- 165 # reset
ggsave(gsurv_ht_PerfSpell, file=paste0(genFigPath,"FULL SET/gsurv_ht_PerfSpell-", mainEventName,"_Hazard-KaplanMeier-SpellLevel-MultiSpell-LatentComp-incLeftTrunc_Correct.png"),
       width=1200/dpi, height=1000/dpi,dpi=dpi, bg="white")

# Housekeeping
rm(dat, datCredit_PWPST, gsurv_Ft_PerfSpells, gsurv_ht_PerfSpell, haz_dat, km_PerfSpells, km_PerfSpells_tableSummary)




