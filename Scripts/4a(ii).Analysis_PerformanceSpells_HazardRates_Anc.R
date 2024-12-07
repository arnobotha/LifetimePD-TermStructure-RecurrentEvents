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

# ----------------- 1. Kaplan-Meier Analysis

# -------- 1.1 Kaplan-Meier analysis on first performance spell when left-truncation is incorporated

# Initialize a dataset with only the first performance spell
dat <- datCredit_TFD[PerfSpell_Num==1,]

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


# --- Graphing survival and related quantities from fitted KM-model | S(t), h(t)

# -- Graphing parameters
vCol <- brewer.pal(8, "Set1")[c(1)] # for F(t)
#vCol2 <- brewer.pal(10, "Paired")[c(10,9)] # for h(t)
#sSpan <- 0.1; # span for LOESS-smoother in h(t)
#vlabel <- paste0("Loess-smoothed hazard [span: ", sSpan, "]") # for h(t)
medLife <- summary(km_TFD)$table["median"]
median_survival <- data.frame(x = medLife, y = 0.5)
mainEventName <- "Default"
chosenFont <- "Cambria"
dpi <- 200

# -- Cumulative event/lifetime probability: F(t)=1-S(t)
gsurv_TFD <- ggsurvplot(km_TFD, fun="event", conf.int=T, surv.scale = "percent", legend="none", 
                        break.time.by=round(max(km_TFD$time)/8), palette=vCol,
                        xlab = bquote(Discrete~time~italic(t)*" (months) in spell: first-spell"),
                        ylab = bquote(CLD~"["*.(mainEventName)*"]"*~italic(F(t))*": Kaplan-Meier (first spell inc. left-truncation)"),
                        xlim=c(0, max(km_TFD$time)+1), censor=F, 
                        ggtheme = theme_bw(base_family=chosenFont), tables.theme = theme_cleantable(),
                        tables.height=0.10, tables.y.text=F, tables.y.text.col=T, risk.table = "abs_pct", risk.table.pos = "out",
                        cumevents=T, cumevents.title="Cumulative number of events", 
                        cumcensor=T, cumcensor.title="Cumulative number of censored observations (incl. competing risks)",
                        risk.table.title = "Number in (% of) sample at risk of main event",
                        font.family=chosenFont, fontsize=2.5,
                        surv.median.line = "hv")
# Add median annotation for illustrative purposes
(gsurv_TFD$plot <- gsurv_TFD$plot  +    
                  annotate("text", x = medLife, y = 0.5,
                  label = paste0("50th Percentile: ", medLife, " months"),
                  hjust = -0.1, color = "black", size = 3, family = chosenFont))
# Add Custom percentile line segment for illustrative purposes   
(gsurv_TFD$plot <- gsurv_TFD$plot  +    
  geom_segment(x = 0, xend=percTimepoint, y=1-survPercentile, yend = 1-survPercentile, linetype = "dashed", color = "black") +
  geom_segment(x = percTimepoint, xend=percTimepoint, y=0, yend = 1-survPercentile, linetype = "dashed", color = "black") +
  annotate("text", x = percTimepoint, y = 1 - survPercentile, label = paste0(comma((1-survPercentile)*100), "th Percentile: ", round(percTimepoint, 2), " months"),
           hjust = -0.1, color = "black", size = 3, family = chosenFont))

# - Save graph and object
ggsave(gsurv_TFD$plot, file=paste0(genFigPath, "FULL SET/KM_TFD_incLT.png"), width=1200/dpi, height=1000/dpi, dpi=dpi, bg="white")
pack.ffdf(paste0(genObjPath,"KM_TFD_RiskTable"), km_TFD_tableSummary)

# Housekeeping
rm(dat,gsurv_TFD, km_TFD, km_TFD_tableSummary, median_survival)



# -------- 1.2 Kaplan-Meier analysis on first performance spell when left-truncation is not incorporated

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
#vCol2 <- brewer.pal(10, "Paired")[c(10,9)] # for h(t)
#sSpan <- 0.1; # span for LOESS-smoother in h(t)
#vlabel <- paste0("Loess-smoothed hazard [span: ", sSpan, "]") # for h(t)
maxProb <- min(km_TFD_noLeftTrunc_tableSummary$surv)
mainEventName <- "Default"
chosenFont <- "Cambria"
dpi <- 200

# -- Cumulative event/lifetime probability: F(t)=1-S(t)
(gsurv_TFD_noLeftTrunc <- ggsurvplot(km_TFD_noLeftTrunc, fun="event", conf.int=T, surv.scale = "percent", legend="none", 
                         break.time.by=round(max(km_TFD_noLeftTrunc$time)/8), palette=vCol,
                         xlab = bquote(Discrete~time~italic(t)*" (months) in spell: first-spell"),
                         ylab = bquote(CLD~"["*.(mainEventName)*"]"*~italic(F(t))*": Kaplan-Meier (first spell exc. left-truncation)"),
                         xlim=c(0, max(km_TFD_noLeftTrunc$time)+1), censor=F, 
                         ggtheme = theme_bw(base_family=chosenFont), tables.theme = theme_cleantable(),
                         tables.height=0.10, tables.y.text=F, tables.y.text.col=T, risk.table = "abs_pct", risk.table.pos = "out",
                         cumevents=T, cumevents.title="Cumulative number of events", 
                         cumcensor=T, cumcensor.title="Cumulative number of censored observations (incl. competing risks)",
                         risk.table.title = "Number in (% of) sample at risk of main event",
                         font.family=chosenFont, fontsize=2.5))
# Add Custom percentile line segment for illustrative purposes   
(gsurv_TFD_noLeftTrunc$plot <- gsurv_TFD_noLeftTrunc$plot  +    
    geom_segment(x = 0, xend=192, y=1-maxProb, yend = 1-maxProb, linetype = "dashed", color = "black") +
    geom_segment(x = 192, xend=192, y=0, yend = 1-maxProb, linetype = "dashed", color = "black") +
    annotate("text", x = 192, y = 1 - maxProb, label = paste0(comma((1-maxProb)*100), "th Percentile: 192 months"),
             vjust=-0.5, hjust=0.5, color = "black", size = 3, family = chosenFont))

# - Save graph and object
ggsave(gsurv_TFD_noLeftTrunc$plot, file=paste0(genFigPath, "FULL SET/KM_TFD_excLT.png"), width=1200/dpi, height=1000/dpi, dpi=dpi, bg="white")
pack.ffdf(paste0(genObjPath,"KM_TFD_noLeftTrunc_RiskTable"), km_TFD_noLeftTrunc_tableSummary)

# Housekeeping
rm(dat,gsurv_TFD_noLeftTrunc, km_TFD_noLeftTrunc, km_TFD_noLeftTrunc_tableSummary)





# -------- 2 Kaplan-Meier analysis on different performance spells
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
dat[, PerfSpell_grp := fifelse(PerfSpell_Num==1,"Spells 1",fifelse(PerfSpell_Num <= 3, "Spells 2-3", fifelse(PerfSpell_Num <= 5, "Spells 4-5", "Spells 6-8")))]
# Extract median survival times from the survival fit
medLife <- summary(km_PerfSpells)$table[, "median"]
median_survival <- data.frame("Spell" = names(medLife),
                              "Median Life" = medLife,
                              x = medLife,
                              y = rep(0,length(medLife)),
                              label = paste0("Median: ", round(medLife,1)),
                              Spell_grp=c("1","1","1","2","2","3","3","3"),
                              vjust=c(0,1,0,1,0,1,0,2))
mainEventName <- "Default"
chosenFont <- "Cambria"
dpi <- 200
# -- Cumulative event/lifetime probability: F(t)=1-S(t)
(gsurv_PerfSpells <- ggsurvplot(km_PerfSpells, data=dat, fun="event", conf.int=T, surv.scale="percent",
                                palette=vCol, xlab = bquote(Discrete~time~italic(t)*" (months) in spells"),
                                ylab = bquote(CLD~"["*.(mainEventName)*"]"*~italic(F(t))*": Kaplan-Meier (By Performance Spells)"),
                                legend.title="Performance Spells", legend="bottom",
                                legend.labs=c("Spell 1", "Spell 2","Spell 3", "Spell 4",
                                              "Spell 5","Spell 6","Spell 7", "Spell 8"),
                                censor=F, ggtheme = theme_bw(base_family=chosenFont),scales="free_x",
                                facet.by="PerfSpell_grp", short.panel.labs=T))

# - Save graph and object
ggsave(gsurv_PerfSpells$plot, file=paste0(genFigPath, "FULL SET/KM_PerfSpells.png"), width=1200/dpi, height=1000/dpi, dpi=dpi, bg="white")
pack.ffdf(paste0(genObjPath,"KM_PerfSpells_RiskTable"), km_PerfSpells_tableSummary)

# -- Cumulative event/lifetime probability: F(t)=1-S(t)
(gsurv_PerfSpells <- ggsurvplot(km_PerfSpells, data=dat, fun="event", conf.int=T, surv.scale = "percent",
                                break.time.by=round(max(km_PerfSpells$time)/8), palette=vCol,
                                xlab = bquote(Discrete~time~italic(t)*" (months) in spells"),
                                ylab = bquote(CLD~"["*.(mainEventName)*"]"*~italic(F(t))*": Kaplan-Meier (By Performance Spells)"),
                                xlim=c(0, max(km_PerfSpells$time)+1), censor=F, 
                                ggtheme = theme_bw(base_family=chosenFont), facet.by="PerfSpell_Num"))
# Add Custom percentile line segment for illustrative purposes   
(gsurv_TFD_noLeftTrunc$plot <- gsurv_TFD_noLeftTrunc$plot  +    
    geom_segment(x = 0, xend=192, y=1-maxProb, yend = 1-maxProb, linetype = "dashed", color = "black") +
    geom_segment(x = 192, xend=192, y=0, yend = 1-maxProb, linetype = "dashed", color = "black") +
    annotate("text", x = 192, y = 1 - maxProb, label = paste0(comma((1-maxProb)*100), "th Percentile: 192 months"),
             vjust=-0.5, hjust=0.5, color = "black", size = 3, family = chosenFont))

# - Save graph and object
ggsave(gsurv_TFD_noLeftTrunc$plot, file=paste0(genFigPath, "FULL SET/KM_TFD_excLT.png"), width=1200/dpi, height=1000/dpi, dpi=dpi, bg="white")
pack.ffdf(paste0(genObjPath,"KM_TFD_noLeftTrunc_RiskTable"), km_TFD_noLeftTrunc_tableSummary)








# Investigate the Kaplan-Meier curves to decide on spell groupings

# - Perform analysis on the PWP ST time definition - Results are similar to that of AG/PWP TT except for the origin of each graph
if (!exists('datCredit_PWPST')) unpack.ffdf(paste0(genPath,"creditdata_final_PWPST"), tempPath)

# Create a graphing data set for the Kaplan-Meier analysis,
# 1) Remove the 9th Performance spell since only a single loan reached this point.
# 2) Remove left-truncated loans since its unclear in which loan they are when the study period started.
datAggr3 <- datCredit_PWPST[PerfSpell_Num!=9,list(Date,LoanID,Start,End,Default_Ind,PerfSpell_Num)]

# - Aesthetic engineering
chosenFont <- "Cambria"; dpi <- 200; colPalette <- "BrBG"
# Grouping spells to create a facted graph
datAggr3[ ,Spell_grp := fifelse(PerfSpell_Num <= 3, "1", fifelse(PerfSpell_Num <= 5, "2", "3"))]
# Fit survival function for the analysis
cox <- survfit(Surv(Start, End, Default_Ind) ~ PerfSpell_Num, data = datAggr3);gc()
# Extract median survival times from the survival fit
medLife <- summary(cox)$table[, "median"]
median_survival <- data.frame("Spell" = names(medLife),
                              "Median Life" = medLife,
                              x = medLife,
                              y = rep(0,length(medLife)),
                              label = paste0("Median: ", round(medLife,1)),
                              Spell_grp=c("1","1","1","2","2","3","3","3"),
                              vjust=c(0,1,0,1,0,1,0,2))

(g3 <- ggsurvplot(cox,data=datAggr3,palette = brewer.pal(8,"Dark2"),
                  censor=FALSE,ggtheme = theme_minimal() +
                    theme(text=element_text(family=chosenFont), strip_text=element_blank(), strip.background = element_blank()),
                  xlab = "Months", legend="bottom",conf.int=TRUE,
                  legend.title="Performance Spells",
                  legend.labs=c("Spell 1", "Spell 2","Spell 3", "Spell 4",
                                "Spell 5","Spell 6","Spell 7", "Spell 8"),
                  facet.by="Spell_grp", scales="free_x", surv.median.line = "v") +
    geom_label(data=median_survival, aes(x = x, y = y,label = label,
                                         group=Spell_grp),
               vjust=median_survival$vjust, alpha=0.7) +
    scale_y_continuous(label = percent_format()))
# - Create risk table object
riskTable <- g3$data.survtable

# - Save graph and object
ggsave(g3$plot, file=paste0(genFigPath, "FULL SET/Kaplan-Meier Analysis.png"), width=1200/dpi, height=1000/dpi, dpi=dpi, bg="white")
pack.ffdf(paste0(genObjPath,"Kaplan-MeierAnalysis_RiskTable"), riskTable)

# - House keeping
rm(cox, datAggr3, g3, riskTable)






























