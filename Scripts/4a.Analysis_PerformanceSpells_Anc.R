# ============================================ GENRIC ANALYSES =========================================
# High-level exploratory analysis on certain diverse aspects
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

# -- Inputs:
#   - datCredit_real | Prepared from script 2f.
#
# -- Outputs:
#   - Bar chart of LoanIDs by performance spells
#   - Bar chart of risk events by performance spells
#   - Cumulative baseline hazard rates by performance spells (2 different groupings)
#   - Kaplan-Meyer analysis by performance spells
# ------------------------------------------------------------------------------------------------------


# ----------------- 1. Max spell number ---
# - Confirm prepared datasets are loaded into memory
if (!exists('datCredit_real')) unpack.ffdf(paste0(genPath,"creditdata_final4"), tempPath)

# lookup
lookup <- subset(datCredit_real, LoanID == datCredit_real[PerfSpell_Num == 8 & PerfSpell_Counter == 1, LoanID][1] )

# - Aggregation to account-level
# NOTE: Assign max conditionally since there are loans that are forever in default and hence will have no 
# information on performing spells
datAggr <- datCredit_real[, list(MaxPerfNum = ifelse(all(is.na(PerfSpell_Num)), 0, 
                                   max(PerfSpell_Num, na.rm=T)) ), by=list(LoanID)]

# - Analysis on Maximum performing spell number
describe(datAggr$MaxPerfNum); hist(datAggr$MaxPerfNum)
### RESULTS: Mean of 1.091 max spells (median: 1), with 5%-95% at [1, 2]. Large outliers of up to 10 spells

# - Rebin
datAggr[, MaxPerfNum_Binned := ifelse(MaxPerfNum >= 5, 5, MaxPerfNum)]
describe(datAggr$MaxPerfNum_Binned); hist(datAggr$MaxPerfNum_Binned)
### REUSLTS: Mean of 1.09 max spells. 1 spell: 92.6%; 2 spells: 5%; 3 spells: 1.2%; 4 Spells: 0.4%; 5+ spells: 0.2%
# Note: rebinning at 4 cap also tried, though this resulted in a mean of 1.89, which is too far from 'true' mean of 0.92,
# at least anecdotally. However, when subsampling, this binning scheme may need to be revisited to allow feasible sample sizes

# - Aesthetic engineering
chosenFont <- "Cambria"; dpi <- 200; colPalette <- "BrBG"
datAggr[MaxPerfNum_Binned > 0, MaxPerfNum_Binned_Total := .N]
totFreq <- datAggr[MaxPerfNum_Binned > 0, .N]
datAggr2 <- unique(datAggr[MaxPerfNum_Binned > 0, list(MaxPerfNum_Binned_Pc = .N / MaxPerfNum_Binned_Total,
                                  MaxPerfNum_Binned_Freq = .N), by=list(MaxPerfNum_Binned)])
datAggr2[, MaxPerfNum_Binned_Pc_labelY := MaxPerfNum_Binned_Freq + totFreq*0.005]
vLabelX <- c("1"="1", "2"="2", "3"="3", "4"="4", "5"="5+")
vBreaks <- 1:length(vLabelX)
vCol <- brewer.pal(10, "Paired")

# - Graph
(g1 <- ggplot(datAggr[MaxPerfNum_Binned > 0,], aes(x=MaxPerfNum_Binned)) + theme_minimal() + 
  theme(text=element_text(family=chosenFont)) + 
  labs(y="Frequency", x="Maximum Number of Performance Spells (pre-binned)") + 
  geom_bar(fill = vCol[2]) +
  geom_label(data=datAggr2, aes(y=MaxPerfNum_Binned_Pc_labelY, label=paste(percent(MaxPerfNum_Binned_Pc, accuracy=0.1))), 
             family=chosenFont, fill = vCol[1]) +
  scale_y_continuous(label=comma) +
  scale_x_continuous(labels=vLabelX, breaks=vBreaks) )

# - Save graph
ggsave(g1, file=paste0(genFigPath, "FULL SET/MaxPerfSpellNum_hist.png"), width=1200/dpi, height=1000/dpi, dpi=dpi, bg="white")

# House keeping
rm(datAggr,datAggr2,lookup,g1);gc()





# ----------------- 2. Risk Events according to different performance spells
# Evaluate the rate at which loans move between different performance spells over time

# Create temporary dataset for graphing that excludes the last month's data. A large proportion of loans are censored in the last month and would then distort the graphs.
datAggr2 <- datCredit_real %>% subset(!is.na(PerfSpell_Key) & Date != as.Date("2022-12-31"),
                                      select=c("Date","LoanID","PerfSpell_Key","PerfSpell_Num",
                                               "DefaultStatus1","PerfSpell_Censored", "WOff_Ind",
                                               "EarlySettle_Ind"))

# Create a new variable to indicate the risk event that occurred.
datAggr2 <- datAggr2[,RiskEvent := 
                       ifelse(DefaultStatus1 == 1 & !is.na(PerfSpell_Key),"Defaulted",
                              ifelse(WOff_Ind == 1 | EarlySettle_Ind == 1,"Competed",
                                     ifelse(PerfSpell_Censored==1, "Censored",NA)))]

# Calculate the number of risk events that occurred for each spell.
datAggr2[!is.na(RiskEvent), RiskEvent_Num := .N,
         by=list(Date,PerfSpell_Num,RiskEvent)]

datAggr2 <- datAggr2[complete.cases(datAggr2),list(Date,PerfSpell_Num,RiskEvent,
                                                   RiskEvent_Num)] %>% unique()

# Sanity Check
sum(datAggr2[PerfSpell_Num == 1 & Date == as.Date("2007-02-28"),RiskEvent_Num]) ==
  datCredit_real[PerfSpell_Num == 1 & Date == "2007-02-28",sum(DefaultStatus1,PerfSpell_Censored)] # Should be TRUE

# Calculate the total and pecentage risk events for each event type and performance spell
datAggr2[,RiskEvent_Total := sum(RiskEvent_Num), by=list(PerfSpell_Num)]
datAggr2[,RiskEventPerc := sum(RiskEvent_Num)/RiskEvent_Total, by=list(PerfSpell_Num,RiskEvent)]

# Sanity Check
sum(datAggr2[PerfSpell_Num == 1, unique(RiskEventPerc)]) == 1 # Should be TRUE
sum(unique(datAggr2[,RiskEvent_Total])) ==
  datCredit_real[!is.na(PerfSpell_Key) & Date!=as.Date("2022-12-31"),
                 sum(PerfSpell_Censored,DefaultStatus1)] # Should be TRUE

# - Aesthetic engineering
chosenFont <- "Cambria"; dpi <- 200; vCol <- brewer.pal(8,"Set1")[c(3,2,1)]
vLabel <- c("Censored","Competed","Defaulted")
maxSpell <- max(datAggr2[Date!=as.Date("2022-12-31"),PerfSpell_Num])

# Create a custom label for each facet
label <- pivot_wider(unique(datAggr2[,list(PerfSpell_Num,RiskEvent,RiskEventPerc)]),
                     names_from=RiskEvent, values_from=RiskEventPerc) %>%
          replace_na(list(Competed  = 0, Defaulted  = 0, Censored = 0)) %>%
          mutate(Label:=paste0("Defaulted: ",percent(Defaulted,accuracy=0.1),
                               ", Competed: ",percent(Competed,accuracy=0.1),
                               ", Censored: ", percent(Censored,accuracy=0.1))) %>%
          subset(select=c(PerfSpell_Num,Label))
datAggr2 <- merge(datAggr2,label,by="PerfSpell_Num")
datAggr2[, PerfSpell_Num := paste0("Spell: ", PerfSpell_Num)]

# - Graph
(g2 <- ggplot(datAggr2, aes(x=Date,y=RiskEvent_Num,fill=factor(RiskEvent))) +
    theme_minimal() + 
    theme(text=element_text(family=chosenFont), axis.text.x=element_text(angle=90),
          strip.background = element_rect(fill="snow2", colour="snow2"),
          legend.position = "inside",legend.position.inside=c(0.1,0.15),
          legend.background=element_rect(fill="snow2"),
          strip.text=element_text(colour="black")) + 
    labs(y="Number of Risk Events", x=bquote("Cohort Date ("*italic(mmm)*" "*italic(ccyy)*")")) + 
    geom_col() + geom_text(aes(x=as.Date("2008-01-31"),family=chosenFont,y=Inf,label=Label),
                           vjust=1, hjust=-0.35, size=2.6, color="snow4") +
    facet_wrap(PerfSpell_Num ~ ., scales="free_y",nrow=maxSpell, strip.position="right") +
    scale_fill_manual(name="Risk Type", values=vCol, labels=vLabel) +
    scale_y_continuous(breaks=pretty_breaks(n=4)) +#, labels=label_number(accuracy=1)) +
    scale_x_date(date_labels = "%b %y", date_breaks = "6 months",
                 limits = c(as.Date("2007-01-31"),as.Date("2022-12-31"))))
### RESULTS:  1)  The default behavior of the first performance spell is different from that of the rest, since it 
#                 significant competing risks are present.
#             2)  As the Performance Spell number increase the sample size decreases. This poses a problem if models are
#                 to be built on the performance spells. This can be remedied by binning some of the later performance spells together.

# - Save graph
ggsave(g2, file=paste0(genFigPath, "/FULL SET/Performance Spell Risk Events.png"), width=1200/dpi, height=1000/dpi, dpi=dpi, bg="white")

# House keeping
rm(datAggr2,censored,g);gc()




# ----------------- 3. Kaplan-Meier Analysis
# Investigate the Kaplan-Meier curves to decide on spell groupings

# Create appropriate variables (according to PWP ST to ensure a common origination point)
datCredit_real <- datCredit_real %>% mutate(Start = ifelse(is.na(PerfSpell_Counter),NA,PerfSpell_Counter-1), # Records the start of time interval
                                            End = ifelse(is.na(PerfSpell_Counter),NA,PerfSpell_Counter), # Records the end of time interval
                                            Default_Ind = ifelse(!is.na(PerfSpell_Num) & DefaultStatus1==1,1,
                                                                 ifelse(!is.na(PerfSpell_Num),0,NA))) %>%
                                      filter(!is.na(PerfSpell_Num));gc()

# Create a graphing data set for the Kaplan-Meier analysis and without the 9th performance spell since it contains a singel spell that was censored
datAggr3 <- datCredit_real[PerfSpell_Num!=9,list(Date,LoanID,Start,End,Default_Ind,PerfSpell_Num)]

# - Aesthetic engineering
chosenFont <- "Cambria"; dpi <- 200; colPalette <- "BrBG"

cox <- survfit(Surv(Start, End, Default_Ind) ~ PerfSpell_Num, data = datAggr3);gc()
(g3 <- ggsurvplot(cox,data=datCredit_real,palette = brewer.pal(8,"Dark2"),
                  censor=FALSE,ggtheme = theme_minimal() +
                  theme(text=element_text(family=chosenFont),
                  plot.title=element_text(hjust=0.5,size=15, face="bold")),
                 xlab="Time (Months)", legend="bottom",conf.int=TRUE,
                 legend.title="Performance Spells",
                 legend.labs=c("Spell 1", "Spell 2","Spell 3", "Spell 4",
                               "Spell 5","Spell 6","Spell 7", "Spell 8")))
# - Create risk table object
riskTable <- g3$data.survtable

# - Save graph and object
ggsave(g3$plot, file=paste0(genFigPath, "FULL SET/Kaplan-Meier Analysis.png"), width=1200/dpi, height=1000/dpi, dpi=dpi, bg="white")
pack.ffdf(paste0(genObjPath,"Kaplan-MeierAnalysis_RiskTable"), riskTable)

# - House keeping
rm(cox, datAggr3, g3, riskTable)




# ----------------- 4. Baseline hazard rates for different performance spells
# Display the baseline hazard rate for different performance spells

describe(datCredit_real[DefaultStatus1 == 1 & !is.na(PerfSpell_Num) | PerfSpell_Censored == 1,PerfSpell_Num])
hist(datCredit_real[DefaultStatus1 == 1 & !is.na(PerfSpell_Num) | PerfSpell_Censored == 1,PerfSpell_Num])
### RESULTS: Mean of 1.127 performance spell numbers; 1 spell: 91.1%, 2 spells: 6.3%, 3 spells: 1.7%, 4 spells: .6%, 5 spells: .2%, 6+ spells .1%.

# Attempt to group performance spells greater than 4 together
# Create variable containing the grouping
datCredit_real[,PerfSpell_Group_5 := ifelse(PerfSpell_Num < 5, PerfSpell_Num, 5)]
describe(datCredit_real[DefaultStatus1 == 1 & !is.na(PerfSpell_Num) | PerfSpell_Censored == 1,PerfSpell_Group_5])
hist(datCredit_real[DefaultStatus1 == 1 & !is.na(PerfSpell_Num) | PerfSpell_Censored == 1,PerfSpell_Group_5])
### RESULTS: Mean of 1.126 performance spell numbers; 1 spell: 91.1%, 2 spells: 6.3%, 3 spells: 1.7%, 4 spells: .6%. and 5 spells: .3%.
### RESULTS: Negligible change in the mean value 

basehaz <- basehaz(coxph(Surv(Start,End,Default_Ind)~1,datCredit_real[PerfSpell_Num == 1])) %>%
  data.table() %>% subset(select=c("time","hazard"))
colnames(basehaz) <- c("Time","cumHazard_1")
for (i in 2:5){
  temp_basehaz <- basehaz(coxph(Surv(Start,End,Default_Ind)~1,
                                datCredit_real[PerfSpell_Group_5 == i])) %>%
    data.table()
  colnames(temp_basehaz) <- c(paste0("cumHazard_",i),"Time")
  basehaz <- left_join(basehaz,temp_basehaz,by="Time");gc()
};

# - House keeping
rm(temp_basehaz); gc()

# - Aesthetic engineering
chosenFont <- "Cambria"; dpi <- 200; colPalette <- "BrBG"
cols1 <- c("cumHazard_1","cumHazard_2","cumHazard_3","cumHazard_4", "cumHazard_5")
datAggr4 <- basehaz[,c("Time",..cols1)]
datAggr4 <- pivot_longer(datAggr4,cols=all_of(cols1), names_to="Type",values_to="cumHazard")


# - Graph
(g4 <- ggplot(datAggr4, aes(x = Time, y = cumHazard, color = Type)) + 
    theme_minimal() +
    theme(text = element_text(family = chosenFont),
          axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
          legend.position=c(0.15,0.8), legend.background=element_rect(fill="snow2"),
          strip.text=element_text(size=8, colour="black")) +
    labs(y = "Cumulative Baseline Hazard Rate (5+ Grouping)",
         x = "Time (months)") + 
    geom_line(size = 1) +
    scale_color_manual(values = brewer.pal(5, "Set1"), name="Performance Spell",
                       labels=c("Spell 1", "Spell 2", "Spell 3", "Spell 4", "Spell 5+")))  # Use Set1 palette for color


# - Save graph
ggsave(g4, file=paste0(genFigPath, paste0("FULL SET/Cumulative Baseline Hazards (5+).png")), width=1200/dpi, height=1000/dpi, dpi=dpi, bg="white")

# - House keeping
rm(basehaz, datAggr4, g4)

# Attempt to group performance spells greater than 3 together
# Create variable containing the grouping
datCredit_real[,PerfSpell_Group_4 := ifelse(PerfSpell_Num < 4, PerfSpell_Num, 4)]
describe(datCredit_real[DefaultStatus1 == 1 & !is.na(PerfSpell_Num) | PerfSpell_Censored == 1,PerfSpell_Group_4])
hist(datCredit_real[DefaultStatus1 == 1 & !is.na(PerfSpell_Num) | PerfSpell_Censored == 1,PerfSpell_Group_5])
### RESULTS: Mean of 1.123 performance spell numbers; 1 spell: 91.1%, 2 spells: 6.3%, 3 spells: 1.7%, 4 spells: .9%
### RESULTS: Minor change in the mean value 

basehaz <- basehaz(coxph(Surv(Start,End,Default_Ind)~1,datCredit_real[PerfSpell_Num == 1])) %>%
  data.table() %>% subset(select=c("time","hazard"))
colnames(basehaz) <- c("Time","cumHazard_1")
for (i in 2:4){
  temp_basehaz <- basehaz(coxph(Surv(Start,End,Default_Ind)~1,
                                datCredit_real[PerfSpell_Group_5 == i])) %>%
    data.table()
  colnames(temp_basehaz) <- c(paste0("cumHazard_",i),"Time")
  basehaz <- left_join(basehaz,temp_basehaz,by="Time");gc()
};

# - House keeping
rm(temp_basehaz); gc()

# - Aesthetic engineering
chosenFont <- "Cambria"; dpi <- 200; colPalette <- "BrBG"
cols1 <- c("cumHazard_1","cumHazard_2","cumHazard_3","cumHazard_4")
datAggr5 <- basehaz[,c("Time",..cols1)]
datAggr5 <- pivot_longer(datAggr5,cols=all_of(cols1), names_to="Type",values_to="cumHazard")


# - Graph
(g5 <- ggplot(datAggr5, aes(x = Time, y = cumHazard, color = Type)) + 
    theme_minimal() +
    theme(text = element_text(family = chosenFont),
          axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
          legend.position=c(0.15,0.8), legend.background=element_rect(fill="snow2"),
          strip.text=element_text(size=8, colour="black")) +
    labs(y = "Cumulative Baseline Hazard Rate (4+ Grouping)",
         x = "Time (months)") + 
    geom_line(size = 1) +
    scale_color_manual(values = brewer.pal(5, "Set1"), name="Performance Spell",
                       labels=c("Spell 1", "Spell 2", "Spell 3", "Spell 4+")))  # Use Set1 palette for color


# - Save graph
ggsave(g5, file=paste0(genFigPath, paste0("FULL SET/Cumulative Baseline Hazards (4+).png")), width=1200/dpi, height=1000/dpi, dpi=dpi, bg="white")

# - House keeping
rm(basehaz, datAggr5, g5)
