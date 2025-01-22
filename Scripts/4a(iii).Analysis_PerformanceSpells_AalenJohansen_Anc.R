# ========================++========= AALEN-JOHANSEN ESTIMATES =========================================
# Estimation and graphing of Aalen-Johansen estimator of the cause-specific failure probability, ie.,
# the cumulative incidence function for cause z.
# ------------------------------------------------------------------------------------------------------
# PROJECT TITLE: Default survival modelling
# SCRIPT AUTHOR(S): Dr Arno Botha
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
#   - <analytics>
# ------------------------------------------------------------------------------------------------------




# -------- 1 Aalen-Johansen analysis is on first performance spell | TFD-definition

# - Confirm prepared datasets are loaded into memory
if (!exists('datCredit_real')) unpack.ffdf(paste0(genPath,"creditdata_final4a"), tempPath)

# --- Initial data preparation for analytical purposes

# Initialize a dataset with only the first performance spell (Time to first default definition)
datTFD <- subset(datCredit_real, !is.na(PerfSpell_Num) & PerfSpell_Num==1, 
                 select=c("LoanID", "Age_Adj","PerfSpell_Key","PerfSpell_Num","PerfSpell_Counter", 
                          "TimeInPerfSpell", "PerfSpell_Age", "PerfSpell_Event","PerfSpell_Censored", "PerfSpellResol_Type_Hist"))

# - Create variables for spell-level analysis
# Create time-dependent factor variable from resolution type for period-level analysis
datTFD[, PerfSpell_Resol_factor := factor(as.character(PerfSpellResol_Type_Hist), levels=c("Censored", "Defaulted", "Paid-up", "Settled", "Written-off"))]
levels(datTFD$PerfSpell_Resol_factor)
describe(datTFD[TimeInPerfSpell==PerfSpell_Age, PerfSpell_Resol_factor]) #last spell record only

# create time of entry
datTFD[, From := TimeInPerfSpell[1]-1, by=list(PerfSpell_Key)]


# - Create variables for period-level analysis
# Create time-dependent factor variable from resolution type for period-level analysis
datTFD[, PerfSpell_Resol_time := ifelse(TimeInPerfSpell != PerfSpell_Age, "Censored", PerfSpellResol_Type_Hist)]
datTFD[, PerfSpell_Resol_time_factor := factor(as.character(PerfSpell_Resol_time), levels=c("Censored", "Defaulted", "Paid-up", "Settled", "Written-off"))]
levels(datTFD$PerfSpell_Resol_time_factor)
describe(datTFD$PerfSpell_Resol_time_factor)

# lookup
lookup <- subset(datTFD, LoanID == unique(datTFD[PerfSpell_Num == 1, LoanID])[1])
lookup <- subset(datTFD, LoanID == unique(datTFD[is.na(TimeInPerfSpell), LoanID])[1]) # should be empty



# --- Fit Aalen-Johansen model | TFD-definition | Spell-level analysis
# Analysis informed by https://cran.r-project.org/web/packages/survivalVignettes/vignettes/tutorial.html#aalen-johansen-curves
# Compute failure probability estimates using the Aalen-Johansen estimator of the cumulative incidence function per cause-specific hazard
aj_fit_spell <- survfit(Surv(time=From, time2=PerfSpell_Age, event=PerfSpell_Resol_factor,type="counting") ~ 1, 
                      id=PerfSpell_Key, data=datTFD[TimeInPerfSpell==PerfSpell_Age,])
summary(aj_fit_spell)
aj_fit_spell$transitions



# --- Graph the cumulative failure probability (cumulative incidence) per cause (resolution type)

# - Create main data object from fitted model estimates
datTransprob <- data.table(Time=aj_fit_spell$time, aj_fit_spell$pstate)
# Extract lower & upper probability bands for the 95% confidence interval
datTransprob_lower <- data.table(Time=aj_fit_spell$time, aj_fit_spell$lower)
colnames(datTransprob_lower)= c("Time", "(s0)_lower", "Defaulted_lower", "Paid-up_lower", "Settled_lower", "Written-off_lower")
datTransprob_upper <- data.table(Time=aj_fit_spell$time, aj_fit_spell$upper)
colnames(datTransprob_upper)= c("Time", "(s0)_upper", "Defaulted_upper", "Paid-up_upper", "Settled_upper", "Written-off_upper")
# merge objects
datTransprob_uplow <- merge(datTransprob_lower, datTransprob_upper, by="Time")
datTransprob <- merge(datTransprob, datTransprob_uplow, by="Time")
# cleanup
rm(datTransprob_lower, datTransprob_upper, datTransprob_uplow)


# - Pivot data object towards creating graphing object
datGraph <- pivot_longer(datTransprob, cols=all_of(c("Defaulted", "Paid-up", "Settled", "Written-off",
                                                    "(s0)_lower", "Defaulted_lower", "Paid-up_lower", "Settled_lower", "Written-off_lower",
                                                    "(s0)_upper", "Defaulted_upper", "Paid-up_upper", "Settled_upper", "Written-off_upper")),
                        names_to="Type", values_to="Value") %>% as.data.table()
# Create primary vs secondary indicator depending on whether the row pertains to main value or its error bands
datGraph[, Main := ifelse(grepl("_lower", Type, fixed=T) | grepl("_upper", Type, fixed=T), F, T )]
# Group by main type by erasing everything before "_" for secondary rows, if found
datGraph[, Group := ifelse(!Main, sub("_.*","", Type), Type)]
# Indicate lower or upper values for secondary rows
datGraph[, Value_lower := ifelse(grepl("_lower", Type, fixed=T), Value , NA)]
datGraph[, Value_upper := ifelse(grepl("_upper", Type, fixed=T), Value , NA)]
# Coalesce lower and upper values across Group by abusing an aggregation function
datGraph[, Value_lower := ifelse(!Main, mean(Value_lower, na.rm=T), NA), by=list(Time, Group)]
datGraph[, Value_upper := ifelse(!Main, mean(Value_upper, na.rm=T), NA), by=list(Time, Group)]
# Now that we have the information in one row for both upper and lower bounds, we can remove either for graphing purposes
datGraph_err <- datGraph[grepl("_lower", Type, fixed=T),]


# - Aesthetic engineering
vCol <- brewer.pal(10,"Paired")[c(8, 4, 2, 6)]
vCol2 <- brewer.pal(10,"Paired")[c(7, 3, 1, 5)]
chosenFont <- "Cambria"

# - Create main graph
(gAJ_spell <- ggplot(datGraph[Main==T & Time <= 360,], aes(x=Time,y=Value,group=Group)) + theme_minimal() +
  theme(text=element_text(family=chosenFont), legend.position="bottom") + 
  labs(x="Time spent in performing spell (months)", y="Cumulative probability (%): Aalen-Johansen") + 
  # Main graphs
  geom_line(aes(x=Time,y=Value, colour=Type, linetype=Type), linewidth=0.5) + 
  #geom_point(aes(x=Time,y=Value, colour=Type, shape=Type), size=1.5) + 
  geom_ribbon(aes(x=Time, ymin=Value_lower, ymax=Value_upper, fill=Group, colour=Group)) + 
  # Scales & options
  scale_colour_manual(name="Resolution type", values=vCol) + 
  scale_fill_manual(name="Resolution type", values=vCol2) + 
  scale_linetype_discrete(name="Resolution type") + 
  scale_shape_discrete(name="Resolution type") + 
  scale_y_continuous(labels=percent))

# - Save graph
dpi <- 200
ggsave(gAJ_spell, file=paste0(genFigPath, "FULL SET/CumulProbFail_Aalen-Johansen.png"), width=1200/dpi, height=1000/dpi, dpi=dpi, bg="white")


