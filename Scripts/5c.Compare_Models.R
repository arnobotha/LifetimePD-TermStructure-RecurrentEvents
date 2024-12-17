# Compare univariate models ======================================================
# Unpack tables
if (!exists('Table_TFD')) unpack.ffdf(paste0(genObjPath,"TFD_Univariate_Models"), tempPath);gc()
if (!exists('Table_PWPST')) unpack.ffdf(paste0(genObjPath,"PWPST_Univariate_Models"), tempPath);gc()

# Set up tables
Table_TFD <- subset(Table_TFD, select=c("Variable", "B_Statistic", "Concordance"))
Table_PWPST <- subset(Table_PWPST, select=c("Variable", "B_Statistic", "Concordance"))

# Merge datasets
Table_Comp <- merge(Table_TFD,Table_PWPST,by="Variable", all=TRUE,
                    suffixes = c("TFD", "PWPST")) %>% data.table()

# Melt dataset
datAggr <- melt(Table_Comp[,c(1,2,4)], id.vars = "Variable", 
                    variable.name = "Model", value.name = "B_Statistic")

# -- Graphing parameters
vCol <- brewer.pal(8, "Set1")[1:2]
vlabel <- c("TFD", "PWP-ST")
chosenFont <- "Cambria"
dpi <- 180

# B_Statistic Comparison
# Plot
(g_comp_TFD <- ggplot(datAggr, aes(x = Variable, y = B_Statistic, fill=Model)) +
                    theme_minimal() + geom_bar(stat = "identity",
                                               position = position_dodge(width = 0.8)) +  # Adjust bar spacing
                    scale_fill_manual(values = vCol, labels=vlabel) +  # Define colors for models
                    coord_flip() +  # Flip the axes for horizontal bars
                    labs(x = "Variables",y = "B Statistic", fill = "Time Definition") +
                    theme(text=element_text(family=chosenFont),legend.position="bottom") +
                    scale_y_continuous(breaks=breaks_pretty(n=8), label=percent, limits=c(0,1)))

ggsave(g_comp_TFD, file=paste0(genFigPath,"FULL SET/GoF_Comparison.png"),
       width=1200/dpi, height=1000/dpi,dpi=dpi, bg="white")

# Melt dataset
datAggr <- melt(Table_Comp[,c(1,3,5)], id.vars = "Variable", 
                variable.name = "Model", value.name = "Concordance")

# -- Graphing parameters
vCol <- brewer.pal(8, "Set1")[1:2]
vlabel <- c("TFD", "PWP-ST")
chosenFont <- "Cambria"
dpi <- 180

# B_Statistic Comparison
# Plot
(g_comp_Acc <- ggplot(datAggr, aes(x = Variable, y = Concordance, fill=Model)) +
    theme_minimal() + geom_bar(stat = "identity",
                               position = position_dodge(width = 0.8)) +  # Adjust bar spacing
    scale_fill_manual(values = vCol, labels=vlabel) +  # Define colors for models
    coord_flip() +  # Flip the axes for horizontal bars
    labs(x = "Variables",y = "Harrell's c", fill = "Time Definition") +
    theme(text=element_text(family=chosenFont),legend.position="bottom") +
    scale_y_continuous(breaks=breaks_pretty(n=8), label=percent, limits=c(0,1)))

ggsave(g_comp_Acc, file=paste0(genFigPath,"FULL SET/Acc_Comparison.png"),
       width=1200/dpi, height=1000/dpi,dpi=dpi, bg="white")

#=================================================================================
# Actual vs expected defatult rate
if (!exists('datCredit_train_TFD')) unpack.ffdf(paste0(genPath,"creditdata_train_TFD"), tempPath);gc()
if (!exists('datCredit_valid_TFD')) unpack.ffdf(paste0(genPath,"creditdata_valid_TFD"), tempPath);gc()
if (!exists('datCredit_valid_PWPST')) unpack.ffdf(paste0(genPath,"creditdata_valid_PWPST"), tempPath);gc()
if (!exists('datCredit_train_PWPST')) unpack.ffdf(paste0(genPath,"creditdata_train_PWPST"), tempPath);gc()

### BS: Move to 3c
datCredit_train_TFD[, slc_acct_arr_dir_3_Change_Ind := ifelse(slc_acct_arr_dir_3 != "SAME", 1,0)]
datCredit_valid_TFD[, slc_acct_arr_dir_3_Change_Ind := ifelse(slc_acct_arr_dir_3 != "SAME", 1,0)]
datCredit_train_PWPST[, slc_acct_arr_dir_3_Change_Ind := ifelse(slc_acct_arr_dir_3 != "SAME", 1,0)]
datCredit_valid_PWPST[, slc_acct_arr_dir_3_Change_Ind := ifelse(slc_acct_arr_dir_3 != "SAME", 1,0)]

if (!exists('cox_TFD')) unpack.ffdf(paste0(genObjPath,"TFD_Cox_Model"), tempPath);gc()
if (!exists('cox_PWPST')) unpack.ffdf(paste0(genObjPath,"PWPST_Cox_Model"), tempPath);gc()

# Term Structure======================================================

# Kaplan-Meier survival fit and hazard data preparation
km_TFD <- survfit(Surv(Start, End, Default_Ind == 1, type = "counting") ~ 1, 
                  id = PerfSpell_Key, data = datCredit_train_TFD)
haz_dat <- data.table(Time = km_TFD$time, Actual_Hazard = km_TFD$n.event / km_TFD$n.risk,
                      Surv_KM = km_TFD$surv) # Actual hazard
describe(dat$Actual_Hazard); hist(dat$Actual_Hazard, breaks="FD")

haz_dat[, Actual_Hazard_Spline := spline_estimation(Time, Actual_Hazard, 10, 3)] # Actual hazard spline
describe(dat$Actual_Hazard_Spline); hist(dat$Actual_Hazard_Spline, breaks="FD")

haz_dat[, Date := seq(as.Date(min(datCredit_train_TFD$Date)), by="month", length.out=.N)]

# Define cohort and filter data
Cohort <- "2007-11-30"
PerfSpell_Cohort <- datCredit_valid_TFD[Date == as.Date(Cohort) & Start == 0, PerfSpell_Key]
dat <- datCredit_valid_TFD[PerfSpell_Key %in% PerfSpell_Cohort]
dat <- merge(dat,haz_dat,by.x="End", by.y="Time", all.x=T)

check1 <- subset(dat, PerfSpell_Key %in% unique(dat[,PerfSpell_Key])[2])
datCredit_valid_TFD[PerfSpell_Counter==1,.N]

datCredit_valid_TFD[Start==0,.N]
dat[Start==0,.N]

datCredit_valid_TFD[Start==0 & Date==as.Date(Cohort),.N]

# Score data using models
dat[, Risk_Score := predict(cox_TFD, newdata = dat, type = "risk", id=PerfSpell_Key)]
dat[, Predicted_Hazard := Actual_Hazard_Spline * Risk_Score]
#describe(dat$Predicted_Hazard); hist(dat[Predicted_Hazard <= 0.03,Predicted_Hazard],breaks="FD")

#dat[, Predicted_cumHazard := cumsum(Predicted_Hazard), by=PerfSpell_Key]
#dat[, Predicted_SurvProb := round(exp(-Predicted_cumHazard),10)]
dat[, Predicted_SurvProb := Surv_KM^Risk_Score]
#describe(dat$Predicted_SurvProb); hist(dat$Predicted_SurvProb, breaks="FD")
#plot(dat$Predicted_SurvProb)

dat[, Predicted_SurvProb_1 := shift(Predicted_SurvProb, n=1, fill=1), by=PerfSpell_Key]
dat[, DefProb := Actual_Hazard_Spline*Predicted_SurvProb_1]
dat[, DefProb_Spline := spline_estimation(End, DefProb, 5, 3)]
#describe(dat$DefProb); hist(dat$DefProb, breaks="FD")

dat[, DefActual := Actual_Hazard*shift(Surv_KM,n=1,fill=1)]
dat[, DefActual_Spline := spline_estimation(End, DefActual, 10, 3)]
#describe(dat$DefActual); hist(dat$DefActual, breaks="FD")

#dat[, Actual_SurvProb := exp(-Actual_cumHazard)]
dat[, Survival_Prob := exp(-predict(cox_TFD, newdata=dat, type="expected", id=PerfSpel_Key))]
dat[, Default_Prob := Actual_Hazard_Spline*Survival_Prob]
dat[, Predicted_Hazard_Spline := spline_estimation(End, Predicted_Hazard, 10, 3)]

# Graphing dataset
datShow <- dat[, list(End, PerfSpell_Key, Actual_Hazard, Actual_Hazard_Spline,
                      Predicted_Hazard, Predicted_SurvProb,
                      Predicted_SurvProb_1, DefProb, DefProb_Spline, DefActual, DefActual_Spline, Risk_Score, Predicted_Hazard,
                      Survival_Prob, Default_Prob)]
Spline_MAE <- mean(abs(dat$DefActual_Spline - dat$DefProb_Spline), na.rm=T)
datAggr <- rbind(datShow[,.(End, "Probability"=DefProb, "Spline"=DefProb_Spline,
                            Type="B_Predicted")], datShow[,.(End, "Probability"=DefActual,
                                                           "Spline"=DefActual_Spline,
                                                           Type="A_Actual")])
datAggr[, ColourVar1 := paste0(Type, "_Line")]
datAggr[, ColourVar2 := paste0(Type, "_Point")]

# -- Graphing parameters
vCol_Point <- brewer.pal(8, "Pastel1")[c(1,2)]
vCol_Line <- brewer.pal(8, "Set1")[c(1,2)]
mainEventName <- "Default"
chosenFont <- "Cambria"
vLabel <- c("Actual Probability", "Actual Probability Cubic Spline")
vLabel2 <- c("A_Actual_Line"="Spline: Actual", "B_Predicted_Line"="Spline: Scored",
             "A_Actual_Point"="Case: Actual", "B_Predicted_Point"="Case: Scored")
vValues <- c("A_Actual_Line"=vCol_Line[1], "B_Predicted_Line"=vCol_Line[2],
             "A_Actual_Point"=vCol_Point[1], "B_Predicted_Point"=vCol_Point[2])
ann <- bquote(plain("Spline "*~italic(MAE)*": "*percent(Spline_MAE)))

# - Graph object for shorter time, informed by previous graphs
(Term_Structure_TFD <- ggplot(datAggr, aes(x=End, group=Type)) + theme_minimal() +
                      geom_point(aes(y=Probability, color=ColourVar2)) + 
                      geom_line(aes(y=Spline, color=ColourVar1),
                                linetype="solid") +
                      labs(y=bquote(plain(~italic(P(T~"="~t~"|"~T>=t-1))*" ["*.(mainEventName)*"]: Term Structure")), 
                           x=bquote(plain(Discrete~time~italic(t)*" (months) in spell starting "*.(Cohort)))) + 
                      annotate("text", x=80, y=0.004, family=chosenFont, label=paste0("Spline MAE: ", sprintf("%.2f", Spline_MAE*100), "%")) + 
                      theme(text=element_text(family=chosenFont),legend.position="bottom") + 
                      scale_color_manual(name = "",values = vValues, labels=vLabel2) +
                      scale_y_continuous(breaks=breaks_pretty(), label=percent) + 
                      scale_x_continuous(breaks=breaks_pretty(n=8), label=comma))
dpi <- 180 # reset
ggsave(Term_Structure_TFD, file=paste0(genFigPath,"FULL SET/Term_Structure_TFD.png"),
       width=1200/dpi, height=1000/dpi,dpi=dpi, bg="white")
#=====================================================================================
# Kaplan-Meier survival fit and hazard data preparation
km_PWPST <- survfit(Surv(Start, End, Default_Ind == 1, type = "counting") ~ 1, 
                  id = PerfSpell_Key, data = datCredit_train_PWPST)
haz_dat <- data.table(Time = km_PWPST$time, Actual_Hazard = km_PWPST$n.event / km_PWPST$n.risk,
                      Surv_KM = km_PWPST$surv) # Actual hazard
describe(dat$Actual_Hazard); hist(dat$Actual_Hazard, breaks="FD")

haz_dat[, Actual_Hazard_Spline := spline_estimation(Time, Actual_Hazard, 10, 3)] # Actual hazard spline
describe(dat$Actual_Hazard_Spline); hist(dat$Actual_Hazard_Spline, breaks="FD")

haz_dat[, Date := seq(as.Date(min(datCredit_train_PWPST$Date)), by="month", length.out=.N)]

# Define cohort and filter data
Cohort <- "2007-11-30"
PerfSpell_Cohort <- datCredit_valid_PWPST[Date == as.Date(Cohort) & Start == 0, PerfSpell_Key]
dat <- datCredit_valid_PWPST[PerfSpell_Key %in% PerfSpell_Cohort]
dat <- merge(dat,haz_dat,by.x="End", by.y="Time", all.x=T)

check1 <- subset(dat, PerfSpell_Key %in% unique(dat[,PerfSpell_Key])[2])
datCredit_valid_PWPST[PerfSpell_Counter==1,.N]

datCredit_valid_PWPST[Start==0,.N]
dat[Start==0,.N]

datCredit_valid_PWPST[Start==0 & Date==as.Date(Cohort),.N]

# Score data using models
dat[, Risk_Score := predict(cox_PWPST, newdata = dat, type = "risk", id=PerfSpell_Key)]
dat[, Predicted_Hazard := Actual_Hazard_Spline * Risk_Score]
#describe(dat$Predicted_Hazard); hist(dat[Predicted_Hazard <= 0.03,Predicted_Hazard],breaks="FD")

#dat[, Predicted_cumHazard := cumsum(Predicted_Hazard), by=PerfSpell_Key]
#dat[, Predicted_SurvProb := round(exp(-Predicted_cumHazard),10)]
dat[, Predicted_SurvProb := Surv_KM^Risk_Score]
#describe(dat$Predicted_SurvProb); hist(dat$Predicted_SurvProb, breaks="FD")
#plot(dat$Predicted_SurvProb)

dat[, Predicted_SurvProb_1 := shift(Predicted_SurvProb, n=1, fill=1), by=PerfSpell_Key]
dat[, DefProb := Actual_Hazard_Spline*Predicted_SurvProb_1]
dat[, DefProb_Spline := spline_estimation(End, DefProb, 10, 3)]
#describe(dat$DefProb); hist(dat$DefProb, breaks="FD")

dat[, DefActual := Actual_Hazard*shift(Surv_KM,n=1,fill=1)]
dat[, DefActual_Spline := spline_estimation(End, DefActual, 10, 3)]
#describe(dat$DefActual); hist(dat$DefActual, breaks="FD")

#dat[, Actual_SurvProb := exp(-Actual_cumHazard)]
#dat[, Survival_Prob := exp(-predict(cox_PWPST, newdata=dat, type="expected", id=PerfSpel_Key))]
#dat[, Default_Prob := Actual_Hazard_Spline*Survival_Prob]
#dat[, Predicted_Hazard_Spline := spline_estimation(End, Predicted_Hazard, 10, 3)]

# Graphing dataset
datShow <- dat[, list(End, PerfSpell_Key, Actual_Hazard, Actual_Hazard_Spline,
                      Predicted_Hazard, Predicted_SurvProb,
                      Predicted_SurvProb_1, DefProb, DefProb_Spline, DefActual, DefActual_Spline, Risk_Score, Predicted_Hazard)]
Spline_MAE <- mean(abs(dat$DefActual_Spline - dat$DefProb_Spline), na.rm=T)
datAggr <- rbind(datShow[,.(End, "Probability"=DefProb, "Spline"=DefProb_Spline,
                            Type="B_Predicted")], datShow[,.(End, "Probability"=DefActual,
                                                             "Spline"=DefActual_Spline,
                                                             Type="A_Actual")])
datAggr[, ColourVar1 := paste0(Type, "_Line")]
datAggr[, ColourVar2 := paste0(Type, "_Point")]

# -- Graphing parameters
vCol_Point <- brewer.pal(8, "Pastel1")[c(1,2)]
vCol_Line <- brewer.pal(8, "Set1")[c(1,2)]
mainEventName <- "Default"
chosenFont <- "Cambria"
vLabel <- c("Actual Probability", "Actual Probability Cubic Spline")
vLabel2 <- c("A_Actual_Line"="Spline: Actual", "B_Predicted_Line"="Spline: Scored",
             "A_Actual_Point"="Case: Actual", "B_Predicted_Point"="Case: Scored")
vValues <- c("A_Actual_Line"=vCol_Line[1], "B_Predicted_Line"=vCol_Line[2],
             "A_Actual_Point"=vCol_Point[1], "B_Predicted_Point"=vCol_Point[2])
ann <- bquote(plain("Spline "*~italic(MAE)*": "*percent(Spline_MAE)))

# - Graph object for shorter time, informed by previous graphs
(Term_Structure_PWPST <- ggplot(datAggr, aes(x=End, group=Type)) + theme_minimal() +
    geom_point(aes(y=Probability, color=ColourVar2)) + 
    geom_line(aes(y=Spline, color=ColourVar1),
              linetype="solid") +
    labs(y=bquote(plain(~italic(P(T~"="~t~"|"~T>=t-1))*" ["*.(mainEventName)*"]: Term Structure")), 
         x=bquote(plain(Discrete~time~italic(t)*" (months) in spell starting "*.(Cohort)))) + 
    annotate("text", x=80, y=0.004, family=chosenFont, label=paste0("Spline MAE: ", sprintf("%.2f", Spline_MAE*100), "%")) + 
    theme(text=element_text(family=chosenFont),legend.position="bottom") + 
    scale_color_manual(name = "",values = vValues, labels=vLabel2) +
    scale_y_continuous(breaks=breaks_pretty(), label=percent) + 
    scale_x_continuous(breaks=breaks_pretty(n=8), label=comma))
dpi <- 180 # reset
ggsave(Term_Structure_PWPST, file=paste0(genFigPath,"FULL SET/Term_Structure_PWPST.png"),
       width=1200/dpi, height=1000/dpi,dpi=dpi, bg="white")

#Actual vs Predicted ======================================================

# Aggregate to monthly level and observe up to given point
datAggr <- unique(datCredit_valid_TFD[DefaultStatus1==0, .(Actual = mean(DefaultStatus1_lead_12_max, na.rm=T)), by=list(Date)])

# Quick plot for visual inspection
plot(datAggr,type="l",xlab="Date", ylab="Probability", main="Default incidence rate")



# Assign risk score period level
# for loop for 1:12

# if (!exists('datMV')) unpack.ffdf(paste0(genPath,"datMV"), tempPath)
# 
# # - Lags of all MVs
# # Specifying the lags (monthly) that should be applied
# lags <- c(1,2,3,6,9,12)
# # Creating a dataset with which to check if the lags are applied correctly to the macroeconomic variables
# datMV_Check1 <- data.table(Variable = NULL, Check = NULL)
# # Getting the column names with which to apply the lags
# ColNames <- colnames(datMV)[-1]
# # Looping over the specified lags and applying each to each of the specified columns
# for (i in seq_along(lags)){
#   for (j in seq_along(ColNames)){
#     datMV[, (paste0(ColNames[j],"_",lags[i])) := fcoalesce(shift(get(ColNames[j]), n=lags[i], type="lag"), get(ColNames[j]))]
#   }
# }

dat_Future <- copy(dat) %>% mutate(Start=Start+12, End=End+12)
dat[, RiskScore := predict(cox_TFD, newdata=dat_Future, type="risk")]

datAggr <- unique(datCredit_valid_TFD[DefaultStatus1==0, .(Actual = mean(DefaultStatus1_lead_12_max, na.rm=T)), by=list(Date)])

















