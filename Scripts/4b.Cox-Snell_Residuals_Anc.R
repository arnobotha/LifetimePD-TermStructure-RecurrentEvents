# =============================== Cox-Snell method Exploration =========================================
# Exploratory analysis on Cox-Snell residual methods, particularly how it relates to the high censoring
# events. The script compares different methods to evaluate Cox-Snell residuals graphically by using the
# same underlying cox proportional hazard model fitted on a single variable.
# ------------------------------------------------------------------------------------------------------
# PROJECT TITLE: Classifier Diagnostics
# SCRIPT AUTHOR(S): Bernard Scheepers
# ------------------------------------------------------------------------------------------------------
# -- Script dependencies:
#   - 0.Setup.R
#   - 1.Data_Import.R
#   - 2a.Data_Prepare_Macro.R
#   - 2b.Data_Prepare_Credit.R
#   - 2c.Data_Enrich.R
#   - 3a.Data_Transform.R
#   - 3c.Data_Fusion_TFD.R

# -- Inputs:
#   - datCredit_train_TFD | Prepared from script 3c.
#
# -- Outputs:
#   - General analysis
# ------------------------------------------------------------------------------------------------------

# - Confirm prepared datasets are loaded into memory
if (!exists('datCredit_train_TFD')) unpack.ffdf(paste0(genPath,"creditdata_train_TFD"), tempPath);gc()

# Subset the data set to select specific columns
dat <- datCredit_train_TFD %>% subset(select=c(Date, LoanID, PerfSpell_Max_Date, Start, End, Default_Ind, BalanceToTerm))
dat[,PerfSpell_Exit_Ind := ifelse(Date==PerfSpell_Max_Date,1,0)] # Create an indicator vaiable for detecting the end of a performance spell

# Initialize variables
valid <- dat$PerfSpell_Exit_Ind==1 # Binary vector to index the end observation of each performance spell.
cox <- coxph(Surv(Start, End, Default_Ind) ~ BalanceToTerm, data = datCredit_train_TFD) # Build a cox model based on a single variable
rexp <- rexp(nrow(dat),1);hist(rexp,breaks="FD") # Evaluate the unit exponential distribution of Cox-Snell residuals





# ----------------- 1. Trivial Cox-Snell residuals calculation ---

# Cox-Snell residuals
cs.1 <- get_csvec(cox); hist(cs.1,breaks="FD") # Evaluate distribution of Cox-Snell residuals
mean(cs.1)
### Result: 0.002201133
### Conclusion: Extremely low mean, i.e. it should tend to 1. This is likely due to the influence of censored Cox-Snell residuals.

# Cox-Snell hazard rate
datGraph.1 <- survfit(coxph(Surv(cs.1, cox$y[, "status"]) ~ 1, method = "breslow"), type = "aalen") %>% tidy() %>%
              mutate(cumu_hazard = -log(.data$estimate)) %>% rename(coxsnell = "time", survival = "estimate") %>% subset(select=c("coxsnell","cumu_hazard"))

# KS statistic
KS.1 <- round(1 - ks.test(datGraph.1$coxsnell,rexp)$statistic,2)

# Graph Cox-Snell residuals
(g.1 <- ggplot(datGraph.1,aes(x=coxsnell, y=cumu_hazard )) + geom_point() + geom_step() + xlab("Cox-Snell Residual") + ylab("Cumulative Hazard Function") + 
        geom_abline(slope=1,intercept=0,color='red',linetype="dashed", linewidth=1) + geom_point() + geom_step() + theme_minimal() + theme(text = element_text(family="Cambria")) +
        annotate("label",label=paste0("KS statistic: ",KS.1),fill="lightgrey",
                 x=min(datGraph.1$coxsnell) + 0.1 * (max(datGraph.1$coxsnell) - min(datGraph.1$coxsnell)),
                 y=max(datGraph.1$cumu_hazard) - 0.1 * (max(datGraph.1$cumu_hazard) - min(datGraph.1$cumu_hazard))))
### RESULTS:  The curve seems appropriate with it meandering around the 45 degree line, however, the Cox-Snell residuals are 
#             extremely low with a max at about 0.006. It seems unlikely that the residuals follow a unit exponential
#             distribution, which is emphasized with the KS statistic of 0.01. This is likely due to the high censoring rate.

# House keeping
rm(datGraph.1,g.1,cs.1,KS.1); gc()





# ----------------- 2. Reduced Cox-Snell residuals calculation ---
# This method tries to reduce the amount of censoring by only limiting the Cox-Snell residuals to
# observations that contain risk events, i.e. end of each performance spell.

# Cox-Snell residuals
cs.2 <- dat[PerfSpell_Exit_Ind==1,Default_Ind] - residuals(cox,type="martingale",collapse=dat$LoanID) # Collapse residuals to only the observations at the end of each performance spell.
mean(cs.2)
### Result: 0.1474528
### Conclusion: Slightly improved mean, i.e. it is higher than previous. Therefore influence of censored Cox-Snell residuals is still severe.

# Cox-Snell hazard rate
datGraph.2 <- survfit(coxph(Surv(cs.2, cox$y[, "status"][valid]) ~ 1, method = "breslow"), type = "aalen") %>% tidy() %>%
              mutate(cumu_hazard = -log(.data$estimate)) %>% rename(coxsnell = "time", survival = "estimate") %>% subset(select=c("coxsnell","cumu_hazard"))

# KS statistic
KS.2 <- round(1 - ks.test(datGraph.2$coxsnell,rexp)$statistic,2)

# Graph Cox-Snell residuals
(g.2 <- ggplot(datGraph.2,aes(x=coxsnell, y=cumu_hazard )) + geom_point() + geom_step() + xlab("Cox-Snell Residual") + ylab("Cumulative Hazard Function") + 
        geom_abline(slope=1,intercept=0,color='red',linetype="dashed", linewidth=1) + geom_point() + geom_step() + theme_minimal() + theme(text = element_text(family="Cambria")) +
        annotate("label",label=paste0("KS statistic: ",KS.2),fill="lightgrey",
                 x=min(datGraph.2$coxsnell) + 0.2 * (max(datGraph.2$coxsnell) - min(datGraph.2$coxsnell)),
                 y=max(datGraph.2$cumu_hazard) - 0.1 * (max(datGraph.2$cumu_hazard) - min(datGraph.2$cumu_hazard))))
### RESULTS:  The curve seems less appropriate with it meandering less around the 45 degree line, however, the Cox-Snell residuals are 
#             higher with many values exceeding 1. It seems likely that the residuals follow a unit exponential
#             distribution, which is emphasized with the KS statistic of 0.45. However the is still a high censoring rate of 85%.

# House keeping
rm(datGraph.2,g.2,cs.2,KS.2); gc()




# ----------------- 3. Adjusted Cox-Snell residuals calculation ---
# This method tries to mitigate the amount of censoring by adjusting the Cox-Snell residuals of
# observations that contain an actual value of 0. This is done by adding a term these Cox-Snell residuals.
# It is recommended to add log(2) since it represents the median value of a unit exponential distribution.

# Cox-Snell residuals
cs.3 <- get_csvec(cox) + log(2)*(1 - dat[,Default_Ind]) # Add log(2) to all observations that have a 0.
mean(cs.3)
### Result: 0.6938226
### Conclusion: Reduced censoring rate, although it is still close to the threshold of 30%.

# Cox-Snell hazard rate
datGraph.3 <- survfit(coxph(Surv(cs.3, cox$y[, "status"]) ~ 1, method = "breslow"), type = "aalen") %>% tidy() %>%
              mutate(cumu_hazard = -log(.data$estimate)) %>% rename(coxsnell = "time", survival = "estimate") %>% subset(select=c("coxsnell","cumu_hazard"))

# KS statistic
KS.3 <- round(1 - ks.test(datGraph.3$coxsnell,rexp)$statistic,2)

# Graph Cox-Snell residuals
(g.3 <- ggplot(datGraph.3,aes(x=coxsnell, y=cumu_hazard )) + geom_point() + geom_step() + xlab("Cox-Snell Residual") + ylab("Cumulative Hazard Function") + 
        geom_abline(slope=1,intercept=0,color='red',linetype="dashed", linewidth=1) + geom_point() + geom_step() + theme_minimal() + theme(text = element_text(family="Cambria")) +
        annotate("label",label=paste0("KS statistic: ",KS.3),fill="lightgrey",
                x=min(datGraph.3$coxsnell) + 0.2 * (max(datGraph.3$coxsnell) - min(datGraph.3$coxsnell)),
                y=max(datGraph.3$cumu_hazard) - 0.1 * (max(datGraph.3$cumu_hazard) - min(datGraph.3$cumu_hazard))))
### RESULTS:  The curve seems appropriate with it meandering around the 45 degree line, although it seems to have some extreme values.
#             The Cox-Snell residuals are also lower with none exceeding 0.8. It seems less likely that the residuals follow a unit exponential
#             distribution, although the KS statistic is 0.5.

# House keeping
rm(datGraph.3,g.3,cs.3,KS.3); gc()





# ----------------- 4. Adjusted reduced Cox-Snell residuals calculation ---
# This method tries to mitigate the amount of censoring further by adjusting and reducing the Cox-Snell residuals
# of observations that contain an actual value of 0. This is done by adding a term these Cox-Snell residuals.
# It only includes observations that contain risk events, i.e. end of each performance spell.

# Cox-Snell residuals
cs.4 <- dat[PerfSpell_Exit_Ind==1,Default_Ind] - residuals(cox,type="martingale",collapse=dat$LoanID) + log(2)*(1 - dat[PerfSpell_Exit_Ind==1,Default_Ind])
mean(cs.4)
### Result: 0.7383935
### Conclusion: Reduced censoring rate considerably.

# Cox-Snell hazard rate
datGraph.4 <- survfit(coxph(Surv(cs.4, cox$y[, "status"][valid]) ~ 1, method = "breslow"), type = "aalen") %>% tidy() %>%
              mutate(cumu_hazard = -log(.data$estimate)) %>% rename(coxsnell = "time", survival = "estimate") %>% subset(select=c("coxsnell","cumu_hazard"))

# KS statistic
KS.4 <- round(1 - ks.test(datGraph.4$coxsnell,rexp)$statistic,2)

# Graph Cox-Snell residuals
(g.4 <- ggplot(datGraph.4,aes(x=coxsnell, y=cumu_hazard )) + geom_point() + geom_step() + xlab("Cox-Snell Residual") + ylab("Cumulative Hazard Function") + 
        geom_abline(slope=1,intercept=0,color='red',linetype="dashed", linewidth=1) + geom_point() + geom_step() + theme_minimal() + theme(text = element_text(family="Cambria")) +
        annotate("label",label=paste0("KS statistic: ",KS.4),fill="lightgrey",
                x=min(datGraph.4$coxsnell) + 0.2 * (max(datGraph.4$coxsnell) - min(datGraph.4$coxsnell)),
                y=max(datGraph.4$cumu_hazard) - 0.1 * (max(datGraph.4$cumu_hazard) - min(datGraph.4$cumu_hazard))))
### RESULTS:  The curve seems inappropriate for the 45 degree line, although it seems to have some extreme negative values.
#             This is in contrast with the KS statistic that do indeed to have improved from previous graphs, although the 
#             extreme negative values do seem to be peculiar.


# House keeping
rm(datGraph.4,g.4,cs.4,KS.4); gc()





# ----------------- 5. Censoring rate for various methods---

# Censoring rate for all Cox-Snell residuals
sum(datCredit_train_TFD$Default_Ind)/length(datCredit_train_TFD$Default_Ind)
# Results: 0.002201133 -> a default proportion of about 0.2%

# Censoring rate for Cox-Snell residuals conatining an event occurence
sum(datCredit_train_TFD$Default_Ind[valid])/length(datCredit_train_TFD$Default_Ind[valid])
# Results: 0.1474528 -> a default proportion of about 15%

# House keeping
rm(datCredit_train_TFD,dat,cox,rexp,valid); gc()