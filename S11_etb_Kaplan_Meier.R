library(survival)
library(ggplot2)
library(survminer)
library(openxlsx)
#load required packages

setwd("C:/Users/clarkb/OneDrive - Nexus365/Documents/GitBC/place_proj") #set working drive
etb_pd_df <- read.xlsx("dummy_pd.xlsx") #call PD file

 
etb_pd_df[] <- lapply(etb_pd_df, as.numeric) #conevrt to numeric data for handling by Surv


sum(etb_pd_df$CXL)
#22 patients with a low Cmax (CXL = 1) 
#CXL = 1 means that Cmax < 2 mg/L

#Create a Surv object fo KM fitting
surv_obj <- Surv(etb_pd_df$TIME, etb_pd_df$DV)

#Fit and create a Kaplan-Meier plot
fit <- survfit(surv_obj ~ CXL, data = etb_pd_df)

#Toggle the associated varible here to stratify by other factors
#For HIV 
#fit <- survfit(surv_obj ~ HIV, data = etb_pd_df)

#Plot using ggsurvplot
ggsurvplot(fit, data = etb_pd_df, conf.int=TRUE, ggtheme=theme_minimal(),
           xlab = "Time (hours)", ylab = "Probability of Bacteriological Clearance", 
           title = "Kaplan-Meier Survival Curves Stratified by CXL",
           legend.title = "CXL",
           legend.labs = c("CXL = 0", "CXL = 1"))



