#logistic regression
using NCA
using NCAUtilities
#load required NCA packages
using GLM
#load GLM for regression building
using CSV

etb_day0 = DataFrame(XLSX.readtable(raw"dummy_data_day0_final.xlsx",1, infer_eltypes=true))
#load date for day0
etb_week6 = DataFrame(XLSX.readtable(raw"dummy_data_week6_final",1, infer_eltypes=true))
#load date for week6
etb_pd_df = DataFrame(XLSX.readtable(raw"dummy_pd.xlsx",1, infer_eltypes=true))
#load in PD data


#add in the rate column 
@rtransform!(etb_day0, :RATE = :EVID == 1 ? -2 : missing)
@rtransform!(etb_week6, :RATE = :EVID == 1 ? -2 : missing)
#create Rate data for all data frames

#define columns that need to be reclassed
columns_to_float = [:TIME, :CREAT, :ALT, :AMT, :WT, :BILRB, :DV]
columns_to_int = [:ID, :EVID, :CMT, :SEX, :HIV, :RATE]

#convert these columns to float
for col in columns_to_float
    replace!(etb_day0[!, col], missing => 0)
    replace!(etb_week6[!, col], missing => 0)
#replace missing values with 0 for conversion
    etb_day0[!, col] = parse.(Float64, string.(etb_day0[!, col]))
    etb_week6[!, col] = parse.(Float64, string.(etb_week6[!, col]))
end

# Convert these columns to integer
for col in columns_to_int
    replace!(etb_day0[!, col], missing => 0)
    replace!(etb_week6[!, col], missing => 0)
    etb_day0[!, col] = parse.(Int, string.(etb_day0[!, col]))
    etb_week6[!, col] = parse.(Int, string.(etb_week6[!, col]))
end



#create day0 NCA population
etb0_nca = read_nca(
    etb_day0;
    id = :ID, #int
    time = :TIME, #numeric
    amt = :AMT, #float
    observations = :DV, #float
    #evid = :EVID,
    #group = [:Dose], #we do not have groups in just the parent dataset
    route = :ROUTE, #route of administration
    blq = :BLQ #1 or 0
    #covariates = [:WT, :AGE, :SEX, :HIV, :CREAT, :ALT]
)

day0_auc_df = NCA.auc(etb0_nca, interval = (0, 24)) #calculate AUC0-24 for day 0
day0_cmax_df = NCA.cmax(etb0_nca) #find Cmax for day 0
day0_thalf_df = NCA.thalf(etb0_nca) #find half life for day 0


#missing data issue solved
week6_nca = read_nca(
    etb_week6;
    id = :ID, #int
    time = :TIME, #numeric
    amt = :AMT, #float
    observations = :DV, #float
    #evid = :EVID,
    #group = [:Dose], #we do not have groups in just the parent dataset
    route = :ROUTE, #route of administration
    #blq = :BLQ #1 or 0
    #covariates = [:WT, :AGE, :SEX, :HIV, :CREAT, :ALT]
)

week6_auc_df = NCA.auc(week6_nca, interval = (0, 24)) #calculate AUC0-24 for week 6
week6_cmax_df = NCA.cmax(week6_nca) #find Cmax for week 6
week6_thalf_df = NCA.thalf(week6_nca) #find half life for week 6

expose_df = outerjoin(day0_auc_df, week6_auc_df, day0_cmax_df, week6_cmax_df, on = :id, makeunique = true)
#merge AUC and cmax data

rename!(expose_df, Dict(:auc0_24 => :AUC_d0, :auc0_24_1 => :AUC_w6, :cmax => :Cmax_d0, :cmax_1 => :Cmax_w6))
#rename cols for clarity

#calculate averages (if only one value exists use this instead)
expose_df.Average_AUC = [mean(skipmissing([row.AUC_d0, row.AUC_w6])) for row in eachrow(expose_df)]
expose_df.Average_Cmax = [mean(skipmissing([row.Cmax_d0, row.Cmax_w6])) for row in eachrow(expose_df)]


##################### form PD data set ##########################################
etb_pd_df.TIME = (etb_pd_df.Bact_clearance)*24*7
#adjust time to hours
replace!(etb_pd_df.TIME, missing => 0);
#missing time here is only for those who did not get cured (death/left study)

etb_pd_df.DV = ifelse.(ismissing.(etb_pd_df.Bact_clearance), 0, 1)
#create DV column, 1 for successful clearance, 0 for not
#0 here means the event does not happen, missing values will throw an error 


etb_pd_df.EVID = ifelse.(isinteger.(etb_pd_df.ID), 0, 1)
#Create EVID column

@rtransform! etb_pd_df :AMT = 0; 
#create AMT column to allow for read in 


### add in AUC for each ID ###

#order both data frames by ID
sort!(etb_pd_df, :ID)
sort!(expose_df, :id)

#add in AUC and cmax values to pd df
etb_pd_df.AUC = expose_df.Average_AUC
etb_pd_df.Cmax = expose_df.Average_Cmax


etb_pd_df.AUC = etb_pd_df.AUC ./ 1000 #convert AUC into mg/L
etb_pd_df.Cmax = etb_pd_df.Cmax ./ 1000 #convert AUC into mg/L

# Fit the logistic regression model for each factor seperately
#HIV
formula_HIV = @formula(DelayBC ~ HIV)
logit_model_HIV = glm(formula_HIV, etb_pd_df, Binomial(), LogitLink())
coeftable(logit_model_HIV) |> c -> c.cols[c.pvalcol][c.rownms .== "HIV"]
#p value: 0.666
bic(logit_model_HIV)
aic(logit_model_HIV)

#AGE
formula_AGE = @formula(DelayBC ~ AGE)
logit_model_AGE = glm(formula_AGE, etb_pd_df, Binomial(), LogitLink())
coeftable(logit_model_AGE) |> c -> c.cols[c.pvalcol][c.rownms .== "AGE"]
#p value: 0.189
bic(logit_model_AGE)
aic(logit_model_AGE)

#AUC
formula_AUC = @formula(DelayBC ~ AUC)
logit_model_AUC = glm(formula_AUC, etb_pd_df, Binomial(), LogitLink())
coeftable(logit_model_AUC) |> c -> c.cols[c.pvalcol][c.rownms .== "AUC"]
#p value: 0.746
bic(logit_model_AUC)
aic(logit_model_AUC)

#Cmax
formula_Cmax = @formula(DelayBC ~ Cmax)
logit_model_Cmax = glm(formula_Cmax, etb_pd_df, Binomial(), LogitLink())
coeftable(logit_model_Cmax) |> c -> c.cols[c.pvalcol][c.rownms .== "Cmax"]
#p value: 0.1437
bic(logit_model_Cmax)
aic(logit_model_Cmax)

#Initial weight 
formula_WtM0 = @formula(DelayBC ~ WtM0)
logit_model_WtM0 = glm(formula_WtM0, etb_pd_df, Binomial(), LogitLink())
coeftable(logit_model_WtM0) |> c -> c.cols[c.pvalcol][c.rownms .== "WtM0"]
#p value: 0.167
bic(logit_model_WtM0)
aic(logit_model_WtM0)

#CXL
formula_CXL = @formula(DelayBC ~ CXL)
logit_model_CXL = glm(formula_CXL, etb_pd_df, Binomial(), LogitLink())
coeftable(logit_model_CXL) |> c -> c.cols[c.pvalcol][c.rownms .== "CXL"]
#p value: 0.04744
bic(logit_model_Cmax)
aic(logit_model_Cmax)

#coefficient of 1.22
exp(1.22) #3.39

#Patients with a low Cmax have 3.39 times the odds of delayed BC compared to a Cmax > 2 mg/L 

#At p < 0.05 only CXL is signifcant


