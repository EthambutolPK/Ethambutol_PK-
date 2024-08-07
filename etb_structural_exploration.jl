#Region >>Load all required packages<<
using Pumas
using PumasUtilities
using NCA
using NCAUtilities
#load these packages for the majority of PUMAS functions
using GLM: lm, @formula
using Random
using CSV
using DataFramesMeta
using CairoMakie
#load other useful packages
using PharmaDatasets
using XLSX


#Add your own file path to load data
etb_day0 = DataFrame(XLSX.readtable(raw"dummy_data_day0_final.xlsx",1, infer_eltypes=true))
#load date for day0
etb_week6 = DataFrame(XLSX.readtable(raw"dummy_data_week6_final.xlsx",1, infer_eltypes=true))
#load date for week 6
etb_all = DataFrame(XLSX.readtable(raw"dummy_data_all_final.xlsx",1, infer_eltypes=true))
#load all data with dosing events for modelling

#add in the MDV column
@rtransform!(etb_day0, :MDV = ismissing(:DV) ? 1 : 0)
@rtransform!(etb_week6, :MDV = ismissing(:DV) ? 1 : 0)
@rtransform!(etb_all, :MDV = ismissing(:DV) ? 1 : 0)
#create MDV data for all data frames

#add in the rate column 
@rtransform!(etb_day0, :RATE = :EVID == 1 ? -2 : missing)
@rtransform!(etb_week6, :RATE = :EVID == 1 ? -2 : missing)
@rtransform!(etb_all, :RATE = :EVID == 1 ? -2 : missing)
#create Rate data for all data frames

#Adjust dose from mg to ng 
etb_day0.AMT = (etb_day0.AMT)*1000000
etb_week6.AMT = (etb_week6.AMT)*1000000
etb_all.AMT = (etb_all.AMT)*1000000


#Add in TFLAG 
#First time of week 6 events is 960.1 hours, 950 used as cutoff
#0 for less than --> 1 used to mark events over 6 weeks
etb_day0.TFLAG = ifelse.(etb_day0.TIME .< 950, 0, 1)
etb_week6.TFLAG = ifelse.(etb_week6.TIME .< 950, 0, 1)
etb_all.TFLAG = ifelse.(etb_all.TIME .< 950, 0, 1)


#define columns that need to be reclassed
columns_to_float = [:TIME, :CREAT, :ALT, :AMT, :WT, :BILRB]
columns_to_int = [:ID, :EVID, :MDV, :CMT, :SEX, :HIV, :RATE]

for col in columns_to_float
    replace!(etb_day0[!, col], missing => 0)
    replace!(etb_week6[!, col], missing => 0)
    replace!(etb_all[!, col], missing => 0)
#replace missing values with 0 for conversion
    etb_day0[!, col] = parse.(Float64, string.(etb_day0[!, col]))
    etb_week6[!, col] = parse.(Float64, string.(etb_week6[!, col]))
    etb_all[!, col] = parse.(Float64, string.(etb_all[!, col]))
end

# Convert these columns to integer
for col in columns_to_int
    replace!(etb_day0[!, col], missing => 0)
    replace!(etb_week6[!, col], missing => 0)
    replace!(etb_all[!, col], missing => 0)
    etb_day0[!, col] = parse.(Int, string.(etb_day0[!, col]))
    etb_week6[!, col] = parse.(Int, string.(etb_week6[!, col]))
    etb_all[!, col] = parse.(Int, string.(etb_all[!, col]))
end


#Begin modelling
#modelling will use LNDV as observations

#day 0 puma pop
etb_d0_pumapop  = read_pumas(
    etb_day0;
    id = :ID, #int
    time = :TIME, #numeric
    amt = :AMT, #float
    observations = [:LNDV], #float
    evid = :EVID,
    cmt = :CMT,
    #group = [:Dose], #we do not have groups in just the parent dataset
    route = :ROUTE, #route of administration
    #blq = :BLQ, #1 or 0
    covariates = [:WT, :AGE, :SEX, :HIV, :CREAT, :ALT, :TFLAG]
)

#create week 6 puma pop
etb_w6_pumapop  = read_pumas(
    etb_week6;
    id = :ID, #int
    time = :TIME, #numeric
    amt = :AMT, #float
    observations = [:LNDV], #float
    evid = :EVID,
    cmt = :CMT,
    #group = [:Dose], #we do not have groups in just the parent dataset
    route = :ROUTE, #route of administration
    #blq = :BLQ, #1 or 0
    covariates = [:WT, :AGE, :SEX, :HIV, :CREAT, :ALT, :TFLAG]
)



#create modelling population using all data
etb_all_pumapop  = read_pumas(
    etb_all;
    id = :ID, #int
    time = :TIME, #numeric
    amt = :AMT, #float
    observations = [:LNDV], #float
    evid = :EVID,
    cmt = :CMT,
    #group = [:Dose], #we do not have groups in just the parent dataset
    route = :ROUTE, #route of administration
    #blq = :BLQ, #1 or 0
    covariates = [:WT, :AGE, :SEX, :HIV, :CREAT, :ALT, :BILRB, :TFLAG]
)

#successful creation of pumas population



#### create 1 comp model ####
etb_1cmp = @model begin

    @metadata begin
        desc = "One Compartment Model"  #describe model in metadata
        timeu = u"hr" #describe time scale
    end

    @param begin
        """
        Clearance (L/hr)
        """
        tvcl ∈ RealDomain(; lower = 0) #bound and estimate clearance
        """
        Volume (L)
        """
        tvv ∈ RealDomain(; lower = 0) #bound and estimate volume
        """
        Absorption rate constant (h-1)
        """
        tvka ∈ RealDomain(; lower = 0) #bind and estimate absorption

        """
          - ΩCL
          - ΩVc
          - ΩKa
        """
        Ω ∈ PDiagDomain(init = [0.04, 0.04, 0.04, 0.04]) #initial matrix estimates
        """
        Proportional RUV
        """
        σ_p ∈ RealDomain(; lower = 0.0001) #proportional error, initially estimated as 20%
    end

    @random begin #define the nature of parameter variability 
        η ~ MvNormal(Ω) #define all parameters as normally distributed
    end

    @dosecontrol begin
        bioav = (; Depot1 = 1*exp(η[4]), Central = 1 / t)
        """
        Relative bioavailability (F) %
       """
    end

    #no covariates block included just yet, this is where covairate names will be loaded 

    @pre begin  #use the pre block to define parameters for the model
        CL = tvcl * exp(η[1]) #all of these use an exponential variability
        Vc = tvv * exp(η[2]) #define all parameters from their typical values
        Ka = tvka * exp(η[3])
    end

    @dynamics begin
        Depot1' = -Ka*Depot1 #define the absoprtion compartment dynamics
        Central' = Ka*Depot1 - (CL/Vc)*Central #define the central compartment dynamics
    end


    @derived begin
        cp := @. log(0.0001+ (Central / Vc)/1000) #add in a small value to ensure no 0 values (will cause error)
        #/1000 to appropriately get ng/ml from  ng/L
        """ 
        ETB Concentration (ng/mL)
        """
        LNDV ~ @. Normal(cp, σ_p) #add in error and estimate log conc (LNDV)
    end

end

#use week 0 NCA report estimates
iparms_1cmp = (;
tvcl = 1.6, tvv = 122, tvka = 5,
Ω = Diagonal([0.05, 0.05, 0.05, 0.05]),
 σ_p = 0.3)

 iparms_1cmp2 = (;
tvcl = 2, tvv = 140, tvka = 10,
Ω = Diagonal([0.05, 0.05, 0.05, 0.05]),
 σ_p = 0.3)



#run fits for all populations
fit_1cmp_day0 = fit(etb_1cmp, etb_d0_pumapop, iparms_1cmp, FOCE()) #TRUE
fit_1cmp_week6 = fit(etb_1cmp, etb_w6_pumapop, iparms_1cmp, FOCE()) #TRUE
fit_1cmp_all = fit(etb_1cmp, etb_all_pumapop, iparms_1cmp2, FOCE()) #TRUE


#day 0 VPC
vpc_1cmp_day0 = vpc(fit_1cmp_day0; samples = 100)
vpc_plot(vpc_1cmp_day0)

#week 6 vpc
vpc_1cmp_week6 = vpc(fit_1cmp_week6; samples = 100)
vpc_plot(vpc_1cmp_week6)
#VPC not very useful to assess fit over the full time period#
#using 1 day of data makes VPC far more interpretable 
#1 comp does a poor job at capturing data shape

infer_day0_1cmp = infer(fit_1cmp_day0) #get inference and parameter ests
infer_week6_1cmp = infer(fit_1cmp_week6) #get inference and parameter ests
infer_all_1cmp = infer(fit_1cmp_all) #get inference and parameter ests


inspect_day0_1cmp_df = DataFrame(inspect(fit_1cmp_day0)) #inspect data and get as df
inspect_day0_1cmp = inspect(fit_1cmp_day0) #inspect data

bic(fit_1cmp_day0) #get the bayesian (Schwarz) information criterion
bic(fit_1cmp_week6)
bic(fit_1cmp_all)





#### create 2 comp model ####
etb_2cmp = @model begin

    @metadata begin
        desc = "Two Compartment Model"  #describe model in metadata
        timeu = u"hr" #describe time scale
    end

    @param begin
        """
        Clearance (L/hr)
        """
        tvcl ∈ RealDomain(; lower = 0) #bound and estimate clearance
        """
        Volume (L)
        """
        tvv ∈ RealDomain(; lower = 0) #bound and estimate volume
        """
        Absorption rate constant (h-1)
        """
        tvka ∈ RealDomain(; lower = 0) #bind and estimate absorption

        """
        Typical volume (L) of peripheral 1 compartment
        """
        tvp1 ∈ RealDomain(; lower = 0) 
        """
        Transfer rate from Central to peripheral 1 (L/h)
        """
        tvq1  ∈ RealDomain(; lower = 0) 

        """
          - ΩCL
          - ΩVc
          - ΩKa
          - ΩVp1
          - ΩQ1
        """
        Ω ∈ PDiagDomain(init = [0.04, 0.04, 0.04, 0.04, 0.04, 0.04]) #initial matrix estimates
        """
        Proportional RUV
        """
        σ_p ∈ RealDomain(; lower = 0.0001) #proportional error, initially estimated as 20%
    end

    @random begin #define the nature of parameter variability 
        η ~ MvNormal(Ω) #define all parameters as normally distributed
    end

    @dosecontrol begin
        bioav = (; Depot1 = 1*exp(η[6]), Central = 1 / t) #initialise an eta on bioavailability
        """
        Relative bioavailability (F) %
       """
    end

    #no covariates block included just yet, this is where covairate names will be loaded 

    @pre begin  #use the pre block to define parameters for the model
        CL = tvcl * exp(η[1]) #all of these use an exponential variability
        Vc = tvv * exp(η[2]) #define all parameters from their typical values
        Ka = tvka * exp(η[3])

        Vp1 = tvp1 * exp(η[4]) 
        Q1 = tvq1 * exp(η[5])
    end

    @dynamics begin
        Depot1' = -Ka*Depot1 #define the absoprtion compartment dynamics
        Central' = Ka*Depot1 - (CL/Vc)*Central - (Q1/Vc)*Central + (Q1/Vp1)*Peripheral1 #define the central compartment dynamics
        Peripheral1' = (Q1/Vc)*Central - (Q1/Vp1)*Peripheral1
    end


    @derived begin
        cp := @. log(0.0001+ (Central / Vc)/1000) #add in a small value to ensure no 0 values (will cause error)
        #/1000 to appropriately get ng/ml from  ng/L
        """ 
        ETB Concentration (ng/mL)
        """
        LNDV ~ @. Normal(cp, σ_p) #add in error and estimate log conc (LNDV)
    end

end


#initial estimates for 2 compartment
iparms_2cmp = (;
tvcl = 60, tvv = 300, tvka = 0.4, tvp1 = 300, tvq1 = 20,
Ω = Diagonal([0.05, 0.05, 0.05, 0.05, 0.05, 0.05]),
 σ_p = 0.3)

iparms_2cmp2 = (;
tvcl = 100, tvv = 500, tvka = 0.4, tvp1 = 500, tvq1 = 20, tvp2 = 500, tvq2 = 20,
Ω = Diagonal([0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05]),
 σ_p = 0.3)


#run fits for all populations
fit_2cmp_day0 = fit(etb_2cmp, etb_d0_pumapop, iparms_2cmp, FOCE()) #TRUE
fit_2cmp_week6 = fit(etb_2cmp, etb_w6_pumapop, iparms_2cmp, FOCE())  #TRUE
fit_2cmp_all = fit(etb_2cmp, etb_all_pumapop, iparms_2cmp, FOCE()) #TRUE

#inference and parameter ests
infer_2cmp_day0 = infer(fit_2cmp_day0) 
infer_2cmp_week6 = infer(fit_2cmp_week6) 
infer_2cmp_all = infer(fit_2cmp_all) 
#parameters are massivley smaller than initial estimates

#bayesian (Schwarz) information criterion
bic(fit_2cmp_day0) 
bic(fit_2cmp_week6) 
bic(fit_2cmp_all)


#VPC
#Day 0 VPC
vpc_2cmp_day0 = vpc(fit_2cmp_day0; samples = 100) #domain error, params2 allows for plot
vpc_plot(vpc_2cmp_day0) #using params2 captures data pretty well
#week 6 VPC
vpc_2cmp_week6 = vpc(fit_2cmp_week6; samples = 100) #domain error
vpc_plot(vpc_2cmp_week6)











#### create 3 comp model ####
etb_3cmp = @model begin

    @metadata begin
        desc = "Three Compartment Model"  #describe model in metadata
        timeu = u"hr" #describe time scale
    end

    @param begin
        """
        Clearance (L/hr)
        """
        tvcl ∈ RealDomain(; lower = 0) #bound and estimate clearance
        """
        Volume (L)
        """
        tvv ∈ RealDomain(; lower = 0) #bound and estimate volume
        """
        Absorption rate constant (h-1)
        """
        tvka ∈ RealDomain(; lower = 0) #bind and estimate absorption

        """
        Typical volume (L) of peripheral 1 compartment
        """
        tvp1 ∈ RealDomain(; lower = 0) 
        """
        Transfer rate from Central to peripheral 1 (L/h)
        """
        tvq1  ∈ RealDomain(; lower = 0) 

        """
        Typical volume (L) of peripheral 2 compartment
        """
        tvp2 ∈ RealDomain(; lower = 0) 
        """
        Transfer rate from Central to peripheral 2 (L/h)
        """
        tvq2  ∈ RealDomain(; lower = 0) 

        """
          - ΩCL
          - ΩVc
          - ΩKa
          - ΩVp1
          - ΩQ1
          - ΩVp2
          - ΩQ2
        """
        Ω ∈ PDiagDomain(init = [0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04]) #initial matrix estimates
        """
        Proportional RUV
        """
        σ_p ∈ RealDomain(; lower = 0.0001) #proportional error, initially estimated as 20%
    end

    @random begin #define the nature of parameter variability 
        η ~ MvNormal(Ω) #define all parameters as normally distributed
    end

    @dosecontrol begin
        bioav = (; Depot1 = 1*exp(η[8]), Central = 1 / t) #initialise an eta on bioavailability
        """
        Relative bioavailability (F) %
       """
    end

    #no covariates block included just yet, this is where covairate names will be loaded 

    @pre begin  #use the pre block to define parameters for the model
        CL = tvcl * exp(η[1]) #all of these use an exponential variability
        Vc = tvv * exp(η[2]) #define all parameters from their typical values
        Ka = tvka * exp(η[3])

        Vp1 = tvp1 * exp(η[4]) 
        Q1 = tvq1 * exp(η[5])

        Vp2 = tvp2 * exp(η[6]) #define new peripheral (2) compartmnet variables
        Q2 = tvq2 * exp(η[7])
    end

    @dynamics begin
        Depot1' = -Ka*Depot1 #define the absoprtion compartment dynamics
        Central' = Ka*Depot1 - (CL/Vc)*Central - (Q1/Vc)*Central + (Q1/Vp1)*Peripheral1 - (Q2/Vc)*Central + (Q2/Vp2)*Peripheral2 #define the central compartment dynamics
        Peripheral1' = (Q1/Vc)*Central - (Q1/Vp1)*Peripheral1
        Peripheral2' = (Q2/Vc)*Central - (Q2/Vp2)*Peripheral2
    end


    @derived begin
        cp := @. log(0.0001+ (Central / Vc)/1000) #add in a small value to ensure no 0 values (will cause error)
        #/1000 to appropriately get ng/ml from  ng/L
        """ 
        ETB Concentration (ng/mL)
        """
        LNDV ~ @. Normal(cp, σ_p) #add in error and estimate log conc (LNDV)
    end

end


#initial estimates for 3 compartment
iparms_3cmp = (;
tvcl = 60, tvv = 300, tvka = 0.4, tvp1 = 300, tvq1 = 20, tvp2 = 300, tvq2 = 20,
Ω = Diagonal([0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05]),
 σ_p = 0.3)

iparms_3cmp2 = (;
tvcl = 100, tvv = 500, tvka = 0.4, tvp1 = 500, tvq1 = 20, tvp2 = 500, tvq2 = 20,
Ω = Diagonal([0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05]),
 σ_p = 0.3)


#run fits for all populations
fit_3cmp_day0 = fit(etb_3cmp, etb_d0_pumapop, iparms_3cmp2, FOCE()) #TRUE with 2
fit_3cmp_week6 = fit(etb_3cmp, etb_w6_pumapop, iparms_3cmp2, FOCE())  #TRUE with 2
fit_3cmp_all = fit(etb_3cmp, etb_all_pumapop, iparms_3cmp, FOCE()) #TRUE with 2

#inference and parameter ests
infer_3cmp_day0 = infer(fit_3cmp_day0) 
infer_3cmp_week6 = infer(fit_3cmp_week6) 
infer_3cmp_all = infer(fit_3cmp_all) 
#parameters are massivley smaller than initial estimates

#bayesian (Schwarz) information criterion
bic(fit_3cmp_day0) 
bic(fit_3cmp_week6) 
bic(fit_3cmp_all)


#VPC
#Day 0 VPC
vpc_3cmp_day0 = vpc(fit_3cmp_day0; samples = 100) #domain error, params2 allows for plot
vpc_plot(vpc_3cmp_day0) #using params2 captures data pretty well
#week 6 VPC
vpc_3cmp_week6 = vpc(fit_3cmp_week6; samples = 100) #domain error
vpc_plot(vpc_3cmp_week6)



#### 2 compartment is best fit ####
