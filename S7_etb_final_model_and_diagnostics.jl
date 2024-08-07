#### Load packages ####
using Pumas
using PumasUtilities
using NCA
using NCAUtilities
#load these packages for the majority of PUMAS functions
using GLM: lm, @formula
using Loess
using Random
using CSV
using DataFramesMeta
using CairoMakie
#load other useful packages
using PharmaDatasets
#Endregion
using XLSX
#### End ####

#### Load data set ####
etb_all = DataFrame(XLSX.readtable(raw"dummy_data_all_final.xlsx",1, infer_eltypes=true))
#load all data with dosing events for modelling
etb_data = DataFrame(XLSX.readtable(raw"dummy_medians.xlsx",1, infer_eltypes=true))

#### End ####

#### Apply transformations ####

#add in the MDV column
@rtransform!(etb_all, :MDV = ismissing(:DV) ? 1 : 0)

#add in the rate column 
@rtransform!(etb_all, :RATE = :EVID == 1 ? -2 : missing)

#Adjust dose from mg to ng 
etb_all.AMT = (etb_all.AMT)*1000000

#### End ####

#### Re-class objects ####

#define columns that need to be reclassed
columns_to_float = [:TIME, :CREAT, :ALT, :AMT, :WT, :BILRB]
columns_to_int = [:ID, :EVID, :MDV, :CMT, :SEX, :HIV, :RATE]

for col in columns_to_float
    replace!(etb_all[!, col], missing => 0)
    etb_all[!, col] = parse.(Float64, string.(etb_all[!, col]))
end

# Convert these columns to integer
for col in columns_to_int
    replace!(etb_all[!, col], missing => 0) #replace missing value with 0 to allow for transfrormation
    etb_all[!, col] = parse.(Int, string.(etb_all[!, col]))
end

#### End ####

#### Create Pumas population ####
etb_all_pumapop  = read_pumas(
    etb_all;
    id = :ID, #int
    time = :TIME, #numeric
    amt = :AMT, #float
    observations = [:LNDV], #float, Log(ng/ml) units
    evid = :EVID, #int
    cmt = :CMT, #numeric
    route = :ROUTE, #route of administration
    covariates = [:WT, :HIV] #add in covariates
)

#### End ####


#### Load in final model ####

#WT + MTT
#renamed as PK_final
etb_pk_final = @model begin

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
        tvmtt ∈ RealDomain(; lower = 0) #bind and estimate absorption
        """
        Typical volume (L) of peripheral 1 compartment
        """
        tvp1 ∈ RealDomain(; lower = 0) 
        """
        Transfer rate from Central to peripheral 1 (L/h)
        """
        tvq1  ∈ RealDomain(; lower = 0) 

        """
        effect of HIV on Volume hivV
        """
        hivMTT  ∈ RealDomain(;) #should we have a lower bound

        """
          - ΩCL
          - ΩVc
          - Ωmtt
          - ΩVp1
          - ΩQ1
        """
        Ω ∈ PDiagDomain(init = [0.04, 0.04, 0.04, 0.04, 0.04]) #initial matrix estimates
        """
        Proportional RUV
        """
        σ_p ∈ RealDomain(; lower = 0.0001) #proportional error, initially estimated as 20%
    end


    @random begin #define the nature of parameter variability 
        η ~ MvNormal(Ω) #define all parameters as normally distributed
    end

    @covariates begin
        WT #add in weight as a covairate
        HIV #add in HIV status as a covariate
    end

    #no covariates block included just yet, this is where covairate names will be loaded 

    @pre begin  #use the pre block to define parameters for the model
        CL = tvcl * exp(η[1]) * (WT / median(etb_data.WT))^0.75 #add in the effect of WT on CL/Q
        Vc = tvv * (WT / median(etb_data.WT))  #add in the effect of weight on volume

        #no IIV on central volume (Vc)

        MTT = tvmtt * exp(η[2]) * (1 + (hivMTT*HIV))#define mean transit time
        ktr = 2/MTT #define rate of Transfer

        Vp1 = tvp1 * exp(η[3]) * (WT / median(etb_data.WT))  #add in the effect of weight on volume
        Q1 = tvq1 * exp(η[4]) * (WT / median(etb_data.WT))^0.75 #add in the effect of WT on CL/Q
    end

    @dosecontrol begin
        bioav = (; Depot1 = 1 * exp(η[5])) #add an eta into bioavailbility

        """
        Relative bioavailability (F) %
        All explored implementations ad no impact on fit
       """
    end


    @dynamics begin
        Depot1' = -ktr*Depot1 #define the absoprtion compartment dynamics
        Transit1' = ktr*Depot1 - ktr*Transit1
        Central' = ktr*Transit1 - (CL/Vc)*Central - (Q1/Vc)*Central + (Q1/Vp1)*Peripheral1 #define the central compartment dynamics
        Peripheral1' = (Q1/Vc)*Central - (Q1/Vp1)*Peripheral1
    end


    @derived begin
        cp := @. log(0.0001+ (Central / Vc)/1000) #add in a small value to ensure no 0 values (will cause error)
        #/1000 to appropriately get ng/ml from ng/L
        """ 
        ETB Concentration (ng/mL)
        """
        LNDV ~ @. Normal(cp, σ_p) #add in error and estimate log conc (LNDV)
    end

end


iparms_pk_final = (; #initial parameters as model parameter estimates
tvcl = 53.547, tvv = 286.32, tvmtt = 1.6521, tvp1 = 755.29, tvq1 = 43.645, hivMTT = 0.44871,
Ω = Diagonal([0.024833, 0.13914, 0.52548, 0.06841, 0.048107]), #using zeroed etas 
σ_p = 0.38293)

#fit model
fit_pk_final_all = fit(etb_pk_final, etb_all_pumapop, iparms_pk_final, FOCE())#True minimization


#bayesian (Schwarz) information criterion
bic(fit_pk_final_all)

#TAD VPC
vpc_pk_final_all = vpc(fit_pk_final_all; samples = 2000, covariates = [:tad])
vpc_plot(vpc_pk_final_all)


#inference and parameter ests
infer_pk_final_all = infer(fit_pk_final_all) 
#parameters are massivley smaller than initial estimates


#inspect fit element for GOF plots
etb_pk_final_inspect = inspect(fit_pk_final_all)


goodness_of_fit(etb_pk_final_inspect; figure = (; fontsize = 12))
#GOF plots

########## Create custom TAD GOF plot ##################################

etb_final_inspect_df = DataFrame(etb_pk_final_inspect)
#get inpsect element as dataframe

describe(etb_final_inspect_df[!,(12:32)])
#final data frame contains columns: LNDV_pred, LNDV_ipred, LNDV_wres, LNDV_iwres, and tad

#extract desired oclumns for custom plotting
gof_df = select(etb_final_inspect_df, [:LNDV_pred, :LNDV_ipred, :LNDV_wres, :LNDV_iwres, :tad])
#lots of missing values 
#all prediction values appear to be missing in the same place, all at tad = 0 (dosing events)

gof_df2 = dropmissing(gof_df)
#get rid of missing predictions

any(gof_df2.tad == 0)
#no dosing events

describe(gof_df2)
#inspect columns

gof_df2 = mapcols(x -> Float32.(x), gof_df2)
#map to Float32 for scatter plot

#ols requires float 64 type
df_ols = mapcols(x -> Float64.(x), gof_df2)

#create OLS
tad_wres_ols_model = lm(@formula(LNDV_wres ~ tad), df_ols)
tad_wres_ols_fit = predict(tad_wres_ols_model)

#create a new df for Loess to ensure data is ordered by :tad
df_loess = mapcols(x -> Float32.(x), gof_df2)
df_loess_sorted = sort(df_loess, :tad)

#create Loess line
x_vals = LinRange(minimum(df_loess_sorted.tad), maximum(df_loess_sorted.tad), 100) #use a range for a smoother fit
tad_wres_loess_model = loess(df_loess_sorted.tad, df_loess_sorted.LNDV_wres)
tad_wres_loess_fit = predict(tad_wres_loess_model, x_vals) #get Loess predctions

#creat figure for plot
tad_fig = Figure(resolution = (500, 300), title = "LNDV") #initialise a new figure
ax = Axis(tad_fig[1, 1], title = "ETB concentration (ng/ml)", xlabel = "Time after dose (h)", ylabel = "Weighted population residuals") #create an axis

#add scatter points of WRES and TAD
scatter!(ax, gof_df2.tad, gof_df2.LNDV_wres, color = :black)

#OLS fit line in red
lines!(ax, gof_df2.tad, tad_wres_ols_fit, color = :green)

#LOESS fit line in green
lines!(ax, x_vals, tad_wres_loess_fit, color = :red)

display(tad_fig) #call plot


####### SIR ######################

etb_sir = infer(fit_pk_final_all, Pumas.SIR(samples=500, resamples = 100))
#run Sampling importance resampling 

#get results in a tabulated form
etb_sir_coef = coeftable(etb_sir)

etb_sir_coef
#call tabulated results


###### Simulated observations ######






