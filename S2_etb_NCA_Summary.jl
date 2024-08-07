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


etb_day0 = DataFrame(XLSX.readtable(raw"ummy_data_day0_final.xlsx",1, infer_eltypes=true))
#load date for day0
etb_week6 = DataFrame(XLSX.readtable(raw"dummy_data_week6_final.xlsx",1, infer_eltypes=true))
#load date for week 6

#add in the MDV column
@rtransform!(etb_day0, :MDV = ismissing(:DV) ? 1 : 0)
@rtransform!(etb_week6, :MDV = ismissing(:DV) ? 1 : 0)
#create MDV data for all data frames

#add in the rate column 
@rtransform!(etb_day0, :RATE = :EVID == 1 ? -2 : missing)
@rtransform!(etb_week6, :RATE = :EVID == 1 ? -2 : missing)
#create Rate data for all data frame

#define columns that need to be reclassed
columns_to_float = [:TIME, :CREAT, :ALT, :AMT, :WT, :BILRB]
columns_to_int = [:ID, :EVID, :MDV, :CMT, :SEX, :HIV, :RATE]

#convert columns to float
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


#df  already loaded with columns :ID and :TIME

# Function to check monotonicity
function check_monotonicity(df)
    grouped = groupby(df, :ID)
    non_monotonic_ids = Int[]
    problematic_times = Tuple{Int, Float64, Float64}[]

    for group in grouped
        times = group.TIME
        for i in 2:length(times)
            if times[i] <= times[i-1]
                push!(non_monotonic_ids, group.ID[1])  # Capture the ID with issues
                push!(problematic_times, (group.ID[1], times[i-1], times[i]))  # Capture the problematic times
                break  # Stop checking further as we found a non-monotonic case
            end
        end
    end

    return non_monotonic_ids, problematic_times
end

# Check monotonicity in the dataset
non_monotonic_ids, problematic_times = check_monotonicity(etb_week6)


###########################################################################################################

#create NCA population for day 0
etb0_nca = read_nca(
    etb_day0;
    id = :ID, #int
    time = :TIME, #numeric
    amt = :AMT, #float
    observations = :LNDV, #float
    #evid = :EVID,
    route = :ROUTE, #route of administration
    blq = :BLQ #1 or 0
    #covariates = [:WT, :AGE, :SEX, :HIV, :CREAT, :ALT]
)


#check observations over time
pd0 = observations_vs_time(
    etb0_nca;
    paginate = true,
    axis = (; xlabel = "Time (hr)", ylabel = "ETB Conc (ng/mL)"),
    facet = (; combinelabels = true)
)
pd0[1]


summary_observations_vs_time(
    etb0_nca,
    figure = (; fontsize = 22, resolution = (800,1000)),
    linewidth = 3,
    axis = (; xlabel = "Time (hr)", ylabel = "Log ETB conc (ng/ml)")
)
#plot population data for each dose

sfit0 = subject_fits(
    etb0_nca,
    paginate = true, 
    axis = (; xlabel = "Time (hr)", ylabel = "Log ETB conc (ng/ml)"),
    facet= (; combinelabels = true, linkaxes = true)
)
#plot individual fits
sfit0[2]



#use the NCA population to create NCA (create NCA pop through read_NCA)
run_etb0_nca = run_nca(etb0_nca;
studyid = "Ethambutol_Study",
studytitle = "Day_0_Evaluations",
conclabel = "ng/mL",
author=[("Anon", "MORU")],
sigdigits=3,
)



#create report this requires an NCA run
D0report = report(
    run_etb0_nca,
    output = raw"",
    header = "Pumas-AI",
    footer = "Sensitive",
    plot_fontsize = 12,
    plot_resolution = (800, 400),
)


pk_params = [:cmax, :tmax, :aucinf_obs, :half_life] #choose parameters to be seen
stats = [minimum, maximum, mean, median,std] # choose statistics to be viewed
day0_etb_pk = summarize(run_etb0_nca.reportdf; parameters = pk_params, stats = stats)
#get summary

#create NCA population for week6
week6_nca = read_nca(
    etb_week6;
    id = :ID, #int
    time = :TIME, #numeric
    amt = :AMT, #float
    observations = :LNDV, #float
    #evid = :EVID,
    route = :ROUTE, #route of administration
    #blq = :BLQ #1 or 0
    #covariates = [:WT, :AGE, :SEX, :HIV, :CREAT, :ALT]
)


#check observations over time
pd6 = observations_vs_time(
    week6_nca;
    paginate = true,
    axis = (; xlabel = "Time (hr)", ylabel = "ETB Conc (ng/mL)"),
    facet = (; combinelabels = true)
)
pd6[1]


#check observations against time
summary_observations_vs_time(
    week6_nca,
    figure = (; fontsize = 22, resolution = (800,1000)),
    linewidth = 3,
    axis = (; xlabel = "Time (hr)", ylabel = "Log ETB conc (ng/ml)")
)
#plot population data for each dose

sfit0 = subject_fits(
    etb0_nca,
    paginate = true, 
    axis = (; xlabel = "Time (hr)", ylabel = "Log ETB conc (ng/ml)"),
    facet= (; combinelabels = true, linkaxes = true)
)
#plot individual fits
sfit0[2]


#use the NCA population to create NCA (create NCA pop through read_NCA)
run_week6_nca = run_nca(week6_nca;
studyid = "Ethambutol_Study",
studytitle = "Week_6_Evaluations",
conclabel = "ng/mL",
author=[("Anon", "MORU")],
sigdigits=3,
)



#create report this requires an NCA run ideally
W6report = report(
    run_week6_nca,
    output = raw"",
    header = "Pumas-AI",
    footer = "Sensitive",
    plot_fontsize = 12,
    plot_resolution = (800, 400),
)
#ID 60 appears to be causing an issue
#ID 60 Timepoint W6 H24 appears anomolous 


week6_etb_pk = summarize(run_week6_nca.reportdf; parameters = pk_params, stats = stats)
#get summary
#appear to be missing AUC and half_life calculations
#any way to gert these back?

#try removing ID 60
etb_week6_2 = filter(row -> row.ID != 60, etb_week6)

#re run week 6 NCA
week6_nca_2 = read_nca(
    etb_week6_2;
    id = :ID, #int
    time = :TIME, #numeric
    amt = :AMT, #float
    observations = :LNDV, #float
    #evid = :EVID,
    #group = [:Dose], #we do not have groups in just the parent dataset
    route = :ROUTE, #route of administration
    #blq = :BLQ #1 or 0
    #covariates = [:WT, :AGE, :SEX, :HIV, :CREAT, :ALT]
)



#use the NCA population to create NCA (create NCA pop through read_NCA)
run_week6_nca_2 = run_nca(week6_nca_2;
studyid = "Ethambutol_Study",
studytitle = "Week_6_Evaluations",
conclabel = "ng/mL",
author=[("Anon", "MORU")],
sigdigits=3,
)


#create report this requires an NCA run ideally
W6report_2 = report(
    run_week6_nca_2,
    output = raw"",
    header = "Pumas-AI",
    footer = "Sensitive",
    plot_fontsize = 12,
    plot_resolution = (800, 400),
)

week6_etb_pk_2 = summarize(run_week6_nca_2.reportdf; parameters = pk_params, stats = stats)
#removal of ID 60 allows for successful calculation of NCA parameters
