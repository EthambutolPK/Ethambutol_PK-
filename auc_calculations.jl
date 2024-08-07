using HypothesisTests
using Random
using Plots
using StatsPlots
using DataStructures
using Pumas
using PumasUtilities
using NCA
using NCAUtilities
#load these packages for the majority of PUMAS functions
using GLM: lm, @formula
using CSV
using DataFramesMeta
using CairoMakie
#load packages

#Set seed 
Random.seed!(2024);

#CALCULATING AUC PREDICTIONS
    weight_bands = [
        (21, 34, 550),
        (35, 39, 687.5),
        (40, 54, 825),
        (55, 70, 1100), #lower, upper, and associated dose
        (71, 100, 1375) # Assuming an upper limit of 100 kg for simplicity
    ]
    
    # Initialize a structure to hold AUC results
    auc_results = Dict{String, Vector{Float64}}() #create dict with string (for weight band label) and numerical vector (for auc data)
    
    #Generic function to calculate AUC
    #AUC formula: AUC = (C1 + C2) / 2 * (t2 - t1) 
    #where C1 and C2 are the corrosponding concentrations for time t1 and t2
    function calculate_auc(sim_data) #function to be applied to simulation data
        auc = 0.0 #initialise auc as 0
        time_points = sim_data.time #get time
        concentrations = sim_data.observations.LNDV #get concentration data
        for i in 1:(length(time_points) - 1) #initialise auc calculating loop
            dt = time_points[i+1] - time_points[i] #calculate dt
            auc += (exp(concentrations[i]) + exp(concentrations[i+1])) / (2*1000) * dt #use += to update auc value, x += 1 is the same as x = x + 1
            #simulated data returns LN(mg/L) here, remove /1000 to get in ng/ml
        end
        return auc
    end
    



#Run simulations for each weight band and HIV status
for (lower, upper, dose) in weight_bands #run for all weight bands
    for hiv_status in ["positive", "negative"] #run for each HIV status option
        weight_band_label = "$lower-$upper kg (HIV $hiv_status)"
        auc_results[weight_band_label] = Float64[]
        
        for _ in 1:1000
            # Generate a random weight within the current weight band
            weight = rand(lower:upper)
            
            #Define the dosage regimen
            regimen = DosageRegimen(dose * 1e6, cmt = 1, time = 0, ii = 24, addl = 1)  #multiply dose by 1e6 to convert from mg to ng
            
            #Create a subject with the given weight and corrosponding dosage and HIV status
            hiv = hiv_status == "positive" ? 1 : 0  #HIV status pos =1 neg = 0
            subject = Subject(id=1, events=regimen, covariates=(WT=weight, HIV=hiv))
            
            #Use model to simulate for current subject
            sim = simobs(etb_pk_final, subject, iparms_pk_final, obstimes=0:1:24)
            
            #Call AUC function to calculate AUC for the simulated data
            auc = calculate_auc(sim)
    
            push!(auc_results[weight_band_label], auc) #push! used to add each auc as an element itself for the corrosponding weight band
        end
    end
end


#Define a mapping function to handle label replacement
function map_labels(labels)
    return replace.(labels, r"( \(HIV .+\))" => "") #regex to handle HIV status
end

#Sort keys to ensure the box plots are in the correct order
sorted_keys = sort(collect(keys(auc_results)), by = key -> parse(Int, split(key, "-")[1]))

#Create a new dictionary with sorted keys and corresponding values
auc_results_place = OrderedDict{String, Vector{Float64}}()
for key in sorted_keys
    auc_results_place[key] = auc_results[key]
end
auc_results = auc_results_place #reassign ordered auc results

#Collect labels and corresponding AUC data
weight_bands_labels = collect(keys(auc_results))
weight_bands = unique(map_labels(weight_bands_labels))
auc_data = [auc_results[band] for band in weight_bands_labels]

#Flatten the auc_data and create corresponding labels
flattened_auc_data = vcat(auc_data...)
simplified_labels = map_labels(weight_bands_labels)

#Convert labels to integer indices for simplified x-axis (only one weight label per weight band)
simplified_label_indices = Dict(weight_bands .=> 1:length(weight_bands))
flattened_simplified_labels = vcat([fill(simplified_label_indices[replace(label, r"( \(HIV .+\))" => "")], length(auc_results[label])) for label in weight_bands_labels]...)

#Create a color mapping for HIV status
colors = [if contains(label, "HIV positive") "red" else "blue" end for label in weight_bands_labels]
flattened_colors = vcat([fill(colors[i], length(auc_results[weight_bands_labels[i]])) for i in 1:length(weight_bands_labels)]...)

#Create a dodge mapping for HIV status
dodges = [if contains(label, "HIV positive") 0.2 else -0.2 end for label in weight_bands_labels]
flattened_dodges = vcat([fill(dodges[i], length(auc_results[weight_bands_labels[i]])) for i in 1:length(weight_bands_labels)]...)

#Convert flattened_auc_data to Float32 for use by boxplot!
flattened_auc_data_float32 = Float32.(flattened_auc_data)

#Create the grouped box plot with weight bands and HIV status on the x-axis
fig = Figure(resolution = (800, 600))
ax = Axis(fig[1, 1], title = "AUC for Different Weight Bands and HIV Status", 
          xlabel = "Weight Band", ylabel = "AUC0-24 (mg/L)", 
          xticks = (1:length(weight_bands), weight_bands)) #removed  

#Plot boxplots with simplified labels and dodging
for (i, label) in enumerate(weight_bands_labels)
    label_index = simplified_label_indices[replace(label, r" \(HIV .+\)" => "")]
    dodge_value = if contains(label, "HIV positive") 0.2 else -0.2 end
    color_value = if contains(label, "HIV positive") "red" else "blue" end
    boxplot!(ax, fill(label_index, length(auc_results[label])) .+ dodge_value, Float32.(auc_results[label]), color = color_value, transparency = true, width = 0.5)
end

#Add custom colour legend to show HIV status
text!(ax, "HIV Positive", position = (1.6, 0.9), color = "red", fontsize = 14, halign = :left)
text!(ax, "HIV Negative", position = (0.8, 0.9), color = "blue", fontsize = 14, halign = :left)
scatter!(ax, [1.55], [0.9], color = "red", markersize = 10, marker = :circle)
scatter!(ax, [0.75], [0.9], color = "blue", markersize = 10, marker = :circle)

#Rotate x-axis labels for better readability
ax.xticklabelrotation = pi/4

fig #call plot

################################## MannWhitney U test #######################################

# Prepare to store results
significance_results = Dict{String, Tuple{Float64, Bool}}()

# Perform the Wilcoxon/ Mann-Whiteny U test rank-sum test for each weight band
for weight_band in weight_bands
    hiv_positive_data = Float64[] #create empty vectors for each weight band loop
    hiv_negative_data = Float64[]
    for label in weight_bands_labels
        if startswith(label, weight_band)
            if contains(label, "HIV positive") #check current label
                append!(hiv_positive_data, auc_results[label]) #add corresponding data for that label
            elseif contains(label, "HIV negative")
                append!(hiv_negative_data, auc_results[label])
            end
        end
    end
    
    # Perform Wilcoxon rank-sum test/ Mann-Whiteny U test
    if length(hiv_positive_data) > 0 && length(hiv_negative_data) > 0 #ensure that there is something in the wieghtbands before testing
        p_value =  HypothesisTests.pvalue(MannWhitneyUTest(hiv_positive_data, hiv_negative_data)) #perform test and extract p value
        significance = p_value < 0.05 #test if significant
        significance_results[weight_band] = (p_value, significance) #add in significance data
    else
        significance_results[weight_band] = (NaN, false) #if no data
    end
end

# Display the results
for (weight_band, (p_value, significance)) in significance_results
    println("Weight Band: $weight_band, p-value: $p_value, Significant: $significance")
end


########### Prepare data for each weight band #####################

grouped_data = Dict{String, Vector{Float64}}() #initialse a vector for data


#group together HIV positive and negative data
for weight_band in weight_bands
    combined_data = Float64[]
    for label in weight_bands_labels
        if startswith(label, weight_band)
            append!(combined_data, auc_results[label])
        end
    end
    grouped_data[weight_band] = combined_data
end

#Perform the Kruskal-Wallis test
kw_data = [grouped_data[band] for band in weight_bands]
kw_result = KruskalWallisTest(kw_data...) #inpsect result
p_value = HypothesisTests.pvalue(kw_result) #extract  p value

println("Kruskal-Wallis test p-value: $p_value")

# Display the results
if p_value < 0.05
    println("There is a significant difference between weight bands.")
else
    println("There is no significant difference between weight bands.")
end




################## Separate data by HIV status and perform KW test ##################################

#initialse results vectors
grouped_data_hiv_positive = Dict{String, Vector{Float64}}()
grouped_data_hiv_negative = Dict{String, Vector{Float64}}()

#add the weight band labels for each results vector
for label in weight_bands_labels
    weight_band = map_labels([label])[1]
    if contains(label, "HIV positive")
        if !haskey(grouped_data_hiv_positive, weight_band)
            grouped_data_hiv_positive[weight_band] = Float64[]
        end
        append!(grouped_data_hiv_positive[weight_band], auc_results[label])
    elseif contains(label, "HIV negative")
        if !haskey(grouped_data_hiv_negative, weight_band)
            grouped_data_hiv_negative[weight_band] = Float64[]
        end
        append!(grouped_data_hiv_negative[weight_band], auc_results[label])
    end
end

#Perform the Kruskal-Wallis test for HIV positive
if length(grouped_data_hiv_positive) > 1
    data_hiv_positive = [grouped_data_hiv_positive[band] for band in keys(grouped_data_hiv_positive)] #get all HIV positive data
    kw_result_hiv_positive = KruskalWallisTest(data_hiv_positive...) #perform KW test 
    p_value_hiv_positive = HypothesisTests.pvalue(kw_result_hiv_positive) #extract p value
    println("Kruskal-Wallis test p-value for HIV Positive: $p_value_hiv_positive")
    if p_value_hiv_positive < 0.05
        println("There is a significant difference between weight bands for HIV Positive individuals.")
    else
        println("There is no significant difference between weight bands for HIV Positive individuals.")
    end
else
    println("Insufficient data to perform Kruskal-Wallis test for HIV Positive group.")
end

#Perform the Kruskal-Wallis test for HIV negative
if length(grouped_data_hiv_negative) > 1
    data_hiv_negative = [grouped_data_hiv_negative[band] for band in keys(grouped_data_hiv_negative)] #get all HIV negative data
    kw_result_hiv_negative = KruskalWallisTest(data_hiv_negative...) #perform KW test
    p_value_hiv_negative = HypothesisTests.pvalue(kw_result_hiv_negative) #extract p value
    println("Kruskal-Wallis test p-value for HIV Negative: $p_value_hiv_negative")
    if p_value_hiv_negative < 0.05
        println("There is a significant difference between weight bands for HIV Negative individuals.")
    else
        println("There is no significant difference between weight bands for HIV Negative individuals.")
    end
else
    println("Insufficient data to perform Kruskal-Wallis test for HIV Negative group.")
end


##################### Get proportion of each categpry above a chosen AUC threshold ##################################

#set an AUC threshold
threshold = 10

#Calculate proportion of AUC values above the threshold for HIV positive patients
proportion_above_threshold_hiv_positive = Dict{String, Float64}() #initialise vector to store results
for (weight_band, auc_values) in grouped_data_hiv_positive
    count_above_threshold = count(x -> x > threshold, auc_values) #assess whether a value is above the threshold
    proportion = count_above_threshold / length(auc_values) #get proportion
    proportion_above_threshold_hiv_positive[weight_band] = proportion #store proportion for that weight band
end

# Calculate proportion of AUC values above the threshold for HIV negative
proportion_above_threshold_hiv_negative = Dict{String, Float64}() #initialise vector to store results
for (weight_band, auc_values) in grouped_data_hiv_negative
    count_above_threshold = count(x -> x > threshold, auc_values) #assess whether a value is above the threshold
    proportion = count_above_threshold / length(auc_values) #get proportion
    proportion_above_threshold_hiv_negative[weight_band] = proportion #store proportion for that weight band
end

#Display the results
println("Proportion of AUC values above $threshold for HIV Positive:")
for (weight_band, proportion) in proportion_above_threshold_hiv_positive
    println("  Weight Band: $weight_band, Proportion: $proportion")
end

println("\nProportion of AUC values above $threshold for HIV Negative:")
for (weight_band, proportion) in proportion_above_threshold_hiv_negative
    println("  Weight Band: $weight_band, Proportion: $proportion")
end