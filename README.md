Welcome to the Ethambutol Pharmacokinetics Repo! Please read me to understand how to interact with all of the files here.
The prefixed with “dummy” contain dummy data to run each of these files. They are loaded in a few of the Julia and R files here, please make sure to adjust to your own file path before running!
etb_structural_exploration.jl : This file explores the number of distribution compartments for ethambutol and also contains the code to create Pumas populations! These are used in other files, as such it is recommended to run this file first.
Etb_NCA_summary.jl: Perform Non-compartmental Analysis for observations on Day 0 and week 6.
Etb_WT_HIV_exploration.jl: Explore allometric weight and HIV covariates (stepwise forward step). REQUIRES POPULATIONS FROM etb_structural_exploration.
Etb_forward_backwards_selection.jl: Forward and backwards selection of covariates.  REQUIRES POPULATIONS FROM etb_structural_exploration.jl AND DATA FROM etb_WT_HIV_exploration.jl
Etb_IIV_stepwsie.jl: Stepwise removal of Inter-individual variability. REQUIRES POPULATIONS FROM etb_structural_exploration.jl AND DATA FROM etb_WT_HIV_exploration.jl
Etb_final_model_and_diagnostics.jl: Final PK model and accompanying diagnostics.
Etb_logistic_regression: Logistic regression of different factors for delayed bacteriological clearance.
Auc_calculations.jl: Simulate AUC using the model and test for any significant differences.
Cmax_predictions.jl: Simulate Cmax using the model and test for any significant differences.
Etb_kaplan_meier: Kaplan Meier curve, stratifying by CXL.
