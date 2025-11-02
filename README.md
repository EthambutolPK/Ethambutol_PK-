Welcome to the Ethambutol Pharmacokinetics Repo! Please read me to understand how to interact with all of the files here.
The prefixed with “dummy” contain dummy data to run each of these files. They are loaded in a few of the Julia and R files here, please make sure to adjust to your own file path before running!

--------------------------Dummy data files removed whilst upgrading these for better performance-------------------------

Please also note that the models here are sensitive to intial parameters: should a model fail to fit initial parameters may need to be adjusted.

S1_etb_structural_exploration.jl : This file explores the number of distribution compartments for ethambutol and also contains the code to create Pumas populations! These are used to fit the models in other files, as such it is recommended to run this file first.

S2_Etb_NCA_summary.jl: Perform Non-compartmental Analysis for observations on Day 0 and week 6.

S3_Etb_WT_HIV_exploration.jl: Explore allometric weight and HIV covariates (stepwise forward step). REQUIRES POPULATIONS FROM etb_structural_exploration.

S4_etb_backwards_selection.jl: Backwards elimination of covariates. REQUIRES POPULATIONS FROM etb_structural_exploration.jl AND DATA FROM etb_WT_HIV_exploration.jl

S5_etb_forward_selection.jl: Forward selection of covariates.REQUIRES POPULATIONS FROM etb_structural_exploration.jl AND DATA FROM etb_WT_HIV_exploration.jl

S6_etb_IIV_stepwsie.jl: Stepwise removal of Inter-individual variability. REQUIRES POPULATIONS FROM etb_structural_exploration.jl AND DATA FROM etb_WT_HIV_exploration.jl

S7_etb_final_model_and_diagnostics.jl: Final PK model and accompanying diagnostics.

S8_etb_logistic_regression: Logistic regression of different factors for delayed bacteriological clearance.

S9_auc_calculations.jl: Simulate AUC using the model and test for any significant differences.

S10_cmax_predictions.jl: Simulate Cmax using the model and test for any significant differences.

S11_etb_kaplan_meier: Kaplan Meier curve, stratifying by CXL.
