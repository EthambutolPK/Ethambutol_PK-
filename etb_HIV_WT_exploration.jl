etb_data = DataFrame(XLSX.readtable(raw"dummy_medians.xlsx",1, infer_eltypes=true))
#load ethambutol data to be able to extract median values of covariates


#### create 2 comp model ####
etb_covariate_WT = @model begin

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
        tvmtt ∈ RealDomain(; lower = 0) # bound and estimate absorption
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
          - Ωmtt
          - ΩVp1
          - ΩQ1
          - ΩF
        """
        Ω ∈ PDiagDomain(init = [0.04, 0.04, 0.04, 0.04, 0.04, 0.04]) #initial matrix estimates
        """
        Proportional RUV
        """
        σ_p ∈ RealDomain(; lower = 0.0001) #proportional error 
    end


    @random begin #define the nature of parameter variability 
        η ~ MvNormal(Ω) #define all parameters as normally distributed
    end

    @covariates begin
        WT #add in weight as a covairate
    end

    #no covariates block included just yet, this is where covairate names will be loaded 

    @pre begin  #use the pre block to define parameters for the model
        CL = tvcl * exp(η[1]) * (WT / median(etb_data.WT))^0.75 #add in the effect of WT on CL/Q
        Vc = tvv * exp(η[2]) * (WT / median(etb_data.WT)) #add in the effect of weight on volume

        MTT = tvmtt * exp(η[3]) #define mean transit time
        ktr = 2/MTT #define rate of Transfer

        Vp1 = tvp1 * exp(η[4]) * (WT / median(etb_data.WT)) #add in the effect of weight on volume
        Q1 = tvq1 * exp(η[5]) * (WT / median(etb_data.WT))^0.75 #add in the effect of WT on CL/Q
    end

    @dosecontrol begin
        bioav = (; Depot1 = 1 * exp(η[6])) #add an eta into bioavailbility

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
        #/1000 to appropriately get ng/ml from  ng/L
        """ 
        ETB Concentration (ng/mL)
        """
        LNDV ~ @. Normal(cp, σ_p) #add in error and estimate log conc (LNDV)
    end

end


#initial parameter estimates
iparms_covariate_WT = (;
tvcl = 50, tvv = 250, tvmtt = 1.5, tvp1 = 800, tvq1 = 40,
Ω = Diagonal([0.05, 0.05, 0.05, 0.05, 0.05, 0.05]),
σ_p = 0.3)


iparms_covariate_WT2 = (;
tvcl = 65, tvv = 285, tvmtt = 1.79, tvp1 = 345, tvq1 = 40,
Ω = Diagonal([0.05, 0.05, 0.05, 0.05, 0.05, 0.05]),
σ_p = 0.3)

#run fits for all populations
fit_covariate_WT_day0 = fit(etb_covariate_WT, etb_d0_pumapop, iparms_covariate_WT2, FOCE()) #true with 2
fit_covariate_WT_week6 = fit(etb_covariate_WT, etb_w6_pumapop, iparms_covariate_WT, FOCE()) #true
fit_covariate_WT_all = fit(etb_covariate_WT, etb_all_pumapop, iparms_covariate_WT, FOCE())
fit_covariate_WT_all2 = fit(etb_covariate_WT, etb_all_pumapop, iparms_covariate_WT2, FOCE()) #2 is a better BIC
#inference and parameter ests
infer_covariate_WT_day0 = infer(fit_covariate_WT_day0) 
infer_covariate_WT_week6 = infer(fit_covariate_WT_week6)  
infer_covariate_WT_all = infer(fit_covariate_WT_all2) 
#parameters are massivley smaller than initial estimates

#bayesian (Schwarz) information criterion
bic(fit_covariate_WT_day0) 
bic(fit_covariate_WT_week6) 
bic(fit_covariate_WT_all2)

#Day 0 VPC
vpc_covariate_WT_day0 = vpc(fit_covariate_WT_day0; samples = 100)
vpc_plot(vpc_covariate_WT_day0)
#week 6 VPC
vpc_covariate_WT_week6 = vpc(fit_covariate_WT_week6; samples = 100)
vpc_plot(vpc_covariate_WT_week6)




#### create 2 comp model ####
etb_WT_HIV_V = @model begin

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
        tvmtt ∈ RealDomain(; lower = 0) # bound and estimate absorption
        """
        Typical volume (L) of peripheral 1 compartment
        """
        tvp1 ∈ RealDomain(; lower = 0) 
        """
        Transfer rate from Central to peripheral 1 (L/h)
        """
        tvq1  ∈ RealDomain(; lower = 0) 

        """
        effect of HIV on central Volume
        """
        hivV  ∈ RealDomain(;) #should we have a lower bound

        """
          - ΩCL
          - ΩVc
          - Ωmtt
          - ΩVp1
          - ΩQ1
          - ΩF
        """
        Ω ∈ PDiagDomain(init = [0.04, 0.04, 0.04, 0.04, 0.04, 0.04]) #initial matrix estimates
        """
        Proportional RUV
        """
        σ_p ∈ RealDomain(; lower = 0.0001) #proportional error 
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
        Vc = tvv * exp(η[2]) * (WT / median(etb_data.WT)) * (1 + (hivV*HIV)) #add in the effect of weight on volume

        MTT = tvmtt * exp(η[3]) #define mean transit time
        ktr = 2/MTT #define rate of Transfer

        Vp1 = tvp1 * exp(η[4]) * (WT / median(etb_data.WT)) #add in the effect of weight on volume
        Q1 = tvq1 * exp(η[5]) * (WT / median(etb_data.WT))^0.75 #add in the effect of WT on CL/Q
    end

    @dosecontrol begin
        bioav = (; Depot1 = 1 * exp(η[6])) #add an eta into bioavailbility

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
        #/1000 to appropriately get ng/ml from  ng/L
        """ 
        ETB Concentration (ng/mL)
        """
        LNDV ~ @. Normal(cp, σ_p) #add in error and estimate log conc (LNDV)
    end

end

iparms_WT_HIV_V = (;
tvcl = 50, tvv = 250, tvmtt = 1.5, tvp1 = 800, tvq1 = 40, hivV = 0.4,
Ω = Diagonal([0.05, 0.05, 0.05, 0.05, 0.05, 0.05]), #add in eta for bioavailibility
 σ_p = 0.3)



#run fits for all populations
fit_WT_HIV_V_day0 = fit(etb_WT_HIV_V, etb_d0_pumapop, iparms_WT_HIV_V, FOCE())
fit_WT_HIV_V_week6 = fit(etb_WT_HIV_V, etb_w6_pumapop, iparms_WT_HIV_V, FOCE()) 
fit_WT_HIV_V_all = fit(etb_WT_HIV_V, etb_all_pumapop, iparms_WT_HIV_V, FOCE()) #true minimization

#inference and parameter ests
infer_WT_HIV_V_day0 = infer(fit_WT_HIV_V_day0) 
infer_WT_HIV_V_week6 = infer(fit_WT_HIV_V_week6) 
infer_WT_HIV_V_all = infer(fit_WT_HIV_V_all) 
#parameters are massivley smaller than initial estimates

#bayesian (Schwarz) information criterion
bic(fit_WT_HIV_V_day0) 
bic(fit_WT_HIV_V_week6) 
bic(fit_WT_HIV_V_all)

#Day 0 VPC
vpc_WT_HIV_V_day0 = vpc(fit_WT_HIV_V_day0; samples = 100)
vpc_plot(vpc_WT_HIV_V_day0)
#week 6 VPC
vpc_WT_HIV_V_week6 = vpc(fit_WT_HIV_V_week6; samples = 100)
vpc_plot(vpc_WT_HIV_V_week6)







#### create 2 comp model ####
etb_WT_HIV_CL = @model begin

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
        tvmtt ∈ RealDomain(; lower = 0) # bound and estimate absorption
        """
        Typical volume (L) of peripheral 1 compartment
        """
        tvp1 ∈ RealDomain(; lower = 0) 
        """
        Transfer rate from Central to peripheral 1 (L/h)
        """
        tvq1  ∈ RealDomain(; lower = 0) 

        """
        effect of HIV on clearance
        """
        hivCL  ∈ RealDomain(;) #should we have a lower bound

        """
          - ΩCL
          - ΩVc
          - Ωmtt
          - ΩVp1
          - ΩQ1
          - ΩF
        """
        Ω ∈ PDiagDomain(init = [0.04, 0.04, 0.04, 0.04, 0.04, 0.04]) #initial matrix estimates
        """
        Proportional RUV
        """
        σ_p ∈ RealDomain(; lower = 0.0001) #proportional error 
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
        CL = tvcl * exp(η[1]) * (WT / median(etb_data.WT))^0.75 * (1 + (hivCL*HIV)) #add in the effect of WT on CL/Q
        Vc = tvv * exp(η[2]) * (WT / median(etb_data.WT))  #add in the effect of weight on volume

        MTT = tvmtt * exp(η[3]) #define mean transit time
        ktr = 2/MTT #define rate of Transfer

        Vp1 = tvp1 * exp(η[4]) * (WT / median(etb_data.WT))  #add in the effect of weight on volume
        Q1 = tvq1 * exp(η[5]) * (WT / median(etb_data.WT))^0.75 #add in the effect of WT on CL/Q
    end

    @dosecontrol begin
        bioav = (; Depot1 = 1 * exp(η[6])) #add an eta into bioavailbility

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
        #/1000 to appropriately get ng/ml from  ng/L
        """ 
        ETB Concentration (ng/mL)
        """
        LNDV ~ @. Normal(cp, σ_p) #add in error and estimate log conc (LNDV)
    end

end

iparms_WT_HIV_CL = (;
tvcl = 50, tvv = 250, tvmtt = 1.5, tvp1 = 800, tvq1 = 40, hivCL = 0.4,
Ω = Diagonal([0.05, 0.05, 0.05, 0.05, 0.05, 0.05]), #add in eta for bioavailibility
 σ_p = 0.3)


iparms_WT_HIV_CL2 = (;
 tvcl = 54, tvv = 289, tvmtt = 1.88, tvp1 = 754, tvq1 = 44, hivCL = -0.023,
 Ω = Diagonal([0.05, 0.05, 0.05, 0.05, 0.05, 0.05]), #add in eta for bioavailibility
  σ_p = 0.3)
 
 #run fits for all populations
fit_WT_HIV_CL_day0 = fit(etb_WT_HIV_CL, etb_d0_pumapop, iparms_WT_HIV_CL, FOCE())
fit_WT_HIV_CL_week6 = fit(etb_WT_HIV_CL, etb_w6_pumapop, iparms_WT_HIV_CL, FOCE()) 
fit_WT_HIV_CL_all = fit(etb_WT_HIV_CL, etb_all_pumapop, iparms_WT_HIV_CL, FOCE()) #better fit with 1

#inference and parameter ests
infer_WT_HIV_CL_day0 = infer(fit_WT_HIV_CL_day0) 
infer_WT_HIV_CL_week6 = infer(fit_WT_HIV_CL_week6) 
infer_WT_HIV_CL_all = infer(fit_WT_HIV_CL_all) 
#parameters are massivley smaller than initial estimates

#bayesian (Schwarz) information criterion
bic(fit_WT_HIV_CL_day0) 
bic(fit_WT_HIV_CL_week6) 
bic(fit_WT_HIV_CL_all) 

#Day 0 VPC
vpc_WT_HIV_CL_day0 = vpc(fit_WT_HIV_CL_day0; samples = 100)
vpc_plot(vpc_WT_HIV_CL_day0)
#week 6 VPC
vpc_WT_HIV_CL_week6 = vpc(fit_WT_HIV_CL_week6; samples = 100)
vpc_plot(vpc_WT_HIV_CL_week6)

















#### create 2 comp model ####
etb_WT_HIV_VP1 = @model begin

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
        tvmtt ∈ RealDomain(; lower = 0) # bound and estimate absorption
        """
        Typical volume (L) of peripheral 1 compartment
        """
        tvp1 ∈ RealDomain(; lower = 0) 
        """
        Transfer rate from Central to peripheral 1 (L/h)
        """
        tvq1  ∈ RealDomain(; lower = 0) 

        """
        effect of HIV on peripheral volume
        """
        hivVP1  ∈ RealDomain(;) #should we have a lower bound

        """
          - ΩCL
          - ΩVc
          - Ωmtt
          - ΩVp1
          - ΩQ1
          - ΩF
        """
        Ω ∈ PDiagDomain(init = [0.04, 0.04, 0.04, 0.04, 0.04, 0.04]) #initial matrix estimates
        """
        Proportional RUV
        """
        σ_p ∈ RealDomain(; lower = 0.0001) #proportional error 
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
        Vc = tvv * exp(η[2]) * (WT / median(etb_data.WT)) #add in the effect of weight on volume

        MTT = tvmtt * exp(η[3]) #define mean transit time
        ktr = 2/MTT #define rate of Transfer

        Vp1 = tvp1 * exp(η[4]) * (WT / median(etb_data.WT)) * (1 + (hivVP1*HIV))#add in the effect of weight on volume
        Q1 = tvq1 * exp(η[5]) * (WT / median(etb_data.WT))^0.75 #add in the effect of WT on CL/Q
    end

    @dosecontrol begin
        bioav = (; Depot1 = 1 * exp(η[6])) #add an eta into bioavailbility

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
        #/1000 to appropriately get ng/ml from  ng/L
        """ 
        ETB Concentration (ng/mL)
        """
        LNDV ~ @. Normal(cp, σ_p) #add in error and estimate log conc (LNDV)
    end

end

iparms_WT_HIV_VP1 = (;
tvcl = 50, tvv = 250, tvmtt = 1.5, tvp1 = 800, tvq1 = 40, hivVP1 = 0.4,
Ω = Diagonal([0.05, 0.05, 0.05, 0.05, 0.05, 0.05]), #add in eta for bioavailibility
 σ_p = 0.3)

iparms_WT_HIV_VP12 = (;
 tvcl = 54, tvv = 295, tvmtt = 1.87, tvp1 = 1020, tvq1 = 44, hivVP1 = -0.51,
 Ω = Diagonal([0.05, 0.05, 0.05, 0.05, 0.05, 0.05]), #add in eta for bioavailibility
  σ_p = 0.3)

 #run fits for all populations
fit_WT_HIV_VP1_day0 = fit(etb_WT_HIV_VP1, etb_d0_pumapop, iparms_WT_HIV_VP1, FOCE())
fit_WT_HIV_VP1_week6 = fit(etb_WT_HIV_VP1, etb_w6_pumapop, iparms_WT_HIV_VP1, FOCE()) 
fit_WT_HIV_VP1_all = fit(etb_WT_HIV_VP1, etb_all_pumapop, iparms_WT_HIV_VP1, FOCE()) #BIC: 1865.603
fit_WT_HIV_VP1_all2 = fit(etb_WT_HIV_VP1, etb_all_pumapop, iparms_WT_HIV_VP12, FOCE()) #BIC 1865.621
#inference and parameter ests
infer_WT_HIV_VP1_day0 = infer(fit_WT_HIV_VP1_day0) 
infer_WT_HIV_VP1_week6 = infer(fit_WT_HIV_VP1_week6) 
infer_WT_HIV_VP1_all = infer(fit_WT_HIV_VP1_all) 
#parameters are massivley smaller than initial estimates

#bayesian (Schwarz) information criterion
bic(fit_WT_HIV_VP1_day0) 
bic(fit_WT_HIV_VP1_week6) 
bic(fit_WT_HIV_VP1_all)
bic(fit_WT_HIV_VP1_all2)
#Day 0 VPC
vpc_WT_HIV_VP1_day0 = vpc(fit_WT_HIV_VP1_day0; samples = 100)
vpc_plot(vpc_WT_HIV_VP1_day0)
#week 6 VPC
vpc_WT_HIV_VP1_week6 = vpc(fit_WT_HIV_VP1_week6; samples = 100)
vpc_plot(vpc_WT_HIV_VP1_week6)







#### create 2 comp model ####
etb_WT_HIV_Q = @model begin

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
        tvmtt ∈ RealDomain(; lower = 0) # bound and estimate absorption
        """
        Typical volume (L) of peripheral 1 compartment
        """
        tvp1 ∈ RealDomain(; lower = 0) 
        """
        Transfer rate from Central to peripheral 1 (L/h)
        """
        tvq1  ∈ RealDomain(; lower = 0) 

        """
        effect of HIV on intercompartmental clearance
        """
        hivQ  ∈ RealDomain(;) #should we have a lower bound

        """
          - ΩCL
          - ΩVc
          - Ωmtt
          - ΩVp1
          - ΩQ1
          - ΩF
        """
        Ω ∈ PDiagDomain(init = [0.04, 0.04, 0.04, 0.04, 0.04, 0.04]) #initial matrix estimates
        """
        Proportional RUV
        """
        σ_p ∈ RealDomain(; lower = 0.0001) #proportional error 
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
        Vc = tvv * exp(η[2]) * (WT / median(etb_data.WT))  #add in the effect of weight on volume

        MTT = tvmtt * exp(η[3]) #define mean transit time
        ktr = 2/MTT #define rate of Transfer

        Vp1 = tvp1 * exp(η[4]) * (WT / median(etb_data.WT))  #add in the effect of weight on volume
        Q1 = tvq1 * exp(η[5]) * (WT / median(etb_data.WT))^0.75 * (1 + (hivQ*HIV)) #add in the effect of WT on CL/Q
    end

    @dosecontrol begin
        bioav = (; Depot1 = 1 * exp(η[6])) #add an eta into bioavailbility

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
        #/1000 to appropriately get ng/ml from  ng/L
        """ 
        ETB Concentration (ng/mL)
        """
        LNDV ~ @. Normal(cp, σ_p) #add in error and estimate log conc (LNDV)
    end

end

iparms_WT_HIV_Q = (;
tvcl = 50, tvv = 250, tvmtt = 1.5, tvp1 = 800, tvq1 = 40, hivQ = 0.4,
Ω = Diagonal([0.05, 0.05, 0.05, 0.05, 0.05, 0.05]), #add in eta for bioavailibility
 σ_p = 0.3)

 iparms_WT_HIV_Q2 = (;
 tvcl = 54, tvv = 286, tvmtt = 1.9, tvp1 = 754, tvq1 = 41, hivQ = 0.25,
 Ω = Diagonal([0.05, 0.05, 0.05, 0.05, 0.05, 0.05]), #add in eta for bioavailibility
  σ_p = 0.3)
 #run fits for all populations
fit_WT_HIV_Q_day0 = fit(etb_WT_HIV_Q, etb_d0_pumapop, iparms_WT_HIV_Q, FOCE())
fit_WT_HIV_Q_week6 = fit(etb_WT_HIV_Q, etb_w6_pumapop, iparms_WT_HIV_Q, FOCE()) 
fit_WT_HIV_Q_all = fit(etb_WT_HIV_Q, etb_all_pumapop, iparms_WT_HIV_Q, FOCE()) #BIC 1872.851
fit_WT_HIV_Q_all2 = fit(etb_WT_HIV_Q, etb_all_pumapop, iparms_WT_HIV_Q2, FOCE()) #BIC 1872.857
#inference and parameter ests
infer_WT_HIV_Q_day0 = infer(fit_WT_HIV_Q_day0) 
infer_WT_HIV_Q_week6 = infer(fit_WT_HIV_Q_week6) 
infer_WT_HIV_Q_all = infer(fit_WT_HIV_Q_all) 
#parameters are massivley smaller than initial estimates

#bayesian (Schwarz) information criterion
bic(fit_WT_HIV_Q_day0) 
bic(fit_WT_HIV_Q_week6) 
bic(fit_WT_HIV_Q_all) 

#Day 0 VPC
vpc_WT_HIV_Q_day0 = vpc(fit_WT_HIV_Q_day0; samples = 100)
vpc_plot(vpc_WT_HIV_Q_day0)
#week 6 VPC
vpc_WT_HIV_Q_week6 = vpc(fit_WT_HIV_Q_week6; samples = 100)
vpc_plot(vpc_WT_HIV_Q_week6)





#### create 2 comp model ####
etb_WT_HIV_F = @model begin

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
        tvmtt ∈ RealDomain(; lower = 0) # bound and estimate absorption
        """
        Typical volume (L) of peripheral 1 compartment
        """
        tvp1 ∈ RealDomain(; lower = 0) 
        """
        Transfer rate from Central to peripheral 1 (L/h)
        """
        tvq1  ∈ RealDomain(; lower = 0) 

        """
        effect of HIV on bioavailability
        """
        hivF  ∈ RealDomain(;) #should we have a lower bound

        """
          - ΩCL
          - ΩVc
          - Ωmtt
          - ΩVp1
          - ΩQ1
          - ΩF
        """
        Ω ∈ PDiagDomain(init = [0.04, 0.04, 0.04, 0.04, 0.04, 0.04]) #initial matrix estimates
        """
        Proportional RUV
        """
        σ_p ∈ RealDomain(; lower = 0.0001) #proportional error 
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
        Vc = tvv * exp(η[2]) * (WT / median(etb_data.WT))  #add in the effect of weight on volume

        MTT = tvmtt * exp(η[3]) #define mean transit time
        ktr = 2/MTT #define rate of Transfer

        Vp1 = tvp1 * exp(η[4]) * (WT / median(etb_data.WT))  #add in the effect of weight on volume
        Q1 = tvq1 * exp(η[5]) * (WT / median(etb_data.WT))^0.75 #add in the effect of WT on CL/Q
    end

    @dosecontrol begin
        bioav = (; Depot1 = 1 * exp(η[6]) * (1 + (hivF*HIV))) #add an eta into bioavailbility

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
        #/1000 to appropriately get ng/ml from  ng/L
        """ 
        ETB Concentration (ng/mL)
        """
        LNDV ~ @. Normal(cp, σ_p) #add in error and estimate log conc (LNDV)
    end

end

iparms_WT_HIV_F = (;
tvcl = 50, tvv = 250, tvmtt = 1.5, tvp1 = 800, tvq1 = 40, hivF = 0.4,
Ω = Diagonal([0.05, 0.05, 0.05, 0.05, 0.05, 0.05]), #add in eta for bioavailibility
 σ_p = 0.3)


iparms_WT_HIV_F2 = (;
 tvcl = 51, tvv = 277, tvmtt = 1.89, tvp1 = 722, tvq1 = 41, hivF = -0.1,
 Ω = Diagonal([0.05, 0.05, 0.05, 0.05, 0.05, 0.05]), #add in eta for bioavailibility
  σ_p = 0.3)
 
 #run fits for all populations
fit_WT_HIV_F_day0 = fit(etb_WT_HIV_F, etb_d0_pumapop, iparms_WT_HIV_F2, FOCE())
fit_WT_HIV_F_week6 = fit(etb_WT_HIV_F, etb_w6_pumapop, iparms_WT_HIV_F2, FOCE()) 
fit_WT_HIV_F_all = fit(etb_WT_HIV_F, etb_all_pumapop, iparms_WT_HIV_F, FOCE())
fit_WT_HIV_F_all2 = fit(etb_WT_HIV_F, etb_all_pumapop, iparms_WT_HIV_F2, FOCE()) #true with 2

#inference and parameter ests
infer_WT_HIV_F_day0 = infer(fit_WT_HIV_F_day0) 
infer_WT_HIV_F_week6 = infer(fit_WT_HIV_F_week6) 
infer_WT_HIV_F_all = infer(fit_WT_HIV_F_all) 
#parameters are massivley smaller than initial estimates

#bayesian (Schwarz) information criterion
bic(fit_WT_HIV_F_day0) 
bic(fit_WT_HIV_F_week6) 
bic(fit_WT_HIV_F_all)

#Day 0 VPC
vpc_WT_HIV_F_day0 = vpc(fit_WT_HIV_F_day0; samples = 100)
vpc_plot(vpc_WT_HIV_F_day0)
#week 6 VPC
vpc_WT_HIV_F_week6 = vpc(fit_WT_HIV_F_week6; samples = 100)
vpc_plot(vpc_WT_HIV_F_week6)





#### create 2 comp model ####
etb_WT_HIV_MTT = @model begin

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
        tvmtt ∈ RealDomain(; lower = 0) # bound and estimate absorption
        """
        Typical volume (L) of peripheral 1 compartment
        """
        tvp1 ∈ RealDomain(; lower = 0) 
        """
        Transfer rate from Central to peripheral 1 (L/h)
        """
        tvq1  ∈ RealDomain(; lower = 0) 

        """
        effect of HIV on MTT
        """
        hivMTT  ∈ RealDomain(;) #should we have a lower bound

        """
          - ΩCL
          - ΩVc
          - Ωmtt
          - ΩVp1
          - ΩQ1
          - ΩF
        """
        Ω ∈ PDiagDomain(init = [0.04, 0.04, 0.04, 0.04, 0.04, 0.04]) #initial matrix estimates
        """
        Proportional RUV
        """
        σ_p ∈ RealDomain(; lower = 0.0001) #proportional error 
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
        Vc = tvv * exp(η[2]) * (WT / median(etb_data.WT))  #add in the effect of weight on volume

        MTT = tvmtt * exp(η[3]) * (1 + (hivMTT*HIV))#define mean transit time
        ktr = 2/MTT #define rate of Transfer

        Vp1 = tvp1 * exp(η[4]) * (WT / median(etb_data.WT))  #add in the effect of weight on volume
        Q1 = tvq1 * exp(η[5]) * (WT / median(etb_data.WT))^0.75 #add in the effect of WT on CL/Q
    end

    @dosecontrol begin
        bioav = (; Depot1 = 1 * exp(η[6])) #add an eta into bioavailbility

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
        #/1000 to appropriately get ng/ml from  ng/L
        """ 
        ETB Concentration (ng/mL)
        """
        LNDV ~ @. Normal(cp, σ_p) #add in error and estimate log conc (LNDV)
    end

end

iparms_WT_HIV_MTT = (;
tvcl = 50, tvv = 250, tvmtt = 1.5, tvp1 = 800, tvq1 = 40, hivMTT = 0.4,
Ω = Diagonal([0.05, 0.05, 0.05, 0.05, 0.05, 0.05]), #add in eta for bioavailibility
 σ_p = 0.3)


iparms_WT_HIV_MTT2 = (;
 tvcl = 53, tvv = 286, tvmtt = 1.65, tvp1 = 755, tvq1 = 43, hivMTT = 0.45,
 Ω = Diagonal([0.05, 0.05, 0.05, 0.05, 0.05, 0.05]), #add in eta for bioavailibility
  σ_p = 0.3)
 
 #run fits for all populations
fit_WT_HIV_MTT_day0 = fit(etb_WT_HIV_MTT, etb_d0_pumapop, iparms_WT_HIV_MTT, FOCE())
fit_WT_HIV_MTT_week6 = fit(etb_WT_HIV_MTT, etb_w6_pumapop, iparms_WT_HIV_MTT, FOCE()) 
fit_WT_HIV_MTT_all = fit(etb_WT_HIV_MTT, etb_all_pumapop, iparms_WT_HIV_MTT, FOCE()) #BIC 1864.94012
fit_WT_HIV_MTT_all2 = fit(etb_WT_HIV_MTT, etb_all_pumapop, iparms_WT_HIV_MTT2, FOCE()) #BIC 1864.94019

#inference and parameter ests
infer_WT_HIV_MTT_day0 = infer(fit_WT_HIV_MTT_day0) 
infer_WT_HIV_MTT_week6 = infer(fit_WT_HIV_MTT_week6) 
infer_WT_HIV_MTT_all = infer(fit_WT_HIV_MTT_all) 
#parameters are massivley smaller than initial estimates

#bayesian (Schwarz) information criterion
bic(fit_WT_HIV_MTT_day0) 
bic(fit_WT_HIV_MTT_week6) 
bic(fit_WT_HIV_MTT_all)

#Day 0 VPC
vpc_WT_HIV_MTT_day0 = vpc(fit_WT_HIV_MTT_day0; samples = 100)
vpc_plot(vpc_WT_HIV_MTT_day0)
#week 6 VPC
vpc_WT_HIV_MTT_week6 = vpc(fit_WT_HIV_MTT_week6; samples = 100)
vpc_plot(vpc_WT_HIV_MTT_week6)


#fit improvement on MTT
#positive value for hivMTT
#HIV positive increases MTT
#Ktr = (n+1)/mtt
#increase in MTT results in slower absorption of Ethambutol
#HIV pos patients absorb ethambutol at a slower rate

















###################### ROUND 2 ###############################################################################


#### MTT + VP1 #####
#### create 2 comp model ####
etb_WT_HIV_MTT_VP1 = @model begin

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
        tvmtt ∈ RealDomain(; lower = 0) # bound and estimate absorption
        """
        Typical volume (L) of peripheral 1 compartment
        """
        tvp1 ∈ RealDomain(; lower = 0) 
        """
        Transfer rate from Central to peripheral 1 (L/h)
        """
        tvq1  ∈ RealDomain(; lower = 0) 

        """
        effect of HIV on MTT
        """
        hivMTT  ∈ RealDomain(;) 
        """
        effect of HIV on Peripheral volume
        """
        hivVP1  ∈ RealDomain(;) 

        """
          - ΩCL
          - ΩVc
          - Ωmtt
          - ΩVp1
          - ΩQ1
          - ΩF
        """
        Ω ∈ PDiagDomain(init = [0.04, 0.04, 0.04, 0.04, 0.04, 0.04]) #initial matrix estimates
        """
        Proportional RUV
        """
        σ_p ∈ RealDomain(; lower = 0.0001) #proportional error 
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
        Vc = tvv * exp(η[2]) * (WT / median(etb_data.WT))  #add in the effect of weight on volume

        MTT = tvmtt * exp(η[3]) * (1 + (hivMTT*HIV))#define mean transit time
        ktr = 2/MTT #define rate of Transfer

        Vp1 = tvp1 * exp(η[4]) * (WT / median(etb_data.WT)) * (1+(hivVP1*HIV)) #add in the effect of weight on volume
        Q1 = tvq1 * exp(η[5]) * (WT / median(etb_data.WT))^0.75 #add in the effect of WT on CL/Q
    end

    @dosecontrol begin
        bioav = (; Depot1 = 1 * exp(η[6])) #add an eta into bioavailbility

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
        #/1000 to appropriately get ng/ml from  ng/L
        """ 
        ETB Concentration (ng/mL)
        """
        LNDV ~ @. Normal(cp, σ_p) #add in error and estimate log conc (LNDV)
    end

end

iparms_WT_HIV_MTT_VP1 = (;
tvcl = 50, tvv = 250, tvmtt = 1.5, tvp1 = 800, tvq1 = 40, hivMTT = 0.4, hivVP1 = 0.4,
Ω = Diagonal([0.05, 0.05, 0.05, 0.05, 0.05, 0.05]), #add in eta for bioavailibility
 σ_p = 0.3)

iparms_WT_HIV_MTT_VP12 = (;
tvcl = 54, tvv = 292, tvmtt = 1.65, tvp1 = 997, tvq1 = 43.5, hivMTT = 0.42, hivVP1 = -0.48,
Ω = Diagonal([0.05, 0.05, 0.05, 0.05, 0.05, 0.05]), #add in eta for bioavailibility
 σ_p = 0.3)

 #run fits for all populations
fit_WT_HIV_MTT_VP1_day0 = fit(etb_WT_HIV_MTT_VP1, etb_d0_pumapop, iparms_WT_HIV_MTT_VP12, FOCE())
fit_WT_HIV_MTT_VP1_week6 = fit(etb_WT_HIV_MTT_VP1, etb_w6_pumapop, iparms_WT_HIV_MTT_VP12, FOCE()) #True
fit_WT_HIV_MTT_VP1_all = fit(etb_WT_HIV_MTT_VP1, etb_all_pumapop, iparms_WT_HIV_MTT_VP1, FOCE()) #BIC 1863.3560
fit_WT_HIV_MTT_VP1_all2 = fit(etb_WT_HIV_MTT_VP1, etb_all_pumapop, iparms_WT_HIV_MTT_VP12, FOCE()) #BIC 1863.3047


#inference and parameter ests
infer_WT_HIV_MTT_VP1_day0 = infer(fit_WT_HIV_MTT_VP1_day0) 
infer_WT_HIV_MTT_VP1_week6 = infer(fit_WT_HIV_MTT_VP1_week6) 
infer_WT_HIV_MTT_VP1_all = infer(fit_WT_HIV_MTT_VP1_all) 
#parameters are massivley smaller than initial estimates

#bayesian (Schwarz) information criterion
bic(fit_WT_HIV_MTT_VP1_day0) 
bic(fit_WT_HIV_MTT_VP1_week6) 
bic(fit_WT_HIV_MTT_VP1_all)
bic(fit_WT_HIV_MTT_VP1_all2)
#Day 0 VPC
vpc_WT_HIV_MTT_VP1_day0 = vpc(fit_WT_HIV_MTT_VP1_day0; samples = 100)
vpc_plot(vpc_WT_HIV_MTT_VP1_day0)
#week 6 VPC
vpc_WT_HIV_MTT_VP1_week6 = vpc(fit_WT_HIV_MTT_VP1_week6; samples = 100)
vpc_plot(vpc_WT_HIV_MTT_VP1_week6)

#TAD VPC
vpc_WT_HIV_MTT_VP1_all = vpc(fit_WT_HIV_MTT_VP1_all; samples = 100, covariates = [:tad])
vpc_plot(vpc_WT_HIV_MTT_VP1_all)


####HIV MTT + V ####
etb_WT_HIV_MTT_V = @model begin

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
        tvmtt ∈ RealDomain(; lower = 0) # bound and estimate absorption
        """
        Typical volume (L) of peripheral 1 compartment
        """
        tvp1 ∈ RealDomain(; lower = 0) 
        """
        Transfer rate from Central to peripheral 1 (L/h)
        """
        tvq1  ∈ RealDomain(; lower = 0) 

        """
        effect of HIV on MTT
        """
        hivMTT  ∈ RealDomain(;) 
        """
        effect of HIV on central volume
        """
        hivV  ∈ RealDomain(;) 

        """
          - ΩCL
          - ΩVc
          - Ωmtt
          - ΩVp1
          - ΩQ1
          - ΩF
        """
        Ω ∈ PDiagDomain(init = [0.04, 0.04, 0.04, 0.04, 0.04, 0.04]) #initial matrix estimates
        """
        Proportional RUV
        """
        σ_p ∈ RealDomain(; lower = 0.0001) #proportional error 
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
        Vc = tvv * exp(η[2]) * (WT / median(etb_data.WT)) * (1+(hivV*HIV))  #add in the effect of weight on volume

        MTT = tvmtt * exp(η[3]) * (1 + (hivMTT*HIV))#define mean transit time
        ktr = 2/MTT #define rate of Transfer

        Vp1 = tvp1 * exp(η[4]) * (WT / median(etb_data.WT)) #add in the effect of weight on volume
        Q1 = tvq1 * exp(η[5]) * (WT / median(etb_data.WT))^0.75 #add in the effect of WT on CL/Q
    end

    @dosecontrol begin
        bioav = (; Depot1 = 1 * exp(η[6])) #add an eta into bioavailbility

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
        #/1000 to appropriately get ng/ml from  ng/L
        """ 
        ETB Concentration (ng/mL)
        """
        LNDV ~ @. Normal(cp, σ_p) #add in error and estimate log conc (LNDV)
    end

end

iparms_WT_HIV_MTT_V = (;
tvcl = 50, tvv = 250, tvmtt = 1.5, tvp1 = 800, tvq1 = 40, hivMTT = 0.4, hivV = 0.4,
Ω = Diagonal([0.05, 0.05, 0.05, 0.05, 0.05, 0.05]), #add in eta for bioavailibility
σ_p = 0.3)

iparms_WT_HIV_MTT_V2 = (;
tvcl = 53.5, tvv = 277.4, tvmtt = 1.68, tvp1 = 748, tvq1 = 43, hivMTT = 0.35, hivV = 0.15,
Ω = Diagonal([0.05, 0.05, 0.05, 0.05, 0.05, 0.05]), #add in eta for bioavailibility
σ_p = 0.3)

 #run fits for all populations
fit_WT_HIV_MTT_V_day0 = fit(etb_WT_HIV_MTT_V, etb_d0_pumapop, iparms_WT_HIV_MTT_V, FOCE())
fit_WT_HIV_MTT_V_week6 = fit(etb_WT_HIV_MTT_V, etb_w6_pumapop, iparms_WT_HIV_MTT_V, FOCE()) #True
fit_WT_HIV_MTT_V_all = fit(etb_WT_HIV_MTT_V, etb_all_pumapop, iparms_WT_HIV_MTT_V, FOCE()) #BIC 1870.966
fit_WT_HIV_MTT_V_all2 = fit(etb_WT_HIV_MTT_V, etb_all_pumapop, iparms_WT_HIV_MTT_V2, FOCE()) #BIC 1871.042


#inference and parameter ests
infer_WT_HIV_MTT_V_day0 = infer(fit_WT_HIV_MTT_V_day0) 
infer_WT_HIV_MTT_V_week6 = infer(fit_WT_HIV_MTT_V_week6) 
infer_WT_HIV_MTT_V_all = infer(fit_WT_HIV_MTT_V_all) 
#parameters are massivley smaller than initial estimates

#bayesian (Schwarz) information criterion
bic(fit_WT_HIV_MTT_V_day0) 
bic(fit_WT_HIV_MTT_V_week6) 
bic(fit_WT_HIV_MTT_V_all)

#Day 0 VPC
vpc_WT_HIV_MTT_V_day0 = vpc(fit_WT_HIV_MTT_V_day0; samples = 100)
vpc_plot(vpc_WT_HIV_MTT_V_day0)
#week 6 VPC
vpc_WT_HIV_MTT_V_week6 = vpc(fit_WT_HIV_MTT_V_week6; samples = 100)
vpc_plot(vpc_WT_HIV_MTT_V_week6)




##### HIV MTT + CL #####
etb_WT_HIV_MTT_CL = @model begin

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
        tvmtt ∈ RealDomain(; lower = 0) # bound and estimate absorption
        """
        Typical volume (L) of peripheral 1 compartment
        """
        tvp1 ∈ RealDomain(; lower = 0) 
        """
        Transfer rate from Central to peripheral 1 (L/h)
        """
        tvq1  ∈ RealDomain(; lower = 0) 

        """
        effect of HIV on MTT
        """
        hivMTT  ∈ RealDomain(;) 
        """
        effect of HIV on clearance
        """
        hivCL  ∈ RealDomain(;) 

        """
          - ΩCL
          - ΩVc
          - Ωmtt
          - ΩVp1
          - ΩQ1
          - ΩF
        """
        Ω ∈ PDiagDomain(init = [0.04, 0.04, 0.04, 0.04, 0.04, 0.04]) #initial matrix estimates
        """
        Proportional RUV
        """
        σ_p ∈ RealDomain(; lower = 0.0001) #proportional error 
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
        CL = tvcl * exp(η[1]) * (WT / median(etb_data.WT))^0.75 * (1+(hivCL*HIV))  #add in the effect of WT on CL/Q
        Vc = tvv * exp(η[2]) * (WT / median(etb_data.WT))  #add in the effect of weight on volume

        MTT = tvmtt * exp(η[3]) * (1 + (hivMTT*HIV))#define mean transit time
        ktr = 2/MTT #define rate of Transfer

        Vp1 = tvp1 * exp(η[4]) * (WT / median(etb_data.WT)) #add in the effect of weight on volume
        Q1 = tvq1 * exp(η[5]) * (WT / median(etb_data.WT))^0.75 #add in the effect of WT on CL/Q
    end

    @dosecontrol begin
        bioav = (; Depot1 = 1 * exp(η[6])) #add an eta into bioavailbility

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
        #/1000 to appropriately get ng/ml from  ng/L
        """ 
        ETB Concentration (ng/mL)
        """
        LNDV ~ @. Normal(cp, σ_p) #add in error and estimate log conc (LNDV)
    end

end

iparms_WT_HIV_MTT_CL = (;
tvcl = 50, tvv = 250, tvmtt = 1.5, tvp1 = 800, tvq1 = 40, hivMTT = 0.4, hivCL = 0.4,
Ω = Diagonal([0.05, 0.05, 0.05, 0.05, 0.05, 0.05]), #add in eta for bioavailibility
 σ_p = 0.3)

iparms_WT_HIV_MTT_CL2 = (;
tvcl = 53, tvv = 286, tvmtt = 1.6, tvp1 = 755, tvq1 = 43, hivMTT = 0.45, hivCL = 0.016,
Ω = Diagonal([0.05, 0.05, 0.05, 0.05, 0.05, 0.05]), #add in eta for bioavailibility
σ_p = 0.3)
 
 #run fits for all populations
fit_WT_HIV_MTT_CL_day0 = fit(etb_WT_HIV_MTT_CL, etb_d0_pumapop, iparms_WT_HIV_MTT_CL2, FOCE())
fit_WT_HIV_MTT_CL_week6 = fit(etb_WT_HIV_MTT_CL, etb_w6_pumapop, iparms_WT_HIV_MTT_CL2, FOCE()) 
fit_WT_HIV_MTT_CL_all = fit(etb_WT_HIV_MTT_CL, etb_all_pumapop, iparms_WT_HIV_MTT_CL, FOCE()) #BIC 1872.605
fit_WT_HIV_MTT_CL_all2 = fit(etb_WT_HIV_MTT_CL, etb_all_pumapop, iparms_WT_HIV_MTT_CL2, FOCE()) #TRUE BIC 1872.198 

#inference and parameter ests
infer_WT_HIV_MTT_CL_day0 = infer(fit_WT_HIV_MTT_CL_day0) 
infer_WT_HIV_MTT_CL_week6 = infer(fit_WT_HIV_MTT_CL_week6) 
infer_WT_HIV_MTT_CL_all = infer(fit_WT_HIV_MTT_CL_all2) 
#parameters are massivley smaller than initial estimates

#bayesian (Schwarz) information criterion
bic(fit_WT_HIV_MTT_CL_day0) 
bic(fit_WT_HIV_MTT_CL_week6) 
bic(fit_WT_HIV_MTT_CL_all2)

#Day 0 VPC
vpc_WT_HIV_MTT_CL_day0 = vpc(fit_WT_HIV_MTT_CL_day0; samples = 100)
vpc_plot(vpc_WT_HIV_MTT_CL_day0)
#week 6 VPC
vpc_WT_HIV_MTT_CL_week6 = vpc(fit_WT_HIV_MTT_CL_week6; samples = 100)
vpc_plot(vpc_WT_HIV_MTT_CL_week6)







##### HIV MTT + Q #####
etb_WT_HIV_MTT_Q = @model begin

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
        tvmtt ∈ RealDomain(; lower = 0) # bound and estimate absorption
        """
        Typical volume (L) of peripheral 1 compartment
        """
        tvp1 ∈ RealDomain(; lower = 0) 
        """
        Transfer rate from Central to peripheral 1 (L/h)
        """
        tvq1  ∈ RealDomain(; lower = 0) 

        """
        effect of HIV on MTT
        """
        hivMTT  ∈ RealDomain(;) 
        """
        effect of HIV on intercompartmental clearance
        """
        hivQ  ∈ RealDomain(;) 

        """
          - ΩCL
          - ΩVc
          - Ωmtt
          - ΩVp1
          - ΩQ1
          - ΩF
        """
        Ω ∈ PDiagDomain(init = [0.04, 0.04, 0.04, 0.04, 0.04, 0.04]) #initial matrix estimates
        """
        Proportional RUV
        """
        σ_p ∈ RealDomain(; lower = 0.0001) #proportional error 
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
        CL = tvcl * exp(η[1]) * (WT / median(etb_data.WT))^0.75  #add in the effect of WT on CL/Q
        Vc = tvv * exp(η[2]) * (WT / median(etb_data.WT))  #add in the effect of weight on volume

        MTT = tvmtt * exp(η[3]) * (1 + (hivMTT*HIV))#define mean transit time
        ktr = 2/MTT #define rate of Transfer

        Vp1 = tvp1 * exp(η[4]) * (WT / median(etb_data.WT)) #add in the effect of weight on volume
        Q1 = tvq1 * exp(η[5]) * (WT / median(etb_data.WT))^0.75 * (1+(hivQ*HIV))#add in the effect of WT on CL/Q
    end

    @dosecontrol begin
        bioav = (; Depot1 = 1 * exp(η[6])) #add an eta into bioavailbility

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
        #/1000 to appropriately get ng/ml from  ng/L
        """ 
        ETB Concentration (ng/mL)
        """
        LNDV ~ @. Normal(cp, σ_p) #add in error and estimate log conc (LNDV)
    end

end

iparms_WT_HIV_MTT_Q = (;
tvcl = 40, tvv = 200, tvmtt = 1.3, tvp1 = 700, tvq1 = 30, hivMTT = 0.4, hivQ = 0.4,
Ω = Diagonal([0.05, 0.05, 0.05, 0.05, 0.05, 0.05]), #add in eta for bioavailibility
 σ_p = 0.3)

iparms_WT_HIV_MTT_Q2 = (;
tvcl = 53.5, tvv = 293, tvmtt = 1.65, tvp1 = 754, tvq1 = 40, hivMTT = 0.47, hivQ = 0.3,
Ω = Diagonal([0.05, 0.05, 0.05, 0.05, 0.05, 0.05]), #add in eta for bioavailibility
σ_p = 0.3)

 #run fits for all populations
fit_WT_HIV_MTT_Q_day0 = fit(etb_WT_HIV_MTT_Q, etb_d0_pumapop, iparms_WT_HIV_MTT_Q, FOCE())
fit_WT_HIV_MTT_Q_week6 = fit(etb_WT_HIV_MTT_Q, etb_w6_pumapop, iparms_WT_HIV_MTT_Q, FOCE()) 
fit_WT_HIV_MTT_Q_all = fit(etb_WT_HIV_MTT_Q, etb_all_pumapop, iparms_WT_HIV_MTT_Q, FOCE()) #BIC 1868.00800
fit_WT_HIV_MTT_Q_all2 = fit(etb_WT_HIV_MTT_Q, etb_all_pumapop, iparms_WT_HIV_MTT_Q2, FOCE()) #BIC 1868.00829

#inference and parameter ests
infer_WT_HIV_MTT_Q_day0 = infer(fit_WT_HIV_MTT_Q_day0) 
infer_WT_HIV_MTT_Q_week6 = infer(fit_WT_HIV_MTT_Q_week6) 
infer_WT_HIV_MTT_Q_all = infer(fit_WT_HIV_MTT_Q_all) 
#parameters are massivley smaller than initial estimates

#bayesian (Schwarz) information criterion
bic(fit_WT_HIV_MTT_Q_day0) 
bic(fit_WT_HIV_MTT_Q_week6) 
bic(fit_WT_HIV_MTT_Q_all2)

#Day 0 VPC
vpc_WT_HIV_MTT_Q_day0 = vpc(fit_WT_HIV_MTT_Q_day0; samples = 100)
vpc_plot(vpc_WT_HIV_MTT_Q_day0)
#week 6 VPC
vpc_WT_HIV_MTT_Q_week6 = vpc(fit_WT_HIV_MTT_Q_week6; samples = 100)
vpc_plot(vpc_WT_HIV_MTT_Q_week6)




##### HIV MTT + F #####
etb_WT_HIV_MTT_F = @model begin

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
        tvmtt ∈ RealDomain(; lower = 0) # bound and estimate absorption
        """
        Typical volume (L) of peripheral 1 compartment
        """
        tvp1 ∈ RealDomain(; lower = 0) 
        """
        Transfer rate from Central to peripheral 1 (L/h)
        """
        tvq1  ∈ RealDomain(; lower = 0) 

        """
        effect of HIV on MTT
        """
        hivMTT  ∈ RealDomain(;) 
        """
        effect of HIV on bioavailability
        """
        hivF  ∈ RealDomain(;) 

        """
          - ΩCL
          - ΩVc
          - Ωmtt
          - ΩVp1
          - ΩQ1
          - ΩF
        """
        Ω ∈ PDiagDomain(init = [0.04, 0.04, 0.04, 0.04, 0.04, 0.04]) #initial matrix estimates
        """
        Proportional RUV
        """
        σ_p ∈ RealDomain(; lower = 0.0001) #proportional error 
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
        CL = tvcl * exp(η[1]) * (WT / median(etb_data.WT))^0.75  #add in the effect of WT on CL/Q
        Vc = tvv * exp(η[2]) * (WT / median(etb_data.WT))  #add in the effect of weight on volume

        MTT = tvmtt * exp(η[3]) * (1 + (hivMTT*HIV))#define mean transit time
        ktr = 2/MTT #define rate of Transfer

        Vp1 = tvp1 * exp(η[4]) * (WT / median(etb_data.WT)) #add in the effect of weight on volume
        Q1 = tvq1 * exp(η[5]) * (WT / median(etb_data.WT))^0.75#add in the effect of WT on CL/Q
    end

    @dosecontrol begin
        bioav = (; Depot1 = 1  * (1+(hivF*HIV)) * exp(η[6])) #add an eta into bioavailbility

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
        #/1000 to appropriately get ng/ml from  ng/L
        """ 
        ETB Concentration (ng/mL)
        """
        LNDV ~ @. Normal(cp, σ_p) #add in error and estimate log conc (LNDV)
    end

end

iparms_WT_HIV_MTT_F = (;
tvcl = 50, tvv = 250, tvmtt = 1.5, tvp1 = 800, tvq1 = 40, hivMTT = 0.4, hivF = 0.4,
Ω = Diagonal([0.05, 0.05, 0.05, 0.05, 0.05, 0.05]), #add in eta for bioavailibility
 σ_p = 0.3)

 iparms_WT_HIV_MTT_F2 = (;
 tvcl = 52, tvv = 278, tvmtt = 1.65, tvp1 = 730, tvq1 = 42, hivMTT = 0.42, hivF = -0.1,
 Ω = Diagonal([0.05, 0.05, 0.05, 0.05, 0.05, 0.05]), #add in eta for bioavailibility
  σ_p = 0.3)
 #run fits for all populations
fit_WT_HIV_MTT_F_day0 = fit(etb_WT_HIV_MTT_F, etb_d0_pumapop, iparms_WT_HIV_MTT_F, FOCE())
fit_WT_HIV_MTT_F_week6 = fit(etb_WT_HIV_MTT_F, etb_w6_pumapop, iparms_WT_HIV_MTT_F, FOCE()) 
fit_WT_HIV_MTT_F_all = fit(etb_WT_HIV_MTT_F, etb_all_pumapop, iparms_WT_HIV_MTT_F, FOCE()) #BIC 1871.225
fit_WT_HIV_MTT_F_all2 = fit(etb_WT_HIV_MTT_F, etb_all_pumapop, iparms_WT_HIV_MTT_F2, FOCE()) #BIC 1871.230

#inference and parameter ests
infer_WT_HIV_MTT_F_day0 = infer(fit_WT_HIV_MTT_F_day0) 
infer_WT_HIV_MTT_F_week6 = infer(fit_WT_HIV_MTT_F_week6) 
infer_WT_HIV_MTT_F_all = infer(fit_WT_HIV_MTT_F_all) 
#parameters are massivley smaller than initial estimates

#bayesian (Schwarz) information criterion
bic(fit_WT_HIV_MTT_F_day0) 
bic(fit_WT_HIV_MTT_F_week6) 
bic(fit_WT_HIV_MTT_F_all2)

#Day 0 VPC
vpc_WT_HIV_MTT_F_day0 = vpc(fit_WT_HIV_MTT_F_day0; samples = 100)
vpc_plot(vpc_WT_HIV_MTT_F_day0)
#week 6 VPC
vpc_WT_HIV_MTT_F_week6 = vpc(fit_WT_HIV_MTT_F_week6; samples = 100)
vpc_plot(vpc_WT_HIV_MTT_F_week6)







##################################### WT + MTT + VP1 #############################################################

#### WT + MTT + VP1  == WMV #####
#### add V ####
etb_WMV_V = @model begin

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
        tvmtt ∈ RealDomain(; lower = 0) # bound and estimate absorption
        """
        Typical volume (L) of peripheral 1 compartment
        """
        tvp1 ∈ RealDomain(; lower = 0) 
        """
        Transfer rate from Central to peripheral 1 (L/h)
        """
        tvq1  ∈ RealDomain(; lower = 0) 

        """
        effect of HIV on MTT
        """
        hivMTT  ∈ RealDomain(;) 
        """
        effect of HIV on Peripheral volume
        """
        hivVP1  ∈ RealDomain(;) 

        """
        effect of HIV on Central volume
        """
        hivV ∈ RealDomain(;) 

        """
          - ΩCL
          - ΩVc
          - Ωmtt
          - ΩVp1
          - ΩQ1
          - ΩF
        """
        Ω ∈ PDiagDomain(init = [0.04, 0.04, 0.04, 0.04, 0.04, 0.04]) #initial matrix estimates
        """
        Proportional RUV
        """
        σ_p ∈ RealDomain(; lower = 0.0001) #proportional error 
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
        Vc = tvv * exp(η[2]) * (WT / median(etb_data.WT)) * (1+(hivV*HIV))  #add in the effect of weight on volume

        MTT = tvmtt * exp(η[3]) * (1 + (hivMTT*HIV))#define mean transit time
        ktr = 2/MTT #define rate of Transfer

        Vp1 = tvp1 * exp(η[4]) * (WT / median(etb_data.WT)) * (1+(hivVP1*HIV)) #add in the effect of weight on volume
        Q1 = tvq1 * exp(η[5]) * (WT / median(etb_data.WT))^0.75 #add in the effect of WT on CL/Q
    end

    @dosecontrol begin
        bioav = (; Depot1 = 1 * exp(η[6])) #add an eta into bioavailbility

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
        #/1000 to appropriately get ng/ml from  ng/L
        """ 
        ETB Concentration (ng/mL)
        """
        LNDV ~ @. Normal(cp, σ_p) #add in error and estimate log conc (LNDV)
    end

end

iparms_WMV_V = (;
tvcl = 50, tvv = 250, tvmtt = 1.5, tvp1 = 800, tvq1 = 40, hivMTT = 0.4, hivVP1 = 0.4, hivV = 0.4,
Ω = Diagonal([0.05, 0.05, 0.05, 0.05, 0.05, 0.05]), #add in eta for bioavailibility
 σ_p = 0.3)

 
iparms_WMV_V2 = (;
tvcl = 53, tvv = 286, tvmtt = 1.65, tvp1 = 981, tvq1 = 43, hivMTT = 0.4, hivVP1 = 0.4, hivV = 0.4,
Ω = Diagonal([0.05, 0.05, 0.05, 0.05, 0.05, 0.05]), #add in eta for bioavailibility
 σ_p = 0.3)

 #run fits for all populations
fit_WMV_V_day0 = fit(etb_WMV_V, etb_d0_pumapop, iparms_WMV_V2, FOCE())
fit_WMV_V_week6 = fit(etb_WMV_V, etb_w6_pumapop, iparms_WMV_V2, FOCE()) 
fit_WMV_V_all = fit(etb_WMV_V, etb_all_pumapop, iparms_WMV_V, FOCE())
fit_WMV_V_all2 = fit(etb_WMV_V, etb_all_pumapop, iparms_WMV_V2, FOCE()) #True, use 2

#inference and parameter ests
infer_WMV_V_day0 = infer(fit_WMV_V_day0) 
infer_WMV_V_week6 = infer(fit_WMV_V_week6) 
infer_WMV_V_all = infer(fit_WMV_V_all) 
#parameters are massivley smaller than initial estimates

#bayesian (Schwarz) information criterion
bic(fit_WMV_V_day0) 
bic(fit_WMV_V_week6) 
bic(fit_WMV_V_all)

#Day 0 VPC
vpc_WMV_V_day0 = vpc(fit_WMV_V_day0; samples = 100)
vpc_plot(vpc_WMV_V_day0)
#week 6 VPC
vpc_WMV_V_week6 = vpc(fit_WMV_V_week6; samples = 100)
vpc_plot(vpc_WMV_V_week6)





#### add CL ####
etb_WMV_CL = @model begin

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
        tvmtt ∈ RealDomain(; lower = 0) # bound and estimate absorption
        """
        Typical volume (L) of peripheral 1 compartment
        """
        tvp1 ∈ RealDomain(; lower = 0) 
        """
        Transfer rate from Central to peripheral 1 (L/h)
        """
        tvq1  ∈ RealDomain(; lower = 0) 

        """
        effect of HIV on MTT
        """
        hivMTT  ∈ RealDomain(;) 
        """
        effect of HIV on Peripheral volume
        """
        hivVP1  ∈ RealDomain(;) 

        """
        effect of HIV on clearance
        """
        hivCL ∈ RealDomain(;) 

        """
          - ΩCL
          - ΩVc
          - Ωmtt
          - ΩVp1
          - ΩQ1
          - ΩF
        """
        Ω ∈ PDiagDomain(init = [0.04, 0.04, 0.04, 0.04, 0.04, 0.04]) #initial matrix estimates
        """
        Proportional RUV
        """
        σ_p ∈ RealDomain(; lower = 0.0001) #proportional error 
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
        CL = tvcl * exp(η[1]) * (WT / median(etb_data.WT))^0.75 * (1+(hivCL*HIV)) #add in the effect of WT on CL/Q
        Vc = tvv * exp(η[2]) * (WT / median(etb_data.WT))  #add in the effect of weight on volume

        MTT = tvmtt * exp(η[3]) * (1 + (hivMTT*HIV))#define mean transit time
        ktr = 2/MTT #define rate of Transfer

        Vp1 = tvp1 * exp(η[4]) * (WT / median(etb_data.WT)) * (1+(hivVP1*HIV)) #add in the effect of weight on volume
        Q1 = tvq1 * exp(η[5]) * (WT / median(etb_data.WT))^0.75 #add in the effect of WT on CL/Q
    end

    @dosecontrol begin
        bioav = (; Depot1 = 1 * exp(η[6])) #add an eta into bioavailbility

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
        #/1000 to appropriately get ng/ml from  ng/L
        """ 
        ETB Concentration (ng/mL)
        """
        LNDV ~ @. Normal(cp, σ_p) #add in error and estimate log conc (LNDV)
    end

end

iparms_WMV_CL = (;
tvcl = 50, tvv = 250, tvmtt = 1.5, tvp1 = 800, tvq1 = 40, hivMTT = 0.4, hivVP1 = 0.4, hivCL = 0.4,
Ω = Diagonal([0.05, 0.05, 0.05, 0.05, 0.05, 0.05]), #add in eta for bioavailibility
 σ_p = 0.3)

iparms_WMV_CL2 = (;
tvcl = 53, tvv = 292, tvmtt = 1.65, tvp1 = 998, tvq1 = 43, hivMTT = 0.4, hivVP1 = -0.5, hivCL = 0.1,
Ω = Diagonal([0.05, 0.05, 0.05, 0.05, 0.05, 0.05]), #add in eta for bioavailibility
σ_p = 0.3)
 

 #run fits for all populations
fit_WMV_CL_day0 = fit(etb_WMV_CL, etb_d0_pumapop, iparms_WMV_CL2, FOCE())
fit_WMV_CL_week6 = fit(etb_WMV_CL, etb_w6_pumapop, iparms_WMV_CL2, FOCE()) 
fit_WMV_CL_all = fit(etb_WMV_CL, etb_all_pumapop, iparms_WMV_CL, FOCE()) #BIC 1870.637
fit_WMV_CL_all2 = fit(etb_WMV_CL, etb_all_pumapop, iparms_WMV_CL2, FOCE()) #BIC 1870.582

#inference and parameter ests
infer_WMV_CL_day0 = infer(fit_WMV_CL_day0) 
infer_WMV_CL_week6 = infer(fit_WMV_CL_week6) 
infer_WMV_CL_all = infer(fit_WMV_CL_all) 
#parameters are massivley smaller than initial estimates

#bayesian (Schwarz) information criterion
bic(fit_WMV_CL_day0) 
bic(fit_WMV_CL_week6) 
bic(fit_WMV_CL_all2)

#Day 0 VPC
vpc_WMV_CL_day0 = vpc(fit_WMV_CL_day0; samples = 100)
vpc_plot(vpc_WMV_CL_day0)
#week 6 VPC
vpc_WMV_CL_week6 = vpc(fit_WMV_CL_week6; samples = 100)
vpc_plot(vpc_WMV_CL_week6)




#### add Q ####
etb_WMV_Q = @model begin

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
        tvmtt ∈ RealDomain(; lower = 0) # bound and estimate absorption
        """
        Typical volume (L) of peripheral 1 compartment
        """
        tvp1 ∈ RealDomain(; lower = 0) 
        """
        Transfer rate from Central to peripheral 1 (L/h)
        """
        tvq1  ∈ RealDomain(; lower = 0) 

        """
        effect of HIV on MTT
        """
        hivMTT  ∈ RealDomain(;) 
        """
        effect of HIV on Peripheral volume
        """
        hivVP1  ∈ RealDomain(;) 

        """
        effect of HIV on intercompartmental clearance
        """
        hivQ ∈ RealDomain(;) 

        """
          - ΩCL
          - ΩVc
          - Ωmtt
          - ΩVp1
          - ΩQ1
          - ΩF
        """
        Ω ∈ PDiagDomain(init = [0.04, 0.04, 0.04, 0.04, 0.04, 0.04]) #initial matrix estimates
        """
        Proportional RUV
        """
        σ_p ∈ RealDomain(; lower = 0.0001) #proportional error 
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
        Vc = tvv * exp(η[2]) * (WT / median(etb_data.WT))  #add in the effect of weight on volume

        MTT = tvmtt * exp(η[3]) * (1 + (hivMTT*HIV))#define mean transit time
        ktr = 2/MTT #define rate of Transfer

        Vp1 = tvp1 * exp(η[4]) * (WT / median(etb_data.WT)) * (1+(hivVP1*HIV)) #add in the effect of weight on volume
        Q1 = tvq1 * exp(η[5]) * (WT / median(etb_data.WT))^0.75 * (1+(hivQ*HIV)) #add in the effect of WT on CL/Q
    end

    @dosecontrol begin
        bioav = (; Depot1 = 1 * exp(η[6])) #add an eta into bioavailbility

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
        #/1000 to appropriately get ng/ml from  ng/L
        """ 
        ETB Concentration (ng/mL)
        """
        LNDV ~ @. Normal(cp, σ_p) #add in error and estimate log conc (LNDV)
    end

end

iparms_WMV_Q = (;
tvcl = 50, tvv = 250, tvmtt = 1.5, tvp1 = 800, tvq1 = 40, hivMTT = 0.4, hivVP1 = 0.4, hivQ = 0.4,
Ω = Diagonal([0.05, 0.05, 0.05, 0.05, 0.05, 0.05]), #add in eta for bioavailibility
 σ_p = 0.3)

iparms_WMV_Q2 = (;
tvcl = 53, tvv = 289, tvmtt = 1.65, tvp1 = 973, tvq1 = 40, hivMTT = 0.4, hivVP1 = -0.5, hivQ = 0.4,
Ω = Diagonal([0.05, 0.05, 0.05, 0.05, 0.05, 0.05]), #add in eta for bioavailibility
σ_p = 0.3)

 #run fits for all populations
fit_WMV_Q_day0 = fit(etb_WMV_Q, etb_d0_pumapop, iparms_WMV_Q2, FOCE())
fit_WMV_Q_week6 = fit(etb_WMV_Q, etb_w6_pumapop, iparms_WMV_Q2, FOCE()) 
fit_WMV_Q_all = fit(etb_WMV_Q, etb_all_pumapop, iparms_WMV_Q, FOCE()) #BIC 1867.971
fit_WMV_Q_all2 = fit(etb_WMV_Q, etb_all_pumapop, iparms_WMV_Q2, FOCE()) #BIC 1867.968

#inference and parameter ests
infer_WMV_Q_day0 = infer(fit_WMV_Q_day0) 
infer_WMV_Q_week6 = infer(fit_WMV_Q_week6) 
infer_WMV_Q_all = infer(fit_WMV_Q_all) 
#parameters are massivley smaller than initial estimates

#bayesian (Schwarz) information criterion
bic(fit_WMV_Q_day0) 
bic(fit_WMV_Q_week6) 
bic(fit_WMV_Q_all2)

#Day 0 VPC
vpc_WMV_Q_day0 = vpc(fit_WMV_Q_day0; samples = 100)
vpc_plot(vpc_WMV_Q_day0)
#week 6 VPC
vpc_WMV_Q_week6 = vpc(fit_WMV_Q_week6; samples = 100)
vpc_plot(vpc_WMV_Q_week6)



#### add F ####
etb_WMV_F = @model begin

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
        tvmtt ∈ RealDomain(; lower = 0) # bound and estimate absorption
        """
        Typical volume (L) of peripheral 1 compartment
        """
        tvp1 ∈ RealDomain(; lower = 0) 
        """
        Transfer rate from Central to peripheral 1 (L/h)
        """
        tvq1  ∈ RealDomain(; lower = 0) 

        """
        effect of HIV on MTT
        """
        hivMTT  ∈ RealDomain(;) 
        """
        effect of HIV on Peripheral volume
        """
        hivVP1  ∈ RealDomain(;) 

        """
        effect of HIV on bioavailability
        """
        hivF ∈ RealDomain(;) 

        """
          - ΩCL
          - ΩVc
          - Ωmtt
          - ΩVp1
          - ΩQ1
          - ΩF
        """
        Ω ∈ PDiagDomain(init = [0.04, 0.04, 0.04, 0.04, 0.04, 0.04]) #initial matrix estimates
        """
        Proportional RUV
        """
        σ_p ∈ RealDomain(; lower = 0.0001) #proportional error 
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
        Vc = tvv * exp(η[2]) * (WT / median(etb_data.WT))  #add in the effect of weight on volume

        MTT = tvmtt * exp(η[3]) * (1 + (hivMTT*HIV))#define mean transit time
        ktr = 2/MTT #define rate of Transfer

        Vp1 = tvp1 * exp(η[4]) * (WT / median(etb_data.WT)) * (1+(hivVP1*HIV)) #add in the effect of weight on volume
        Q1 = tvq1 * exp(η[5]) * (WT / median(etb_data.WT))^0.75 #add in the effect of WT on CL/Q
    end

    @dosecontrol begin
        bioav = (; Depot1 = 1  * (1+(hivF*HIV)) * exp(η[6])) #add an eta into bioavailbility

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
        #/1000 to appropriately get ng/ml from  ng/L
        """ 
        ETB Concentration (ng/mL)
        """
        LNDV ~ @. Normal(cp, σ_p) #add in error and estimate log conc (LNDV)
    end

end

iparms_WMV_F = (;
tvcl = 50, tvv = 250, tvmtt = 1.5, tvp1 = 800, tvq1 = 40, hivMTT = 0.4, hivVP1 = 0.4, hivF = 0.4,
Ω = Diagonal([0.05, 0.05, 0.05, 0.05, 0.05, 0.05]), #add in eta for bioavailibility
 σ_p = 0.3)


iparms_WMV_F2 = (;
 tvcl = 51, tvv = 280, tvmtt = 1.65, tvp1 = 971, tvq1 = 42, hivMTT = 0.4, hivVP1 = -0.5, hivF = 0.1,
 Ω = Diagonal([0.05, 0.05, 0.05, 0.05, 0.05, 0.05]), #add in eta for bioavailibility
  σ_p = 0.3)
 
 #run fits for all populations
fit_WMV_F_day0 = fit(etb_WMV_F, etb_d0_pumapop, iparms_WMV_F2, FOCE())
fit_WMV_F_week6 = fit(etb_WMV_F, etb_w6_pumapop, iparms_WMV_F2, FOCE()) 
fit_WMV_F_all = fit(etb_WMV_F, etb_all_pumapop, iparms_WMV_F, FOCE()) #BIC 1868.373
fit_WMV_F_all2 = fit(etb_WMV_F, etb_all_pumapop, iparms_WMV_F2, FOCE()) #BIC 1868.360

#inference and parameter ests
infer_WMV_F_day0 = infer(fit_WMV_F_day0) 
infer_WMV_F_week6 = infer(fit_WMV_F_week6) 
infer_WMV_F_all = infer(fit_WMV_F_all) 
#parameters are massivley smaller than initial estimates

#bayesian (Schwarz) information criterion
bic(fit_WMV_F_day0) 
bic(fit_WMV_F_week6) 
bic(fit_WMV_F_all2)

#Day 0 VPC
vpc_WMV_F_day0 = vpc(fit_WMV_F_day0; samples = 100)
vpc_plot(vpc_WMV_F_day0)
#week 6 VPC
vpc_WMV_F_week6 = vpc(fit_WMV_F_week6; samples = 100)
vpc_plot(vpc_WMV_F_week6)