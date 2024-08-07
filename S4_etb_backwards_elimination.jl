
############################################## Backwards select #######################################################

#Full PK model: WT + HIV MTT

#Remove MTT

### Just WT covariate included #####
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
          - ΩCL
          - ΩVc
          - Ωmtt
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



#initial parameters
iparms_WT_only = (;
tvcl = 65, tvv = 285, tvmtt = 1.79, tvp1 = 345, tvq1 = 40,
Ω = Diagonal([0.05, 0.05, 0.05, 0.05, 0.05, 0.05]),
σ_p = 0.3)


#run fit
fit_WT_only_all = fit(etb_covariate_WT, etb_all_pumapop, iparms_WT_only, FOCE())
#log likeliood = -890.523

#bayesian (Schwarz) information criterion
bic(fit_WT_only_all) #1868.603

#TAD vpc
vpc_WT_HIV_MTT_VP1_all = vpc(fit_WT_HIV_MTT_VP1_all; samples = 100, covariates = [:tad])
vpc_plot(vpc_WT_HIV_MTT_VP1_all)



###################### HIV mtt only ##############################################


### Just HIV covariate on MTT included #####
etb_HIV_MTT = @model begin

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
        Ω ∈ PDiagDomain(init = [0.04, 0.04, 0.04, 0.04, 0.04, 0.04]) #initial matrix estimates
        """
        Proportional RUV
        """
        σ_p ∈ RealDomain(; lower = 0.0001) #proportional error, initially estimated as 20%
    end


    @random begin #define the nature of parameter variability 
        η ~ MvNormal(Ω) #define all parameters as normally distributed
    end

    @covariates begin
        HIV #add in HIV status as a covariate
    end

    #no covariates block included just yet, this is where covairate names will be loaded 

    @pre begin  #use the pre block to define parameters for the model
        CL = tvcl * exp(η[1])
        Vc = tvv * exp(η[2])

        MTT = tvmtt * exp(η[3]) * (1 + (hivMTT*HIV))#define mean transit time
        ktr = 2/MTT #define rate of Transfer

        Vp1 = tvp1 * exp(η[4])
        Q1 = tvq1 * exp(η[5])
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
        #/1000 to appropriately get ng/ml from ng/L
        """ 
        ETB Concentration (ng/mL)
        """
        LNDV ~ @. Normal(cp, σ_p) #add in error and estimate log conc (LNDV)
    end

end



#initial parameters
iparms_MTT_only = (;
tvcl = 65, tvv = 285, tvmtt = 1.79, tvp1 = 345, tvq1 = 40, hivMTT = 0.2,
Ω = Diagonal([0.05, 0.05, 0.05, 0.05, 0.05, 0.05]),
σ_p = 0.3)

#run fit
fit_MTT_only_all = fit(etb_HIV_MTT, etb_all_pumapop, iparms_MTT_only, FOCE())
#Log likeliood = -891.011

#bayesian (Schwarz) information criterion
bic(fit_MTT_only_all) #1876

#TAD vpc
vpc_MTT_only_all = vpc(fit_MTT_only_all; samples = 100, covariates = [:tad])
vpc_plot(vpc_MTT_only_all)