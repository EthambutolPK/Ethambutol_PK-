#Final PK covairate model

############### NOTE: the terms IIV and eta are used interchangeably in this document ############


#rename to reflect which IIVs are being included in this model
etb_IIV_all = @model begin

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
        ktr = 2/MTT #define rate of transfer

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
        #/1000 to appropriately get ng/ml from mg/L
        """ 
        ETB Concentration (ng/mL)
        """
        LNDV ~ @. Normal(cp, σ_p) #add in error and estimate log conc (LNDV)
    end

end


iparms_IIV_all = (; #initial parameters as model parameter estimates
tvcl = 53.547, tvv = 286.32, tvmtt = 1.6521, tvp1 = 755.29, tvq1 = 43.645, hivMTT = 0.44871,
Ω = Diagonal([0.024833, 0.000001744, 0.13914, 0.52548, 0.06841, 0.048107]), #using zeroed etas 
σ_p = 0.38293)

fit_etb_IIV_all = fit(etb_IIV_all, etb_all_pumapop, iparms_IIV_all, FOCE()) #BIC 1864.940
bic(fit_etb_IIV_all)

############################################ NO IIV on VOLUME ################################################################################################################

etb_IIV_no_V = @model begin

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
        Vc = tvv * (WT / median(etb_data.WT))  #add in the effect of weight on volume

        #IIV on Central volume (Vc) is removed

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
        #/1000 to appropriately get ng/ml from mg/L
        """ 
        ETB Concentration (ng/mL)
        """
        LNDV ~ @. Normal(cp, σ_p) #add in error and estimate log conc (LNDV)
    end

end


iparms_IIV_no_V = (; #initial parameters as model parameter estimates
tvcl = 53.547, tvv = 286.32, tvmtt = 1.6521, tvp1 = 755.29, tvq1 = 43.645, hivMTT = 0.44871,
Ω = Diagonal([0.024833, 0.13914, 0.52548, 0.06841, 0.048107]), #using zeroed etas 
σ_p = 0.38293)

fit_etb_IVV_no_V = fit(etb_IIV_no_V, etb_all_pumapop, iparms_IIV_no_V, FOCE())
bic(fit_etb_IVV_no_V) #1857.643
############################################ NO IIV on VOLUME or Clearance ################################################################################################################

etb_IIV_no_V_CL = @model begin

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
        Ω ∈ PDiagDomain(init = [0.04, 0.04, 0.04, 0.04]) #initial matrix estimates
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
        CL = tvcl * (WT / median(etb_data.WT))^0.75 #add in the effect of WT on CL/Q
        Vc = tvv * (WT / median(etb_data.WT))  #add in the effect of weight on volume

        #IIV removed on both central volume and clearance

        MTT = tvmtt * exp(η[1]) * (1 + (hivMTT*HIV))#define mean transit time
        ktr = 2/MTT #define rate of Transfer

        Vp1 = tvp1 * exp(η[2]) * (WT / median(etb_data.WT))  #add in the effect of weight on volume
        Q1 = tvq1 * exp(η[3]) * (WT / median(etb_data.WT))^0.75 #add in the effect of WT on CL/Q
    end

    @dosecontrol begin
        bioav = (; Depot1 = 1 * exp(η[4])) #add an eta into bioavailbility

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
        #/1000 to appropriately get ng/ml from mg/L
        """ 
        ETB Concentration (ng/mL)
        """
        LNDV ~ @. Normal(cp, σ_p) #add in error and estimate log conc (LNDV)
    end

end


iparms_IIV_no_V_CL = (; #initial parameters as model parameter estimates
tvcl = 53.547, tvv = 286.32, tvmtt = 1.6521, tvp1 = 755.29, tvq1 = 43.645, hivMTT = 0.44871,
Ω = Diagonal([0.13914, 0.52548, 0.06841, 0.048107]), #using zeroed etas 
σ_p = 0.38293)

fit_etb_IVV_no_V_CL = fit(etb_IIV_no_V_CL, etb_all_pumapop, iparms_IIV_no_V_CL, FOCE())
bic(fit_etb_IVV_no_V_CL) #1856.562

############################################ NO IIV on VOLUME or Clearance or bioavailability ################################################################################################################

etb_IIV_no_V_CL_F = @model begin

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
        Ω ∈ PDiagDomain(init = [0.04, 0.04, 0.04]) #initial matrix estimates
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
        CL = tvcl * (WT / median(etb_data.WT))^0.75 #add in the effect of WT on CL/Q
        Vc = tvv * (WT / median(etb_data.WT))  #add in the effect of weight on volume

        MTT = tvmtt * exp(η[1]) * (1 + (hivMTT*HIV))#define mean transit time
        ktr = 2/MTT #define rate of Transfer

        Vp1 = tvp1 * exp(η[2]) * (WT / median(etb_data.WT))  #add in the effect of weight on volume
        Q1 = tvq1 * exp(η[3]) * (WT / median(etb_data.WT))^0.75 #add in the effect of WT on CL/Q
    end

    @dosecontrol begin
        bioav = (; Depot1 = 1) #remove IIV on bioavailbility

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
        #/1000 to appropriately get ng/ml from mg/L
        """ 
        ETB Concentration (ng/mL)
        """
        LNDV ~ @. Normal(cp, σ_p) #add in error and estimate log conc (LNDV)
    end

end


iparms_IIV_no_V_CL_F = (; #initial parameters as model parameter estimates
tvcl = 53.547, tvv = 286.32, tvmtt = 1.6521, tvp1 = 755.29, tvq1 = 43.645, hivMTT = 0.44871,
Ω = Diagonal([0.13914, 0.52548, 0.06841]), #using zeroed etas 
σ_p = 0.38293)

fit_etb_IVV_no_V_CL_F = fit(etb_IIV_no_V_CL_F, etb_all_pumapop, iparms_IIV_no_V_CL_F, FOCE())
bic(fit_etb_IVV_no_V_CL_F) #1857.643
