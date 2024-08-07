#### create 2 comp model ####
etb_cov_forward = @model begin

    @metadata begin
        desc = "Two Compartment Model"
        timeu = u"hr"
    end

    @param begin
        """
        Clearance (L/hr)
        """
        tvcl ∈ RealDomain(; lower = 0)
        """
        Volume (L)
        """
        tvv ∈ RealDomain(; lower = 0) 
        """
        Absorption rate constant (h-1)
        """
        tvmtt ∈ RealDomain(; lower = 0) 
        """
        Typical volume (L) of peripheral 1 compartment
        """
        tvp1 ∈ RealDomain(; lower = 0) 
        """
        Transfer rate from Central to peripheral 1 (L/h)
        """
        tvq1 ∈ RealDomain(; lower = 0) 

        """
        HIV on mtt
        """
        hivMTT ∈ RealDomain(;)

        """
        Sex effects
        """
        sexV ∈ RealDomain(;)
        sexCL ∈ RealDomain(;) 
        sexF ∈ RealDomain(;) 
        sexMTT ∈ RealDomain(;)
        sexVP1 ∈ RealDomain(;)
        sexQ ∈ RealDomain(;)    

        """
        Age effects
        """
        ageV ∈ RealDomain(;)
        ageCL ∈ RealDomain(;) 
        ageF ∈ RealDomain(;) 
        ageMTT ∈ RealDomain(;)
        ageVP1 ∈ RealDomain(;)
        ageQ ∈ RealDomain(;)  

        """
        CREAT effects
        """
        creatCL ∈ RealDomain(;) 

        """
        ALT effects
        """
        altCL ∈ RealDomain(;) 

        """
        BILRB effects
        """
        bilrbCL ∈ RealDomain(;) 

        """
        TIME FLAG effect on clearence
        """
        tflagCL ∈ RealDomain(;) 
        """
        TIME FLAG effect on volume
        """
        tflagV ∈ RealDomain(;) 
        """
          - ΩCL
          - ΩVc
          - Ωmtt
          - ΩVp1
          - ΩQ1
        """
        Ω ∈ PDiagDomain(init = [0.04, 0.04, 0.04, 0.04, 0.04, 0.04])
        """
        Proportional RUV
        """
        σ_p ∈ RealDomain(; lower = 0.0001)
    end


    @random begin
        η ~ MvNormal(Ω) 
    end

    @covariates begin
        """
        Define all covariate names
        """
        WT 
        HIV 
        SEX
        AGE
        CREAT
        ALT
        BILRB
        TFLAG
    end

    @pre begin
        
        """
        Define all model parameters with every possible covariate reltionship 

        This allows the forward selection to freely experiment with different covariate permutations

        """
        
        CL = tvcl * (WT / median(etb_data.WT))^0.75 * (1 + sexCL*SEX) * (1+ageCL*(AGE - median(etb_data.AGE))) * (1+creatCL*(CREAT - median(etb_data.CREAT))) * (1+altCL*(ALT - median(etb_data.ALT))) * (1+bilrbCL*(BILRB - median(etb_data.BILRB))) * (1 + tflagCL*TFLAG)* exp(η[1])
        Vc = tvv * (WT / median(etb_data.WT)) * (1 + sexV*SEX) * (1+ageV*(AGE - median(etb_data.AGE))) * (1 + tflagV*TFLAG) * exp(η[2])

        MTT = tvmtt * (1 + (hivMTT*HIV)) * (1 + sexMTT*SEX) * (1+ageMTT*(AGE - median(etb_data.AGE))) * exp(η[3])
        ktr = 2/MTT 

        Vp1 = tvp1 * (WT / median(etb_data.WT)) * (1 + sexVP1*SEX) * (1+ageVP1*(AGE - median(etb_data.AGE))) * exp(η[4])
        Q1 = tvq1 * (WT / median(etb_data.WT))^0.75 * (1 + sexQ*SEX) * (1+ageQ*(AGE - median(etb_data.AGE))) * exp(η[5])
    
        """
        Predefine rates between compartments to help declutter model equations
        """
        kel = CL/Vc
        kc1 = Q1/Vc
        k1c = Q1/Vp1

         end

    @dosecontrol begin
        bioav = (;Depot1 = 1 * (1 + sexF*SEX) * (1 + ageF*(AGE - median(etb_data.AGE))) * exp(η[6]))
        """
        Relative bioavailability (F) %
       """
    end

    @dynamics begin
        Depot1' = -ktr*Depot1
        Transit1' = ktr*Depot1 - ktr*Transit1
        Central' = ktr*Transit1 - kel*Central - kc1*Central + k1c*Peripheral1
        Peripheral1' = kc1*Central - k1c*Peripheral1
    end


    @derived begin
        cp := @. log(0.0001+ (Central / Vc)/1000) 
        """ 
        ETB Concentration (ng/mL)
        """
        LNDV ~ @. Normal(cp, σ_p) 
    end

end

#set up initial parameter estimates using HIV MTT model [arameter estimates
iparms_cov_forward_explore = (;
tvcl = 53.394, tvv = 291.900, tvmtt = 1.643, tvp1 = 997.31, tvq1 = 43.368, hivMTT = 0.4191,
sexV = 0, sexCL = 0, sexF = 0, sexMTT = 0, sexVP1 = 0, sexQ= 0, ageV = 0, ageCL = 0, ageF = 0, ageMTT = 0, ageVP1 = 0, ageQ = 0, creatCL = 0, altCL = 0,  bilrbCL = 0, tflagCL =0, tflagV,
Ω = Diagonal([0.0259, 0.000317, 0.141, 0.431, 0.0578, 0.0501]), #add in eta for bioavailibility
 σ_p = 0.383)


#setting initial estimates of covariates to 0 allows for initial paramter check to be passed
fit_all_cov = fit(etb_cov_forward, etb_all_pumapop, iparms_cov_forward_explore, FOCE())
#call fit just to check that intiial neg log likelihood and its gradient are finite, allowing fitting to occur

#Call and run the forward selection
etb_forward_selection = covariate_select(
    etb_cov_forward,
    etb_all_pumapop,
    iparms_cov_forward_explore,
    FOCE();
    method = CovariateSelection.Forward, #run a forward selection
    criterion = bic, #use Schwarz criterion (BIC) to discriminate between models
    control_param = (:sexV, :sexCL, :sexF, :sexMTT, :sexVP1, :sexQ, :ageV, :ageCL, :ageF, :ageMTT, :ageVP1, :ageQ, :creatCL, :altCL, :bilrbCL, :tflagCL
    )
)
#note that control parameters do not include Weight and hivMTT, thus these will always be included in whatever model explored

forward_df = DataFrame(etb_forward_selection.fits)
#convert the fit data into a dataframe
forward_df.fitted_model
#call the fitted models

etb_forward_selection.best_model
#best model has no further additions
bic(etb_forward_selection.best_model)
#very similar BIC/parameters
infer(etb_forward_selection.best_model)



