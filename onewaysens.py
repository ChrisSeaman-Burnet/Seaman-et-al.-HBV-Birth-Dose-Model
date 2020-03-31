#One way sensitivity analyses
from numpy import array

if region==1:
    if runs==1:
        
        #To reset paramaters after each sensitivity analysis
        globz_copy=array(globz.iloc[:,1])
        hbv_pars_copy=array(hbv_pars.iloc[:,:])
        
#%% HBsAg Prevalence
        hbsag_sens=zeros((3,9)) 
        
        hbv_pars.iloc[6,8]=0.02     #2% HBsAg prevalence
        
        fac_cov_weighted, com_cov_weighted = facility_weight(runs, weight, hbv_pars, regions_perturbed)
        saturation, sat_costing=coverage_saturation(fac_cov_weighted, com_cov_weighted, hbv_pars, regions_perturbed,runs)
        for scen in range(0,5,4): 
            fac_cc, fac_ctc, com_cc, com_ctc, com_cpad= model_scenarios(runs, scen, hbv_pars, fac_cov_weighted, com_cov_weighted, saturation)
            
            if scen==0:
                model_init=model_initializer(runs, globz, hbv_pars, globz_perturbed, regions_perturbed, fac_cc, fac_ctc, com_cc, com_ctc, com_cpad, t_steps,ve_warm)
                sensitivity_vcost, sensitivity_coverage=vaccine_costs(runs, hbv_pars, fac_cc, fac_ctc, com_cc, com_ctc, com_cpad,sat_costing, saturation, regions_perturbed)
                for reg in range(len(hbv_pars)):
                    model_init[reg,0,-2]=sensitivity_vcost[reg,5]
                sensivity_model=hbv_model(runs,t_steps,dt,model_init, hbv_pars, regions_perturbed, globz, globz_perturbed,discount)
                
                hbsag_sens[0,0]=sensivity_model[6,len(t_steps)-1,13]      #DALYs
                hbsag_sens[0,1]=sensivity_model[6,len(t_steps)-1,11]      #Disease Cost
                hbsag_sens[0,2]=sensivity_model[6,0,12]                   #Vaccine Cost
                hbsag_sens[0,3]=sensivity_model[6,len(t_steps)-1,12]      #Total Cost (for ICER)
            
            if scen==4:
                model_init=model_initializer(runs, globz, hbv_pars, globz_perturbed, regions_perturbed, fac_cc, fac_ctc, com_cc, com_ctc, com_cpad, t_steps,ve_warm)
                sensitivity_vcost, sensitivity_coverage=vaccine_costs(runs, hbv_pars, fac_cc, fac_ctc, com_cc, com_ctc, com_cpad,sat_costing, saturation, regions_perturbed)
                for reg in range(len(hbv_pars)):
                    model_init[reg,0,-2]=sensitivity_vcost[reg,5]
                sensitivity_model=hbv_model(runs,t_steps,dt,model_init, hbv_pars, regions_perturbed, globz, globz_perturbed,discount)
                
                hbsag_sens[0,4]=sensitivity_model[6,len(t_steps)-1,13]      #DALYs
                hbsag_sens[0,5]=sensitivity_model[6,len(t_steps)-1,11]      #Disease Cost
                hbsag_sens[0,6]=sensitivity_model[6,0,12]                   #Vaccine Cost
                hbsag_sens[0,7]=sensitivity_model[6,len(t_steps)-1,12]      #Total Cost (for ICER)
            
        
        hbv_pars.iloc[6,8]=0.05
        fac_cov_weighted, com_cov_weighted = facility_weight(runs, weight, hbv_pars, regions_perturbed)
        saturation, sat_costing=coverage_saturation(fac_cov_weighted, com_cov_weighted, hbv_pars, regions_perturbed,runs)
        for scen in range(0,5,4):       #baseline and scenario 4
            fac_cc, fac_ctc, com_cc, com_ctc, com_cpad= model_scenarios(runs, scen, hbv_pars, fac_cov_weighted, com_cov_weighted, saturation)
            
            if scen==0:
                model_init=model_initializer(runs, globz, hbv_pars, globz_perturbed, regions_perturbed, fac_cc, fac_ctc, com_cc, com_ctc, com_cpad, t_steps,ve_warm)
                sensitivity_vcost, sensitivity_coverage=vaccine_costs(runs, hbv_pars, fac_cc, fac_ctc, com_cc, com_ctc, com_cpad,sat_costing, saturation, regions_perturbed)
                for reg in range(len(hbv_pars)):
                    model_init[reg,0,-2]=sensitivity_vcost[reg,5]
                sensivity_model=hbv_model(runs,t_steps,dt,model_init, hbv_pars, regions_perturbed, globz, globz_perturbed,discount)
                
                hbsag_sens[1,0]=sensivity_model[6,len(t_steps)-1,13]      #DALYs
                hbsag_sens[1,1]=sensivity_model[6,len(t_steps)-1,11]      #Disease Cost
                hbsag_sens[1,2]=sensivity_model[6,0,12]                   #Vaccine Cost
                hbsag_sens[1,3]=sensivity_model[6,len(t_steps)-1,12]      #Total Cost (for ICER)
            
            if scen==4:
                model_init=model_initializer(runs, globz, hbv_pars, globz_perturbed, regions_perturbed, fac_cc, fac_ctc, com_cc, com_ctc, com_cpad, t_steps,ve_warm)
                sensitivity_vcost, sensitivity_coverage=vaccine_costs(runs, hbv_pars, fac_cc, fac_ctc, com_cc, com_ctc, com_cpad,sat_costing, saturation, regions_perturbed)
                for reg in range(len(hbv_pars)):
                    model_init[reg,0,-2]=sensitivity_vcost[reg,5]
                sensitivity_model=hbv_model(runs,t_steps,dt,model_init, hbv_pars, regions_perturbed, globz, globz_perturbed,discount)
                
                hbsag_sens[1,4]=sensitivity_model[6,len(t_steps)-1,13]      #DALYs
                hbsag_sens[1,5]=sensitivity_model[6,len(t_steps)-1,11]      #Disease Cost
                hbsag_sens[1,6]=sensitivity_model[6,0,12]                   #Vaccine Cost
                hbsag_sens[1,7]=sensitivity_model[6,len(t_steps)-1,12]      #Total Cost (for ICER)
            
        
        hbv_pars.iloc[6,8]=0.1
        fac_cov_weighted, com_cov_weighted = facility_weight(runs, weight, hbv_pars, regions_perturbed)
        saturation, sat_costing=coverage_saturation(fac_cov_weighted, com_cov_weighted, hbv_pars, regions_perturbed,runs)
        for scen in range(0,5,4):       #baseline and scenario 4
            fac_cc, fac_ctc, com_cc, com_ctc, com_cpad= model_scenarios(runs, scen, hbv_pars, fac_cov_weighted, com_cov_weighted, saturation)
            
            if scen==0:
                model_init=model_initializer(runs, globz, hbv_pars, globz_perturbed, regions_perturbed, fac_cc, fac_ctc, com_cc, com_ctc, com_cpad, t_steps,ve_warm)
                sensitivity_vcost, sensitivity_coverage=vaccine_costs(runs, hbv_pars, fac_cc, fac_ctc, com_cc, com_ctc, com_cpad,sat_costing, saturation, regions_perturbed)
                for reg in range(len(hbv_pars)):
                    model_init[reg,0,-2]=sensitivity_vcost[reg,5]
                sensivity_model=hbv_model(runs,t_steps,dt,model_init, hbv_pars, regions_perturbed, globz, globz_perturbed,discount)
                
                hbsag_sens[2,0]=sensivity_model[6,len(t_steps)-1,13]      #DALYs
                hbsag_sens[2,1]=sensivity_model[6,len(t_steps)-1,11]      #Disease Cost
                hbsag_sens[2,2]=sensivity_model[6,0,12]                   #Vaccine Cost
                hbsag_sens[2,3]=sensivity_model[6,len(t_steps)-1,12]      #Total Cost (for ICER)
            
            if scen==4:
                model_init=model_initializer(runs, globz, hbv_pars, globz_perturbed, regions_perturbed, fac_cc, fac_ctc, com_cc, com_ctc, com_cpad, t_steps,ve_warm)
                sensitivity_vcost, sensitivity_coverage=vaccine_costs(runs, hbv_pars, fac_cc, fac_ctc, com_cc, com_ctc, com_cpad,sat_costing, saturation, regions_perturbed)
                for reg in range(len(hbv_pars)):
                    model_init[reg,0,-2]=sensitivity_vcost[reg,5]
                sensitivity_model=hbv_model(runs,t_steps,dt,model_init, hbv_pars, regions_perturbed, globz, globz_perturbed,discount)
                
                hbsag_sens[2,4]=sensitivity_model[6,len(t_steps)-1,13]      #DALYs
                hbsag_sens[2,5]=sensitivity_model[6,len(t_steps)-1,11]      #Disease Cost
                hbsag_sens[2,6]=sensitivity_model[6,0,12]                   #Vaccine Cost
                hbsag_sens[2,7]=sensitivity_model[6,len(t_steps)-1,12]      #Total Cost (for ICER)
        
        hbsag_sens[:,8]=-((hbsag_sens[:,7]-hbsag_sens[:,3])/(hbsag_sens[:,4]-hbsag_sens[:,0]))
        
        
#%% HBeAg Prevalence
        hbe_sens=zeros((2,9))
        hbv_pars.iloc[6,8]=hbv_pars_copy[6,8]
        #5% HBeAg prevalence
        hbv_pars.iloc[6,48]=0.05
        
        fac_cov_weighted, com_cov_weighted = facility_weight(runs, weight, hbv_pars, regions_perturbed)
        saturation, sat_costing=coverage_saturation(fac_cov_weighted, com_cov_weighted, hbv_pars, regions_perturbed,runs)
        for scen in range(0,5,4):       #baseline and scenario 4
            fac_cc, fac_ctc, com_cc, com_ctc, com_cpad= model_scenarios(runs, scen, hbv_pars, fac_cov_weighted, com_cov_weighted, saturation)
            
            if scen==0:
                model_init=model_initializer(runs, globz, hbv_pars, globz_perturbed, regions_perturbed, fac_cc, fac_ctc, com_cc, com_ctc, com_cpad, t_steps,ve_warm)
                sensitivity_vcost, sensitivity_coverage=vaccine_costs(runs, hbv_pars, fac_cc, fac_ctc, com_cc, com_ctc, com_cpad,sat_costing, saturation, regions_perturbed)
                for reg in range(len(hbv_pars)):
                    model_init[reg,0,-2]=sensitivity_vcost[reg,5]
                sensivity_model=hbv_model(runs,t_steps,dt,model_init, hbv_pars, regions_perturbed, globz, globz_perturbed,discount)
                
                hbe_sens[0,0]=sensivity_model[6,len(t_steps)-1,13]      #DALYs
                hbe_sens[0,1]=sensivity_model[6,len(t_steps)-1,11]      #Disease Cost
                hbe_sens[0,2]=sensivity_model[6,0,12]                   #Vaccine Cost
                hbe_sens[0,3]=sensivity_model[6,len(t_steps)-1,12]      #Total Cost (for ICER)
            
            if scen==4:
                model_init=model_initializer(runs, globz, hbv_pars, globz_perturbed, regions_perturbed, fac_cc, fac_ctc, com_cc, com_ctc, com_cpad, t_steps,ve_warm)
                sensitivity_vcost, sensitivity_coverage=vaccine_costs(runs, hbv_pars, fac_cc, fac_ctc, com_cc, com_ctc, com_cpad,sat_costing, saturation, regions_perturbed)
                for reg in range(len(hbv_pars)):
                    model_init[reg,0,-2]=sensitivity_vcost[reg,5]
                sensitivity_model=hbv_model(runs,t_steps,dt,model_init, hbv_pars, regions_perturbed, globz, globz_perturbed,discount)
                
                hbe_sens[0,4]=sensitivity_model[6,len(t_steps)-1,13]      #DALYs
                hbe_sens[0,5]=sensitivity_model[6,len(t_steps)-1,11]      #Disease Cost
                hbe_sens[0,6]=sensitivity_model[6,0,12]                   #Vaccine Cost
                hbe_sens[0,7]=sensitivity_model[6,len(t_steps)-1,12]      #Total Cost (for ICER)
        
        #50% HBeAg prevalence
        hbv_pars.iloc[6,48]=0.5
        
        fac_cov_weighted, com_cov_weighted = facility_weight(runs, weight, hbv_pars, regions_perturbed)
        saturation, sat_costing=coverage_saturation(fac_cov_weighted, com_cov_weighted, hbv_pars, regions_perturbed,runs)
        for scen in range(0,5,4):       #baseline and scenario 4
            fac_cc, fac_ctc, com_cc, com_ctc, com_cpad= model_scenarios(runs, scen, hbv_pars, fac_cov_weighted, com_cov_weighted, saturation)
            
            if scen==0:
                model_init=model_initializer(runs, globz, hbv_pars, globz_perturbed, regions_perturbed, fac_cc, fac_ctc, com_cc, com_ctc, com_cpad, t_steps,ve_warm)
                sensitivity_vcost, sensitivity_coverage=vaccine_costs(runs, hbv_pars, fac_cc, fac_ctc, com_cc, com_ctc, com_cpad,sat_costing, saturation, regions_perturbed)
                for reg in range(len(hbv_pars)):
                    model_init[reg,0,-2]=sensitivity_vcost[reg,5]
                sensivity_model=hbv_model(runs,t_steps,dt,model_init, hbv_pars, regions_perturbed, globz, globz_perturbed,discount)
                
                hbe_sens[1,0]=sensivity_model[6,len(t_steps)-1,13]      #DALYs
                hbe_sens[1,1]=sensivity_model[6,len(t_steps)-1,11]      #Disease Cost
                hbe_sens[1,2]=sensivity_model[6,0,12]                   #Vaccine Cost
                hbe_sens[1,3]=sensivity_model[6,len(t_steps)-1,12]      #Total Cost (for ICER)
            
            if scen==4:
                model_init=model_initializer(runs, globz, hbv_pars, globz_perturbed, regions_perturbed, fac_cc, fac_ctc, com_cc, com_ctc, com_cpad, t_steps,ve_warm)
                sensitivity_vcost, sensitivity_coverage=vaccine_costs(runs, hbv_pars, fac_cc, fac_ctc, com_cc, com_ctc, com_cpad,sat_costing, saturation, regions_perturbed)
                for reg in range(len(hbv_pars)):
                    model_init[reg,0,-2]=sensitivity_vcost[reg,5]
                sensitivity_model=hbv_model(runs,t_steps,dt,model_init, hbv_pars, regions_perturbed, globz, globz_perturbed,discount)
                
                hbe_sens[1,4]=sensitivity_model[6,len(t_steps)-1,13]      #DALYs
                hbe_sens[1,5]=sensitivity_model[6,len(t_steps)-1,11]      #Disease Cost
                hbe_sens[1,6]=sensitivity_model[6,0,12]                   #Vaccine Cost
                hbe_sens[1,7]=sensitivity_model[6,len(t_steps)-1,12]      #Total Cost (for ICER)
            
            hbe_sens[:,8]=-((hbe_sens[:,7]-hbe_sens[:,3])/(hbe_sens[:,4]-hbe_sens[:,0]))

#%% Facility Births Proportion
        facb_sens=zeros((2,9))
        hbv_pars.iloc[6,48]=hbv_pars_copy[6,48]
        
        #10% Facility Births
        hbv_pars.iloc[6,14]=0.10
        
        fac_cov_weighted, com_cov_weighted = facility_weight(runs, weight, hbv_pars, regions_perturbed)
        saturation, sat_costing=coverage_saturation(fac_cov_weighted, com_cov_weighted, hbv_pars, regions_perturbed,runs)
        for scen in range(0,5,4):       #baseline and scenario 4
            fac_cc, fac_ctc, com_cc, com_ctc, com_cpad= model_scenarios(runs, scen, hbv_pars, fac_cov_weighted, com_cov_weighted, saturation)
            
            if scen==0:
                model_init=model_initializer(runs, globz, hbv_pars, globz_perturbed, regions_perturbed, fac_cc, fac_ctc, com_cc, com_ctc, com_cpad, t_steps,ve_warm)
                sensitivity_vcost, sensitivity_coverage=vaccine_costs(runs, hbv_pars, fac_cc, fac_ctc, com_cc, com_ctc, com_cpad,sat_costing, saturation, regions_perturbed)
                for reg in range(len(hbv_pars)):
                    model_init[reg,0,-2]=sensitivity_vcost[reg,5]
                sensivity_model=hbv_model(runs,t_steps,dt,model_init, hbv_pars, regions_perturbed, globz, globz_perturbed,discount)
                
                facb_sens[0,0]=sensivity_model[6,len(t_steps)-1,13]      #DALYs
                facb_sens[0,1]=sensivity_model[6,len(t_steps)-1,11]      #Disease Cost
                facb_sens[0,2]=sensivity_model[6,0,12]                   #Vaccine Cost
                facb_sens[0,3]=sensivity_model[6,len(t_steps)-1,12]      #Total Cost (for ICER)
            
            if scen==4:
                model_init=model_initializer(runs, globz, hbv_pars, globz_perturbed, regions_perturbed, fac_cc, fac_ctc, com_cc, com_ctc, com_cpad, t_steps,ve_warm)
                sensitivity_vcost, sensitivity_coverage=vaccine_costs(runs, hbv_pars, fac_cc, fac_ctc, com_cc, com_ctc, com_cpad,sat_costing, saturation, regions_perturbed)
                for reg in range(len(hbv_pars)):
                    model_init[reg,0,-2]=sensitivity_vcost[reg,5]
                sensitivity_model=hbv_model(runs,t_steps,dt,model_init, hbv_pars, regions_perturbed, globz, globz_perturbed,discount)
                
                facb_sens[0,4]=sensitivity_model[6,len(t_steps)-1,13]      #DALYs
                facb_sens[0,5]=sensitivity_model[6,len(t_steps)-1,11]      #Disease Cost
                facb_sens[0,6]=sensitivity_model[6,0,12]                   #Vaccine Cost
                facb_sens[0,7]=sensitivity_model[6,len(t_steps)-1,12]      #Total Cost (for ICER)
        
        #99% Facility Births
        hbv_pars.iloc[6,14]=0.99
        fac_cov_weighted, com_cov_weighted = facility_weight(runs, weight, hbv_pars, regions_perturbed)
        saturation, sat_costing=coverage_saturation(fac_cov_weighted, com_cov_weighted, hbv_pars, regions_perturbed,runs)
        for scen in range(0,5,4):       #baseline and scenario 4
            fac_cc, fac_ctc, com_cc, com_ctc, com_cpad= model_scenarios(runs, scen, hbv_pars, fac_cov_weighted, com_cov_weighted, saturation)
            
            if scen==0:
                model_init=model_initializer(runs, globz, hbv_pars, globz_perturbed, regions_perturbed, fac_cc, fac_ctc, com_cc, com_ctc, com_cpad, t_steps,ve_warm)
                sensitivity_vcost, sensitivity_coverage=vaccine_costs(runs, hbv_pars, fac_cc, fac_ctc, com_cc, com_ctc, com_cpad,sat_costing, saturation, regions_perturbed)
                for reg in range(len(hbv_pars)):
                    model_init[reg,0,-2]=sensitivity_vcost[reg,5]
                sensivity_model=hbv_model(runs,t_steps,dt,model_init, hbv_pars, regions_perturbed, globz, globz_perturbed,discount)
                
                facb_sens[1,0]=sensivity_model[6,len(t_steps)-1,13]      #DALYs
                facb_sens[1,1]=sensivity_model[6,len(t_steps)-1,11]      #Disease Cost
                facb_sens[1,2]=sensivity_model[6,0,12]                   #Vaccine Cost
                facb_sens[1,3]=sensivity_model[6,len(t_steps)-1,12]      #Total Cost (for ICER)
            
            if scen==4:
                model_init=model_initializer(runs, globz, hbv_pars, globz_perturbed, regions_perturbed, fac_cc, fac_ctc, com_cc, com_ctc, com_cpad, t_steps,ve_warm)
                sensitivity_vcost, sensitivity_coverage=vaccine_costs(runs, hbv_pars, fac_cc, fac_ctc, com_cc, com_ctc, com_cpad,sat_costing, saturation, regions_perturbed)
                for reg in range(len(hbv_pars)):
                    model_init[reg,0,-2]=sensitivity_vcost[reg,5]
                sensitivity_model=hbv_model(runs,t_steps,dt,model_init, hbv_pars, regions_perturbed, globz, globz_perturbed,discount)
                
                facb_sens[1,4]=sensitivity_model[6,len(t_steps)-1,13]      #DALYs
                facb_sens[1,5]=sensitivity_model[6,len(t_steps)-1,11]      #Disease Cost
                facb_sens[1,6]=sensitivity_model[6,0,12]                   #Vaccine Cost
                facb_sens[1,7]=sensitivity_model[6,len(t_steps)-1,12]      #Total Cost (for ICER)
        
            facb_sens[:,8]=-((facb_sens[:,7]-facb_sens[:,3])/(facb_sens[:,4]-facb_sens[:,0]))
#%% Attendance at Community Births
        comsba_sens=zeros((2,9))
        hbv_pars.iloc[6,14]=hbv_pars_copy[6,14]
        
        #1% Attendance at Community Births
        hbv_pars.iloc[6,47]=0.01
        fac_cov_weighted, com_cov_weighted = facility_weight(runs, weight, hbv_pars, regions_perturbed)
        saturation, sat_costing=coverage_saturation(fac_cov_weighted, com_cov_weighted, hbv_pars, regions_perturbed,runs)
        for scen in range(0,5,4):       #baseline and scenario 4
            fac_cc, fac_ctc, com_cc, com_ctc, com_cpad= model_scenarios(runs, scen, hbv_pars, fac_cov_weighted, com_cov_weighted, saturation)
            
            if scen==0:
                model_init=model_initializer(runs, globz, hbv_pars, globz_perturbed, regions_perturbed, fac_cc, fac_ctc, com_cc, com_ctc, com_cpad, t_steps,ve_warm)
                sensitivity_vcost, sensitivity_coverage=vaccine_costs(runs, hbv_pars, fac_cc, fac_ctc, com_cc, com_ctc, com_cpad,sat_costing, saturation, regions_perturbed)
                for reg in range(len(hbv_pars)):
                    model_init[reg,0,-2]=sensitivity_vcost[reg,5]
                sensivity_model=hbv_model(runs,t_steps,dt,model_init, hbv_pars, regions_perturbed, globz, globz_perturbed,discount)
                
                comsba_sens[0,0]=sensivity_model[6,len(t_steps)-1,13]      #DALYs
                comsba_sens[0,1]=sensivity_model[6,len(t_steps)-1,11]      #Disease Cost
                comsba_sens[0,2]=sensivity_model[6,0,12]                   #Vaccine Cost
                comsba_sens[0,3]=sensivity_model[6,len(t_steps)-1,12]      #Total Cost (for ICER)
            
            if scen==4:
                model_init=model_initializer(runs, globz, hbv_pars, globz_perturbed, regions_perturbed, fac_cc, fac_ctc, com_cc, com_ctc, com_cpad, t_steps,ve_warm)
                sensitivity_vcost, sensitivity_coverage=vaccine_costs(runs, hbv_pars, fac_cc, fac_ctc, com_cc, com_ctc, com_cpad,sat_costing, saturation, regions_perturbed)
                for reg in range(len(hbv_pars)):
                    model_init[reg,0,-2]=sensitivity_vcost[reg,5]
                sensitivity_model=hbv_model(runs,t_steps,dt,model_init, hbv_pars, regions_perturbed, globz, globz_perturbed,discount)
                
                comsba_sens[0,4]=sensitivity_model[6,len(t_steps)-1,13]      #DALYs
                comsba_sens[0,5]=sensitivity_model[6,len(t_steps)-1,11]      #Disease Cost
                comsba_sens[0,6]=sensitivity_model[6,0,12]                   #Vaccine Cost
                comsba_sens[0,7]=sensitivity_model[6,len(t_steps)-1,12]      #Total Cost (for ICER)
        
        #90% Attendance at Community Births
        hbv_pars.iloc[6,47]=0.90
        fac_cov_weighted, com_cov_weighted = facility_weight(runs, weight, hbv_pars, regions_perturbed)
        saturation, sat_costing=coverage_saturation(fac_cov_weighted, com_cov_weighted, hbv_pars, regions_perturbed,runs)
        for scen in range(0,5,4):       #baseline and scenario 4
            fac_cc, fac_ctc, com_cc, com_ctc, com_cpad= model_scenarios(runs, scen, hbv_pars, fac_cov_weighted, com_cov_weighted, saturation)
            
            if scen==0:
                model_init=model_initializer(runs, globz, hbv_pars, globz_perturbed, regions_perturbed, fac_cc, fac_ctc, com_cc, com_ctc, com_cpad, t_steps,ve_warm)
                sensitivity_vcost, sensitivity_coverage=vaccine_costs(runs, hbv_pars, fac_cc, fac_ctc, com_cc, com_ctc, com_cpad,sat_costing, saturation, regions_perturbed)
                for reg in range(len(hbv_pars)):
                    model_init[reg,0,-2]=sensitivity_vcost[reg,5]
                sensivity_model=hbv_model(runs,t_steps,dt,model_init, hbv_pars, regions_perturbed, globz, globz_perturbed,discount)
                
                comsba_sens[1,0]=sensivity_model[6,len(t_steps)-1,13]      #DALYs
                comsba_sens[1,1]=sensivity_model[6,len(t_steps)-1,11]      #Disease Cost
                comsba_sens[1,2]=sensivity_model[6,0,12]                   #Vaccine Cost
                comsba_sens[1,3]=sensivity_model[6,len(t_steps)-1,12]      #Total Cost (for ICER)
            
            if scen==4:
                model_init=model_initializer(runs, globz, hbv_pars, globz_perturbed, regions_perturbed, fac_cc, fac_ctc, com_cc, com_ctc, com_cpad, t_steps,ve_warm)
                sensitivity_vcost, sensitivity_coverage=vaccine_costs(runs, hbv_pars, fac_cc, fac_ctc, com_cc, com_ctc, com_cpad,sat_costing, saturation, regions_perturbed)
                for reg in range(len(hbv_pars)):
                    model_init[reg,0,-2]=sensitivity_vcost[reg,5]
                sensitivity_model=hbv_model(runs,t_steps,dt,model_init, hbv_pars, regions_perturbed, globz, globz_perturbed,discount)
                
                comsba_sens[1,4]=sensitivity_model[6,len(t_steps)-1,13]      #DALYs
                comsba_sens[1,5]=sensitivity_model[6,len(t_steps)-1,11]      #Disease Cost
                comsba_sens[1,6]=sensitivity_model[6,0,12]                   #Vaccine Cost
                comsba_sens[1,7]=sensitivity_model[6,len(t_steps)-1,12]      #Total Cost (for ICER)
                
            comsba_sens[:,8]=-((comsba_sens[:,7]-comsba_sens[:,3])/(comsba_sens[:,4]-comsba_sens[:,0]))
#%%  Baseline HBV Birth Dose Coverage 
        hbv_pars.iloc[6,47]=hbv_pars_copy[6,47]
        blcov_sens=zeros((2,9))
        
        #10% Baseline Coverage
        hbv_pars.iloc[6,6]=0.10
        fac_cov_weighted, com_cov_weighted = facility_weight(runs, weight, hbv_pars, regions_perturbed)
        saturation, sat_costing=coverage_saturation(fac_cov_weighted, com_cov_weighted, hbv_pars, regions_perturbed,runs)
        for scen in range(0,5,4):       #baseline and scenario 4
            fac_cc, fac_ctc, com_cc, com_ctc, com_cpad= model_scenarios(runs, scen, hbv_pars, fac_cov_weighted, com_cov_weighted, saturation)
            
            if scen==0:
                model_init=model_initializer(runs, globz, hbv_pars, globz_perturbed, regions_perturbed, fac_cc, fac_ctc, com_cc, com_ctc, com_cpad, t_steps,ve_warm)
                sensitivity_vcost, sensitivity_coverage=vaccine_costs(runs, hbv_pars, fac_cc, fac_ctc, com_cc, com_ctc, com_cpad,sat_costing, saturation, regions_perturbed)
                for reg in range(len(hbv_pars)):
                    model_init[reg,0,-2]=sensitivity_vcost[reg,5]
                sensivity_model=hbv_model(runs,t_steps,dt,model_init, hbv_pars, regions_perturbed, globz, globz_perturbed,discount)
                
                blcov_sens[0,0]=sensivity_model[6,len(t_steps)-1,13]      #DALYs
                blcov_sens[0,1]=sensivity_model[6,len(t_steps)-1,11]      #Disease Cost
                blcov_sens[0,2]=sensivity_model[6,0,12]                   #Vaccine Cost
                blcov_sens[0,3]=sensivity_model[6,len(t_steps)-1,12]      #Total Cost (for ICER)
            
            if scen==4:
                model_init=model_initializer(runs, globz, hbv_pars, globz_perturbed, regions_perturbed, fac_cc, fac_ctc, com_cc, com_ctc, com_cpad, t_steps,ve_warm)
                sensitivity_vcost, sensitivity_coverage=vaccine_costs(runs, hbv_pars, fac_cc, fac_ctc, com_cc, com_ctc, com_cpad,sat_costing, saturation, regions_perturbed)
                for reg in range(len(hbv_pars)):
                    model_init[reg,0,-2]=sensitivity_vcost[reg,5]
                sensitivity_model=hbv_model(runs,t_steps,dt,model_init, hbv_pars, regions_perturbed, globz, globz_perturbed,discount)
                
                blcov_sens[0,4]=sensitivity_model[6,len(t_steps)-1,13]      #DALYs
                blcov_sens[0,5]=sensitivity_model[6,len(t_steps)-1,11]      #Disease Cost
                blcov_sens[0,6]=sensitivity_model[6,0,12]                   #Vaccine Cost
                blcov_sens[0,7]=sensitivity_model[6,len(t_steps)-1,12]      #Total Cost (for ICER)
        
        #85% Baseline Coverage
        hbv_pars.iloc[6,6]=0.85
        fac_cov_weighted, com_cov_weighted = facility_weight(runs, weight, hbv_pars, regions_perturbed)
        saturation, sat_costing=coverage_saturation(fac_cov_weighted, com_cov_weighted, hbv_pars, regions_perturbed,runs)
        for scen in range(0,5,4):       #baseline and scenario 4
            fac_cc, fac_ctc, com_cc, com_ctc, com_cpad= model_scenarios(runs, scen, hbv_pars, fac_cov_weighted, com_cov_weighted, saturation)
            
            if scen==0:
                model_init=model_initializer(runs, globz, hbv_pars, globz_perturbed, regions_perturbed, fac_cc, fac_ctc, com_cc, com_ctc, com_cpad, t_steps,ve_warm)
                sensitivity_vcost, sensitivity_coverage=vaccine_costs(runs, hbv_pars, fac_cc, fac_ctc, com_cc, com_ctc, com_cpad,sat_costing, saturation, regions_perturbed)
                for reg in range(len(hbv_pars)):
                    model_init[reg,0,-2]=sensitivity_vcost[reg,5]
                sensivity_model=hbv_model(runs,t_steps,dt,model_init, hbv_pars, regions_perturbed, globz, globz_perturbed,discount)
                
                blcov_sens[1,0]=sensivity_model[6,len(t_steps)-1,13]      #DALYs
                blcov_sens[1,1]=sensivity_model[6,len(t_steps)-1,11]      #Disease Cost
                blcov_sens[1,2]=sensivity_model[6,0,12]                   #Vaccine Cost
                blcov_sens[1,3]=sensivity_model[6,len(t_steps)-1,12]      #Total Cost (for ICER)
            
            if scen==4:
                model_init=model_initializer(runs, globz, hbv_pars, globz_perturbed, regions_perturbed, fac_cc, fac_ctc, com_cc, com_ctc, com_cpad, t_steps,ve_warm)
                sensitivity_vcost, sensitivity_coverage=vaccine_costs(runs, hbv_pars, fac_cc, fac_ctc, com_cc, com_ctc, com_cpad,sat_costing, saturation, regions_perturbed)
                for reg in range(len(hbv_pars)):
                    model_init[reg,0,-2]=sensitivity_vcost[reg,5]
                sensitivity_model=hbv_model(runs,t_steps,dt,model_init, hbv_pars, regions_perturbed, globz, globz_perturbed,discount)
                
                blcov_sens[1,4]=sensitivity_model[6,len(t_steps)-1,13]      #DALYs
                blcov_sens[1,5]=sensitivity_model[6,len(t_steps)-1,11]      #Disease Cost
                blcov_sens[1,6]=sensitivity_model[6,0,12]                   #Vaccine Cost
                blcov_sens[1,7]=sensitivity_model[6,len(t_steps)-1,12]      #Total Cost (for ICER)
                
            blcov_sens[:,8]=-((blcov_sens[:,7]-blcov_sens[:,3])/(blcov_sens[:,4]-blcov_sens[:,0]))
        
#%% VE of Ambiently Stored Vaccines
        hbv_pars.iloc[6,6]=hbv_pars_copy[6,6]
        ve_sens=zeros((2,9))        #include baseline of 100%
        
        #20% VE
        ve_warm=0.2
        fac_cov_weighted, com_cov_weighted = facility_weight(runs, weight, hbv_pars, regions_perturbed)
        saturation, sat_costing=coverage_saturation(fac_cov_weighted, com_cov_weighted, hbv_pars, regions_perturbed,runs)
        for scen in range(0,5,4):
            fac_cc, fac_ctc, com_cc, com_ctc, com_cpad= model_scenarios(runs, scen, hbv_pars, fac_cov_weighted, com_cov_weighted, saturation)
            
            if scen==0:
                model_init=model_initializer(runs, globz, hbv_pars, globz_perturbed, regions_perturbed, fac_cc, fac_ctc, com_cc, com_ctc, com_cpad, t_steps,ve_warm)
                sensitivity_vcost, sensitivity_coverage=vaccine_costs(runs, hbv_pars, fac_cc, fac_ctc, com_cc, com_ctc, com_cpad,sat_costing, saturation, regions_perturbed)
                for reg in range(len(hbv_pars)):
                    model_init[reg,0,-2]=sensitivity_vcost[reg,5]
                sensivity_model=hbv_model(runs,t_steps,dt,model_init, hbv_pars, regions_perturbed, globz, globz_perturbed,discount)
                
                ve_sens[0,0]=sensivity_model[6,len(t_steps)-1,13]      #DALYs
                ve_sens[0,1]=sensivity_model[6,len(t_steps)-1,11]      #Disease Cost
                ve_sens[0,2]=sensivity_model[6,0,12]                   #Vaccine Cost
                ve_sens[0,3]=sensivity_model[6,len(t_steps)-1,12]      #Total Cost (for ICER)
            
            if scen==4:
                model_init=model_initializer(runs, globz, hbv_pars, globz_perturbed, regions_perturbed, fac_cc, fac_ctc, com_cc, com_ctc, com_cpad, t_steps,ve_warm)
                sensitivity_vcost, sensitivity_coverage=vaccine_costs(runs, hbv_pars, fac_cc, fac_ctc, com_cc, com_ctc, com_cpad,sat_costing, saturation, regions_perturbed)
                for reg in range(len(hbv_pars)):
                    model_init[reg,0,-2]=sensitivity_vcost[reg,5]
                sensitivity_model=hbv_model(runs,t_steps,dt,model_init, hbv_pars, regions_perturbed, globz, globz_perturbed,discount)
                
                ve_sens[0,4]=sensitivity_model[6,len(t_steps)-1,13]      #DALYs
                ve_sens[0,5]=sensitivity_model[6,len(t_steps)-1,11]      #Disease Cost
                ve_sens[0,6]=sensitivity_model[6,0,12]                   #Vaccine Cost
                ve_sens[0,7]=sensitivity_model[6,len(t_steps)-1,12]      #Total Cost (for ICER)
        
        
        #100% VE (Model Assumption)
        ve_sens[1,0]=baseline[6, len(t_steps)-1, 13]
        ve_sens[1,1]=baseline[6, len(t_steps)-1, 11]
        ve_sens[1,2]=baseline[6,0,12]
        ve_sens[1,3]=baseline[6, len(t_steps)-1, 12]
        
        ve_sens[1,4]=cpad_community[6, len(t_steps)-1, 13]
        ve_sens[1,5]=cpad_community[6, len(t_steps)-1, 11]
        ve_sens[1,6]=cpad_community[6, 0, 12]
        ve_sens[1,7]=cpad_community[6, len(t_steps)-1, 12]
        
        ve_sens[:,8]=-((ve_sens[:,7]-ve_sens[:,3])/(ve_sens[:,4]-ve_sens[:,0]))
#%% Cost of CTC licensure 
        ve_warm=1
        
        ctclic_sens=zeros((2,9))
        #Cost Equal to Cold Chain Vaccines
        hbv_pars.iloc[6,43]=hbv_pars_copy[6,43]-0.58 #factors in wastage
        hbv_pars.iloc[6,45]=hbv_pars_copy[6,45]-0.58
        
        fac_cov_weighted, com_cov_weighted = facility_weight(runs, weight, hbv_pars, regions_perturbed)
        saturation, sat_costing=coverage_saturation(fac_cov_weighted, com_cov_weighted, hbv_pars, regions_perturbed,runs)
        for scen in range(0,5,4):
            fac_cc, fac_ctc, com_cc, com_ctc, com_cpad= model_scenarios(runs, scen, hbv_pars, fac_cov_weighted, com_cov_weighted, saturation)
            
            if scen==0:
                model_init=model_initializer(runs, globz, hbv_pars, globz_perturbed, regions_perturbed, fac_cc, fac_ctc, com_cc, com_ctc, com_cpad, t_steps,ve_warm)
                sensitivity_vcost, sensitivity_coverage=vaccine_costs(runs, hbv_pars, fac_cc, fac_ctc, com_cc, com_ctc, com_cpad,sat_costing, saturation, regions_perturbed)
                for reg in range(len(hbv_pars)):
                    model_init[reg,0,-2]=sensitivity_vcost[reg,5]
                sensivity_model=hbv_model(runs,t_steps,dt,model_init, hbv_pars, regions_perturbed, globz, globz_perturbed,discount)
                
                ctclic_sens[0,0]=sensivity_model[6,len(t_steps)-1,13]      #DALYs
                ctclic_sens[0,1]=sensivity_model[6,len(t_steps)-1,11]      #Disease Cost
                ctclic_sens[0,2]=sensivity_model[6,0,12]                   #Vaccine Cost
                ctclic_sens[0,3]=sensivity_model[6,len(t_steps)-1,12]      #Total Cost (for ICER)
            
            if scen==4:
                model_init=model_initializer(runs, globz, hbv_pars, globz_perturbed, regions_perturbed, fac_cc, fac_ctc, com_cc, com_ctc, com_cpad, t_steps,ve_warm)
                sensitivity_vcost, sensitivity_coverage=vaccine_costs(runs, hbv_pars, fac_cc, fac_ctc, com_cc, com_ctc, com_cpad,sat_costing, saturation, regions_perturbed)
                for reg in range(len(hbv_pars)):
                    model_init[reg,0,-2]=sensitivity_vcost[reg,5]
                sensitivity_model=hbv_model(runs,t_steps,dt,model_init, hbv_pars, regions_perturbed, globz, globz_perturbed,discount)
                
                ctclic_sens[0,4]=sensitivity_model[6,len(t_steps)-1,13]      #DALYs
                ctclic_sens[0,5]=sensitivity_model[6,len(t_steps)-1,11]      #Disease Cost
                ctclic_sens[0,6]=sensitivity_model[6,0,12]                   #Vaccine Cost
                ctclic_sens[0,7]=sensitivity_model[6,len(t_steps)-1,12]      #Total Cost (for ICER)
        
        
        #Cost Quadruple Cold Chain Vaccines 
        hbv_pars.iloc[6,43]=hbv_pars_copy[6,43]+1.16 #factors in wastage
        hbv_pars.iloc[6,45]=hbv_pars_copy[6,45]+1.16
        
        fac_cov_weighted, com_cov_weighted = facility_weight(runs, weight, hbv_pars, regions_perturbed)
        saturation, sat_costing=coverage_saturation(fac_cov_weighted, com_cov_weighted, hbv_pars, regions_perturbed,runs)
        for scen in range(0,5,4):
            fac_cc, fac_ctc, com_cc, com_ctc, com_cpad= model_scenarios(runs, scen, hbv_pars, fac_cov_weighted, com_cov_weighted, saturation)
            
            if scen==0:
                model_init=model_initializer(runs, globz, hbv_pars, globz_perturbed, regions_perturbed, fac_cc, fac_ctc, com_cc, com_ctc, com_cpad, t_steps,ve_warm)
                sensitivity_vcost, sensitivity_coverage=vaccine_costs(runs, hbv_pars, fac_cc, fac_ctc, com_cc, com_ctc, com_cpad,sat_costing, saturation, regions_perturbed)
                for reg in range(len(hbv_pars)):
                    model_init[reg,0,-2]=sensitivity_vcost[reg,5]
                sensivity_model=hbv_model(runs,t_steps,dt,model_init, hbv_pars, regions_perturbed, globz, globz_perturbed,discount)
                
                ctclic_sens[1,0]=sensivity_model[6,len(t_steps)-1,13]      #DALYs
                ctclic_sens[1,1]=sensivity_model[6,len(t_steps)-1,11]      #Disease Cost
                ctclic_sens[1,2]=sensivity_model[6,0,12]                   #Vaccine Cost
                ctclic_sens[1,3]=sensivity_model[6,len(t_steps)-1,12]      #Total Cost (for ICER)
            
            if scen==4:
                model_init=model_initializer(runs, globz, hbv_pars, globz_perturbed, regions_perturbed, fac_cc, fac_ctc, com_cc, com_ctc, com_cpad, t_steps,ve_warm)
                sensitivity_vcost, sensitivity_coverage=vaccine_costs(runs, hbv_pars, fac_cc, fac_ctc, com_cc, com_ctc, com_cpad,sat_costing, saturation, regions_perturbed)
                for reg in range(len(hbv_pars)):
                    model_init[reg,0,-2]=sensitivity_vcost[reg,5]
                sensitivity_model=hbv_model(runs,t_steps,dt,model_init, hbv_pars, regions_perturbed, globz, globz_perturbed,discount)
                
                ctclic_sens[1,4]=sensitivity_model[6,len(t_steps)-1,13]      #DALYs
                ctclic_sens[1,5]=sensitivity_model[6,len(t_steps)-1,11]      #Disease Cost
                ctclic_sens[1,6]=sensitivity_model[6,0,12]                   #Vaccine Cost
                ctclic_sens[1,7]=sensitivity_model[6,len(t_steps)-1,12]      #Total Cost (for ICER)
        
        ctclic_sens[:,8]=-((ctclic_sens[:,7]-ctclic_sens[:,3])/(ctclic_sens[:,4]-ctclic_sens[:,0]))
#%% Cost of CPAD vaccine
        hbv_pars.iloc[6,43]=hbv_pars_copy[6,43]
        hbv_pars.iloc[6,45]=hbv_pars_copy[6,45]
        
        cpadc_sens=zeros((2,9))
        #Equal to Cold Chain Vaccines
        hbv_pars.iloc[6,46]=hbv_pars_copy[6,46]-1.10 #Factors in wastage
        
        fac_cov_weighted, com_cov_weighted = facility_weight(runs, weight, hbv_pars, regions_perturbed)
        saturation, sat_costing=coverage_saturation(fac_cov_weighted, com_cov_weighted, hbv_pars, regions_perturbed,runs)
        for scen in range(0,5,4):
            fac_cc, fac_ctc, com_cc, com_ctc, com_cpad= model_scenarios(runs, scen, hbv_pars, fac_cov_weighted, com_cov_weighted, saturation)
            
            if scen==0:
                model_init=model_initializer(runs, globz, hbv_pars, globz_perturbed, regions_perturbed, fac_cc, fac_ctc, com_cc, com_ctc, com_cpad, t_steps,ve_warm)
                sensitivity_vcost, sensitivity_coverage=vaccine_costs(runs, hbv_pars, fac_cc, fac_ctc, com_cc, com_ctc, com_cpad,sat_costing, saturation, regions_perturbed)
                for reg in range(len(hbv_pars)):
                    model_init[reg,0,-2]=sensitivity_vcost[reg,5]
                sensivity_model=hbv_model(runs,t_steps,dt,model_init, hbv_pars, regions_perturbed, globz, globz_perturbed,discount)
                
                cpadc_sens[0,0]=sensivity_model[6,len(t_steps)-1,13]      #DALYs
                cpadc_sens[0,1]=sensivity_model[6,len(t_steps)-1,11]      #Disease Cost
                cpadc_sens[0,2]=sensivity_model[6,0,12]                   #Vaccine Cost
                cpadc_sens[0,3]=sensivity_model[6,len(t_steps)-1,12]      #Total Cost (for ICER)
            
            if scen==4:
                model_init=model_initializer(runs, globz, hbv_pars, globz_perturbed, regions_perturbed, fac_cc, fac_ctc, com_cc, com_ctc, com_cpad, t_steps,ve_warm)
                sensitivity_vcost, sensitivity_coverage=vaccine_costs(runs, hbv_pars, fac_cc, fac_ctc, com_cc, com_ctc, com_cpad,sat_costing, saturation, regions_perturbed)
                for reg in range(len(hbv_pars)):
                    model_init[reg,0,-2]=sensitivity_vcost[reg,5]
                sensitivity_model=hbv_model(runs,t_steps,dt,model_init, hbv_pars, regions_perturbed, globz, globz_perturbed,discount)
                
                cpadc_sens[0,4]=sensitivity_model[6,len(t_steps)-1,13]      #DALYs
                cpadc_sens[0,5]=sensitivity_model[6,len(t_steps)-1,11]      #Disease Cost
                cpadc_sens[0,6]=sensitivity_model[6,0,12]                   #Vaccine Cost
                cpadc_sens[0,7]=sensitivity_model[6,len(t_steps)-1,12]      #Total Cost (for ICER)
        
        
        #Double Estimated
        hbv_pars.iloc[6,46]=hbv_pars_copy[6,46]+1.66 #Factors in wastage
        
        fac_cov_weighted, com_cov_weighted = facility_weight(runs, weight, hbv_pars, regions_perturbed)
        saturation, sat_costing=coverage_saturation(fac_cov_weighted, com_cov_weighted, hbv_pars, regions_perturbed,runs)
        for scen in range(0,5,4):
            fac_cc, fac_ctc, com_cc, com_ctc, com_cpad= model_scenarios(runs, scen, hbv_pars, fac_cov_weighted, com_cov_weighted, saturation)
            
            if scen==0:
                model_init=model_initializer(runs, globz, hbv_pars, globz_perturbed, regions_perturbed, fac_cc, fac_ctc, com_cc, com_ctc, com_cpad, t_steps,ve_warm)
                sensitivity_vcost, sensitivity_coverage=vaccine_costs(runs, hbv_pars, fac_cc, fac_ctc, com_cc, com_ctc, com_cpad,sat_costing, saturation, regions_perturbed)
                for reg in range(len(hbv_pars)):
                    model_init[reg,0,-2]=sensitivity_vcost[reg,5]
                sensivity_model=hbv_model(runs,t_steps,dt,model_init, hbv_pars, regions_perturbed, globz, globz_perturbed,discount)
                
                cpadc_sens[1,0]=sensivity_model[6,len(t_steps)-1,13]      #DALYs
                cpadc_sens[1,1]=sensivity_model[6,len(t_steps)-1,11]      #Disease Cost
                cpadc_sens[1,2]=sensivity_model[6,0,12]                   #Vaccine Cost
                cpadc_sens[1,3]=sensivity_model[6,len(t_steps)-1,12]      #Total Cost (for ICER)
            
            if scen==4:
                model_init=model_initializer(runs, globz, hbv_pars, globz_perturbed, regions_perturbed, fac_cc, fac_ctc, com_cc, com_ctc, com_cpad, t_steps,ve_warm)
                sensitivity_vcost, sensitivity_coverage=vaccine_costs(runs, hbv_pars, fac_cc, fac_ctc, com_cc, com_ctc, com_cpad,sat_costing, saturation, regions_perturbed)
                for reg in range(len(hbv_pars)):
                    model_init[reg,0,-2]=sensitivity_vcost[reg,5]
                sensitivity_model=hbv_model(runs,t_steps,dt,model_init, hbv_pars, regions_perturbed, globz, globz_perturbed,discount)
                
                cpadc_sens[1,4]=sensitivity_model[6,len(t_steps)-1,13]      #DALYs
                cpadc_sens[1,5]=sensitivity_model[6,len(t_steps)-1,11]      #Disease Cost
                cpadc_sens[1,6]=sensitivity_model[6,0,12]                   #Vaccine Cost
                cpadc_sens[1,7]=sensitivity_model[6,len(t_steps)-1,12]      #Total Cost (for ICER)
        
        cpadc_sens[:,8]=-((cpadc_sens[:,7]-cpadc_sens[:,3])/(cpadc_sens[:,4]-cpadc_sens[:,0]))

       
#%%  Cost of community delivery of CPAD   
        hbv_pars.iloc[6,46]=hbv_pars_copy[6,46]
        comdel_sens=zeros((2,9))        #add baseline representative of equal cost
        
        #20% of professional birth attendant
        hbv_pars.iloc[6,46]=hbv_pars_copy[6,46]-0.69
        
        fac_cov_weighted, com_cov_weighted = facility_weight(runs, weight, hbv_pars, regions_perturbed)
        saturation, sat_costing=coverage_saturation(fac_cov_weighted, com_cov_weighted, hbv_pars, regions_perturbed,runs)
        for scen in range(0,5,4):
            fac_cc, fac_ctc, com_cc, com_ctc, com_cpad= model_scenarios(runs, scen, hbv_pars, fac_cov_weighted, com_cov_weighted, saturation)
            
            if scen==0:
                model_init=model_initializer(runs, globz, hbv_pars, globz_perturbed, regions_perturbed, fac_cc, fac_ctc, com_cc, com_ctc, com_cpad, t_steps,ve_warm)
                sensitivity_vcost, sensitivity_coverage=vaccine_costs(runs, hbv_pars, fac_cc, fac_ctc, com_cc, com_ctc, com_cpad,sat_costing, saturation, regions_perturbed)
                for reg in range(len(hbv_pars)):
                    model_init[reg,0,-2]=sensitivity_vcost[reg,5]
                sensivity_model=hbv_model(runs,t_steps,dt,model_init, hbv_pars, regions_perturbed, globz, globz_perturbed,discount)
                
                comdel_sens[0,0]=sensivity_model[6,len(t_steps)-1,13]      #DALYs
                comdel_sens[0,1]=sensivity_model[6,len(t_steps)-1,11]      #Disease Cost
                comdel_sens[0,2]=sensivity_model[6,0,12]                   #Vaccine Cost
                comdel_sens[0,3]=sensivity_model[6,len(t_steps)-1,12]      #Total Cost (for ICER)
            
            if scen==4:
                model_init=model_initializer(runs, globz, hbv_pars, globz_perturbed, regions_perturbed, fac_cc, fac_ctc, com_cc, com_ctc, com_cpad, t_steps,ve_warm)
                sensitivity_vcost, sensitivity_coverage=vaccine_costs(runs, hbv_pars, fac_cc, fac_ctc, com_cc, com_ctc, com_cpad,sat_costing, saturation, regions_perturbed)
                for reg in range(len(hbv_pars)):
                    model_init[reg,0,-2]=sensitivity_vcost[reg,5]
                sensitivity_model=hbv_model(runs,t_steps,dt,model_init, hbv_pars, regions_perturbed, globz, globz_perturbed,discount)
                
                comdel_sens[0,4]=sensitivity_model[6,len(t_steps)-1,13]      #DALYs
                comdel_sens[0,5]=sensitivity_model[6,len(t_steps)-1,11]      #Disease Cost
                comdel_sens[0,6]=sensitivity_model[6,0,12]                   #Vaccine Cost
                comdel_sens[0,7]=sensitivity_model[6,len(t_steps)-1,12]      #Total Cost (for ICER)
        
        #Equal to professional birth attendant
         #100% VE (Model Assumption)
        comdel_sens[1,0]=baseline[6, len(t_steps)-1, 13]
        comdel_sens[1,1]=baseline[6, len(t_steps)-1, 11]
        comdel_sens[1,2]=baseline[6,0,12]
        comdel_sens[1,3]=baseline[6, len(t_steps)-1, 12]
        
        comdel_sens[1,4]=cpad_community[6, len(t_steps)-1, 13]
        comdel_sens[1,5]=cpad_community[6, len(t_steps)-1, 11]
        comdel_sens[1,6]=cpad_community[6, 0, 12]
        comdel_sens[1,7]=cpad_community[6, len(t_steps)-1, 12]
        
        comdel_sens[:,8]=-((comdel_sens[:,7]-comdel_sens[:,3])/(comdel_sens[:,4]-comdel_sens[:,0]))
        
#%% Training for CTC and CPAD costs   #up to here
        hbv_pars.iloc[6,46]=hbv_pars_copy[6,46]
        training_sens=zeros((2,9))
        
        #50% of assumed
        hbv_pars.iloc[6,43]=hbv_pars_copy[6,43]-0.02
        hbv_pars.iloc[6,45]=hbv_pars_copy[6,45]-0.02
        hbv_pars.iloc[6,46]=hbv_pars_copy[6,46]-0.05
        
        fac_cov_weighted, com_cov_weighted = facility_weight(runs, weight, hbv_pars, regions_perturbed)
        saturation, sat_costing=coverage_saturation(fac_cov_weighted, com_cov_weighted, hbv_pars, regions_perturbed,runs)
        for scen in range(0,5,4):
            fac_cc, fac_ctc, com_cc, com_ctc, com_cpad= model_scenarios(runs, scen, hbv_pars, fac_cov_weighted, com_cov_weighted, saturation)
            
            if scen==0:
                model_init=model_initializer(runs, globz, hbv_pars, globz_perturbed, regions_perturbed, fac_cc, fac_ctc, com_cc, com_ctc, com_cpad, t_steps,ve_warm)
                sensitivity_vcost, sensitivity_coverage=vaccine_costs(runs, hbv_pars, fac_cc, fac_ctc, com_cc, com_ctc, com_cpad,sat_costing, saturation, regions_perturbed)
                for reg in range(len(hbv_pars)):
                    model_init[reg,0,-2]=sensitivity_vcost[reg,5]
                sensivity_model=hbv_model(runs,t_steps,dt,model_init, hbv_pars, regions_perturbed, globz, globz_perturbed,discount)
                
                training_sens[0,0]=sensivity_model[6,len(t_steps)-1,13]      #DALYs
                training_sens[0,1]=sensivity_model[6,len(t_steps)-1,11]      #Disease Cost
                training_sens[0,2]=sensivity_model[6,0,12]                   #Vaccine Cost
                training_sens[0,3]=sensivity_model[6,len(t_steps)-1,12]      #Total Cost (for ICER)
            
            if scen==4:
                model_init=model_initializer(runs, globz, hbv_pars, globz_perturbed, regions_perturbed, fac_cc, fac_ctc, com_cc, com_ctc, com_cpad, t_steps,ve_warm)
                sensitivity_vcost, sensitivity_coverage=vaccine_costs(runs, hbv_pars, fac_cc, fac_ctc, com_cc, com_ctc, com_cpad,sat_costing, saturation, regions_perturbed)
                for reg in range(len(hbv_pars)):
                    model_init[reg,0,-2]=sensitivity_vcost[reg,5]
                sensitivity_model=hbv_model(runs,t_steps,dt,model_init, hbv_pars, regions_perturbed, globz, globz_perturbed,discount)
                
                training_sens[0,4]=sensitivity_model[6,len(t_steps)-1,13]      #DALYs
                training_sens[0,5]=sensitivity_model[6,len(t_steps)-1,11]      #Disease Cost
                training_sens[0,6]=sensitivity_model[6,0,12]                   #Vaccine Cost
                training_sens[0,7]=sensitivity_model[6,len(t_steps)-1,12]      #Total Cost (for ICER)
        
        
        #10x assumed cost (0.4USD and 1.0USD)
        hbv_pars.iloc[6,43]=hbv_pars_copy[6,43]+0.36
        hbv_pars.iloc[6,45]=hbv_pars_copy[6,45]+0.36
        hbv_pars.iloc[6,46]=hbv_pars_copy[6,46]+0.90
        
        fac_cov_weighted, com_cov_weighted = facility_weight(runs, weight, hbv_pars, regions_perturbed)
        saturation, sat_costing=coverage_saturation(fac_cov_weighted, com_cov_weighted, hbv_pars, regions_perturbed,runs)
        for scen in range(0,5,4):
            fac_cc, fac_ctc, com_cc, com_ctc, com_cpad= model_scenarios(runs, scen, hbv_pars, fac_cov_weighted, com_cov_weighted, saturation)
            
            if scen==0:
                model_init=model_initializer(runs, globz, hbv_pars, globz_perturbed, regions_perturbed, fac_cc, fac_ctc, com_cc, com_ctc, com_cpad, t_steps,ve_warm)
                sensitivity_vcost, sensitivity_coverage=vaccine_costs(runs, hbv_pars, fac_cc, fac_ctc, com_cc, com_ctc, com_cpad,sat_costing, saturation, regions_perturbed)
                for reg in range(len(hbv_pars)):
                    model_init[reg,0,-2]=sensitivity_vcost[reg,5]
                sensivity_model=hbv_model(runs,t_steps,dt,model_init, hbv_pars, regions_perturbed, globz, globz_perturbed,discount)
                
                training_sens[1,0]=sensivity_model[6,len(t_steps)-1,13]      #DALYs
                training_sens[1,1]=sensivity_model[6,len(t_steps)-1,11]      #Disease Cost
                training_sens[1,2]=sensivity_model[6,0,12]                   #Vaccine Cost
                training_sens[1,3]=sensivity_model[6,len(t_steps)-1,12]      #Total Cost (for ICER)
            
            if scen==4:
                model_init=model_initializer(runs, globz, hbv_pars, globz_perturbed, regions_perturbed, fac_cc, fac_ctc, com_cc, com_ctc, com_cpad, t_steps,ve_warm)
                sensitivity_vcost, sensitivity_coverage=vaccine_costs(runs, hbv_pars, fac_cc, fac_ctc, com_cc, com_ctc, com_cpad,sat_costing, saturation, regions_perturbed)
                for reg in range(len(hbv_pars)):
                    model_init[reg,0,-2]=sensitivity_vcost[reg,5]
                sensitivity_model=hbv_model(runs,t_steps,dt,model_init, hbv_pars, regions_perturbed, globz, globz_perturbed,discount)
                
                training_sens[1,4]=sensitivity_model[6,len(t_steps)-1,13]      #DALYs
                training_sens[1,5]=sensitivity_model[6,len(t_steps)-1,11]      #Disease Cost
                training_sens[1,6]=sensitivity_model[6,0,12]                   #Vaccine Cost
                training_sens[1,7]=sensitivity_model[6,len(t_steps)-1,12]      #Total Cost (for ICER)
        
        training_sens[:,8]=-((training_sens[:,7]-training_sens[:,3])/(training_sens[:,4]-training_sens[:,0]))
#%%Return HBV pars back to base values
        
        hbv_pars.iloc[6,43]=hbv_pars_copy[6,43]
        hbv_pars.iloc[6,45]=hbv_pars_copy[6,45]
        hbv_pars.iloc[6,46]=hbv_pars_copy[6,46]