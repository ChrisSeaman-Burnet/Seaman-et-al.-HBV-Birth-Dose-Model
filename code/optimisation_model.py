from numpy import zeros, array
from scipy.optimize import minimize
import pandas as pd


################################### Global Burden of Disease Regions, with Uncertainity#######################################################
if region==1: 
    if runs > 1: 
           
    #Fixed Parameters for comparison 
        for reg in range(len(hbv_pars)):
            regions_perturbed[reg,4,:]=hbv_pars.iloc[reg,14] #Facility Births
            regions_perturbed[reg,1,:]=hbv_pars.iloc[reg, 6] #Baseline Birth Dose Coverage
            regions_perturbed[reg,37,:]=hbv_pars.iloc[reg,47] #Professional Health Workers at Facility Births
            regions_perturbed[reg,32,:]=hbv_pars.iloc[reg,42] #Cost of cold chain vaccine at facility
            regions_perturbed[reg,33,:]=hbv_pars.iloc[reg,43] #Cost of CTC vaccine at facility
            regions_perturbed[reg,34,:]=hbv_pars.iloc[reg,44] #Cost of cold-chain vaccine in community
            regions_perturbed[reg,35,:]=hbv_pars.iloc[reg,45] #Cost of CTC vaccine in community
            regions_perturbed[reg,36,:]=hbv_pars.iloc[reg,46] #Cost of CPAD vaccine in community
        
        # Remodelling baseline with fixed paramaters for comparison (ICER) - only needed when runs > 1
        fac_cov_weighted, com_cov_weighted = facility_weight(runs, weight, hbv_pars, regions_perturbed)
        saturation, sat_costing=coverage_saturation(fac_cov_weighted, com_cov_weighted, hbv_pars, regions_perturbed,runs)
        scen=0
        fac_cc, fac_ctc, com_cc, com_ctc, com_cpad= model_scenarios(runs, scen, hbv_pars, fac_cov_weighted, com_cov_weighted, saturation)
        model_init=model_initializer(runs, globz, hbv_pars, globz_perturbed, regions_perturbed, fac_cc, fac_ctc, com_cc, com_ctc, com_cpad, t_steps, ve_warm)
        baseline_opt_vcost, baseline_opt_coverage=vaccine_costs(runs, hbv_pars, fac_cc, fac_ctc, com_cc, com_ctc, com_cpad,sat_costing, saturation, regions_perturbed)
        
        for reg in range(len(hbv_pars)):
            for run in range (runs):
                model_init[(runs*reg)+run,0,-2]=baseline_opt_vcost[reg,run,5]
    
        baseline_opt=hbv_model(runs,t_steps,dt,model_init, hbv_pars, regions_perturbed, globz, globz_perturbed,discount)
    
            
        # Collating data required for optimisation (formatting needs to be constant for runs==1 and runs >1; region==0 and region==1 )
        baseline_expenditure=zeros((len(hbv_pars),1))
        for reg in range(len(hbv_pars)):
            baseline_expenditure[reg,0]=baseline_opt_vcost[reg,0,5]
            
        #Series of one-dimensional arrays so can loop through 
        hosp_n=zeros((len(hbv_pars)))
        for reg in range(len(hbv_pars)):
            hosp_n[reg]=hbv_pars.iloc[reg,14]*1000
            
        comm_n=1000-hosp_n
            
        vac_cost_opt=zeros((len(hbv_pars),5))
            
        for reg in range(len(hbv_pars)):
            vac_cost_opt[reg,0]=hbv_pars.iloc[reg,42]
            vac_cost_opt[reg,1]=hbv_pars.iloc[reg,43]
            vac_cost_opt[reg,2]=hbv_pars.iloc[reg,44]
            vac_cost_opt[reg,3]=hbv_pars.iloc[reg,45]
            vac_cost_opt[reg,4]=hbv_pars.iloc[reg,46]
            
            
        sat_opt=zeros((len(hbv_pars),5))
            
        for reg in range(len(hbv_pars)):
            sat_opt[reg,0]=saturation[reg,0,0]
            sat_opt[reg,1]=saturation[reg,0,1]
            sat_opt[reg,2]=saturation[reg,0,2]
            sat_opt[reg,3]=saturation[reg,0,3]
            sat_opt[reg,4]=saturation[reg,0,4]
            
        #Optimisation of Baseline Expenditure
            
        baseline_optimised_coverage, baseline_optimised_expenditure=[], []
        for reg in range(len(hbv_pars)):
            cons=({'type':'eq', 'fun': constraint_baseline_budget, 'args': (reg,)},
                  {'type': 'ineq', 'fun': constraint_fac_cc_coverage, 'args': (sat_opt, reg)},
                  {'type': 'ineq', 'fun': constraint_fac_cc_ctc_coverage, 'args': (sat_opt, reg)},
                  {'type': 'ineq', 'fun': constraint_com_cc_ctc_coverage, 'args': (sat_opt, reg)},
                  {'type': 'ineq', 'fun': constraint_com_cc_ctc_cpad_coverage, 'args': (sat_opt, reg)})
            bnds=((0,None), (0,None), (0,None), (0,None), (0,None))
            expenditure=array([baseline_expenditure[reg,0]/5, baseline_expenditure[reg,0]/5,baseline_expenditure[reg,0]/5,baseline_expenditure[reg,0]/5,baseline_expenditure[reg,0]/5])
            sol=minimize(objective_coverage, expenditure, method='SLSQP', constraints=cons, bounds=bnds, options={'disp':False})
            expenditure=array([sol.x[0], sol.x[1], sol.x[2], sol.x[3], sol.x[4]])
            baseline_optimised_expenditure.append([sol.x[0], sol.x[1], sol.x[2], sol.x[3], sol.x[4],sol.x[0]+ sol.x[1]+sol.x[2]+sol.x[3]+sol.x[4]])
            baseline_optimised_coverage.append(coverage_from_opt(expenditure))
            
        baseline_optimised_expenditure=array(baseline_optimised_expenditure)
        baseline_optimised_coverage=array(baseline_optimised_coverage)
            
        # Run the model for optimised baseline expenditure
            
        fac_cc=baseline_optimised_coverage[:,1].repeat(runs).reshape(len(hbv_pars),runs,1)
        fac_ctc=baseline_optimised_coverage[:,3].repeat(runs).reshape(len(hbv_pars),runs,1)
        com_cc=baseline_optimised_coverage[:,5].repeat(runs).reshape(len(hbv_pars),runs,1)
        com_ctc=baseline_optimised_coverage[:,7].repeat(runs).reshape(len(hbv_pars),runs,1)
        com_cpad=baseline_optimised_coverage[:,9].repeat(runs).reshape(len(hbv_pars),runs,1)
        
        model_init=model_initializer(runs, globz, hbv_pars, globz_perturbed, regions_perturbed, fac_cc, fac_ctc, com_cc, com_ctc, com_cpad, t_steps, ve_warm)
        
        for reg in range(len(hbv_pars)):
            for run in range (runs):
                model_init[(runs*reg)+run,0,-2]=baseline_optimised_expenditure[reg,5]
    
        baseline_optimised=hbv_model(runs,t_steps,dt,model_init, hbv_pars, regions_perturbed, globz, globz_perturbed,discount)
            
            
        #Optimisation to reach 90% coverage
        optimised_90_coverage, optimised_90_expenditure=[], []
        for reg in range(len(hbv_pars)):
            cons=({'type':'eq', 'fun': constraint_90_coverage},
                  {'type': 'ineq', 'fun': constraint_fac_cc_coverage, 'args': (sat_opt, reg)},
                  {'type': 'ineq', 'fun': constraint_fac_cc_ctc_coverage, 'args': (sat_opt, reg)},
                  {'type': 'ineq', 'fun': constraint_com_cc_ctc_coverage, 'args': (sat_opt, reg)},
                  {'type': 'ineq', 'fun': constraint_com_cc_ctc_cpad_coverage, 'args': (sat_opt, reg)})
            bnds=((0,None), (0,None), (0,None), (0,None), (0,None))
            expenditure=array([1,1,1,1,1])
            sol=minimize(objective_90_coverage, expenditure, method='SLSQP', constraints=cons, bounds=bnds, options={'disp':False})
            expenditure=array([sol.x[0], sol.x[1], sol.x[2], sol.x[3], sol.x[4]])
            optimised_90_expenditure.append([sol.x[0], sol.x[1], sol.x[2], sol.x[3], sol.x[4],sol.x[0]+ sol.x[1]+sol.x[2]+sol.x[3]+sol.x[4]])
            optimised_90_coverage.append(coverage_from_opt(expenditure))
            
        optimised_90_expenditure=array(optimised_90_expenditure)
        optimised_90_coverage=array(optimised_90_coverage)
            
        #Run the model for optimised 90% coverage
            
        fac_cc=optimised_90_coverage[:,1].repeat(runs).reshape(len(hbv_pars),runs,1)
        fac_ctc=optimised_90_coverage[:,3].repeat(runs).reshape(len(hbv_pars),runs,1)
        com_cc=optimised_90_coverage[:,5].repeat(runs).reshape(len(hbv_pars),runs,1)
        com_ctc=optimised_90_coverage[:,7].repeat(runs).reshape(len(hbv_pars),runs,1)
        com_cpad=optimised_90_coverage[:,9].repeat(runs).reshape(len(hbv_pars),runs,1)
        
        model_init=model_initializer(runs, globz, hbv_pars, globz_perturbed, regions_perturbed, fac_cc, fac_ctc, com_cc, com_ctc, com_cpad, t_steps, ve_warm)
        
        for reg in range(len(hbv_pars)):
            for run in range (runs):
                model_init[(runs*reg)+run,0,-2]=optimised_90_expenditure[reg,5]
    
        optimised_90=hbv_model(runs,t_steps,dt,model_init, hbv_pars, regions_perturbed, globz, globz_perturbed,discount)
            
###################################### Global Burden of Disease Regions, No Uncertainity########################################################################
    elif runs ==1:
            
        #Preparing the variables:
        baseline_expenditure =baseline_vcost[:,5].reshape(len(hbv_pars),1)
        
        hosp_n=zeros((len(hbv_pars)))
            
        for reg in range(len(hbv_pars)):
            hosp_n[reg]=hbv_pars.iloc[reg,14]*1000
            
        comm_n=1000-hosp_n
        
        vac_cost_opt=zeros((len(hbv_pars),5))
        for reg in range(len(hbv_pars)):
            vac_cost_opt[reg,0]=hbv_pars.iloc[reg,42]
            vac_cost_opt[reg,1]=hbv_pars.iloc[reg,43]
            vac_cost_opt[reg,2]=hbv_pars.iloc[reg,44]
            vac_cost_opt[reg,3]=hbv_pars.iloc[reg,45]
            vac_cost_opt[reg,4]=hbv_pars.iloc[reg,46]

        sat_opt=zeros((len(hbv_pars),5))
        for reg in range(len(hbv_pars)):
            sat_opt[reg,0]=saturation[reg,0]
            sat_opt[reg,1]=saturation[reg,1]
            sat_opt[reg,2]=saturation[reg,2]
            sat_opt[reg,3]=saturation[reg,3]
            sat_opt[reg,4]=saturation[reg,4]
            
            
        #Baseline Optimisation:
        baseline_optimised_coverage, baseline_optimised_expenditure=[], []
        for reg in range(len(hbv_pars)):
            cons=({'type':'eq', 'fun': constraint_baseline_budget, 'args': (reg,)},
                  {'type': 'ineq', 'fun': constraint_fac_cc_coverage, 'args': (sat_opt, reg)},
                  {'type': 'ineq', 'fun': constraint_fac_cc_ctc_coverage, 'args': (sat_opt, reg)},
                  {'type': 'ineq', 'fun': constraint_com_cc_ctc_coverage, 'args': (sat_opt, reg)},
                  {'type': 'ineq', 'fun': constraint_com_cc_ctc_cpad_coverage, 'args': (sat_opt, reg)})
            bnds=((0,None), (0,None), (0,None), (0,None), (0,None))
            expenditure=array([baseline_expenditure[reg,0]/5, baseline_expenditure[reg,0]/5,baseline_expenditure[reg,0]/5,baseline_expenditure[reg,0]/5,baseline_expenditure[reg,0]/5])
            sol=minimize(objective_coverage, expenditure, method='SLSQP', constraints=cons, bounds=bnds, options={'disp':False})
            expenditure=array([sol.x[0], sol.x[1], sol.x[2], sol.x[3], sol.x[4]])
            baseline_optimised_expenditure.append([sol.x[0], sol.x[1], sol.x[2], sol.x[3], sol.x[4],sol.x[0]+ sol.x[1]+sol.x[2]+sol.x[3]+sol.x[4]])
            baseline_optimised_coverage.append(coverage_from_opt(expenditure))
            
        baseline_optimised_expenditure=array(baseline_optimised_expenditure)
        baseline_optimised_coverage=array(baseline_optimised_coverage)
            
        fac_cc=baseline_optimised_coverage[:,1]
        fac_ctc=baseline_optimised_coverage[:,3]
        com_cc=baseline_optimised_coverage[:,5]
        com_ctc=baseline_optimised_coverage[:,7]
        com_cpad=baseline_optimised_coverage[:,9]
            
        model_init=model_initializer(runs, globz, hbv_pars, globz_perturbed, regions_perturbed, fac_cc, fac_ctc, com_cc, com_ctc, com_cpad, t_steps, ve_warm)
            
        for reg in range(len(hbv_pars)):
            model_init[reg,0,-2]=baseline_optimised_expenditure[reg,5]
            
        baseline_optimised=hbv_model(runs,t_steps,dt,model_init, hbv_pars, regions_perturbed, globz, globz_perturbed,discount)
            
            
        #90% Optimisation:
        optimised_90_coverage, optimised_90_expenditure=[], []
        for reg in range(len(hbv_pars)):
            cons=({'type':'eq', 'fun': constraint_90_coverage},
                  {'type': 'ineq', 'fun': constraint_fac_cc_coverage, 'args': (sat_opt, reg)},
                  {'type': 'ineq', 'fun': constraint_fac_cc_ctc_coverage, 'args': (sat_opt, reg)},
                  {'type': 'ineq', 'fun': constraint_com_cc_ctc_coverage, 'args': (sat_opt, reg)},
                  {'type': 'ineq', 'fun': constraint_com_cc_ctc_cpad_coverage, 'args': (sat_opt, reg)})
            bnds=((0,None), (0,None), (0,None), (0,None), (0,None))
            expenditure=array([1,1,1,1,1])
            sol=minimize(objective_90_coverage, expenditure, method='SLSQP', constraints=cons, bounds=bnds, options={'disp':False})
            expenditure=array([sol.x[0], sol.x[1], sol.x[2], sol.x[3], sol.x[4]])
            optimised_90_expenditure.append([sol.x[0], sol.x[1], sol.x[2], sol.x[3], sol.x[4],sol.x[0]+ sol.x[1]+sol.x[2]+sol.x[3]+sol.x[4]])
            optimised_90_coverage.append(coverage_from_opt(expenditure))
            
        optimised_90_expenditure=array(optimised_90_expenditure)
        optimised_90_coverage=array(optimised_90_coverage)
            
        fac_cc=optimised_90_coverage[:,1]
        fac_ctc=optimised_90_coverage[:,3]
        com_cc=optimised_90_coverage[:,5]
        com_ctc=optimised_90_coverage[:,7]
        com_cpad=optimised_90_coverage[:,9]
            
        model_init=model_initializer(runs, globz, hbv_pars, globz_perturbed, regions_perturbed, fac_cc, fac_ctc, com_cc, com_ctc, com_cpad, t_steps, ve_warm)
            
        for reg in range(len(hbv_pars)):
            model_init[reg,0,-2]=optimised_90_expenditure[reg,5]
            
        optimised_90=hbv_model(runs,t_steps,dt,model_init, hbv_pars, regions_perturbed, globz, globz_perturbed,discount)
            

##############################################Individual LMICs, No Uncertainity Analysis################################################
elif region==0:
    if runs ==1:
            #Remove Countries with >= 90% HBV BD Coverage 
        countries,relevant_pars,baseline_expenditure, sat_opt, baseline_rel, baseline_cov_rel =[],[],[],[],[],[]
        for idx,val in enumerate(hbv_pars["birth_dose"]):
            if val < 0.90:
                countries.append(settings[idx])
                relevant_pars.append(hbv_pars.iloc[idx,:])
                baseline_expenditure.append(baseline_vcost[idx,5])
                sat_opt.append(saturation[idx,:])
                baseline_rel.append(baseline[idx,:,:])
                baseline_cov_rel.append(baseline_coverage[idx])
                    
        hbv_pars=pd.DataFrame(relevant_pars)
        baseline_expenditure=array(baseline_expenditure).reshape(len(hbv_pars),1)
        sat_opt=array(sat_opt)
        baseline_comp=array(baseline_rel)
        baseline_cov_comp=array(baseline_cov_rel)
            
    
        vac_cost_opt=zeros((len(hbv_pars),5))
        for reg in range(len(hbv_pars)):
            vac_cost_opt[reg,0]=hbv_pars.iloc[reg,42]
            vac_cost_opt[reg,1]=hbv_pars.iloc[reg,43]
            vac_cost_opt[reg,2]=hbv_pars.iloc[reg,44]
            vac_cost_opt[reg,3]=hbv_pars.iloc[reg,45]
            vac_cost_opt[reg,4]=hbv_pars.iloc[reg,46]

        hosp_n=zeros((len(hbv_pars)))
        for reg in range(len(hbv_pars)):
            hosp_n[reg]=hbv_pars.iloc[reg,14]*1000
            
        comm_n=1000-hosp_n
           
        #Baseline Optimisation:
        baseline_optimised_coverage, baseline_optimised_expenditure=[], []
        for reg in range(len(hbv_pars)):
            cons=({'type':'eq', 'fun': constraint_baseline_budget, 'args': (reg,)},
                  {'type': 'ineq', 'fun': constraint_fac_cc_coverage, 'args': (sat_opt, reg)},
                  {'type': 'ineq', 'fun': constraint_fac_cc_ctc_coverage, 'args': (sat_opt, reg)},
                  {'type': 'ineq', 'fun': constraint_com_cc_ctc_coverage, 'args': (sat_opt, reg)},
                  {'type': 'ineq', 'fun': constraint_com_cc_ctc_cpad_coverage, 'args': (sat_opt, reg)})
            bnds=((0,None), (0,None), (0,None), (0,None), (0,None))
            expenditure=array([baseline_expenditure[reg,0]/5, baseline_expenditure[reg,0]/5,baseline_expenditure[reg,0]/5,baseline_expenditure[reg,0]/5,baseline_expenditure[reg,0]/5])
            sol=minimize(objective_coverage, expenditure, method='SLSQP', constraints=cons, bounds=bnds, options={'disp':False})
            expenditure=array([sol.x[0], sol.x[1], sol.x[2], sol.x[3], sol.x[4]])
            baseline_optimised_expenditure.append([sol.x[0], sol.x[1], sol.x[2], sol.x[3], sol.x[4],sol.x[0]+ sol.x[1]+sol.x[2]+sol.x[3]+sol.x[4]])
            baseline_optimised_coverage.append(coverage_from_opt(expenditure))
            
        baseline_optimised_expenditure=array(baseline_optimised_expenditure)
        baseline_optimised_coverage=array(baseline_optimised_coverage)
        
        fac_cc=baseline_optimised_coverage[:,1]
        fac_ctc=baseline_optimised_coverage[:,3]
        com_cc=baseline_optimised_coverage[:,5]
        com_ctc=baseline_optimised_coverage[:,7]
        com_cpad=baseline_optimised_coverage[:,9]
            
        model_init=model_initializer(runs, globz, hbv_pars, globz_perturbed, regions_perturbed, fac_cc, fac_ctc, com_cc, com_ctc, com_cpad, t_steps, ve_warm)
            
        for reg in range(len(hbv_pars)):
            model_init[reg,0,-2]=baseline_optimised_expenditure[reg,5]
            
        baseline_optimised=hbv_model(runs,t_steps,dt,model_init, hbv_pars, regions_perturbed, globz, globz_perturbed,discount)
            
            
        #90% Optimisation:
        optimised_90_coverage, optimised_90_expenditure=[], []
        for reg in range(len(hbv_pars)):
            cons=({'type':'eq', 'fun': constraint_90_coverage},
                  {'type': 'ineq', 'fun': constraint_fac_cc_coverage, 'args': (sat_opt, reg)},
                  {'type': 'ineq', 'fun': constraint_fac_cc_ctc_coverage, 'args': (sat_opt, reg)},
                  {'type': 'ineq', 'fun': constraint_com_cc_ctc_coverage, 'args': (sat_opt, reg)},
                  {'type': 'ineq', 'fun': constraint_com_cc_ctc_cpad_coverage, 'args': (sat_opt, reg)})
            bnds=((0,None), (0,None), (0,None), (0,None), (0,None))
            expenditure=array([1,1,1,1,1])
            sol=minimize(objective_90_coverage, expenditure, method='SLSQP', constraints=cons, bounds=bnds, options={'disp':False})
            expenditure=array([sol.x[0], sol.x[1], sol.x[2], sol.x[3], sol.x[4]])
            optimised_90_expenditure.append([sol.x[0], sol.x[1], sol.x[2], sol.x[3], sol.x[4],sol.x[0]+ sol.x[1]+sol.x[2]+sol.x[3]+sol.x[4]])
            optimised_90_coverage.append(coverage_from_opt(expenditure))
            
        optimised_90_expenditure=array(optimised_90_expenditure)
        optimised_90_coverage=array(optimised_90_coverage)
            
        fac_cc=optimised_90_coverage[:,1]
        fac_ctc=optimised_90_coverage[:,3]
        com_cc=optimised_90_coverage[:,5]
        com_ctc=optimised_90_coverage[:,7]
        com_cpad=optimised_90_coverage[:,9]
            
        model_init=model_initializer(runs, globz, hbv_pars, globz_perturbed, regions_perturbed, fac_cc, fac_ctc, com_cc, com_ctc, com_cpad, t_steps, ve_warm)
            
        for reg in range(len(hbv_pars)):
            model_init[reg,0,-2]=optimised_90_expenditure[reg,5]
            
        optimised_90=hbv_model(runs,t_steps,dt,model_init, hbv_pars, regions_perturbed, globz, globz_perturbed,discount)
        
    
    elif runs > 1:
        
        regs_pert,hbv_pars_t,countries =[],[],[]
        for idx, val in enumerate(hbv_pars["birth_dose"]):
            if val < 0.9:
                regs_pert.append(regions_perturbed[idx,:,:])
                hbv_pars_t.append(hbv_pars.iloc[idx,:])
                countries.append(settings[idx])
                
        regions_perturbed=array((regs_pert))
        hbv_pars=pd.DataFrame(hbv_pars_t)
        

        for reg in range(len(hbv_pars)):
            regions_perturbed[reg,4,:]=hbv_pars.iloc[reg,14] #Facility Births
            regions_perturbed[reg,1,:]=hbv_pars.iloc[reg, 6] #Baseline Birth Dose Coverage
            regions_perturbed[reg,37,:]=hbv_pars.iloc[reg,47] #Professional Health Workers at Facility Births
            regions_perturbed[reg,32,:]=hbv_pars.iloc[reg,42] #Cost of cold chain vaccine at facility
            regions_perturbed[reg,33,:]=hbv_pars.iloc[reg,43] #Cost of CTC vaccine at facility
            regions_perturbed[reg,34,:]=hbv_pars.iloc[reg,44] #Cost of cold-chain vaccine in community
            regions_perturbed[reg,35,:]=hbv_pars.iloc[reg,45] #Cost of CTC vaccine in community
            regions_perturbed[reg,36,:]=hbv_pars.iloc[reg,46] #Cost of CPAD vaccine in community
        
        fac_cov_weighted, com_cov_weighted = facility_weight(runs, weight, hbv_pars, regions_perturbed)
        saturation, sat_costing=coverage_saturation(fac_cov_weighted, com_cov_weighted, hbv_pars, regions_perturbed,runs)
        scen=0
        fac_cc, fac_ctc, com_cc, com_ctc, com_cpad= model_scenarios(runs, scen, hbv_pars, fac_cov_weighted, com_cov_weighted, saturation)
        model_init=model_initializer(runs, globz, hbv_pars, globz_perturbed, regions_perturbed, fac_cc, fac_ctc, com_cc, com_ctc, com_cpad, t_steps, ve_warm)
        baseline_opt_vcost, baseline_opt_coverage=vaccine_costs(runs, hbv_pars, fac_cc, fac_ctc, com_cc, com_ctc, com_cpad,sat_costing, saturation, regions_perturbed)
        
        for reg in range(len(hbv_pars)):
            for run in range (runs):
                model_init[(runs*reg)+run,0,-2]=baseline_opt_vcost[reg,run,5]
    
        baseline_opt=hbv_model(runs,t_steps,dt,model_init, hbv_pars, regions_perturbed, globz, globz_perturbed,discount)
    
            
        baseline_expenditure=zeros((len(hbv_pars),1))
        for reg in range(len(hbv_pars)):
            baseline_expenditure[reg,0]=baseline_opt_vcost[reg,0,5]
            
        #Series of one-dimensional arrays so can loop through 
        hosp_n=zeros((len(hbv_pars)))
        for reg in range(len(hbv_pars)):
            hosp_n[reg]=hbv_pars.iloc[reg,14]*1000
            
        comm_n=1000-hosp_n
            
        vac_cost_opt=zeros((len(hbv_pars),5))
            
        for reg in range(len(hbv_pars)):
            vac_cost_opt[reg,0]=hbv_pars.iloc[reg,42]
            vac_cost_opt[reg,1]=hbv_pars.iloc[reg,43]
            vac_cost_opt[reg,2]=hbv_pars.iloc[reg,44]
            vac_cost_opt[reg,3]=hbv_pars.iloc[reg,45]
            vac_cost_opt[reg,4]=hbv_pars.iloc[reg,46]
            
            
        sat_opt=zeros((len(hbv_pars),5))
            
        for reg in range(len(hbv_pars)):
            sat_opt[reg,0]=saturation[reg,0,0]
            sat_opt[reg,1]=saturation[reg,0,1]
            sat_opt[reg,2]=saturation[reg,0,2]
            sat_opt[reg,3]=saturation[reg,0,3]
            sat_opt[reg,4]=saturation[reg,0,4]
            
        #Optimisation of Baseline Expenditure
            
        baseline_optimised_coverage, baseline_optimised_expenditure=[], []
        for reg in range(len(hbv_pars)):
            cons=({'type':'eq', 'fun': constraint_baseline_budget, 'args': (reg,)},
                  {'type': 'ineq', 'fun': constraint_fac_cc_coverage, 'args': (sat_opt, reg)},
                  {'type': 'ineq', 'fun': constraint_fac_cc_ctc_coverage, 'args': (sat_opt, reg)},
                  {'type': 'ineq', 'fun': constraint_com_cc_ctc_coverage, 'args': (sat_opt, reg)},
                  {'type': 'ineq', 'fun': constraint_com_cc_ctc_cpad_coverage, 'args': (sat_opt, reg)})
            bnds=((0,None), (0,None), (0,None), (0,None), (0,None))
            expenditure=array([baseline_expenditure[reg,0]/5, baseline_expenditure[reg,0]/5,baseline_expenditure[reg,0]/5,baseline_expenditure[reg,0]/5,baseline_expenditure[reg,0]/5])
            sol=minimize(objective_coverage, expenditure, method='SLSQP', constraints=cons, bounds=bnds, options={'disp':False})
            expenditure=array([sol.x[0], sol.x[1], sol.x[2], sol.x[3], sol.x[4]])
            baseline_optimised_expenditure.append([sol.x[0], sol.x[1], sol.x[2], sol.x[3], sol.x[4],sol.x[0]+ sol.x[1]+sol.x[2]+sol.x[3]+sol.x[4]])
            baseline_optimised_coverage.append(coverage_from_opt(expenditure))
            
        baseline_optimised_expenditure=array(baseline_optimised_expenditure)
        baseline_optimised_coverage=array(baseline_optimised_coverage)
            
        # Run the model for optimised baseline expenditure
            
        fac_cc=baseline_optimised_coverage[:,1].repeat(runs).reshape(len(hbv_pars),runs,1)
        fac_ctc=baseline_optimised_coverage[:,3].repeat(runs).reshape(len(hbv_pars),runs,1)
        com_cc=baseline_optimised_coverage[:,5].repeat(runs).reshape(len(hbv_pars),runs,1)
        com_ctc=baseline_optimised_coverage[:,7].repeat(runs).reshape(len(hbv_pars),runs,1)
        com_cpad=baseline_optimised_coverage[:,9].repeat(runs).reshape(len(hbv_pars),runs,1)
            
        model_init=model_initializer(runs, globz, hbv_pars, globz_perturbed, regions_perturbed, fac_cc, fac_ctc, com_cc, com_ctc, com_cpad, t_steps, ve_warm)
        
        for reg in range(len(hbv_pars)):
            for run in range (runs):
                model_init[(runs*reg)+run,0,-2]=baseline_optimised_expenditure[reg,5]
    
        baseline_optimised=hbv_model(runs,t_steps,dt,model_init, hbv_pars, regions_perturbed, globz, globz_perturbed,discount)
            
            
        #Optimisation to reach 90% coverage
        optimised_90_coverage, optimised_90_expenditure=[], []
        for reg in range(len(hbv_pars)):
            cons=({'type':'eq', 'fun': constraint_90_coverage},
                  {'type': 'ineq', 'fun': constraint_fac_cc_coverage, 'args': (sat_opt, reg)},
                  {'type': 'ineq', 'fun': constraint_fac_cc_ctc_coverage, 'args': (sat_opt, reg)},
                  {'type': 'ineq', 'fun': constraint_com_cc_ctc_coverage, 'args': (sat_opt, reg)},
                  {'type': 'ineq', 'fun': constraint_com_cc_ctc_cpad_coverage, 'args': (sat_opt, reg)})
            bnds=((0,None), (0,None), (0,None), (0,None), (0,None))
            expenditure=array([1,1,1,1,1])
            sol=minimize(objective_90_coverage, expenditure, method='SLSQP', constraints=cons, bounds=bnds, options={'disp':False})
            expenditure=array([sol.x[0], sol.x[1], sol.x[2], sol.x[3], sol.x[4]])
            optimised_90_expenditure.append([sol.x[0], sol.x[1], sol.x[2], sol.x[3], sol.x[4],sol.x[0]+ sol.x[1]+sol.x[2]+sol.x[3]+sol.x[4]])
            optimised_90_coverage.append(coverage_from_opt(expenditure))
            
        optimised_90_expenditure=array(optimised_90_expenditure)
        optimised_90_coverage=array(optimised_90_coverage)
            
        #Run the model for optimised 90% coverage
            
        fac_cc=optimised_90_coverage[:,1].repeat(runs).reshape(len(hbv_pars),runs,1)
        fac_ctc=optimised_90_coverage[:,3].repeat(runs).reshape(len(hbv_pars),runs,1)
        com_cc=optimised_90_coverage[:,5].repeat(runs).reshape(len(hbv_pars),runs,1)
        com_ctc=optimised_90_coverage[:,7].repeat(runs).reshape(len(hbv_pars),runs,1)
        com_cpad=optimised_90_coverage[:,9].repeat(runs).reshape(len(hbv_pars),runs,1)
            
        model_init=model_initializer(runs, globz, hbv_pars, globz_perturbed, regions_perturbed, fac_cc, fac_ctc, com_cc, com_ctc, com_cpad, t_steps, ve_warm)
        
        for reg in range(len(hbv_pars)):
            for run in range (runs):
                model_init[(runs*reg)+run,0,-2]=optimised_90_expenditure[reg,5]
    
        optimised_90=hbv_model(runs,t_steps,dt,model_init, hbv_pars, regions_perturbed, globz, globz_perturbed,discount)
            
        
  
