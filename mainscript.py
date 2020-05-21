"""
Model Code

Use of the controlled-temperature chain and compact prefilled auto-disable devices 
to reach 2030 hepatitis B birth dose vaccination targets in low and middle-income
countries: A modelling and cost-optimization study.

Last Update: 29/03/2020

"""
#%% Run time of model 
from datetime import datetime
start_time=datetime.now()

#%% Set Working Directory
working_dir= """Set your own working directory"""
from os import chdir as set_wd
set_wd(working_dir)

#%% Analysis Level
#Define the analyses you want run
region=1     #1 = Global Burden of Disease Regions, 0 = Countries
runs=1     # >1 = Uncertainity Analysis of N runs (1000 runs takes approximately 20 minutes computing time for GBD regions)
ve_warm=1       #Assumption of equality to be tested in Sensitivity Analysis
#%% 
###Importing the required functions###
from importdata import data_import
from model_setup import time_steps, facility_weight
from uncertainity import perturb_pars 
from model import  coverage_saturation, model_scenarios, model_initializer, vaccine_costs, hbv_model


#%%
if runs==1:
    globz, hbv_pars=data_import(working_dir, region,runs)
elif runs > 1:
    globz, hbv_pars, hbv_pars_lb, hbv_pars_ub=data_import(working_dir, region,runs)

Tin=globz.iloc[1,1]
t_steps,dt=time_steps(Tin,5)

if len(t_steps) < 400:
    globz.loc[14, 1:]=1/dt                
    globz.loc[15, 1:]=1/dt                 
    globz.loc[17, 1:]=1/dt                 
elif len(t_steps) >=400 or len(t_steps)<2400:
    globz.loc[14, 1:]=1/dt                  

#Monte Carlo Simulations for Variables 
if runs > 1:
    P,weight,discount,globz_perturbed,reg_uc_par_names,regions_perturbed=perturb_pars(runs, globz, hbv_pars, hbv_pars_lb, hbv_pars_ub)
elif runs==1:
    P=globz.iloc[0,1]
    weight=globz.iloc[2,1]
    discount=globz.iloc[3,1]

#Weighting of likelihood to recieve hepatitis B birth dose by birth location
if runs==1:
    regions_perturbed=0
    globz_perturbed=0
    fac_cov_weighted, com_cov_weighted = facility_weight(runs, weight, hbv_pars, regions_perturbed)
elif runs > 1:
    fac_cov_weighted, com_cov_weighted = facility_weight(runs, weight, hbv_pars, regions_perturbed)

#Setting saturation bounds for modelling and optimization
settings=list(hbv_pars["Setting"])
saturation, sat_costing=coverage_saturation(fac_cov_weighted, com_cov_weighted, hbv_pars, regions_perturbed,runs)

for scen in range(5):
    
    fac_cc, fac_ctc, com_cc, com_ctc, com_cpad= model_scenarios(runs, scen, hbv_pars, fac_cov_weighted, com_cov_weighted, saturation)
   
    if scen==0:
        
        #Initializes the model
        model_init=model_initializer(runs, globz, hbv_pars, globz_perturbed, regions_perturbed, fac_cc, fac_ctc, com_cc, com_ctc, com_cpad, t_steps,ve_warm)
        
        #Estimates costs using cost-curve functions
        baseline_vcost, baseline_coverage=vaccine_costs(runs, hbv_pars, fac_cc, fac_ctc, com_cc, com_ctc, com_cpad,sat_costing, saturation, regions_perturbed)
        
        #Puts the vaccine cost into disease simulation array to estimate total cost
        if runs ==1:
            for reg in range(len(hbv_pars)):
                model_init[reg,0,-2]=baseline_vcost[reg,5]
        elif runs> 1:
            for reg in range(len(hbv_pars)):
                for run in range (runs):
                    model_init[(runs*reg)+run,0,-2]=baseline_vcost[reg,run,5]
        
        #Runs the model
        baseline=hbv_model(runs,t_steps,dt,model_init, hbv_pars, regions_perturbed, globz, globz_perturbed,discount)
        

    elif scen==1 and region==1:
        
        model_init=model_initializer(runs, globz, hbv_pars, globz_perturbed, regions_perturbed, fac_cc, fac_ctc, com_cc, com_ctc, com_cpad, t_steps,ve_warm)
        max_current_vcost, max_current_coverage=vaccine_costs(runs, hbv_pars, fac_cc, fac_ctc, com_cc, com_ctc, com_cpad,sat_costing, saturation, regions_perturbed)
        
        if runs ==1:
            for reg in range(len(hbv_pars)):
                model_init[reg,0,-2]=max_current_vcost[reg,5]
        elif runs> 1:
            for reg in range(len(hbv_pars)):
                for run in range (runs):
                    model_init[(runs*reg)+run,0,-2]=max_current_vcost[reg,run,5]
        
        max_current=hbv_model(runs,t_steps,dt,model_init, hbv_pars, regions_perturbed, globz, globz_perturbed,discount)
       
        
    elif scen==2 and region==1:
        
        model_init=model_initializer(runs, globz, hbv_pars, globz_perturbed, regions_perturbed, fac_cc, fac_ctc, com_cc, com_ctc, com_cpad, t_steps, ve_warm)
        ctc_facilities_vcost, ctc_facilities_coverage=vaccine_costs(runs, hbv_pars, fac_cc, fac_ctc, com_cc, com_ctc, com_cpad,sat_costing, saturation, regions_perturbed)
        
        if runs ==1:
            for reg in range(len(hbv_pars)):
                model_init[reg,0,-2]=ctc_facilities_vcost[reg,5]
        elif runs> 1:
            for reg in range(len(hbv_pars)):
                for run in range (runs):
                    model_init[(runs*reg)+run,0,-2]=ctc_facilities_vcost[reg,run,5]
        
        ctc_facilities=hbv_model(runs,t_steps,dt,model_init, hbv_pars, regions_perturbed, globz, globz_perturbed,discount)

        
    elif scen==3 and region==1:
        
        model_init=model_initializer(runs, globz, hbv_pars, globz_perturbed, regions_perturbed, fac_cc, fac_ctc, com_cc, com_ctc, com_cpad, t_steps, ve_warm)
        ctc_community_vcost, ctc_community_coverage=vaccine_costs(runs, hbv_pars, fac_cc, fac_ctc, com_cc, com_ctc, com_cpad,sat_costing, saturation, regions_perturbed)
        
        if runs ==1:
            for reg in range(len(hbv_pars)):
                model_init[reg,0,-2]=ctc_community_vcost[reg,5]
        elif runs> 1:
            for reg in range(len(hbv_pars)):
                for run in range (runs):
                    model_init[(runs*reg)+run,0,-2]=ctc_community_vcost[reg,run,5]
        
        ctc_community=hbv_model(runs,t_steps,dt,model_init, hbv_pars, regions_perturbed, globz, globz_perturbed,discount)
        
    elif scen==4 and region==1:
        
        model_init=model_initializer(runs, globz, hbv_pars, globz_perturbed, regions_perturbed, fac_cc, fac_ctc, com_cc, com_ctc, com_cpad, t_steps, ve_warm)
        cpad_community_vcost, cpad_community_coverage=vaccine_costs(runs, hbv_pars, fac_cc, fac_ctc, com_cc, com_ctc, com_cpad,sat_costing, saturation, regions_perturbed)
        
        if runs ==1:
            for reg in range(len(hbv_pars)):
                model_init[reg,0,-2]=cpad_community_vcost[reg,5]
        elif runs> 1:
            for reg in range(len(hbv_pars)):
                for run in range (runs):
                    model_init[(runs*reg)+run,0,-2]=cpad_community_vcost[reg,run,5]
        
        cpad_community=hbv_model(runs,t_steps,dt,model_init, hbv_pars, regions_perturbed, globz, globz_perturbed,discount)

#Optimizations
exec(open("costcoverage.py").read())
exec(open("optimisation_model.py").read())

#Results
exec(open("tablesandfigs.py").read())

#Sensitivity Analysis
exec(open("onewaysens.py").read())


#%%Get runtime of model
print("Time to completion =" + str(datetime.now()-start_time))    

