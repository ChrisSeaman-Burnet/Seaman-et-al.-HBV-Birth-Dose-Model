"""Uncertainity Analysis - only needs to run if runs > 1 in mainscript"""

def perturb_pars(runs, globz, hbv_pars, hbv_pars_lb, hbv_pars_ub):
    
    """Randomly draws parameters from uncertainity bounds, with values stored in array of
    length runs. All parameters drawn from a unform distribution, except prevalence drawn from
    truncated normal distribution"""
    
    import pandas as pd
    from numpy import random,zeros,array,arange, r_
    from scipy import stats
    
    random.RandomState(1234)
    
    #Export constant values (Population size, Tin, weight and discount)
    P=globz.iloc[0,1]
    weight=globz.iloc[2,1]
    discount=globz.iloc[3,1]
    
    #Monte Carlo Simulations for Global Parameters
    glob_uc_par_names=list(globz.loc[4:,"Symbol"])
   
    glob_uc=zeros((len(glob_uc_par_names),runs))
    for i in range(4,29):
        glob_uc[(i-4),:]=random.uniform(globz.iloc[i,2], globz.iloc[i,3], runs)       
   
    globz_perturbed=pd.DataFrame(glob_uc, index=[i for i in glob_uc_par_names])
    
   #Monte Carlo Simulations for Regional/Country Parameters
    
    regions_perturbed=zeros((len(hbv_pars),41,runs)) 
    ####Perturb of prevalence (Truncated Normal Distribution)
    random.seed(1234)
    for i in range(len(hbv_pars)):
        regions_perturbed[i,0,:]=stats.truncnorm.rvs((hbv_pars.iloc[i,9]-hbv_pars.iloc[i,8])/hbv_pars.iloc[i,11], (hbv_pars.iloc[i,10]-hbv_pars.iloc[i,8])/hbv_pars.iloc[i,11], loc=hbv_pars.iloc[i,8], scale=hbv_pars.iloc[i,11], size=runs)
    
    ###Perturb of othe regional variables (Uniform Distribution)
    reg_uniform_lb=array(hbv_pars_lb.iloc[:,r_[6,7,13:51]])
    reg_uniform_ub=array(hbv_pars_ub.iloc[:,r_[6,7,13:51]])
    
    for reg in range(len(reg_uniform_ub)):
        for val in range(40):
            regions_perturbed[reg, val+1, :]=random.uniform(reg_uniform_lb[reg,val], reg_uniform_ub[reg,val], runs)
    
    return P,weight,discount,globz_perturbed,reg_uc_par_names,regions_perturbed
    