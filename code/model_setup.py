###Model Set-Up Stuff###

def time_steps(Tin, ints):
    """Sets the model run time, and length of interval. 
    
    Tin = Number of years model to be run over
    ints = Number of time intervals per year
    """
    from numpy import linspace
    
    t_steps=linspace(0, Tin, int(ints*Tin))
    dt=(len(t_steps)/Tin)**-1
    
    return t_steps, dt



def facility_weight(runs, weight, hbv_pars, regions_perturbed):
    """Weights vaccination so coverage is proportional to facility and community births in each setting, but
    relative likelihood of vaccination in facilities is higher by a factor determined by the weight"""
    
    from numpy import zeros
    
    if runs==1:
        fac_cov_weighted=zeros(len(hbv_pars))
        com_cov_weighted=zeros(len(hbv_pars))
        
        for idx,(bd,facility) in enumerate(zip(hbv_pars["birth_dose"], hbv_pars["Facility"])):
            if bd <= 0.5:
                fac_cov_weighted[idx]=bd+((bd*(weight-1)/(weight+1))*(2*(1-facility)))
                com_cov_weighted[idx]=bd-((bd*(weight-1)/(weight+1))*(2*facility))
            elif bd >0.5:
                fac_cov_weighted[idx]=bd+(((1-bd)*(weight-1)/(weight+1))*(2*(1-facility)))
                com_cov_weighted[idx]=bd-(((1-bd)*(weight-1)/(weight+1))*(2*facility))
    
    elif runs >1:
        fac_cov_weighted=zeros((len(hbv_pars),runs,1))
        com_cov_weighted=zeros((len(hbv_pars),runs,1))
        
        birth_dose=regions_perturbed[:,1,:].reshape((len(hbv_pars),runs,1))
        facility=regions_perturbed[:,4,:].reshape((len(hbv_pars),runs,1))
        
        for idx,i in enumerate(birth_dose[:,:,:]):
            for idy,bd in enumerate(i):
                if bd <= 0.5:
                    fac_cov_weighted[idx,idy,0]=bd+((bd*(weight-1)/(weight+1))*(2*(1-facility[idx,idy,0])))
                    com_cov_weighted[idx,idy,0]=bd-((bd*(weight-1)/(weight+1))*(2*facility[idx,idy,0]))
                elif bd > 0.5:
                    fac_cov_weighted[idx,idy,0]=bd+(((1-bd)*(weight-1)/(weight+1))*(2*(1-facility[idx,idy,0])))
                    com_cov_weighted[idx,idy,0]=bd-(((1-bd)*(weight-1)/(weight+1))*(2*facility[idx,idy,0]))

    return fac_cov_weighted, com_cov_weighted 


