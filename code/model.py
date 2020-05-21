###model

def coverage_saturation(fac_cov_weighted, com_cov_weighted, hbv_pars, regions_perturbed,runs):
    
    """Calculates maximum coverage capacity (saturation) for each intervention, in each setting. Values 
    are dependent upon epidemiological and demographic parameters. Costing values are also calculated, to avoid zero
    errors when using logistic functions."""
    
    from numpy import zeros
    
    if runs==1:
        saturation=zeros((len(hbv_pars), 5))
        sat_costing=zeros((len(hbv_pars),5))
        
        for reg in range(len(hbv_pars)):
            # Facility Vaccination
            saturation[reg,0]=max(0.80, fac_cov_weighted[reg]+(1-fac_cov_weighted[reg])*0.0001)
            saturation[reg,1]=max(0.90, fac_cov_weighted[reg]+(1-fac_cov_weighted[reg])*0.05)
            
            #Community Vaccination
            saturation[reg,2]=com_cov_weighted[reg]
            
            if hbv_pars.iloc[reg,47] > com_cov_weighted[reg]*0.3:
                saturation[reg,3]=min(0.95, hbv_pars.iloc[reg,47]+(0.7*com_cov_weighted[reg])+(0.1*(1-com_cov_weighted[reg])))
            else:
                saturation[reg,3]=min(0.95, com_cov_weighted[reg]+(1-com_cov_weighted[reg])*0.1)
            
            saturation[reg,4]=0.95
           
            if saturation[reg,2]>saturation[reg,3]:  #This can be modified so that additional coverage occurs in high coverage settings
                saturation[reg,3]=saturation[reg,2]+(1-saturation[reg,2])*0.1
                saturation[reg,4]=saturation[reg,3]
            
            
            sat_costing[reg,0]= saturation[reg,0]
            sat_costing[reg,1]= saturation[reg,1]
            sat_costing[reg,2]= saturation[reg,2]+(1- saturation[reg,2])*0.05
            sat_costing[reg,3]= saturation[reg,3]+(1- saturation[reg,3])*0.05
            sat_costing[reg,4]= saturation[reg,4]+(1- saturation[reg,4])*0.05
        
    elif runs > 1:
        saturation=zeros((len(hbv_pars),runs,5))
        sat_costing=zeros((len(hbv_pars),runs,5))
        
        for reg in range(len(hbv_pars)):
            for run in range(runs):
                #Facility Vaccination
                saturation[reg,run,0]=max(0.80, fac_cov_weighted[reg,run,0]+(1-fac_cov_weighted[reg,run,0])*0.0001)
                saturation[reg,run,1]=max(0.90, fac_cov_weighted[reg,run,0]+(1-fac_cov_weighted[reg,run,0])*0.05)
                
                #Community Vaccination
                saturation[reg,run,2]=com_cov_weighted[reg,run,0]
                
                if regions_perturbed[reg,37,run] > com_cov_weighted[reg,run,0]*0.3:
                    saturation[reg,run,3]=min(0.95, regions_perturbed[reg,37,run]+(0.7*com_cov_weighted[reg,run,0])+(0.1*(1-com_cov_weighted[reg,run,0])))
                else:
                    saturation[reg,run,3]=min(0.95, com_cov_weighted[reg,run,0]+(1-com_cov_weighted[reg,run,0])*0.1)
                
                saturation[reg,run,4]=0.95
                
                if saturation[reg,run,2]>saturation[reg,run,3]:         #This can be modified so that additional coverage occurs in high coverage settings
                    saturation[reg,run,3]=saturation[reg,run,2]+(1-saturation[reg,run,2])*0.1
                    saturation[reg,run,4]=saturation[reg,run,3]
                
                sat_costing[reg,run,0]=saturation[reg,run,0]
                sat_costing[reg,run,1]=saturation[reg,run,1]
                sat_costing[reg,run,2]=saturation[reg,run,2]+(1-saturation[reg,run,2])*0.05
                sat_costing[reg,run,3]=saturation[reg,run,3]+(1-saturation[reg,run,3])*0.05
                sat_costing[reg,run,4]=saturation[reg,run,4]+(1-saturation[reg,run,4])*0.05
                
    return saturation, sat_costing


def model_scenarios(runs, scen, hbv_pars, fac_cov_weighted, com_cov_weighted, saturation):
    
    """Returns the coverage specific to each modelled intervention, as an array
    for each modelled scenario"""
    
    from numpy import zeros
    
    if runs ==1:
        
        if scen==0:      #baseline
            
            fac_cc=fac_cov_weighted[:]
            fac_ctc=zeros(len(hbv_pars))
            com_cc=com_cov_weighted[:]
            com_ctc=zeros(len(hbv_pars))
            com_cpad=zeros(len(hbv_pars))
            
            return fac_cc, fac_ctc, com_cc, com_ctc, com_cpad
        
        elif scen==1:       #Maximized Current Practice
            
            fac_cc=saturation[:,0]
            fac_ctc=zeros(len(hbv_pars))
            com_cc=saturation[:,2]
            com_ctc=zeros(len(hbv_pars))
            com_cpad=zeros(len(hbv_pars))
            
            return fac_cc, fac_ctc, com_cc, com_ctc, com_cpad
        
        elif scen==2:       #S1 + CTC in Facilities
           
            fac_cc=saturation[:,0]
            fac_ctc=saturation[:,1]-saturation[:,0]
            com_cc=saturation[:,2]
            com_ctc=zeros(len(hbv_pars))
            com_cpad=zeros(len(hbv_pars))

            return fac_cc, fac_ctc, com_cc, com_ctc, com_cpad
        
        elif scen==3:       #S2+ CTC in the Community
        
            fac_cc=saturation[:,0]
            fac_ctc=saturation[:,1]-saturation[:,0]
            com_cc=saturation[:,2]
            com_ctc=saturation[:,3]-saturation[:,2]
            com_cpad=zeros(len(hbv_pars))
            
            return fac_cc, fac_ctc, com_cc, com_ctc, com_cpad
            
        
        elif scen==4:       #S3 + CPAD in the Community
            
            fac_cc=saturation[:,0]
            fac_ctc=saturation[:,1]-saturation[:,0]
            com_cc=saturation[:,2]
            com_ctc=saturation[:,3]-saturation[:,2]
            com_cpad=saturation[:,4]-saturation[:,3]
            
            return fac_cc, fac_ctc, com_cc, com_ctc, com_cpad
    
    elif runs > 1:
    
        if scen==0:
            fac_cc=fac_cov_weighted[:,:,:]
            fac_ctc=zeros((len(hbv_pars),runs,1))
            com_cc=com_cov_weighted[:,:,:]
            com_ctc=zeros((len(hbv_pars),runs,1))
            com_cpad=zeros((len(hbv_pars),runs,1))
            
            return fac_cc, fac_ctc, com_cc, com_ctc, com_cpad
            
        elif scen==1:
            fac_cc=saturation[:,:,0].reshape(len(hbv_pars),runs,1)
            fac_ctc=zeros((len(hbv_pars),runs,1))
            com_cc=saturation[:,:,2].reshape(len(hbv_pars),runs,1)
            com_ctc=zeros((len(hbv_pars),runs,1))
            com_cpad=zeros((len(hbv_pars),runs,1))
            
            return fac_cc, fac_ctc, com_cc, com_ctc, com_cpad
            
        elif scen ==2:
            fac_cc=saturation[:,:,0].reshape(len(hbv_pars),runs,1)
            fac_ctc=(saturation[:,:,1]-saturation[:,:,0]).reshape(len(hbv_pars),runs,1)
            com_cc=saturation[:,:,2].reshape(len(hbv_pars),runs,1)
            com_ctc=zeros((len(hbv_pars),runs,1))
            com_cpad=zeros((len(hbv_pars),runs,1))
            
            return fac_cc, fac_ctc, com_cc, com_ctc, com_cpad
            
        elif scen==3:
            fac_cc=saturation[:,:,0].reshape(len(hbv_pars),runs,1)
            fac_ctc=(saturation[:,:,1]-saturation[:,:,0]).reshape(len(hbv_pars),runs,1)
            com_cc=saturation[:,:,2].reshape(len(hbv_pars),runs,1)
            com_ctc=(saturation[:,:,3]-saturation[:,:,2]).reshape(len(hbv_pars),runs,1)
            com_cpad=zeros((len(hbv_pars),runs,1))
            
            return fac_cc, fac_ctc, com_cc, com_ctc, com_cpad
            
        elif scen==4:
            fac_cc=saturation[:,:,0].reshape(len(hbv_pars),runs,1)
            fac_ctc=(saturation[:,:,1]-saturation[:,:,0]).reshape(len(hbv_pars),runs,1)
            com_cc=saturation[:,:,2].reshape(len(hbv_pars),runs,1)
            com_ctc=(saturation[:,:,3]-saturation[:,:,2]).reshape(len(hbv_pars),runs,1)
            com_cpad=(saturation[:,:,4]-saturation[:,:,3]).reshape(len(hbv_pars),runs,1)
            
            return fac_cc, fac_ctc, com_cc, com_ctc, com_cpad


def model_initializer(runs, globz, hbv_pars, globz_perturbed, regions_perturbed, fac_cc, fac_ctc, com_cc, com_ctc, com_cpad, t_steps, ve_warm):
    
    """Returns initial compartments of Immune ,Susceptible, and Acute after distrubuting vaccination into 
    time strata, and vaccine effectiveness calculations."""
    from numpy import array, zeros
    
    fac_cc_timing=array([0.8, 0.2/3, 0.2/3, 0.2/3,0])
    fac_ctc_timing=array([0.8+(0.05*(0.2/3)), 0.2/3, 0.2/3, 0.95*(0.2/3),0])
    
    com_cc_timing=array([0.3, 0.4, 0.15, 0.15,0])
    com_ctc_timing=array([1, 0.375/0.66, 0.15/0.66, 0.135/0.66, 0])
    com_cpad_timing=array([1,0,0,0,0])
  
    if runs ==1:
        
        ve=array(globz.iloc[9:14,1])
        com_sba=array(hbv_pars.iloc[:,-4])
       
        
        fac_cc_ve=zeros((len(hbv_pars),5))
        fac_ctc_ve=zeros((len(hbv_pars),5))
        com_cc_ve=zeros((len(hbv_pars),5))
        com_ctc_ve=zeros((len(hbv_pars),5))
        com_cpad_ve=zeros((len(hbv_pars),5))
        
        fac_cc_ve_fail=zeros((len(hbv_pars),5))
        fac_ctc_ve_fail=zeros((len(hbv_pars),5))
        com_cc_ve_fail=zeros((len(hbv_pars),5))
        com_ctc_ve_fail=zeros((len(hbv_pars),5))
        com_cpad_ve_fail=zeros((len(hbv_pars),5))
        
        
        for reg in range(len(hbv_pars)):
            for i in range(4):
                fac_cc_ve[reg, i]=fac_cc[reg]*fac_cc_timing[i]*ve[i]
                fac_cc_ve_fail[reg, i]=fac_cc[reg]*fac_cc_timing[i]*(1-ve[i])       
            
                com_cc_ve[reg,i]=com_cc[reg]*com_cc_timing[i]*ve[i]
                com_cc_ve_fail[reg,i]=com_cc[reg]*com_cc_timing[i]*(1-ve[i])
                
        for reg in range(len(hbv_pars)):
            fac_cc_ve_fail[reg,4]=1-(fac_cc[reg]+fac_ctc[reg])
            com_cc_ve_fail[reg,4]=1-(com_cc[reg]+com_ctc[reg]+com_cpad[reg])
            
        for reg in range(len(hbv_pars)):
            for i in range(5):
                fac_ctc_ve[reg,i]=fac_ctc[reg]*fac_ctc_timing[i]*(ve[i]*ve_warm)
                fac_ctc_ve_fail[reg,i]=fac_ctc[reg]*fac_ctc_timing[i]*(1-(ve[i]*ve_warm))
            
                com_cpad_ve[reg,i]=com_cpad[reg]*com_cpad_timing[i]*(ve[i]*ve_warm)
                com_cpad_ve_fail[reg,i]=com_cpad[reg]*com_cpad_timing[i]*(1-(ve[i]*ve_warm))
    
        for reg in range(len(hbv_pars)):
            for i in range(1,5,1):
            
                if com_ctc[reg]==0:
                    com_ctc_ve[reg,0]=0
                    com_ctc_ve[reg,i]=0
            
                elif com_sba[reg]> 0.3*com_cc[reg]:            #if existing timely coverage is greater than professsional birth attendance
                    com_ctc_ve[reg,0]=(com_sba[reg]-0.3*com_cc[reg])*(ve[0]*ve_warm)
                    com_ctc_ve_fail[reg,0]=(com_sba[reg]-0.3*com_cc[reg])*(1-(ve[0]*ve_warm))
                
                    com_ctc_ve[reg,i]=(com_ctc[reg]-(com_sba[reg]-0.3*com_cc[reg]))*(ve[i]*ve_warm)*com_ctc_timing[i]
                    com_ctc_ve_fail[reg,i]=(com_ctc[reg]-(com_sba[reg]-0.3*com_cc[reg]))*(1-(ve[i]*ve_warm))*com_ctc_timing[i]
            
                else:
                    com_ctc_ve[reg,0]=com_ctc[reg]*0.34*(ve[0]*ve_warm)
                    com_ctc_ve_fail[reg,0]=com_ctc[reg]*0.34*(1-(ve[0]*ve_warm))
                
                    com_ctc_ve[reg,i]=0.66*com_ctc[reg]*com_ctc_timing[i]*(ve[i]*ve_warm)
                    com_ctc_ve_fail[reg,i]=0.66*com_ctc[reg]*com_ctc_timing[i]*(1-(ve[i]*ve_warm))
                
        fac_cov=zeros((len(hbv_pars)))
        fac_ncov=zeros((len(hbv_pars)))
        com_cov=zeros((len(hbv_pars)))
        com_ncov=zeros((len(hbv_pars)))
     
        for reg in range(len(hbv_pars)):
            fac_cov[reg]=sum(fac_cc_ve[reg,:])+sum(fac_ctc_ve[reg,:])
            fac_ncov[reg]=sum(fac_cc_ve_fail[reg,:])+sum(fac_ctc_ve_fail[reg,:])
            com_cov[reg]=sum(com_cc_ve[reg,:])+sum(com_ctc_ve[reg,:])+sum(com_cpad_ve[reg,:])
            com_ncov[reg]=sum(com_cc_ve_fail[reg,:])+sum(com_ctc_ve_fail[reg,:])+sum(com_cpad_ve_fail[reg,:])
    
        model_init=zeros((len(hbv_pars), len(t_steps),14))

        #Data needed to initialize the model
        P=1000
        prev=array(hbv_pars.iloc[:,8])
        hbv3=array(hbv_pars.iloc[:,7])
        hbe_prev=array(hbv_pars.iloc[:,48])
        hbe_pos=array(hbv_pars.iloc[:,49])
        hbe_neg=array(hbv_pars.iloc[:,50])
        hosp=array(hbv_pars.iloc[:,14])
        comm=1-hosp
    
        for reg in range(len(hbv_pars)):
            model_init[reg,0,0]=(P*hosp[reg]*hbv3[reg]*((prev[reg]*(hbe_prev[reg])*(1-hbe_pos[reg]))+(prev[reg]*hbe_prev[reg]*hbe_pos[reg]*(fac_cov[reg]))+(prev[reg]*(1-hbe_prev[reg])*(1-hbe_neg[reg]))+(prev[reg]*(1-hbe_prev[reg])*hbe_neg[reg]*fac_cov[reg])+(1-prev[reg])))+(P*comm[reg]*hbv3[reg]*((prev[reg]*hbe_prev[reg]*(1-hbe_pos[reg]))+(prev[reg]*hbe_prev[reg]*hbe_pos[reg]*com_cov[reg])+(prev[reg]*(1-hbe_prev[reg])*(1-hbe_neg[reg]))+(prev[reg]*(1-hbe_prev[reg])*hbe_neg[reg]*com_cov[reg])+(1-prev[reg]))) #Z
            model_init[reg,0,1]=(P*hosp[reg]*(1-hbv3[reg])*((prev[reg]*(hbe_prev[reg])*(1-hbe_pos[reg]))+(prev[reg]*hbe_prev[reg]*hbe_pos[reg]*(fac_cov[reg]))+(prev[reg]*(1-hbe_prev[reg])*(1-hbe_neg[reg]))+(prev[reg]*(1-hbe_prev[reg])*hbe_neg[reg]*fac_cov[reg])+(1-prev[reg])))+(P*comm[reg]*(1-hbv3[reg])*((prev[reg]*hbe_prev[reg]*(1-hbe_pos[reg]))+(prev[reg]*hbe_prev[reg]*hbe_pos[reg]*com_cov[reg])+(prev[reg]*(1-hbe_prev[reg])*(1-hbe_neg[reg]))+(prev[reg]*(1-hbe_prev[reg])*hbe_neg[reg]*com_cov[reg])+(1-prev[reg])))
            model_init[reg,0,2]=(P*hosp[reg]*((prev[reg]*hbe_prev[reg]*hbe_pos[reg]*fac_ncov[reg])+(prev[reg]*(1-hbe_prev[reg])*hbe_neg[reg]*fac_ncov[reg])))+(P*comm[reg]*((prev[reg]*hbe_prev[reg]*hbe_pos[reg]*com_ncov[reg])+(prev[reg]*(1-hbe_prev[reg])*hbe_neg[reg]*com_ncov[reg])))
    
        return model_init

    elif runs > 1:
        
        ve=array(globz_perturbed.iloc[5:10,:])
        com_sba=regions_perturbed[:,37,:].reshape(len(hbv_pars),runs,1)
       
        
        fac_cc_ve=zeros((len(hbv_pars),runs,5))
        fac_ctc_ve=zeros((len(hbv_pars),runs,5))
        com_cc_ve=zeros((len(hbv_pars),runs, 5))
        com_ctc_ve=zeros((len(hbv_pars),runs, 5))
        com_cpad_ve=zeros((len(hbv_pars),runs, 5))
        
        fac_cc_ve_fail=zeros((len(hbv_pars),runs,5))
        fac_ctc_ve_fail=zeros((len(hbv_pars), runs, 5))
        com_cc_ve_fail=zeros((len(hbv_pars), runs,5))
        com_ctc_ve_fail=zeros((len(hbv_pars), runs, 5))
        com_cpad_ve_fail=zeros((len(hbv_pars), runs,5))
        
        for reg in range(len(hbv_pars)):
            for run in range(runs):
                for i in range(4):
                    fac_cc_ve[reg, run, i]=fac_cc[reg,run,0]*fac_cc_timing[i]*ve[i,run]
                    fac_cc_ve_fail[reg,run,i]=fac_cc[reg,run,0]*fac_cc_timing[i]*(1-ve[i,run])       
            
                    com_cc_ve[reg,run,i]=com_cc[reg,run,0]*com_cc_timing[i]*ve[i,run]
                    com_cc_ve_fail[reg,run,i]=com_cc[reg,run,0]*com_cc_timing[i]*(1-ve[i,run])
        
        for reg in range(len(hbv_pars)):
            for run in range(runs):
                fac_cc_ve_fail[reg,run,4]=1-(fac_cc[reg,run,0]+fac_ctc[reg,run,0])
                com_cc_ve_fail[reg,run,4]=1-(com_cc[reg,run,0]+com_ctc[reg,run,0]+com_cpad[reg,run,0])
        
        
        for reg in range(len(hbv_pars)):
            for run in range(runs):
                for i in range(5):
                    fac_ctc_ve[reg,run,i]=fac_ctc[reg,run,0]*fac_ctc_timing[i]*(ve[i,run]*ve_warm)
                    fac_ctc_ve_fail[reg,run,i]=fac_ctc[reg,run,0]*fac_ctc_timing[i]*(1-(ve[i,run]*ve_warm))
            
                    com_cpad_ve[reg,run,i]=com_cpad[reg,run,0]*com_cpad_timing[i]*(ve[i,run]*ve_warm)
                    com_cpad_ve_fail[reg,run,i]=com_cpad[reg,run,0]*com_cpad_timing[i]*(1-(ve[i,run]*ve_warm))
        
        for reg in range(len(hbv_pars)):
            for run in range(runs):
                for i in range(1,5,1):
            
                    if com_ctc[reg,run,0]==0:
                        com_ctc_ve[reg,run,0]=0
                        com_ctc_ve[reg,run,i]=0
            
                    elif com_sba[reg,run,0]> 0.3*com_cc[reg,run,0]:            #if existing timely coverage is greater than professsional birth attendance
                        com_ctc_ve[reg,run,0]=(com_sba[reg,run,0]-0.3*com_cc[reg,run,0])*(ve[0,run]*ve_warm)
                        com_ctc_ve_fail[reg,run,0]=(com_sba[reg,run,0]-0.3*com_cc[reg,run,0])*(1-(ve[0,run]*ve_warm))
                
                        com_ctc_ve[reg,run,i]=max((com_ctc[reg,run,0]-(com_sba[reg,run,0]-0.3*com_cc[reg,run,0]))*(ve[i,run]*ve_warm)*com_ctc_timing[i],0)
                        com_ctc_ve_fail[reg,run,i]=max((com_ctc[reg,run,0]-(com_sba[reg,run,0]-0.3*com_cc[reg,run,0]))*(1-(ve[i,run]*ve_warm))*com_ctc_timing[i],0)
            
                    else:
                        com_ctc_ve[reg,run,0]=com_ctc[reg,run,0]*0.34*(ve[0,run]*ve_warm)
                        com_ctc_ve_fail[reg,run,0]=com_ctc[reg,run,0]*0.34*(1-(ve[0,run]*ve_warm))
                
                        com_ctc_ve[reg,run,i]=0.66*com_ctc[reg,run,0]*com_ctc_timing[i]*(ve[i,run]*ve_warm)
                        com_ctc_ve_fail[reg,run,i]=0.66*com_ctc[reg,run,0]*com_ctc_timing[i]*(1-(ve[i,run]*ve_warm))
                    
        fac_cov=zeros((len(hbv_pars),runs,1))
        fac_ncov=zeros((len(hbv_pars),runs,1))
        com_cov=zeros((len(hbv_pars), runs, 1))
        com_ncov=zeros((len(hbv_pars), runs, 1))
     
        for reg in range(len(hbv_pars)):
            for run in range(runs): 
                fac_cov[reg,run,0]=sum(fac_cc_ve[reg,run,:])+sum(fac_ctc_ve[reg,run,:])
                fac_ncov[reg,run, 0]=sum(fac_cc_ve_fail[reg, run,:])+sum(fac_ctc_ve_fail[reg,run,:])
                com_cov[reg,run, 0]=sum(com_cc_ve[reg,run,:])+sum(com_ctc_ve[reg,run,:])+sum(com_cpad_ve[reg,run,:])
                com_ncov[reg, run, 0]=sum(com_cc_ve_fail[reg,run,:])+sum(com_ctc_ve_fail[reg,run,:])+sum(com_cpad_ve_fail[reg,run,:])
                
        #Data Needed to Initialize the Model
        P=1000
        
        prev=regions_perturbed[:,0,:].reshape(len(hbv_pars),runs,1)
        hbv3=regions_perturbed[:,2,:].reshape(len(hbv_pars),runs,1)
        hosp=regions_perturbed[:,4,:].reshape(len(hbv_pars),runs,1)
        comm=1-hosp
        
        hbe_prev=regions_perturbed[:,38,:].reshape(len(hbv_pars),runs,1)
        hbe_pos=regions_perturbed[:,39,:].reshape(len(hbv_pars),runs,1)
        hbe_neg=regions_perturbed[:,40,:].reshape(len(hbv_pars),runs,1)
        
        model_init=zeros((len(hbv_pars)*runs, len(t_steps), 14))
        
        #For simplicity, runs are in consecutive order and aggregated. Will split later on for results, tables and figures
        for reg in range(len(hbv_pars)):
            for run in range(runs):
                model_init[(runs*reg)+run, 0,0]=(P*hosp[reg,run,0]*hbv3[reg,run,0]*((prev[reg,run,0]*(hbe_prev[reg,run,0])*(1-hbe_pos[reg,run,0]))+(prev[reg,run,0]*hbe_prev[reg,run,0]*hbe_pos[reg,run,0]*(fac_cov[reg,run,0]))+(prev[reg,run,0]*(1-hbe_prev[reg,run,0])*(1-hbe_neg[reg,run,0]))+(prev[reg,run,0]*(1-hbe_prev[reg,run,0])*hbe_neg[reg,run,0]*fac_cov[reg,run,0])+(1-prev[reg,run,0])))+(P*comm[reg,run,0]*hbv3[reg,run,0]*((prev[reg,run,0]*hbe_prev[reg,run,0]*(1-hbe_pos[reg,run,0]))+(prev[reg,run,0]*hbe_prev[reg,run,0]*hbe_pos[reg,run,0]*com_cov[reg,run,0])+(prev[reg,run,0]*(1-hbe_prev[reg,run,0])*(1-hbe_neg[reg,run,0]))+(prev[reg,run,0]*(1-hbe_prev[reg,run,0])*hbe_neg[reg,run,0]*com_cov[reg,run,0])+(1-prev[reg,run,0]))) #Z
                model_init[(runs*reg)+run, 0,1]=(P*hosp[reg,run,0]*(1-hbv3[reg,run,0])*((prev[reg,run,0]*(hbe_prev[reg,run,0])*(1-hbe_pos[reg,run,0]))+(prev[reg,run,0]*hbe_prev[reg,run,0]*hbe_pos[reg,run,0]*(fac_cov[reg,run,0]))+(prev[reg,run,0]*(1-hbe_prev[reg,run,0])*(1-hbe_neg[reg,run,0]))+(prev[reg,run,0]*(1-hbe_prev[reg,run,0])*hbe_neg[reg,run,0]*fac_cov[reg,run,0])+(1-prev[reg,run,0])))+(P*comm[reg,run,0]*(1-hbv3[reg,run,0])*((prev[reg,run,0]*hbe_prev[reg,run,0]*(1-hbe_pos[reg,run,0]))+(prev[reg,run,0]*hbe_prev[reg,run,0]*hbe_pos[reg,run,0]*com_cov[reg,run,0])+(prev[reg,run,0]*(1-hbe_prev[reg,run,0])*(1-hbe_neg[reg,run,0]))+(prev[reg,run,0]*(1-hbe_prev[reg,run,0])*hbe_neg[reg,run,0]*com_cov[reg,run,0])+(1-prev[reg,run,0]))) #S
                model_init[(runs*reg)+run, 0,2]=(P*hosp[reg,run,0]*((prev[reg,run,0]*hbe_prev[reg,run,0]*hbe_pos[reg,run,0]*fac_ncov[reg,run,0])+(prev[reg,run,0]*(1-hbe_prev[reg,run,0])*hbe_neg[reg,run,0]*fac_ncov[reg,run,0])))+(P*comm[reg,run,0]*((prev[reg,run,0]*hbe_prev[reg,run,0]*hbe_pos[reg,run,0]*com_ncov[reg,run,0])+(prev[reg,run,0]*(1-hbe_prev[reg,run,0])*hbe_neg[reg,run,0]*com_ncov[reg,run,0]))) #L
        
        return model_init


def vaccine_costs(runs, hbv_pars, fac_cc, fac_ctc, com_cc, com_ctc, com_cpad,sat_costing, saturation, regions_perturbed):
    
    """Returns costs of each intervention, plus total cost, for each model scenario. Also returns the total coverage for 
    each scenario"""
    
    from numpy import array,zeros
    from math import log
    
    if runs ==1:
        
        c_fac_cc=array(hbv_pars.iloc[:,-9])             #Facility Cold Chain costs 
        c_fac_ctc=array(hbv_pars.iloc[:,-8])            #Facility CTC costs
        c_com_cc=array(hbv_pars.iloc[:,-7])             #Community Cold Chain costs
        c_com_ctc=array(hbv_pars.iloc[:,-6])            #Community CTC costs
        c_com_cpad=array(hbv_pars.iloc[:,-5])           #Community CPAD costs
        
        hosp_n=array(hbv_pars.iloc[:,14])*1000
        comm_n=1000-hosp_n
        
        vaccine_costs_perk=zeros((len(hbv_pars), 6))                #Modelling Costs, per 1000 births
        
        for reg in range(len(hbv_pars)):
            vaccine_costs_perk[reg,0]=min(hosp_n[reg]*fac_cc[reg]*c_fac_cc[reg], hosp_n[reg]*sat_costing[reg,0]*c_fac_cc[reg])       #Facility Cold Chain
            vaccine_costs_perk[reg,1]=min(hosp_n[reg]*fac_ctc[reg]*c_fac_ctc[reg], hosp_n[reg]*(sat_costing[reg,1]-sat_costing[reg,0])*c_fac_ctc[reg])   #Facility CTC
            vaccine_costs_perk[reg,2]=abs(-0.5*comm_n[reg]*sat_costing[reg,2]*c_com_cc[reg]*log((sat_costing[reg,2]*comm_n[reg]-(comm_n[reg]*com_cc[reg]))/(sat_costing[reg,2]*comm_n[reg]+(comm_n[reg]*com_cc[reg]))))   #Community, Cold-Chain
            vaccine_costs_perk[reg,3]=abs(-0.5*comm_n[reg]*(sat_costing[reg,3]-saturation[reg,2])*c_com_ctc[reg]*log((((sat_costing[reg,3]-saturation[reg,2])*comm_n[reg])-(comm_n[reg]*com_ctc[reg]))/(((sat_costing[reg,3]-saturation[reg,2])*comm_n[reg])+(comm_n[reg]*com_ctc[reg]))))
            vaccine_costs_perk[reg,4]=abs(-0.5*comm_n[reg]*(sat_costing[reg,4]-saturation[reg,3])*c_com_cpad[reg]*log((((sat_costing[reg,4]-saturation[reg,3])*comm_n[reg])-(comm_n[reg]*com_cpad[reg]))/(((sat_costing[reg,4]-saturation[reg,3])*comm_n[reg])+(comm_n[reg]*com_cpad[reg]))))
            
        for reg in range(len(hbv_pars)):
            vaccine_costs_perk[reg,5]=vaccine_costs_perk[reg,0]+vaccine_costs_perk[reg,1]+vaccine_costs_perk[reg,2]+vaccine_costs_perk[reg,3]+vaccine_costs_perk[reg,4]
        
        total_cov_prop=zeros(len(hbv_pars))
        
        for reg in range(len(hbv_pars)):
            total_cov_prop[reg]=(((fac_cc[reg]+fac_ctc[reg])*hosp_n[reg])/1000)+(((com_cc[reg]+com_ctc[reg]+com_cpad[reg])*comm_n[reg])/1000)
        
        return vaccine_costs_perk, total_cov_prop
    
    elif runs > 1:
        
        c_fac_cc=regions_perturbed[:,-9,:].reshape(len(hbv_pars),runs,1)             #Facility Cold Chain costs 
        c_fac_ctc=regions_perturbed[:,-8,:].reshape(len(hbv_pars),runs,1)            #Facility CTC costs
        c_com_cc=regions_perturbed[:,-7,:].reshape(len(hbv_pars),runs,1)             #Community Cold Chain costs
        c_com_ctc=regions_perturbed[:,-6,:].reshape(len(hbv_pars),runs,1)            #Community CTC costs
        c_com_cpad=regions_perturbed[:,-5,:].reshape(len(hbv_pars),runs,1)           #Community CPAD costs
        
        hosp_n=(regions_perturbed[:,4,:]*1000).reshape(len(hbv_pars),runs,1)
        comm_n=1000-hosp_n
        
        vaccine_costs_perk=zeros((len(hbv_pars),runs,6))
        
        for reg in range(len(hbv_pars)):
            for run in range(runs):
                vaccine_costs_perk[reg,run,0]=min(hosp_n[reg,run,0]*fac_cc[reg,run,0]*c_fac_cc[reg,run,0], hosp_n[reg,run,0]*sat_costing[reg,run,0]*c_fac_cc[reg,run,0])       #Facility Cold Chain
                vaccine_costs_perk[reg,run,1]=min(hosp_n[reg,run,0]*fac_ctc[reg,run,0]*c_fac_ctc[reg,run,0], hosp_n[reg,run,0]*(sat_costing[reg,run,1]-sat_costing[reg,run,0])*c_fac_ctc[reg,run,0])   #Facility CTC
                vaccine_costs_perk[reg,run,2]=abs(-0.5*comm_n[reg,run,0]*sat_costing[reg,run,2]*c_com_cc[reg,run,0]*log((sat_costing[reg,run,2]*comm_n[reg,run,0]-(comm_n[reg,run,0]*com_cc[reg,run,0]))/(sat_costing[reg,run,2]*comm_n[reg,run,0]+(comm_n[reg,run,0]*com_cc[reg,run,0]))))   #Community, Cold-Chain
                vaccine_costs_perk[reg,run,3]=abs(-0.5*comm_n[reg,run,0]*(sat_costing[reg,run,3]-saturation[reg,run,2])*c_com_ctc[reg,run,0]*log((((sat_costing[reg,run,3]-saturation[reg,run,2])*comm_n[reg,run,0])-(comm_n[reg,run,0]*com_ctc[reg,run,0]))/(((sat_costing[reg,run,3]-saturation[reg,run,2])*comm_n[reg,run,0])+(comm_n[reg,run,0]*com_ctc[reg,run,0]))))
                vaccine_costs_perk[reg,run,4]=abs(-0.5*comm_n[reg,run,0]*(sat_costing[reg,run,4]-saturation[reg,run,3])*c_com_cpad[reg,run,0]*log((((sat_costing[reg,run,4]-saturation[reg,run,3])*comm_n[reg,run,0])-(comm_n[reg,run,0]*com_cpad[reg,run,0]))/(((sat_costing[reg,run,4]-saturation[reg,run,3])*comm_n[reg,run,0])+(comm_n[reg,run,0]*com_cpad[reg,run,0]))))
        
        for reg in range(len(hbv_pars)):
            for run in range(runs):
                vaccine_costs_perk[reg,run,5]=vaccine_costs_perk[reg,run,0]+vaccine_costs_perk[reg,run,1]+vaccine_costs_perk[reg,run,2]+vaccine_costs_perk[reg,run,3]+vaccine_costs_perk[reg,run,4]
                
        total_cov_prop=zeros((len(hbv_pars),runs,1))
        
        for reg in range (len(hbv_pars)):
            for run in range(runs):
                total_cov_prop[reg,run,0]=(((fac_cc[reg,run,0]+fac_ctc[reg,run,0])*hosp_n[reg,run,0])/1000)+(((com_cc[reg,run,0]+com_ctc[reg,run,0]+com_cpad[reg,run,0])*comm_n[reg,run,0])/1000)
        
        return vaccine_costs_perk, total_cov_prop


def hbv_model(runs,t_steps,dt,model_init, hbv_pars, regions_perturbed, globz, globz_perturbed,discount):
    
    """Runs the disease simulation. Returns YLL, YLD, DALYs, Disease Cost and Total Costs as part of the array"""
   
    from numpy import array, zeros, exp,ceil 
    
   
    if runs==1:
        
        #Age Dependent Acute to Chronic Transition
        p_AC=zeros((len(t_steps),1))
        
        for idx,years in enumerate(t_steps):
            if years < 0.5:
                p_AC[idx,0]=globz.iloc[16,1]
            else:
                p_AC[idx,0]=exp(-0.645*years**(0.455))

        #All Cause Mortality Rates
        acm=zeros((len(hbv_pars),len(t_steps)))
        
        for reg in range(len(hbv_pars)):
            for idx,years in enumerate(t_steps):
                if years <1:
                    acm[reg,idx]=hbv_pars.iloc[reg,15]
                
                elif years >=1 and years <=4:
                    acm[reg,idx]=hbv_pars.iloc[reg,16]
                
                else:
                    acm[reg,idx]=hbv_pars.iloc[reg,16+int(ceil(years)/5)]

       
        #Disease Progression (Global)
        r_LA=globz.iloc[14,1]
        r_AC=globz.iloc[15,1]
        r_AZ=globz.iloc[17,1]
        r_CZ=globz.iloc[18,1]
        r_CCC=globz.iloc[19,1]
        r_CDC=globz.iloc[20,1]
        r_DCHCC=globz.iloc[21,1]
        r_CHCC=globz.iloc[22,1]
        r_CHC=globz.iloc[23,1]
        
        #HBV Mortality (Global)
        mu_A=globz.iloc[24,1]
        mu_C=globz.iloc[25,1]
        mu_CC=globz.iloc[26,1]
        mu_DC=globz.iloc[27,1]
        mu_HCC=globz.iloc[28,1]
        
        #DALYs (Global)
        DALY_A=globz.iloc[4,1]
        DALY_C=globz.iloc[5,1]
        DALY_CC=globz.iloc[6,1]
        DALY_DC=globz.iloc[7,1]
        DALY_HCC=globz.iloc[8,1]
        
        #Disease Cost (Setting Specific) #-7 is cost of HCC, -11 is cost of Acute
        c_A=array(hbv_pars.iloc[:,37])
        c_C=array(hbv_pars.iloc[:,38])
        c_CC=array(hbv_pars.iloc[:,39])
        c_DC=array(hbv_pars.iloc[:,40])
        c_HCC=array(hbv_pars.iloc[:,41])
        
        #Model Run Here
        for reg in range(len(hbv_pars)):
            for t in range(len(t_steps)-1):
                #0. Immune (Z)
                model_init[reg, t+1, 0]=model_init[reg, t,0]+(dt*(1-p_AC[t,0])*r_AZ*model_init[reg,t,3])+(dt*r_CZ*model_init[reg,t,4])-(dt*acm[reg,t]*model_init[reg,t,0])
                #1. Susceptible (S)
                model_init[reg, t+1, 1]=model_init[reg,t,1]-(dt*acm[reg,t]*model_init[reg,t,1])
                #2. Latent (L)
                model_init[reg, t+1, 2]=max(model_init[reg,t,2]-(dt*r_LA*model_init[reg,t,2])-(dt*acm[reg,t]*model_init[reg,t,2]),0)
                #3. Acute (A)
                model_init[reg, t+1, 3]=max(model_init[reg,t,3]+(dt*r_LA*model_init[reg,t,2])-(dt*r_AC*p_AC[t,0]*model_init[reg,t,3])-(dt*r_AZ*(1-p_AC[t,0])*model_init[reg,t,3])-(dt*mu_A*model_init[reg,t,3])-(dt*acm[reg,t]*model_init[reg,t,3]),0)
                #4. Chronic (C)
                model_init[reg, t+1, 4]=model_init[reg,t,4]+(dt*p_AC[t,0]*r_AC*model_init[reg,t,3])-(dt*(r_CCC+r_CZ+r_CHC)*model_init[reg,t,4])-(dt*mu_C*model_init[reg,t,4])-(dt*acm[reg,t]*model_init[reg,t,4])
                #5. Compensated Cirrhosis (CC)
                model_init[reg, t+1, 5]=model_init[reg,t,5]+(dt*r_CCC*model_init[reg,t,4])-(dt*(r_CDC+r_CHCC)*model_init[reg,t,5])-(dt*mu_CC*model_init[reg,t,5])-(dt*acm[reg,t]*model_init[reg,t,5])
                #6. Decompensated Cirrhosis (DC)
                model_init[reg, t+1, 6]=model_init[reg,t,6]+(dt*r_CDC*model_init[reg,t,5])-(dt*r_DCHCC*model_init[reg,t,6])-(dt*mu_DC*model_init[reg,t,6])-(dt*acm[reg,t]*model_init[reg,t,6])
                #7. Hepatocellular Carcinoma (HCC)
                model_init[reg, t+1, 7]=model_init[reg,t,7]+(dt*r_CHC*model_init[reg,t,4])+(dt*r_CHCC*model_init[reg,t,5])+(dt*r_DCHCC*model_init[reg,t,6])-(dt*mu_HCC*model_init[reg,t,7])-(dt*acm[reg,t]*model_init[reg,t,7])
                #8. Deaths due to HBV (DD)
                model_init[reg, t+1, 8]=model_init[reg,t,8]+(dt*mu_A*model_init[reg,t,3])+(dt*mu_C*model_init[reg,t,4])+(dt*mu_CC*model_init[reg,t,5])+(dt*mu_DC*model_init[reg,t,6])+(dt*mu_HCC*model_init[reg,t,7])
                #9. Years Life Lost due to HBV (discounted 3% per annum)
                model_init[reg, t+1, 9]=model_init[reg,t,9]+exp(-discount*(t*dt))*model_init[reg,t,8]*dt
                #10. Years Lost to Disability from HBV (discounted 3% per annum)
                model_init[reg, t+1, 10]=model_init[reg,t,10]+exp(-discount*(t*dt))*(DALY_A*model_init[reg,t,3]+DALY_C*model_init[reg,t,4]+DALY_CC*model_init[reg,t,5]+DALY_DC*model_init[reg,t,6]+DALY_HCC*model_init[reg,t,7])*(2997/len(t_steps))# normalizing factor, based on equivalent ODES as calculation is sensitive to step size
                #11. Cost of Disease (discounted 3% per annum)
                model_init[reg, t+1, 11]=model_init[reg,t,11]+exp(-discount*(t*dt))*(c_A[reg]*model_init[reg,t,3]*dt+c_C[reg]*model_init[reg,t,4]*dt+c_CC[reg]*model_init[reg,t,5]*dt+c_DC[reg]*model_init[reg,t,6]*dt+c_HCC[reg]*model_init[reg,t,7]*dt)
                #12. Total Cost (disease discounted 3% per annum, vaccine one off and not discounted)
                model_init[reg, t+1, 12]=model_init[reg,t,12]+exp(-discount*(t*dt))*(c_A[reg]*model_init[reg,t,3]*dt+c_C[reg]*model_init[reg,t,4]*dt+c_CC[reg]*model_init[reg,t,5]*dt+c_DC[reg]*model_init[reg,t,6]*dt+c_HCC[reg]*model_init[reg,t,7]*dt)
                #13. Disability Adjusted Life Years (DALYS)
                model_init[reg, t+1, 13]=model_init[reg,t,9]+model_init[reg,t,10]
        
        return model_init
    
    elif runs > 1:
        
        #Age Dependent Acute to Chronic Transition
        p_AC=zeros((runs,len(t_steps),1))
        
        for run in range(runs):
            for idx,years in enumerate(t_steps):
                if years < 0.5:
                    p_AC[run,idx,0]=globz_perturbed.iloc[12,run]
                else:
                    p_AC[run,idx,0]=exp(-0.645*years**(0.455))
         
        #All Cause Mortality

        acm=zeros((len(hbv_pars),runs,len(t_steps)))
        
        for reg in range(len(hbv_pars)):
            for run in range(runs):
                for idx,years in enumerate(t_steps):
                    if years < 1:
                        acm[reg,run,idx]=regions_perturbed[reg,5,run]
                    elif years >=1 and years <= 4 :
                        acm[reg,run,idx]=regions_perturbed[reg,6,run]
                    else:
                        acm[reg,run,idx]=regions_perturbed[reg,6+int(ceil(years)/5),run]

         #Disease Progression (Global)
        r_LA=array(globz_perturbed.iloc[10,:])
        r_AC=array(globz_perturbed.iloc[11,:])
        r_AZ=array(globz_perturbed.iloc[13,:])
        r_CZ=array(globz_perturbed.iloc[14,:])
        r_CCC=array(globz_perturbed.iloc[15,:])
        r_CDC=array(globz_perturbed.iloc[16,:])
        r_DCHCC=array(globz_perturbed.iloc[17,:])
        r_CHCC=array(globz_perturbed.iloc[18,:])
        r_CHC=array(globz_perturbed.iloc[19,:])
         
        #HBV Mortality (Global)
        mu_A=array(globz_perturbed.iloc[20,:])
        mu_C=array(globz_perturbed.iloc[21,:])
        mu_CC=array(globz_perturbed.iloc[22,:])
        mu_DC=array(globz_perturbed.iloc[23,:])
        mu_HCC=array(globz_perturbed.iloc[24,:])
        
        #DALYs (Global)
        DALY_A=array(globz_perturbed.iloc[0,:])
        DALY_C=array(globz_perturbed.iloc[1,:])
        DALY_CC=array(globz_perturbed.iloc[2,:])
        DALY_DC=array(globz_perturbed.iloc[3,:])
        DALY_HCC=array(globz_perturbed.iloc[4,:])
        
        #Disease Cost (Setting Specific)
        c_A=regions_perturbed[:,27,:].reshape(len(hbv_pars),runs,1)
        c_C=regions_perturbed[:,28,:].reshape(len(hbv_pars),runs,1)
        c_CC=regions_perturbed[:,29,:].reshape(len(hbv_pars),runs,1)
        c_DC=regions_perturbed[:,30,:].reshape(len(hbv_pars),runs,1)
        c_HCC=regions_perturbed[:,31,:].reshape(len(hbv_pars),runs,1)
        
         #Model Run Here
        for reg in range(len(hbv_pars)):
            for run in range(runs):
                for t in range(len(t_steps)-1):
                    #0. Immune (Z)
                    model_init[(reg*runs)+run, t+1, 0]=model_init[(reg*runs)+run,t,0]+(dt*(1-p_AC[run,t,0])*r_AZ[run]*model_init[(reg*runs)+run,t,3])+(dt*r_CZ[run]*model_init[(reg*runs)+run,t,4])-(dt*acm[reg,run,t]*model_init[(reg*runs)+run,t,0])
                    #1. Susceptible (S)
                    model_init[(reg*runs)+run, t+1, 1]=model_init[(reg*runs)+run,t,1]-(dt*acm[reg,run,t]*model_init[(reg*runs)+run,t,1])
                    #2. Latent (L)
                    model_init[(reg*runs)+run, t+1, 2]=max(model_init[(reg*runs)+run,t,2]-(dt*r_LA[run]*model_init[(reg*runs)+run,t,2])-(dt*acm[reg,run,t]*model_init[(reg*runs)+run,t,2]),0)
                    #3. Acute (A)
                    model_init[(reg*runs)+run, t+1, 3]=max(model_init[(reg*runs)+run,t,3]+(dt*r_LA[run]*model_init[(reg*runs)+run,t,2])-(dt*r_AC[run]*p_AC[run,t,0]*model_init[(reg*runs)+run,t,3])-(dt*r_AZ[run]*(1-p_AC[run,t,0])*model_init[(reg*runs)+run,t,3])-(dt*mu_A[run]*model_init[(reg*runs)+run,t,3])-(dt*acm[reg,run,t]*model_init[(reg*runs)+run,t,3]),0)
                    #4. Chronic (C)
                    model_init[(reg*runs)+run, t+1, 4]=model_init[(reg*runs)+run,t,4]+(dt*p_AC[run,t,0]*r_AC[run]*model_init[(reg*runs)+run,t,3])-(dt*(r_CCC[run]+r_CZ[run]+r_CHC[run])*model_init[(reg*runs)+run,t,4])-(dt*mu_C[run]*model_init[(reg*runs)+run,t,4])-(dt*acm[reg,run,t]*model_init[(reg*runs)+run,t,4])
                    #5. Compensated Cirrhosis (CC)
                    model_init[(reg*runs)+run, t+1, 5]=model_init[(reg*runs)+run,t,5]+(dt*r_CCC[run]*model_init[(reg*runs)+run,t,4])-(dt*(r_CDC[run]+r_CHCC[run])*model_init[(reg*runs)+run,t,5])-(dt*mu_CC[run]*model_init[(reg*runs)+run,t,5])-(dt*acm[reg,run,t]*model_init[(reg*runs)+run,t,5])
                    #6. Decompensated Cirrhosis (DC)
                    model_init[(reg*runs)+run, t+1, 6]=model_init[(reg*runs)+run,t,6]+(dt*r_CDC[run]*model_init[(reg*runs)+run,t,5])-(dt*r_DCHCC[run]*model_init[(reg*runs)+run,t,6])-(dt*mu_DC[run]*model_init[(reg*runs)+run,t,6])-(dt*acm[reg,run,t]*model_init[(reg*runs)+run,t,6])
                    #7. Hepatocellular Carcinoma (HCC)
                    model_init[(reg*runs)+run, t+1, 7]=model_init[(reg*runs)+run,t,7]+(dt*r_CHC[run]*model_init[(reg*runs)+run,t,4])+(dt*r_CHCC[run]*model_init[(reg*runs)+run,t,5])+(dt*r_DCHCC[run]*model_init[(reg*runs)+run,t,6])-(dt*mu_HCC[run]*model_init[(reg*runs)+run,t,7])-(dt*acm[reg,run,t]*model_init[(reg*runs)+run,t,7])
                    #8. Deaths due to HBV (DD)
                    model_init[(reg*runs)+run, t+1, 8]=model_init[(reg*runs)+run,t,8]+(dt*mu_A[run]*model_init[(reg*runs)+run,t,3])+(dt*mu_C[run]*model_init[(reg*runs)+run,t,4])+(dt*mu_CC[run]*model_init[(reg*runs)+run,t,5])+(dt*mu_DC[run]*model_init[(reg*runs)+run,t,6])+(dt*mu_HCC[run]*model_init[(reg*runs)+run,t,7])
                    #9. Years Life Lost due to HBV (discounted 3% per annum)
                    model_init[(reg*runs)+run, t+1, 9]=model_init[(reg*runs)+run,t,9]+exp(-discount*(t*dt))*model_init[(reg*runs)+run,t,8]*dt
                    #10. Years Lost to Disability from HBV (discounted 3% per annum)
                    model_init[(reg*runs)+run, t+1, 10]=model_init[(reg*runs)+run,t,10]+exp(-discount*(t*dt))*(DALY_A[run]*model_init[(reg*runs)+run,t,3]+DALY_C[run]*model_init[(reg*runs)+run,t,4]+DALY_CC[run]*model_init[(reg*runs)+run,t,5]+DALY_DC[run]*model_init[(reg*runs)+run,t,6]+DALY_HCC[run]*model_init[(reg*runs)+run,t,7])*(2997/len(t_steps))# normalizing factor, based on equivalent ODES as calculation is sensitive to step size
                    #11. Cost of Disease (discounted 3% per annum)
                    model_init[(reg*runs)+run, t+1, 11]=model_init[(reg*runs)+run,t,11]+exp(-discount*(t*dt))*(c_A[reg,run,0]*model_init[(reg*runs)+run,t,3]*dt+c_C[reg,run,0]*model_init[(reg*runs)+run,t,4]*dt+c_CC[reg,run,0]*model_init[(reg*runs)+run,t,5]*dt+c_DC[reg,run,0]*model_init[(reg*runs)+run,t,6]*dt+c_HCC[reg,run,0]*model_init[(reg*runs)+run,t,7]*dt)
                    #12. Total Cost (disease discounted 3% per annum, vaccine one off and not discounted)
                    model_init[(reg*runs)+run, t+1, 12]=model_init[(reg*runs)+run,t,12]+exp(-discount*(t*dt))*(c_A[reg,run,0]*model_init[(reg*runs)+run,t,3]*dt+c_C[reg,run,0]*model_init[(reg*runs)+run,t,4]*dt+c_CC[reg,run,0]*model_init[(reg*runs)+run,t,5]*dt+c_DC[reg,run,0]*model_init[(reg*runs)+run,t,6]*dt+c_HCC[reg,run,0]*model_init[(reg*runs)+run,t,7]*dt)
                    #13. Disability Adjusted Life Years (DALYS)
                    model_init[(reg*runs)+run, t+1, 13]=model_init[(reg*runs)+run,t,9]+model_init[(reg*runs)+run,t,10]
        
        return model_init
