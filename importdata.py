###Importing the data from excel###

def data_import (working_dir, region, runs):
   
    """Imports the excel file, and extracts global paramaters, 
    regional/country level data and uncertainity bounds"""
    
    import pandas as pd  
    xls=pd.ExcelFile("hbvbd_model_inputs_revision.xlsx")        #updated 23/03/20
  
    if region==1:
        level="regions"
    elif region==0:
        level="countries"
    
    if runs == 1:
        globz=pd.read_excel(xls, "global", index_col=None).iloc[0:30, 1:3]      #eAg is no longer a global variable for sensitivity analysis
        hbv_pars=pd.DataFrame(pd.read_excel(xls, level))
        return globz, hbv_pars
    elif runs >1: 
       globz=pd.read_excel(xls, "global", index_col=None).iloc[0:30, 1:5]
       hbv_pars=pd.DataFrame(pd.read_excel(xls, level))
       hbv_pars_lb=pd.DataFrame(pd.read_excel(xls, level+"_lb"))                #hbv_pars for countries is now limited to LMICs, and is longer by ~5 columns
       hbv_pars_ub=pd.DataFrame(pd.read_excel(xls, level+"_ub"))
       
       
       
       return globz, hbv_pars, hbv_pars_lb, hbv_pars_ub
        
    


   
    