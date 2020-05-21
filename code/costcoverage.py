##%% Cost Coverage Functions
from numpy import exp

def cost_coverage_fac_cc(expenditure, reg, hosp_n, vac_cost_opt, sat_opt):
    """ Cost Coverage for Facility Cold Chain Coverage (Linear)"""
    n_covered=expenditure[0]/vac_cost_opt[reg,0]
    prop_covered=expenditure[0]/(vac_cost_opt[reg,0]*hosp_n[reg])
                
    return n_covered, prop_covered
            
def cost_coverage_fac_ctc(expenditure, reg, hosp_n, vac_cost_opt, sat_opt):
    """Cost Coverage for Facility CTC (Linear)"""
    n_covered=expenditure[1]/vac_cost_opt[reg,1]
    prop_covered=expenditure[1]/(vac_cost_opt[reg,1]*hosp_n[reg])
                
    return n_covered, prop_covered
            
def cost_coverage_com_cc(expenditure, reg, comm_n, vac_cost_opt, sat_opt):
    """Cost Coverage Curve for Community, Cold Chain (Logistic)"""
    n_covered=max((2*sat_opt[reg,2]/(1+exp(-2*expenditure[2]/(comm_n[reg]*sat_opt[reg,2]*vac_cost_opt[reg,2])))-sat_opt[reg,2])*comm_n[reg],0)
    prop_covered=max((2*sat_opt[reg,2]/(1+exp(-2*expenditure[2]/(comm_n[reg]*sat_opt[reg,2]*vac_cost_opt[reg,2])))-sat_opt[reg,2]),0)
                
    return n_covered, prop_covered
            
def cost_coverage_com_ctc(expenditure, reg, comm_n, vac_cost_opt, sat_opt):
    """Cost Coverage Curve for Community CTC (Logistic)"""
    adj_cc_n, adj_cc_prop=cost_coverage_com_cc(expenditure, reg, comm_n, vac_cost_opt, sat_opt) #Adjusts for Community Cold-Chain Coverage
                
    n_covered=max((2*(sat_opt[reg,3]-adj_cc_prop)/(1+exp(-2*expenditure[3]/(comm_n[reg]*(sat_opt[reg,3]-adj_cc_prop)*vac_cost_opt[reg,3])))-(sat_opt[reg,3]-adj_cc_prop))*comm_n[reg],0)
    prop_covered=max((2*(sat_opt[reg,3]-adj_cc_prop)/(1+exp(-2*expenditure[3]/(comm_n[reg]*(sat_opt[reg,3]-adj_cc_prop)*vac_cost_opt[reg,3])))-(sat_opt[reg,3]-adj_cc_prop)),0)
                
    return n_covered, prop_covered
            
def cost_coverage_com_cpad(expenditure, reg, comm_n, vac_cost_opt, sat_opt):
    """Cost Coverage Curve for Commuity CPAD (Logistic)"""
    adj_cc_n, adj_cc_prop=cost_coverage_com_cc(expenditure, reg, comm_n, vac_cost_opt, sat_opt) #Adjusts for Community Cold Chain Coverage
    adj_ctc_n, adj_ctc_prop=cost_coverage_com_ctc(expenditure, reg, comm_n, vac_cost_opt, sat_opt) # Adjusts for Community CTC Coverage
                
    if sat_opt[reg,3]==sat_opt[reg,4]:
        n_covered=0
        prop_covered=0
    else:
        n_covered=max((2*(sat_opt[reg,4]-adj_cc_prop-adj_ctc_prop)/(1+exp(-2*expenditure[4]/(comm_n[reg]*(sat_opt[reg,4]-adj_cc_prop-adj_ctc_prop)*vac_cost_opt[reg,4])))-(sat_opt[reg,4]-(adj_cc_prop+adj_ctc_prop)))*comm_n[reg],0)
        prop_covered=max((2*(sat_opt[reg,4]-adj_cc_prop-adj_ctc_prop)/(1+exp(-2*expenditure[4]/(comm_n[reg]*(sat_opt[reg,4]-adj_cc_prop-adj_ctc_prop)*vac_cost_opt[reg,4])))-(sat_opt[reg,4]-(adj_cc_prop+adj_ctc_prop))),0)
                
    return n_covered, prop_covered
            
            
#%% Objective Functions
def objective_coverage(expenditure):
    """Objective function for optimisation of baseline coverage. Note that optimisation uses minimisation, 
    so to maximise coverage, need to return the negative sum of coverage"""
    cc_fac_n, cc_fac_prop=cost_coverage_fac_cc(expenditure, reg, hosp_n, vac_cost_opt, sat_opt)
    ctc_fac_n, ctc_fac_prop=cost_coverage_fac_ctc(expenditure, reg, hosp_n, vac_cost_opt, sat_opt)
    cc_com_n, ctc_com_prop=cost_coverage_com_cc(expenditure, reg, comm_n, vac_cost_opt, sat_opt)
    ctc_com_n, ctc_com_prop=cost_coverage_com_ctc(expenditure, reg, comm_n, vac_cost_opt, sat_opt)
    cpad_com_n, cpad_com_prop=cost_coverage_com_cpad(expenditure, reg, comm_n, vac_cost_opt, sat_opt)
                
    return -(cc_fac_n+ctc_fac_n+cc_com_n+ctc_com_n+cpad_com_n) 
            
def objective_90_coverage(expenditure):
    """Objective function to return the smallest total vaccine expenditure to acheieve 90% coverage. As the goal is
    to achieve a value as small as possible, no negative sign is required for this objective function."""
    return sum(expenditure)
            
#%% Constraints
def constraint_90_coverage(expenditure):
    """Equality Constraint, as objective_coverage returns a negative value, coverage can just be added to target"""
    return 900+objective_coverage(expenditure) #As above, this value is negative, so needs to be plus or coverage will just be 0
            
def constraint_baseline_budget(expenditure,reg):
    """Equality Constraint; ensures that optimised expenditure for baseline optimisations do not exceed initial baseline
    expenditure"""
    return baseline_expenditure[reg,0]-sum(expenditure)
            
def constraint_fac_cc_coverage(expenditure, sat_opt,reg):
    """Inequality Constraint; ensures coverage of cold-chain in facilities does not exceed maximum coverage 
    levels"""
    cc_cov_n, cc_cov_prop=cost_coverage_fac_cc(expenditure, reg, hosp_n, vac_cost_opt, sat_opt)
    return sat_opt[reg,0]-cc_cov_prop
                
def constraint_fac_cc_ctc_coverage(expenditure, sat_opt,reg):
    """Inequality constraint; ensures coverage of CTC in facilities cannot exceed maximum coverage levels"""
    ctc_cov_n, ctc_cov_prop=cost_coverage_fac_ctc(expenditure, reg, hosp_n, vac_cost_opt, sat_opt)
    return (sat_opt[reg,1]-sat_opt[reg,0])-ctc_cov_prop
            
def constraint_com_cc_ctc_coverage(expenditure, sat_opt, reg):
    """Inequality constraint; ensures coverage of cold-chain and CTC in the community cannot exceed
    maximum coverage levels of vaccines delivered by professional health workers"""
    cc_cov_n, cc_cov_prop=cost_coverage_com_cc(expenditure, reg, comm_n, vac_cost_opt, sat_opt)
    ctc_cov_n, ctc_cov_prop=cost_coverage_com_ctc(expenditure, reg, comm_n, vac_cost_opt, sat_opt)
                
    return sat_opt[reg,3]-(cc_cov_prop+ctc_cov_prop)
            
def constraint_com_cc_ctc_cpad_coverage(expenditure,sat_opt,reg):
    """Inequality constraint; ensures coverage in the community (Cold-Chain, CTC and CPAD) cannot exceed maximum 
    coverage levels of vaccinations in the community"""
    cc_cov_n, cc_cov_prop=cost_coverage_com_cc(expenditure, reg, comm_n, vac_cost_opt, sat_opt)
    ctc_cov_n, ctc_cov_prop=cost_coverage_com_ctc(expenditure, reg, comm_n, vac_cost_opt, sat_opt)
    cpad_cov_n, cpad_cov_prop=cost_coverage_com_cpad(expenditure, reg, comm_n, vac_cost_opt, sat_opt)
                
    return sat_opt[reg,4]-(cc_cov_prop+ctc_cov_prop+cpad_cov_prop)
            
#%% Coverage from Costs
def coverage_from_opt(expenditure):
    """Returns coverage of each vaccine from optimisation, both as number and proportion covered. Proportion covered is ready for for modelling"""
    cc_fac_n, cc_fac_prop=cost_coverage_fac_cc(expenditure, reg, hosp_n, vac_cost_opt, sat_opt)
    ctc_fac_n, ctc_fac_prop=cost_coverage_fac_ctc(expenditure, reg, hosp_n, vac_cost_opt, sat_opt)
    cc_com_n, cc_com_prop=cost_coverage_com_cc(expenditure, reg, comm_n, vac_cost_opt, sat_opt)
    ctc_com_n, ctc_com_prop=cost_coverage_com_ctc(expenditure, reg, comm_n, vac_cost_opt, sat_opt)
    cpad_com_n, cpad_com_prop=cost_coverage_com_cpad(expenditure, reg, comm_n, vac_cost_opt, sat_opt)
                
    return(cc_fac_n, cc_fac_prop, ctc_fac_n, ctc_fac_prop, cc_com_n, cc_com_prop, ctc_com_n, ctc_com_prop, cpad_com_n, cpad_com_prop)
