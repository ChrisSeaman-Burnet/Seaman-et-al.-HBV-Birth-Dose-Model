#Results Tables and Figures

######################### Main Text, Figures and Supplementary Tables ####################################################
from matplotlib import pyplot as plt
import pandas as pd
from numpy import zeros, array,percentile,arange, ceil

#Figure 2: Cost-Effectiveness Modelling Scatter Plot (GBD Regions)
if region==1:
    if runs > 1:
        
        ce_modelling_results_cost=zeros((len(hbv_pars),runs,8))     #baseline, s1,22,s3,s4,bl_opt, optimised_bl, 90%_opt
        ce_modelling_results_daly=zeros((len(hbv_pars),runs,8))
        
        for reg in range (len(hbv_pars)):
            for run in range(runs):
                ce_modelling_results_cost[reg,run,0]=baseline[(reg*runs)+run, len(t_steps)-1, 12]
                ce_modelling_results_cost[reg,run,1]=max_current[(reg*runs)+run, len(t_steps)-1, 12]
                ce_modelling_results_cost[reg,run,2]=ctc_facilities[(reg*runs)+run, len(t_steps)-1, 12]
                ce_modelling_results_cost[reg,run,3]=ctc_community[(reg*runs)+run, len(t_steps)-1, 12]
                ce_modelling_results_cost[reg,run,4]=cpad_community[(reg*runs)+run, len(t_steps)-1, 12]
                ce_modelling_results_cost[reg,run,5]=baseline_opt[(reg*runs)+run, len(t_steps)-1, 12]
                ce_modelling_results_cost[reg,run,6]=baseline_optimised[(reg*runs)+run, len(t_steps)-1, 12]
                ce_modelling_results_cost[reg,run,7]=optimised_90[(reg*runs)+run, len(t_steps)-1, 12]
                
                ce_modelling_results_daly[reg,run,0]=baseline[(reg*runs)+run, len(t_steps)-1, 13]
                ce_modelling_results_daly[reg,run,1]=max_current[(reg*runs)+run, len(t_steps)-1, 13]
                ce_modelling_results_daly[reg,run,2]=ctc_facilities[(reg*runs)+run, len(t_steps)-1, 13]
                ce_modelling_results_daly[reg,run,3]=ctc_community[(reg*runs)+run, len(t_steps)-1, 13]
                ce_modelling_results_daly[reg,run,4]=cpad_community[(reg*runs)+run, len(t_steps)-1, 13]
                ce_modelling_results_daly[reg,run,5]=baseline_opt[(reg*runs)+run, len(t_steps)-1, 13]
                ce_modelling_results_daly[reg,run,6]=baseline_optimised[(reg*runs)+run, len(t_steps)-1, 13]
                ce_modelling_results_daly[reg,run,7]=optimised_90[(reg*runs)+run, len(t_steps)-1, 13]
        
        for reg in range(len(hbv_pars)):
            for run in range (runs):
                ce_modelling_results_cost[reg,run,1]=ce_modelling_results_cost[reg,run,1]-ce_modelling_results_cost[reg,run,0]
                ce_modelling_results_cost[reg,run,2]=ce_modelling_results_cost[reg,run,2]-ce_modelling_results_cost[reg,run,0]
                ce_modelling_results_cost[reg,run,3]=ce_modelling_results_cost[reg,run,3]-ce_modelling_results_cost[reg,run,0]
                ce_modelling_results_cost[reg,run,4]=ce_modelling_results_cost[reg,run,4]-ce_modelling_results_cost[reg,run,0]
                ce_modelling_results_cost[reg,run,6]=ce_modelling_results_cost[reg,run,6]-ce_modelling_results_cost[reg,run,5]
                ce_modelling_results_cost[reg,run,7]=ce_modelling_results_cost[reg,run,7]-ce_modelling_results_cost[reg,run,5]
                
                ce_modelling_results_daly[reg,run,1]=ce_modelling_results_daly[reg,run,1]-ce_modelling_results_daly[reg,run,0]
                ce_modelling_results_daly[reg,run,2]=ce_modelling_results_daly[reg,run,2]-ce_modelling_results_daly[reg,run,0]
                ce_modelling_results_daly[reg,run,3]=ce_modelling_results_daly[reg,run,3]-ce_modelling_results_daly[reg,run,0]
                ce_modelling_results_daly[reg,run,4]=ce_modelling_results_daly[reg,run,4]-ce_modelling_results_daly[reg,run,0]
                ce_modelling_results_daly[reg,run,6]=ce_modelling_results_daly[reg,run,6]-ce_modelling_results_daly[reg,run,5]
                ce_modelling_results_daly[reg,run,7]=ce_modelling_results_daly[reg,run,7]-ce_modelling_results_daly[reg,run,5]
        
        figure_2=zeros((len(hbv_pars), 6, 6))  # DALY median for scatter (0), error for error bars (1,2); diff_cost median for scatter (3), error for error bars (4,5) 
        
        for reg in range(len(hbv_pars)):
            #Cost-Effectiveness Modelling
            for i in range(4):
                #Difference in total costs
                figure_2[reg,i,0]=percentile(ce_modelling_results_cost[reg,:,(i+1)],50)
                figure_2[reg,i,1]=abs(figure_2[reg,i,0]-percentile(ce_modelling_results_cost[reg,:,(i+1)],25))
                figure_2[reg,i,2]=abs(percentile(ce_modelling_results_cost[reg,:,(i+1)],75)-figure_2[reg,i,0])
                
                #Difference in DALYs
                figure_2[reg,i,3]=percentile(ce_modelling_results_daly[reg,:,(i+1)],50)
                figure_2[reg,i,4]=abs(figure_2[reg,i,3]-percentile(ce_modelling_results_daly[reg,:,(i+1)],25))
                figure_2[reg,i,5]=abs(percentile(ce_modelling_results_daly[reg,:,(i+1)],75)-figure_2[reg,i,3])
           
            #Optimisations
            for i in range(4,6):
                #Difference in Total Costs
                figure_2[reg,i,0]=percentile(ce_modelling_results_cost[reg,:,(i+2)],50)
                figure_2[reg,i,1]=abs(figure_2[reg,i,0]-percentile(ce_modelling_results_cost[reg,:,(i+2)],25))
                figure_2[reg,i,2]=abs(percentile(ce_modelling_results_cost[reg,:,(i+2)],75)-figure_2[reg,i,0])
                
                #Difference in DALYs
                figure_2[reg,i,3]=percentile(ce_modelling_results_daly[reg,:,(i+2)],50)
                figure_2[reg,i,4]=abs(figure_2[reg,i,3]-percentile(ce_modelling_results_daly[reg,:,(i+2)],25))
                figure_2[reg,i,5]=abs(percentile(ce_modelling_results_daly[reg,:,(i+2)],75)-figure_2[reg,i,3])
                
        fig_2=plt.figure(figsize=(19.2, 9.08))
        
        label=["S1: Maximized Current Practice", "S2: S1+Facility CTC", "S3: S2+Community CTC", "S4: S3+Community CPAD", "Optimised Baseline Expenditure", "Optimised 90% Coverage"]
        markers=["D", "D", "D", "D", "s", "s"]
        size=[60,60,60,60,100,100]
        alph=0.5
        
        #Central and Eastern Europe and Central Asia
        ax1=fig_2.add_subplot(241)
        for i in range(6):
            ax1.scatter(figure_2[0,i,0], -figure_2[0,i,3], marker=markers[i], s=size[i], alpha=alph)
            ax1.errorbar(figure_2[0,i,0], -figure_2[0,i,3], xerr=[[figure_2[0,i,1]], [figure_2[0,i,2]]], yerr=[[figure_2[0,i,5]], [figure_2[0,i,4]]])            
            ax1.set(ylabel="DALYs Averted", xlabel="Additional Cost (USD)", title="Central and Eastern Europe and Central Asia")
            ax1.axvline(0, color='k', linestyle=":")

        
        #East Asia and Pacific
        ax2=fig_2.add_subplot(242)
        for i in range(6):
            ax2.scatter(figure_2[1,i,0], -figure_2[1,i,3], marker=markers[i], s=size[i], alpha=alph)
            ax2.errorbar(figure_2[1,i,0], -figure_2[1,i,3], xerr=[[figure_2[1,i,1]], [figure_2[1,i,2]]], yerr=[[figure_2[1,i,5]], [figure_2[1,i,4]]])            
            ax2.set(ylabel="DALYs Averted", xlabel="Additional Cost (USD)", title="East Asia and Pacific")
            ax2.axvline(0, color='k', linestyle=":")


        #Latin America and Caribbean
        ax3=fig_2.add_subplot(243)
        for i in range(6):
            ax3.scatter(figure_2[2,i,0], -figure_2[2,i,3], marker=markers[i], s=size[i], alpha=alph)
            ax3.errorbar(figure_2[2,i,0], -figure_2[2,i,3], xerr=[[figure_2[2,i,1]], [figure_2[2,i,2]]], yerr=[[figure_2[2,i,5]], [figure_2[2,i,4]]])            
            ax3.set(ylabel="DALYs Averted", xlabel="Additional Cost (USD)", title="Latin America and Caribbean")
            ax3.axvline(0, color='k', linestyle=":")


        #North Africa and Middle East
        ax4=fig_2.add_subplot(244)
        for i in range(6):
            ax4.scatter(figure_2[3,i,0], -figure_2[3,i,3], marker=markers[i], s=size[i], alpha=alph)
            ax4.errorbar(figure_2[3,i,0], -figure_2[3,i,3], xerr=[[figure_2[3,i,1]], [figure_2[3,i,2]]], yerr=[[figure_2[3,i,5]], [figure_2[3,i,4]]])            
            ax4.set(ylabel="DALYs Averted", xlabel="Additional Cost (USD)", title="North Africa and Middle East")
            ax4.axvline(0, color='k', linestyle=":")


        #South Asia
        ax5=fig_2.add_subplot(245)
        for i in range(6):
            ax5.scatter(figure_2[4,i,0], -figure_2[4,i,3], marker=markers[i], s=size[i], alpha=alph)
            ax5.errorbar(figure_2[4,i,0], -figure_2[4,i,3], xerr=[[figure_2[4,i,1]], [figure_2[4,i,2]]], yerr=[[figure_2[4,i,5]], [figure_2[4,i,4]]])            
            ax5.set(ylabel="DALYs Averted", xlabel="Additional Cost (USD)", title="South Asia")
            ax5.axvline(0, color='k', linestyle=":")


        #Sub-Saharan Africa
        ax6=fig_2.add_subplot(246)
        for i in range(6):
            ax6.scatter(figure_2[5,i,0], -figure_2[5,i,3], marker=markers[i], s=size[i], alpha=alph)
            ax6.errorbar(figure_2[5,i,0], -figure_2[5,i,3], xerr=[[figure_2[5,i,1]], [figure_2[5,i,2]]], yerr=[[figure_2[5,i,5]], [figure_2[5,i,4]]])            
            ax6.set(ylabel="DALYs Averted", xlabel="Additional Cost (USD)", title="Sub-Saharan Africa")
            ax6.axvline(0, color='k', linestyle=":")

        #All GBD Regions
        ax7=fig_2.add_subplot(247)
        for i in range(6):
            ax7.scatter(figure_2[6,i,0], -figure_2[6,i,3], marker=markers[i], s=size[i], alpha=alph, label=label[i])
            ax7.errorbar(figure_2[6,i,0], -figure_2[6,i,3], xerr=[[figure_2[6,i,1]], [figure_2[6,i,2]]], yerr=[[figure_2[6,i,5]], [figure_2[6,i,4]]])            
            ax7.legend(loc="center left", bbox_to_anchor=(1.1,0.5), shadow=1, fancybox=1, fontsize=12)
            ax7.set(ylabel="DALYs Averted", xlabel="Additional Cost (USD)", title="All GBD Regions")
            ax7.axvline(0, color='k', linestyle=":")
            
        plt.subplots_adjust(top=0.96, bottom=0.08, left=0.050, right=0.959, hspace=0.35, wspace=0.35)
        plt.savefig("Figure 2 Cost Effectiveness Plane.pdf", dpi=300)
        plt.show()
        
        #Suplementary Table:
        
        #Coverage (IQR), DALYs (IQR), Disease Cost (IQR), Vaccine Costs (IQR), ICER (IQR) - Needs to be calculated seperately (use ce_modelling_results)
        icer_ce_model=zeros((len(hbv_pars),6,runs))
        
        for reg in range(len(hbv_pars)):
            for run in range(runs):
                icer_ce_model[reg,0,run]=-(ce_modelling_results_cost[reg,run,1]/ce_modelling_results_daly[reg,run,1])#S1
                icer_ce_model[reg,1,run]=-(ce_modelling_results_cost[reg,run,2]/ce_modelling_results_daly[reg,run,2])#S2
                icer_ce_model[reg,2,run]=-(ce_modelling_results_cost[reg,run,3]/ce_modelling_results_daly[reg,run,3])#S3
                icer_ce_model[reg,3,run]=-(ce_modelling_results_cost[reg,run,4]/ce_modelling_results_daly[reg,run,4])#S4
                icer_ce_model[reg,4,run]=-(ce_modelling_results_cost[reg,run,6]/ce_modelling_results_daly[reg,run,6])#BO
                icer_ce_model[reg,5,run]=-(ce_modelling_results_cost[reg,run,7]/ce_modelling_results_daly[reg,run,7])#90%
        
        st5_comp=zeros((len(hbv_pars), 5, 15))
        
        baseline_results=zeros((len(hbv_pars), runs, 5))
        max_current_results=zeros((len(hbv_pars), runs, 5))
        ctc_facilities_results=zeros((len(hbv_pars), runs, 5))
        ctc_community_results=zeros((len(hbv_pars), runs, 5))
        cpad_community_results=zeros((len(hbv_pars), runs, 5))
        baseline_opt_results=zeros((len(hbv_pars), runs, 5))
        baseline_optimised_results=zeros((len(hbv_pars), runs, 5))
        optimised_90_results=zeros((len(hbv_pars), runs, 5))
        
        for reg in range(len(hbv_pars)):
            for run in range(runs):
                
                #baseline
                baseline_results[reg,run,0]=baseline_coverage[reg,run,0]                                #coverage
                baseline_results[reg,run,1]=baseline[(reg*runs)+run, len(t_steps)-1, 13]                #DALY
                baseline_results[reg,run,2]=baseline[(reg*runs)+run, len(t_steps)-1, 11]                #Disease Costs
                baseline_results[reg,run,3]=baseline[(reg*runs)+run, 0, 12]                             #Vaccine Costs
                baseline_results[reg,run,4]=0                                                           #ICER
                
                #max_current
                max_current_results[reg,run,0]=max_current_coverage[reg,run,0]
                max_current_results[reg,run,1]=max_current[(reg*runs)+run, len(t_steps)-1, 13]
                max_current_results[reg,run,2]=max_current[(reg*runs)+run, len(t_steps)-1, 11]
                max_current_results[reg,run,3]=max_current[(reg*runs)+run, 0, 12]
                max_current_results[reg,run,4]=icer_ce_model[reg,0,run]
                
                #ctc_facilities
                ctc_facilities_results[reg,run,0]=ctc_facilities_coverage[reg,run,0]
                ctc_facilities_results[reg,run,1]=ctc_facilities[(reg*runs)+run, len(t_steps)-1, 13]
                ctc_facilities_results[reg,run,2]=ctc_facilities[(reg*runs)+run, len(t_steps)-1, 11]
                ctc_facilities_results[reg,run,3]=ctc_facilities[(reg*runs)+run, 0, 12]
                ctc_facilities_results[reg,run,4]=icer_ce_model[reg,1,run]
                
                #ctc_community
                ctc_community_results[reg,run,0]=ctc_community_coverage[reg,run,0]
                ctc_community_results[reg,run,1]=ctc_community[(reg*runs)+run, len(t_steps)-1, 13]
                ctc_community_results[reg,run,2]=ctc_community[(reg*runs)+run, len(t_steps)-1, 11]
                ctc_community_results[reg,run,3]=ctc_community[(reg*runs)+run, 0, 12]
                ctc_community_results[reg,run,4]=icer_ce_model[reg,2,run]
                
                #cpad_community
                cpad_community_results[reg,run,0]=cpad_community_coverage[reg,run,0]
                cpad_community_results[reg,run,1]=cpad_community[(reg*runs)+run, len(t_steps)-1, 13]
                cpad_community_results[reg,run,2]=cpad_community[(reg*runs)+run, len(t_steps)-1, 11]
                cpad_community_results[reg,run,3]=cpad_community[(reg*runs)+run, 0, 12]
                cpad_community_results[reg,run,4]=icer_ce_model[reg,3,run]
                
                #opt_baseline
                baseline_opt_results[reg,run,0]=baseline_opt_coverage[reg,run,0]
                baseline_opt_results[reg,run,1]=baseline_opt[(reg*runs)+run, len(t_steps)-1, 13]
                baseline_opt_results[reg,run,2]=baseline_opt[(reg*runs)+run, len(t_steps)-1, 11]
                baseline_opt_results[reg,run,3]=baseline_opt[(reg*runs)+run, 0, 12]
                baseline_opt_results[reg,run,4]=0
                
                #baseline_optimisation
                baseline_optimised_results[reg,run,0]=(baseline_optimised_coverage[reg,0]+baseline_optimised_coverage[reg,2]+baseline_optimised_coverage[reg,4]+baseline_optimised_coverage[reg,6]+baseline_optimised_coverage[reg,8])/1000
                baseline_optimised_results[reg,run,1]=baseline_optimised[(reg*runs)+run, len(t_steps)-1, 13]
                baseline_optimised_results[reg,run,2]=baseline_optimised[(reg*runs)+run, len(t_steps)-1, 11]
                baseline_optimised_results[reg,run,3]=baseline_optimised[(reg*runs)+run, 0, 12]
                baseline_optimised_results[reg,run,4]=icer_ce_model[reg,4,run]
                
                #90% Optimisattion
                optimised_90_results[reg,run,0]=(optimised_90_coverage[reg,0]+optimised_90_coverage[reg,2]+optimised_90_coverage[reg,4]+optimised_90_coverage[reg,6]+optimised_90_coverage[reg,8])/1000
                optimised_90_results[reg,run,1]=optimised_90[(reg*runs)+run, len(t_steps)-1, 13]
                optimised_90_results[reg,run,2]=optimised_90[(reg*runs)+run, len(t_steps)-1, 11]
                optimised_90_results[reg,run,3]=optimised_90[(reg*runs)+run, 0, 12]
                optimised_90_results[reg,run,4]=icer_ce_model[reg,5,run]
                
                
                
        supp_table_5=zeros((len(hbv_pars),5,15))
        
        for reg in range(len(hbv_pars)):
            for i in range(1,4):
                
                #baseline
                supp_table_5[reg,0,i-1]=percentile(baseline_results[reg,:,0], i*25)
                supp_table_5[reg,0,2+i]=percentile(baseline_results[reg,:,1], i*25)
                supp_table_5[reg,0,5+i]=percentile(baseline_results[reg,:,2], i*25)
                supp_table_5[reg,0,8+i]=percentile(baseline_results[reg,:,3], i*25)
                supp_table_5[reg,0,11+i]=percentile(baseline_results[reg,:,4], i*25)
                
                #max_current
                supp_table_5[reg,1,i-1]=percentile(max_current_results[reg,:,0], i*25)
                supp_table_5[reg,1,2+i]=percentile(max_current_results[reg,:,1], i*25)
                supp_table_5[reg,1,5+i]=percentile(max_current_results[reg,:,2], i*25)
                supp_table_5[reg,1,8+i]=percentile(max_current_results[reg,:,3], i*25)
                supp_table_5[reg,1,11+i]=percentile(max_current_results[reg,:,4], i*25)
                
                #ctc_facilities
                supp_table_5[reg,2,i-1]=percentile(ctc_facilities_results[reg,:,0], i*25)
                supp_table_5[reg,2,2+i]=percentile(ctc_facilities_results[reg,:,1], i*25)
                supp_table_5[reg,2,5+i]=percentile(ctc_facilities_results[reg,:,2], i*25)
                supp_table_5[reg,2,8+i]=percentile(ctc_facilities_results[reg,:,3], i*25)
                supp_table_5[reg,2,11+i]=percentile(ctc_facilities_results[reg,:,4], i*25)
               
                #ctc_community
                supp_table_5[reg,3,i-1]=percentile(ctc_community_results[reg,:,0], i*25)
                supp_table_5[reg,3,2+i]=percentile(ctc_community_results[reg,:,1], i*25)
                supp_table_5[reg,3,5+i]=percentile(ctc_community_results[reg,:,2], i*25)
                supp_table_5[reg,3,8+i]=percentile(ctc_community_results[reg,:,3], i*25)
                supp_table_5[reg,3,11+i]=percentile(ctc_community_results[reg,:,4], i*25)
                
                #cpad_community
                supp_table_5[reg,4,i-1]=percentile(cpad_community_results[reg,:,0], i*25)
                supp_table_5[reg,4,2+i]=percentile(cpad_community_results[reg,:,1], i*25)
                supp_table_5[reg,4,5+i]=percentile(cpad_community_results[reg,:,2], i*25)
                supp_table_5[reg,4,8+i]=percentile(cpad_community_results[reg,:,3], i*25)
                supp_table_5[reg,4,11+i]=percentile(cpad_community_results[reg,:,4], i*25)

        supp_table_5=supp_table_5.reshape(len(hbv_pars)*25, 3)
        
        supp_table_5=pd.DataFrame(supp_table_5)
        
        supp_table_5.to_excel("Supplementary Table 6.xlsx")

#%%
#Figure 3: Optimization Bar Plot (GBD Regions); Supplementary Table 7
if region==1:
    if runs > 1:
        
        #Figure 3 - Regional Optimization Bar Plot 
        figure_3=zeros((len(hbv_pars),9,13))             #Include DALYs, ICER, so can export this as is and use for supplement table, and a column for sum cost
        
        for reg in range(len(hbv_pars)):
            
            #covereage
            figure_3[reg,0,0]= percentile(baseline_opt_coverage[reg,:,0],50)*100#Baseline
            figure_3[reg,1,0]= (baseline_optimised_coverage[reg,0]+baseline_optimised_coverage[reg,2]+baseline_optimised_coverage[reg,4]+baseline_optimised_coverage[reg,6]+baseline_optimised_coverage[reg,8])/10#Optimised Baseline
            figure_3[reg,3,0]= percentile(max_current_coverage[reg,:,0],50)*100#S1
            figure_3[reg,4,0]= percentile(ctc_facilities_coverage[reg,:,0],50)*100#S2
            figure_3[reg,5,0]= percentile(ctc_community_coverage[reg,:,0],50)*100#S3
            figure_3[reg,6,0]= percentile(cpad_community_coverage[reg,:,0],50)*100#S4
            figure_3[reg,8,0]=(optimised_90_coverage[reg,0]+optimised_90_coverage[reg,2]+optimised_90_coverage[reg,4]+optimised_90_coverage[reg,6]+optimised_90_coverage[reg,8])/10 #Optimised 90%
            
            #Facility Cold Chain Expenditure (Total, USD Million)
            figure_3[reg,0,2]=(baseline_opt_vcost[reg,0,0]*(hbv_pars.iloc[reg,5]/1000))/1.e6
            figure_3[reg,1,2]=(baseline_optimised_expenditure[reg,0]*(hbv_pars.iloc[reg,5]/1000))/1.e6
            figure_3[reg,3,2]=(percentile(max_current_vcost[reg,:,0],50)*(hbv_pars.iloc[reg,5]/1000))/1.e6
            figure_3[reg,4,2]=(percentile(ctc_facilities_vcost[reg,:,0],50)*(hbv_pars.iloc[reg,5]/1000))/1.e6
            figure_3[reg,5,2]=(percentile(ctc_community_vcost[reg,:,0],50)*(hbv_pars.iloc[reg,5]/1000))/1.e6
            figure_3[reg,6,2]=(percentile(cpad_community_vcost[reg,:,0],50)*(hbv_pars.iloc[reg,5]/1000))/1.e6
            figure_3[reg,8,2]=(optimised_90_expenditure[reg,0]*(hbv_pars.iloc[reg,5]/1000))/1.e6
            
            #Facility CTC Expenditure (Total)
            figure_3[reg,0,3]=(baseline_opt_vcost[reg,0,1]*(hbv_pars.iloc[reg,5]/1000))/1.e6
            figure_3[reg,1,3]=(baseline_optimised_expenditure[reg,1]*(hbv_pars.iloc[reg,5]/1000))/1.e6
            figure_3[reg,3,3]=(percentile(max_current_vcost[reg,:,1],50)*(hbv_pars.iloc[reg,5]/1000))/1.e6
            figure_3[reg,4,3]=(percentile(ctc_facilities_vcost[reg,:,1],50)*(hbv_pars.iloc[reg,5]/1000))/1.e6
            figure_3[reg,5,3]=(percentile(ctc_community_vcost[reg,:,1],50)*(hbv_pars.iloc[reg,5]/1000))/1.e6
            figure_3[reg,6,3]=(percentile(cpad_community_vcost[reg,:,1],50)*(hbv_pars.iloc[reg,5]/1000))/1.e6
            figure_3[reg,8,3]=(optimised_90_expenditure[reg,1]*(hbv_pars.iloc[reg,5]/1000))/1.e6
            
            #Community Cold Chain Expenditure (Total)
            figure_3[reg,0,4]=(baseline_opt_vcost[reg,0,2]*(hbv_pars.iloc[reg,5]/1000))/1.e6
            figure_3[reg,1,4]=(baseline_optimised_expenditure[reg,2]*(hbv_pars.iloc[reg,5]/1000))/1.e6
            figure_3[reg,3,4]=(percentile(max_current_vcost[reg,:,2],50)*(hbv_pars.iloc[reg,5]/1000))/1.e6
            figure_3[reg,4,4]=(percentile(ctc_facilities_vcost[reg,:,2],50)*(hbv_pars.iloc[reg,5]/1000))/1.e6
            figure_3[reg,5,4]=(percentile(ctc_community_vcost[reg,:,2],50)*(hbv_pars.iloc[reg,5]/1000))/1.e6
            figure_3[reg,6,4]=(percentile(cpad_community_vcost[reg,:,2],50)*(hbv_pars.iloc[reg,5]/1000))/1.e6
            figure_3[reg,8,4]=(optimised_90_expenditure[reg,2]*(hbv_pars.iloc[reg,5]/1000))/1.e6
            
            #Community CTC Expenditure (Total)
            figure_3[reg,0,5]=(baseline_opt_vcost[reg,0,3]*(hbv_pars.iloc[reg,5]/1000))/1.e6
            figure_3[reg,1,5]=(baseline_optimised_expenditure[reg,3]*(hbv_pars.iloc[reg,5]/1000))/1.e6
            figure_3[reg,3,5]=(percentile(max_current_vcost[reg,:,3],50)*(hbv_pars.iloc[reg,5]/1000))/1.e6
            figure_3[reg,4,5]=(percentile(ctc_facilities_vcost[reg,:,3],50)*(hbv_pars.iloc[reg,5]/1000))/1.e6
            figure_3[reg,5,5]=(percentile(ctc_community_vcost[reg,:,3],50)*(hbv_pars.iloc[reg,5]/1000))/1.e6
            figure_3[reg,6,5]=(percentile(cpad_community_vcost[reg,:,3],50)*(hbv_pars.iloc[reg,5]/1000))/1.e6
            figure_3[reg,8,5]=(optimised_90_expenditure[reg,3]*(hbv_pars.iloc[reg,5]/1000))/1.e6
           
            #Community CPAD Expenditure (Total)
            figure_3[reg,0,6]=(baseline_opt_vcost[reg,0,4]*(hbv_pars.iloc[reg,5]/1000))/1.e6
            figure_3[reg,1,6]=(baseline_optimised_expenditure[reg,4]*(hbv_pars.iloc[reg,5]/1000))/1.e6
            figure_3[reg,3,6]=(percentile(max_current_vcost[reg,:,4],50)*(hbv_pars.iloc[reg,5]/1000))/1.e6
            figure_3[reg,4,6]=(percentile(ctc_facilities_vcost[reg,:,4],50)*(hbv_pars.iloc[reg,5]/1000))/1.e6
            figure_3[reg,5,6]=(percentile(ctc_community_vcost[reg,:,4],50)*(hbv_pars.iloc[reg,5]/1000))/1.e6
            figure_3[reg,6,6]=(percentile(cpad_community_vcost[reg,:,4],50)*(hbv_pars.iloc[reg,5]/1000))/1.e6
            figure_3[reg,8,6]=(optimised_90_expenditure[reg,4]*(hbv_pars.iloc[reg,5]/1000))/1.e6
           
        #DALYs (IQR)
        for reg in range(len(hbv_pars)):
            for i in range(1,4):
                figure_3[reg,0,6+i]=(percentile(baseline_opt_results[reg,:,1], 25*i)*(hbv_pars.iloc[reg,5]/1000))/1.e4
                figure_3[reg,1,6+i]=(percentile(baseline_optimised_results[reg,:,1], 25*i)*(hbv_pars.iloc[reg,5]/1000))/1.e4
                figure_3[reg,3,6+i]=(percentile(max_current_results[reg,:,1], 25*i)*(hbv_pars.iloc[reg,5]/1000))/1.e4
                figure_3[reg,4,6+i]=(percentile(ctc_facilities_results[reg,:,1], 25*i)*(hbv_pars.iloc[reg,5]/1000))/1.e4
                figure_3[reg,5,6+i]=(percentile(ctc_community_results[reg,:,1], 25*i)*(hbv_pars.iloc[reg,5]/1000))/1.e4
                figure_3[reg,6,6+i]=(percentile(cpad_community_results[reg,:,1], 25*i)*(hbv_pars.iloc[reg,5]/1000))/1.e4
                figure_3[reg,8,6+i]=(percentile(optimised_90_results[reg,:,1], 25*i)*(hbv_pars.iloc[reg,5]/1000))/1.e4
                   
        #ICER (IQR)
        for reg in range(len(hbv_pars)):
            for i in range(1,4):
                figure_3[reg,0,9+i]=percentile(baseline_opt_results[reg,:,4], 25*i)
                figure_3[reg,1,9+i]=percentile(baseline_optimised_results[reg,:,4], 25*i)
                figure_3[reg,3,9+i]=percentile(max_current_results[reg,:,4], 25*i)
                figure_3[reg,4,9+i]=percentile(ctc_facilities_results[reg,:,4], 25*i)
                figure_3[reg,5,9+i]=percentile(ctc_community_results[reg,:,4], 25*i)
                figure_3[reg,6,9+i]=percentile(cpad_community_results[reg,:,4], 25*i)
                figure_3[reg,8,9+i]=percentile(optimised_90_results[reg,:,4], 25*i)
        
        for reg in range(len(hbv_pars)):
            figure_3[reg,2,0]=-200
            figure_3[reg,7,0]=-200
            figure_3[reg,:,1]=figure_3[reg,:,2]+figure_3[reg,:,3]+figure_3[reg,:,4]+figure_3[reg,:,5]+figure_3[reg,:,6]
        
        supp_table_7=figure_3.reshape(7*9, 13)
        supp_table_7=pd.DataFrame(supp_table_7)
        supp_table_7.to_excel("Supplementary Table 7.xlsx")
        
        index=arange(9)
        region_names=["Central and Eastern Europe and Central Asia", "East Asia and Pacific", "Latin America and Caribbean", "North Africa and Middle East", "South Asia", "Sub-Saharan Africa", "All GBD Regions"]
        
        plt.figure(figsize=(19.2,9.08))
        
        for reg in range(len(hbv_pars)):
            plt.subplot(ceil(len(hbv_pars)/4), ceil(len(hbv_pars)/2), reg+1)
            p1=plt.bar(index,figure_3[reg,:,2])
            p2=plt.bar(index,figure_3[reg,:,3],bottom=figure_3[reg,:,2], color='red')
            p3=plt.bar(index,figure_3[reg,:,4],bottom=figure_3[reg,:,2]+figure_3[reg,:,3], color='green')
            p4=plt.bar(index,figure_3[reg,:,5],bottom=figure_3[reg,:,2]+figure_3[reg,:,3]+figure_3[reg,:,4], color='orange')
            p5=plt.bar(index,figure_3[reg,:,6],bottom=figure_3[reg,:,2]+figure_3[reg,:,3]+figure_3[reg,:,4]+figure_3[reg,:,5], color='purple')
            plt.ylim(ymin=0.0)
            plt.xticks(index, ("Baseline", "Opt_BL", "", "S1", "S2", "S3","S4","","90%_Opt" ), rotation=290, fontsize=10)
            plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=True )
            plt.title(str(region_names[reg]))
            
            if reg==0 or reg==4:
                plt.ylabel("Vaccine Expenditure (USD Million)", fontsize=12)

        plt.legend(("Facility, Cold Chain", "Facility, CTC", "Community, Cold Chain", "Community, CTC", "Community, CPAD"), loc='best',bbox_to_anchor=(1.3,0.5), shadow=1, fancybox=1)
        
        for reg in range(len(hbv_pars)):
            plt.subplot(ceil(len(hbv_pars)/4), ceil(len(hbv_pars)/2), reg+1)
            p6=plt.twinx()
            p6=plt.scatter(index, figure_3[reg,:,0], facecolors="r", edgecolors='k')
            pl6=plt.axhline(90,color='grey', linestyle=":")
            plt.ylim(0,100)
            
            if reg==3 or reg==6:
                plt.ylabel("% Birth Dose Coverage", fontsize=12)
        
        plt.subplots_adjust(left=0.125, right=0.9, bottom=0.1, top=0.9, wspace=0.580, hspace=0.345)
        plt.savefig("Figure 3 Regional Opt Bar Plot.png", dpi=300)
        plt.show()
            
        

# #Figure 4: Optimization Bar Plot (LMICs); Supplementary Table 8
if region==0:
    if runs > 1:
        
        #################################################Figure 4##############################################
        cea,eap,lac,nae,soa,ssa=[],[],[],[],[],[]
        cea_c, eap_c, lac_c, nae_c, soa_c, ssa_c = [],[],[],[],[],[]
        
        for idx, val in enumerate(hbv_pars["worldbank_region"]):
            if val=="Central and Eastern Europe and Central Asia":
                cea.append([optimised_90_expenditure[idx,:]])
                cea_c.append(countries[idx])
            if val=="East Asia and Pacific":
                eap.append([optimised_90_expenditure[idx,:]])
                eap_c.append(countries[idx])
            if val=="Latin America and Caribbean":
                lac.append([optimised_90_expenditure[idx,:]])
                lac_c.append(countries[idx])
            if val=="North Africa and Middle East":
                nae.append([optimised_90_expenditure[idx,:]])
                nae_c.append(countries[idx])
            if val=="South Asia":
                soa.append([optimised_90_expenditure[idx,:]])
                soa_c.append(countries[idx])
            if val=="Sub-Saharan Africa":
                ssa.append([optimised_90_expenditure[idx,:]])
                ssa_c.append(countries[idx])
            
        cea=array(cea).reshape(len(cea_c),6)
        eap=array(eap).reshape(len(eap_c),6)
        lac=array(lac).reshape(len(lac_c),6)
        nae=array(nae).reshape(len(nae_c),6)
        soa=array(soa).reshape(len(soa_c),6)
        ssa=array(ssa).reshape(len(ssa_c),6)
            
        fig_4=plt.figure(figsize=(19.2,9.08))
            
        #Central and Eastern Europe and Central Asia
        index=arange(0,len(cea_c))
        
        ax1=fig_4.add_subplot(2,3,1)
        ax1.barh(index, cea[:,0], height=0.1, color="royalblue")
        ax1.barh(index, cea[:,1], height=0.1, left= cea[:,0], color='red')
        ax1.barh(index, cea[:,2], height=0.1, left= cea[:,0]+cea[:,1], color='orange')
        ax1.barh(index, cea[:,3], height=0.1, left= cea[:,0]+cea[:,1]+cea[:,2], color='green')
        ax1.barh(index, cea[:,4], height=0.1, left= cea[:,0]+cea[:,1]+cea[:,2]+cea[:,3], color='purple')
        ax1.set_yticks(index)
        ax1.set_yticklabels(cea_c[:], fontsize=12, rotation=357)
        ax1.set_ylim(ax1.get_ylim()[::-1]) #Flip the Y-Axis so region is on top
        ax1.set_ylim(-0.12, 0.14)
        ax1.set_xlim(-36, 9048)
        ax1.set(title="Central and Eastern Europe and Central Asia")
            
        #East Asia and Pacific
        index_eap=arange(0,len(eap_c))
           
        ax2=fig_4.add_subplot(2,3,2)
        ax2.barh(index_eap, eap[:,0], color="royalblue")
        ax2.barh(index_eap, eap[:,1], left= eap[:,0], color='red')
        ax2.barh(index_eap, eap[:,2], left= eap[:,0]+eap[:,1], color='orange')
        ax2.barh(index_eap, eap[:,3], left= eap[:,0]+eap[:,1]+eap[:,2], color='green')
        ax2.barh(index_eap, eap[:,4], left= eap[:,0]+eap[:,1]+eap[:,2]+eap[:,3], color='purple')
        ax2.set_yticks(index_eap)
        ax2.set_yticklabels(eap_c[:], fontsize=12, rotation=359)
        ax2.set_ylim(ax2.get_ylim()[::-1]) #Flip the Y-Axis so region is on top
        ax2.set(title="East Asia and Pacific")
        
        
            
        #Latin America and Caribbean
        ax3=fig_4.add_subplot(2,3,3)
        index_lac=arange(0,len(lac_c))
        ax3.barh(index_lac, lac[:,0],  color="royalblue")
        ax3.barh(index_lac, lac[:,1],  left= lac[:,0], color='red')
        ax3.barh(index_lac, lac[:,2],  left= lac[:,0]+lac[:,1], color='orange')
        ax3.barh(index_lac, lac[:,3],  left= lac[:,0]+lac[:,1]+lac[:,2], color='green')
        ax3.barh(index_lac, lac[:,4],  left= lac[:,0]+lac[:,1]+lac[:,2]+lac[:,3], color='purple')
        ax3.set_yticks(index_lac)
        ax3.set_yticklabels(lac_c[:], fontsize=12, rotation=359)
        ax3.set_ylim(ax3.get_ylim()[::-1]) 
        ax3.set(title="Latin America and Caribbean")
            
        #North Africa and Middle East
        ax4=fig_4.add_subplot(2,3,4)
        index_nae=arange(0,len(nae_c))
        ax4.barh(index_nae, nae[:,0], color="royalblue")
        ax4.barh(index_nae, nae[:,1], left= nae[:,0], color='red')
        ax4.barh(index_nae, nae[:,2], left= nae[:,0]+nae[:,1], color='orange')
        ax4.barh(index_nae, nae[:,3], left= nae[:,0]+nae[:,1]+nae[:,2], color='green')
        ax4.barh(index_nae, nae[:,4], left= nae[:,0]+nae[:,1]+nae[:,2]+nae[:,3], color='purple')
        ax4.set_yticks(index_nae)
        ax4.set_yticklabels(nae_c[:], fontsize=12, rotation=357)
        ax4.set_ylim(ax4.get_ylim()[::-1])
        ax4.set(title="North Africa and Middle East", xlabel="Vaccine Expenditure (USD, per 1000 births)")
            
        #South Asia
        ax5=fig_4.add_subplot(2,3,5)
        index_soa=arange(0,len(soa_c))
        ax5.barh(index_soa, soa[:,0], color="royalblue", label="Facility, Cold Chain")
        ax5.barh(index_soa, soa[:,1], left= soa[:,0], color='red', label="Facility, CTC")
        ax5.barh(index_soa, soa[:,2], left= soa[:,0]+soa[:,1], color='orange', label="Community, Cold Chain")
        ax5.barh(index_soa, soa[:,3], left= soa[:,0]+soa[:,1]+soa[:,2], color='green', label="Community, CTC")
        ax5.barh(index_soa, soa[:,4], left= soa[:,0]+soa[:,1]+soa[:,2]+soa[:,3], color='purple', label="Community, CPAD")
        ax5.legend(loc='best', shadow=1, fancybox=1, fontsize=10)
        ax5.set_yticks(index_soa)
        ax5.set_yticklabels(soa_c[:], fontsize=12, rotation=357)
        ax5.set_ylim(ax5.get_ylim()[::-1]) 
        ax5.set(title="South Asia", xlabel="Vaccine Expenditure (USD, per 1000 births)")
            
        #Sub-Saharan Africa
        ax6=fig_4.add_subplot(2,3,6)
        index_ssa=arange(0,len(ssa_c))
        ax6.barh(index_ssa, ssa[:,0], color="royalblue")
        ax6.barh(index_ssa, ssa[:,1], left= ssa[:,0], color='red')
        ax6.barh(index_ssa, ssa[:,2], left= ssa[:,0]+ssa[:,1], color='orange')
        ax6.barh(index_ssa, ssa[:,3], left= ssa[:,0]+ssa[:,1]+ssa[:,2], color='green')
        ax6.barh(index_ssa, ssa[:,4], left= ssa[:,0]+ssa[:,1]+ssa[:,2]+ssa[:,3], color='purple')
        ax6.set_yticks(index_ssa)
        ax6.set_yticklabels(ssa_c[:], fontsize=12, rotation=357)
        ax6.set_ylim(ax6.get_ylim()[::-1])
        ax6.set(title="Sub-Saharan Africa", xlabel="Vaccine Expenditure (USD, per 1000 births)")
        
        plt.subplots_adjust(top=0.957, bottom=0.065, left=0.100, right=0.985, hspace=0.222, wspace=0.986)
        plt.savefig("Figure 4 - Country Level 90% Bar Plot.pdf", dpi=300)
        plt.show()
        
        #####################################Supplementary Table 8#################################################################################
        supp_table_8_bl=zeros((len(hbv_pars), 12))      #Total Cost, Each Intervention, DALYs (IQR), ICER (0's for baseline)
        supp_table_8_90=zeros((len(hbv_pars), 12))      #Total Cost, Each Intervention, DALYS (IQR), ICER (Needs to be calculated)
        
        for reg in range(len(hbv_pars)):
            
            supp_table_8_bl[reg,0]=baseline_opt_vcost[reg,0,5]*(hbv_pars.iloc[reg,5]/1000)
            supp_table_8_bl[reg,1]=baseline_opt_vcost[reg,0,0]*(hbv_pars.iloc[reg,5]/1000)
            supp_table_8_bl[reg,2]=baseline_opt_vcost[reg,0,1]*(hbv_pars.iloc[reg,5]/1000)
            supp_table_8_bl[reg,3]=baseline_opt_vcost[reg,0,2]*(hbv_pars.iloc[reg,5]/1000)
            supp_table_8_bl[reg,4]=baseline_opt_vcost[reg,0,3]*(hbv_pars.iloc[reg,5]/1000)
            supp_table_8_bl[reg,5]=baseline_opt_vcost[reg,0,4]*(hbv_pars.iloc[reg,5]/1000)
            
            supp_table_8_90[reg,0]=optimised_90_expenditure[reg,5]*(hbv_pars.iloc[reg,5]/1000)
            supp_table_8_90[reg,1]=optimised_90_expenditure[reg,0]*(hbv_pars.iloc[reg,5]/1000)
            supp_table_8_90[reg,2]=optimised_90_expenditure[reg,1]*(hbv_pars.iloc[reg,5]/1000)
            supp_table_8_90[reg,3]=optimised_90_expenditure[reg,2]*(hbv_pars.iloc[reg,5]/1000)
            supp_table_8_90[reg,4]=optimised_90_expenditure[reg,3]*(hbv_pars.iloc[reg,5]/1000)
            supp_table_8_90[reg,5]=optimised_90_expenditure[reg,4]*(hbv_pars.iloc[reg,5]/1000)
        
        bl_opt_results=zeros((len(hbv_pars), runs, 2))      #Total Costs, DALYs
        opt_90_results=zeros((len(hbv_pars), runs, 2))
        
        for reg in range(len(hbv_pars)):
            for run in range(runs):
                bl_opt_results[reg,run,0]=baseline_opt[(runs*reg)+run, len(t_steps)-1, 12]
                bl_opt_results[reg,run,1]=baseline_opt[(runs*reg)+run, len(t_steps)-1, 13]
                
                opt_90_results[reg,run,0]=optimised_90[(runs*reg)+run, len(t_steps)-1, 12]
                opt_90_results[reg,run,1]=optimised_90[(runs*reg)+run, len(t_steps)-1, 13]
                
        icer_90=zeros((len(hbv_pars),runs,1))
        
        for reg in range(len(hbv_pars)):
            for run in range(runs):
                icer_90[reg,run,0]=-((opt_90_results[reg,run,0]-bl_opt_results[reg,run,0])/(opt_90_results[reg,run,1]-bl_opt_results[reg,run,1]))
                
        for reg in range(len(hbv_pars)):
            for i in range(1,4):
                supp_table_8_bl[reg, 5+i]=percentile(bl_opt_results[reg,:,1], 25*i)*(hbv_pars.iloc[reg,5]/1000)      #DALYs
                
                supp_table_8_90[reg,5+i]=percentile(opt_90_results[reg,:,1], 25*i)*(hbv_pars.iloc[reg,5]/1000)       #DALYs
                supp_table_8_90[reg,8+i]=percentile(icer_90[reg,:,0],25*i)       #ICER
                
        s_t8_index=[]
        for i in range(len(countries)):
            s_t8_index.append(countries[i])
            s_t8_index.append(countries[i])
            
        supp_table_8=zeros((len(hbv_pars)*2, 12))
        
        for bl in range(0,len(hbv_pars)*2,2):
            for i in range(12):
                supp_table_8[bl,i]=supp_table_8_bl[int(bl/2),i]
        
        for nt in range(1, len(hbv_pars)*2,2):
            for i in range(12):
                supp_table_8[nt,i]=supp_table_8_90[int((nt-1)/2),i]
                
         
        supp_table_8=pd.DataFrame(supp_table_8, index=s_t8_index)
        
        supp_table_8.to_excel("Supplementary Table 8.xlsx")
            
        
