# -*- coding: utf-8 -*-
"""
"averagetRNAAbundance.py"
Calculate average tRNA abundance per anticodon type from multiple samples.
Plot them.
Save results to files.

Also do statistically analyse variations in the tRNA pool

Input file: tRNAabNormal_cells.xlsx from the output of sumtRNAAbundances.py,
which has been modified so each tissue type is in a different Excel sheet.
Output: 5 files corresponding to the data in each sheet plus 
the average, standard deviation and percent error corresponding to each anticodon type.

Created on Tue Mar  3 17:22:23 2020
@author: Marianne Buat
"""
import pandas as pd
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt

def saveToFile(df_in,output_file):
    ##Save updated output to file
    writer = pd.ExcelWriter(output_file)
    df_in.to_excel(writer)
    writer.save()
    print('Output written successfully to Excel File.')
    return

def averagetRNAAbundance(df_tRNAs):
    ## Calculate average tRNA values and standard deviations
    acs=list(df_tRNAs.index) #extract anticodons from index
    av=list()
    std=list()
    percent_err=list()
    std_err=list()
    for ac in acs:    
        #Calculate mean, standard deviation and percent error of tRNA abundance in given cell type
        #Here start with normal bladder cells
        av.append(df_tRNAs.loc[ac].mean())
        std.append(df_tRNAs.loc[ac].std())
        percent_err.append(std[-1]/av[-1])
        std_err.append(std[-1]/len(df_tRNAs.columns))
    av_tRNA = pd.DataFrame({"average":av,"std dev":std, "percent error":percent_err}, index = acs)
    return av_tRNA

def heatmap2d(arr: np.ndarray):
    fig,ax=plt.subplots()
    a=ax.imshow(arr, cmap='coolwarm')
    fig.colorbar(a,ax=ax)
    plt.show()
    return fig,ax

# Functions for ANOVA
def SumSqTr(treatment_list, N = None, k=None, nis=None):
    """
    Calculate the Sum of Squares for treatment for ANOVA
    treatment_list is a list of numpy arrays containing the data for each treatment,
    N the total number of samples, k the number of different treatments and 
    nis the number of samples per treatment
    If any of the last 3 parameters are not provided, they are automatically calculated from the data
    
    """
    if k==None: #number of different treatments
        k=len(treatment_list)
    if nis==None: #number of samples in each treatment
        nis=[]
        for t in treatment_list:
            nis.append(len(t))
    if N==None or sum(nis)!=N: #total number of samples
        N=sum(nis)
    
    tr_means=[np.mean(t) for t in treatment_list] #means of each treatment
    grand_mean=sum([nis[tr]*tr_means[tr] for tr in range(k)])/N #grand mean, 
            #equivalent to weighted average of treatment means 
            #or to average of all samples together regardless of treatment 
    SSTr=0
    for tr in range(0,k):
        SSTr += nis[tr]*(tr_means[tr]-grand_mean)**2
    return SSTr
        
    
def SumSqEr(treatment_list, N = None, k=None, nis=None):
    """Calculate the Sum of Squares due to errors for ANOVA
    treatment_list is a list of numpy arrays containing the data for each treatment,
    N the total number of samples, k the number of different treatments and 
    nis the number of samples per treatment
    If any of the last 3 parameters are not provided, they are automatically calculated from the data
    
    """
    if k==None:
        k=len(treatment_list)
    if nis==None:
        nis=[]
        for t in treatment_list:
            nis.append(len(t))
    if N==None or sum(nis)!=N:
        N=sum(nis)
    tr_means=[np.mean(t) for t in treatment_list] #means of each treatment
    SSE=0
    for tr in range(0,k):
        for sample in treatment_list[tr]:
            SSE += (sample-tr_means[tr])**2
    return SSE

def SumSqTot(treatment_list, N = None, k=None, nis=None):
    """
    Calculate the total sum of squares for ANOVA
    treatment_list is a list of numpy arrays containing the data for each treatment,
    N the total number of samples, k the number of different treatments and 
    nis the number of samples per treatment
    If any of the last 3 parameters are not provided, they are automatically calculated from the data
    
    """
    if nis==None:
        nis=[]
        for t in treatment_list:
            nis.append(len(t))
    if N==None or sum(nis)!=N:
        N=sum(nis)
    grand_mean=sum([sum(t) for t in treatment_list])/N
    return sum([sum([(sample-grand_mean)**2 for sample in t]) for t in treatment_list])
        
def Eta_Squared(F_stat,k,N):
    """
    Another way of calculating eta squared directly from the f statistic
    and the degrees of freedom (from the number of samples and treatments used)
    """
    return F_stat*(k-1)/(F_stat*(k-1)+N-k)
###############################################################################
#                               Main running code                             #
###############################################################################

#Import tRNA abundances for all samples in each tissues
tRNA_brain = pd.read_excel("tRNAs/SummedtRNAab - with aa, iMet SeC included.xlsx", sheet_name = 'Normal brain',index_col=0)
tRNA_Bcell = pd.read_excel("tRNAs/SummedtRNAab - with aa, iMet SeC included.xlsx", sheet_name = 'Normal B-cells',index_col=0)
tRNA_prostate = pd.read_excel("tRNAs/SummedtRNAab - with aa, iMet SeC included.xlsx", sheet_name = 'Normal prostate',index_col=0)
tRNA_bladder = pd.read_excel("tRNAs/SummedtRNAab - with aa, iMet SeC included.xlsx", sheet_name = 'Normal bladder',index_col=0)
tRNA_colon = pd.read_excel("tRNAs/SummedtRNAab - with aa, iMet SeC included.xlsx", sheet_name = 'colon Adjacent Normal Mucosa',index_col=0)

#----------------------------------------------------------------------------#
## Make heatmap of relative expressions of tRNA data (raw)
#----------------------------------------------------------------------------#
tRNA_brain2 = tRNA_brain.sort_values(by='amino acid').drop('amino acid', axis=1)
tRNA_Bcell2 = tRNA_Bcell.sort_values(by='amino acid').drop('amino acid', axis=1)
tRNA_prostate2 = tRNA_prostate.sort_values(by='amino acid').drop('amino acid', axis=1)
tRNA_bladder2 = tRNA_bladder.sort_values(by='amino acid').drop('amino acid', axis=1)
tRNA_colon2 = tRNA_colon.sort_values(by='amino acid').drop('amino acid', axis=1)

#Convert pandas DataFrame to numpy arrays
np_bladder=tRNA_bladder2.to_numpy()
np_Bcell=tRNA_Bcell2.to_numpy()
np_brain=tRNA_brain2.to_numpy()
np_colon=tRNA_colon2.to_numpy()
np_prostate=tRNA_prostate2.to_numpy()

bladder_range=np.arange(0,np.shape(np_bladder)[1])
bladder_range=np.arange(0,np.shape(np_bladder)[1])
bladder_range=np.arange(0,np.shape(np_bladder)[1])
bladder_range=np.arange(0,np.shape(np_bladder)[1])
bladder_range=np.arange(0,np.shape(np_bladder)[1])

#Concatenate numpy version of arrays into single array, 
#conserve rows as corresponding to anticodons 
np_tRNA=np.concatenate((np_brain,np_Bcell,np_prostate,np_bladder,np_colon),axis=1) 
#For each tRNA anticodon, calculate average abundance in all tissues
np_av=np.mean(np_tRNA,axis=1)
#Normalise by average tRNA abundance for fold change
#np_tRNA_log=np.log10(np_tRNA/np_av[:,None]) #each row of np_tRNA divided by corresponding element in np_av
np_tRNA_log=np.log2(np_tRNA/(np_tRNA_av.mean(axis=1)[:,None]))
heatfig,heatax=heatmap2d(np_tRNA_log)
#heatax.set_title("tRNA level fold-change in each tissue relative to average(log2 scale)")
y_pos=np.arange(len(tRNA_bladder.index))
heatax.set_yticks(y_pos)
heatax.set_yticklabels(tRNA_brain.index)
heatax.set_ylabel("Anticodon")

heatax.set_xlabel("Tissue samples")
#Get number of samples for each tissue type
brain_samples=np.shape(np_brain)[1]
Bcell_samples=np.shape(np_Bcell)[1]
prostate_samples=np.shape(np_prostate)[1]
bladder_samples=np.shape(np_bladder)[1]
colon_samples=np.shape(np_colon)[1]
#tick positions corresponding to separation between tissue types on graph
edges=[0,
       brain_samples,
       brain_samples+Bcell_samples,
       brain_samples+Bcell_samples+prostate_samples,
       brain_samples+Bcell_samples+prostate_samples+bladder_samples,
       brain_samples+Bcell_samples+prostate_samples+bladder_samples+colon_samples]
x_labelpos=[np.median(np.arange(0,edges[1])),
            np.median(np.arange(edges[1],edges[2])),
            np.median(np.arange(edges[2],edges[3])),
            np.median(np.arange(edges[3],edges[4])),
            np.median(np.arange(edges[4],edges[5]))]     
edges=[e-0.5 for e in edges]
x_pos=[edges[0],x_labelpos[0],
       edges[1],x_labelpos[1],
       edges[2],x_labelpos[2],
       edges[3],x_labelpos[3],
       edges[4],x_labelpos[4],
       edges[5]]
x_labels=[None,'Brain',None,'B cell',None,'Prostate',None,'Bladder',None,'Colon',None]
#Set tick position
heatax.set_xticks(x_pos)
heatax.set_xticklabels(x_labels, rotation=90)

##################### Repeat using brain as reference for fold-change
'''
np_tRNA=np.concatenate((np_brain,np_Bcell,np_prostate,np_bladder,np_colon),axis=1) 
np_av=np.mean(np_brain,axis=1) #average abundance in brain for each anticodon

np_tRNA=np.log2(np_tRNA/np_av[:,None]) #each row of np_tRNA divided by corresponding element in np_av

heatfig,heatax=heatmap2d(np_tRNA)
heatax.set_title("tRNA level fold-change in each tissue relative to average(log2 scale)")
y_pos=np.arange(len(tRNA_bladder.index))
heatax.set_yticks(y_pos)
heatax.set_yticklabels(tRNA_brain.index)

#Get number of samples for each tissue type
brain_samples=np.shape(np_brain)[1]
Bcell_samples=np.shape(np_Bcell)[1]
prostate_samples=np.shape(np_prostate)[1]
bladder_samples=np.shape(np_bladder)[1]
colon_samples=np.shape(np_colon)[1]
#tick positions corresponding to separation between tissue types on graph
edges=[0,
       brain_samples,
       brain_samples+Bcell_samples,
       brain_samples+Bcell_samples+prostate_samples,
       brain_samples+Bcell_samples+prostate_samples+bladder_samples,
       brain_samples+Bcell_samples+prostate_samples+bladder_samples+colon_samples]
x_labelpos=[np.median(np.arange(0,edges[1])),
            np.median(np.arange(edges[1],edges[2])),
            np.median(np.arange(edges[2],edges[3])),
            np.median(np.arange(edges[3],edges[4])),
            np.median(np.arange(edges[4],edges[5]))]     
edges=[e-0.5 for e in edges]
x_pos=[edges[0],x_labelpos[0],
       edges[1],x_labelpos[1],
       edges[2],x_labelpos[2],
       edges[3],x_labelpos[3],
       edges[4],x_labelpos[4],
       edges[5]]
x_labels=[None,'Brain',None,'B cell',None,'Prostate',None,'Bladder',None,'Colon',None]
#Set tick position
heatax.set_xticks(x_pos)
heatax.set_xticklabels(x_labels, rotation=90)
'''
#---------------------------------------------------------------------------#
## Find averages and add new columns to dataframe
#---------------------------------------------------------------------------#
av_tRNA=averagetRNAAbundance(tRNA_brain2)
tRNA_brain = pd.concat([av_tRNA,tRNA_brain2], axis=1)
av_tRNA=averagetRNAAbundance(tRNA_Bcell2)
tRNA_Bcell = pd.concat([av_tRNA,tRNA_Bcell2], axis=1)
av_tRNA=averagetRNAAbundance(tRNA_prostate2)
tRNA_prostate = pd.concat([av_tRNA,tRNA_prostate2], axis=1)
av_tRNA=averagetRNAAbundance(tRNA_bladder2)
tRNA_bladder = pd.concat([av_tRNA,tRNA_bladder2], axis=1)
av_tRNA=averagetRNAAbundance(tRNA_colon2)
tRNA_colon = pd.concat([av_tRNA,tRNA_colon2], axis=1)

#Plot average tRNA abundances and their std dev on a graph
tRNA_list=[tRNA_brain,tRNA_Bcell,tRNA_prostate,tRNA_bladder,tRNA_colon]
x_pos = np.arange(len(tRNA_brain.index))
fig, ax = plt.subplots(nrows=2, ncols=3)
val=0
for row in ax:
    for col in row:
        if val>=len(tRNA_list):
            break
        col.bar(x_pos,tRNA_list[val]['average'],yerr=tRNA_list[val]['std dev'])
        col.set_ylabel('tRNA abundance')
        col.set_xticks(x_pos)
        col.set_xticklabels(tRNA_list[val].index,rotation=90)
        col.set_title('({}) tRNA abundance per anticodon type'.format(val+1))
        col.yaxis.grid(True)
        val+=1
plt.suptitle("Order: Brain, B cells, Prostate, Bladder, Colon")
plt.show()

## Make a similar heatmap but for averages
np_tRNA_av=np.concatenate((tRNA_brain['average'].to_numpy(),
                          tRNA_Bcell['average'].to_numpy(),
                          tRNA_prostate['average'].to_numpy(),
                          tRNA_bladder['average'].to_numpy(),
                          tRNA_colon['average'].to_numpy()))\
                           .reshape(-1,len(tRNA_brain)).transpose()

np_av=np.mean(np_tRNA_av,axis=1)
np_tRNA_av_log=np.log10(np_tRNA_av/np_av[:,None]) #each row of np_tRNA divided by corresponding element in np_av, take log base 2
heatfig_av,heatax_av=heatmap2d(np_tRNA_av_log)
#heatax_av.set_title("tRNA level fold-change in each tissue relative to average(log10 scale)")
y_pos=np.arange(len(tRNA_bladder.index))
heatax_av.set_yticks(y_pos)
heatax_av.set_yticklabels(tRNA_brain.index)
heatax_av.set_xticks(np.arange(np.shape(np_tRNA_av)[1]))
heatax_av.set_xticklabels(['Brain', 'B cells', 'Prostate', 'Bladder', 'Colon'],rotation=90)
#plt.savefig('tRNA average Abundances in tissues.png')


#-----------------------------------------------------------------------------#
#                       statistical analysis - ANOVA
#-----------------------------------------------------------------------------#

### A) 
#Does a tissue consistently exhibit different abundance for all anticodons? 
# Method 1: compare raw averages for all anticodons in a group for each tissue 
av_fold_change=np.mean(np_tRNA_av,axis=0)
ANOVA_av_tissues=stats.f_oneway(np_tRNA_av[:,0],np_tRNA_av[:,1],np_tRNA_av[:,2],np_tRNA_av[:,3],np_tRNA_av[:,4])
#F=1.43
#p value=0.225
# 
#evaluate magnitude
tissue_list=[np_tRNA_av[:,0],np_tRNA_av[:,1],np_tRNA_av[:,2],np_tRNA_av[:,3],np_tRNA_av[:,4]]
SSTr_av=SumSqTr(tissue_list)
SST_av=SumSqTot(tissue_list)
Eta_sq_av=SSTr_av/SST_av
#Here Eta squared is 0.034, small to medium effect

#Method 2:
tRNAs_list=[np_brain, np_Bcell,np_prostate,np_bladder,np_colon]
k=len(tRNAs_list) #number of "treatments" for ANOVA (tissues)
N=(brain_samples+Bcell_samples+prostate_samples+bladder_samples+colon_samples)*len(tRNA_brain) #total number of data points

ANOVA_all=stats.f_oneway(tRNAs_list[0].flatten(),tRNAs_list[1].flatten(),tRNAs_list[2].flatten(),tRNAs_list[3].flatten(),tRNAs_list[4].flatten())
eta_all=Eta_Squared(ANOVA_all[0], k, N)
# Method 3: Compare difference from mean tRNA level of anticodon for all anticodons and samples from tissue
tRNAs_list=[np_brain, np_Bcell,np_prostate,np_bladder,np_colon]

tRNAs_diff=[]
for tissue_tRNA in tRNAs_list:
    
    d=np.empty_like(tissue_tRNA)
    #for each anticodon type, remove average abundance 
    for i in range(np.shape(tissue_tRNA)[1]): #loop through each column (sample)
        d[:,i]=np_tRNA[:,i]#-np_av 
    #Linearise so all measurements from tissue part of a single group
    tRNAs_diff.append(d.flatten())
k=len(tRNAs_list) #number of "treatments" for ANOVA (tissues)
N=(brain_samples+Bcell_samples+prostate_samples+bladder_samples+colon_samples)*len(tRNA_brain) #total number of data points


f_diff,p_diff=stats.f_oneway(tRNAs_diff[0],tRNAs_diff[1],tRNAs_diff[2],tRNAs_diff[3],tRNAs_diff[4])
eta_sq_diff=Eta_Squared(f_diff, k=len(tRNA_list), N=N)
print("Eta squared for difference from mean, all anticodons together: {}".format(eta_sq_diff))

# Method 4: Compare fold-changes
tRNAs_FC=[]
for tissue_tRNA in tRNAs_list:
    
    #for each anticodon type, divide by average tRNA abundance
    fc=tissue_tRNA/np_av[:,None]
    #Linearise so all measurements from tissue part of a single group
    tRNAs_FC.append(fc.flatten())
k=len(tRNAs_list) #number of "treatments" for ANOVA (tissues)
N=(brain_samples+Bcell_samples+prostate_samples+bladder_samples+colon_samples)*len(tRNA_brain) #total number of data points
f_FC,p_FC=stats.f_oneway(tRNAs_FC[0],tRNAs_FC[1],tRNAs_FC[2],tRNAs_FC[3],tRNAs_FC[4])
eta_sq_FC=Eta_Squared(f_FC, k=len(tRNA_list), N=N)
print("Eta squared for fold-change with respect to mean, all anticodons together: {}".format(eta_sq_FC))



###  B)
#Normalise by tissue average and compare values for each anticodon
normal_bladder=np_bladder/(np.mean(np_bladder,axis=1)[:,None])
normal_brain=np_brain/(np.mean(np_brain,axis=1)[:,None])
normal_Bcell=np_Bcell/(np.mean(np_Bcell,axis=1)[:,None])
normal_colon=np_colon/(np.mean(np_colon,axis=1)[:,None])
normal_prostate=np_prostate/(np.mean(np_prostate,axis=1)[:,None])

tRNAs_normal=[normal_brain, normal_Bcell,normal_prostate,normal_bladder,normal_colon]
nis=[brain_samples,
       Bcell_samples,
       prostate_samples,
       bladder_samples,
       colon_samples] #Number of samples in each tissue
k=len(tRNAs_normal) #number of "treatments" for ANOVA (tissues)
N=brain_samples+Bcell_samples+prostate_samples+bladder_samples+colon_samples #total number of tissue samples


ac_Fstat=[]
ac_pval=[]
SSTs=[] #Total sum of squares
SSTrs=[] #Sum of square for treatment
MSTrs=[] #Mean of squares for treatment
SSEs=[] #Sum of squares due to errors
MSEs=[] #Mean of squares for treatment
SSTots=[]
Etas_squared=[]
for ac in range(len(tRNA_brain.index)):    
    f=stats.f_oneway(normal_bladder[ac,:],normal_brain[ac,:],normal_Bcell[ac,:],normal_prostate[ac,:],normal_colon[ac,:])
    ac_Fstat.append(f[0])
    ac_pval.append(f[1])
    
    treatment_list=[tissue[ac,:] for tissue in tRNAs_normal]
    sstr=SumSqTr(treatment_list)
    SSTrs.append(sstr)
    MSTrs.append(sstr/(k-1))
    sse=SumSqEr(treatment_list)
    SSEs.append(sse)
    MSEs.append(sse/(N-k))
    sstot=sstr+sse #total sum of squares, use fundamental ANOVA identity    
    SSTots.append(sstot)
    eta_sq=sstr/sstot
    Etas_squared.append(eta_sq)
    #print("{}: eta_sq={}".format(list(tRNA_brain.index)[ac],eta_sq))


'''
#Save average tRNA abundances to file
saveToFile(tRNA_brain,'tRNA brain av.xlsx')
saveToFile(tRNA_colon,'tRNA colon av.xlsx')
saveToFile(tRNA_bladder,'tRNA bladder av.xlsx')
saveToFile(tRNA_prostate,'tRNA prostate av.xlsx')
saveToFile(tRNA_Bcell,'tRNA B cells av.xlsx')
'''
mean_colon_FC=(np_colon/(np_tRNA_av.mean(axis=1)[:,None])).mean()
mean_brain_FC=(np_brain/(np_tRNA_av.mean(axis=1)[:,None])).mean()
mean_prostate_FC=(np_prostate/(np_tRNA_av.mean(axis=1)[:,None])).mean()
mean_Bcell_FC=(np_Bcell/(np_tRNA_av.mean(axis=1)[:,None])).mean()
mean_bladder_FC=(np_bladder/(np_tRNA_av.mean(axis=1)[:,None])).mean()