# -*- coding: utf-8 -*-
"""
"averagetRNAAbundance.py"
Calculate average tRNA abundance per anticodon type from multiple samples.
Plot them.
Save results to files.

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
    for ac in acs:    
        #Calculate mean, standard deviation and percent error of tRNA abundance in given cell type
        #Here start with normal bladder cells
        av.append(df_tRNAs.loc[ac].mean())
        std.append(df_tRNAs.loc[ac].std())
        percent_err.append(std[-1]/av[-1])
        
    av_tRNA = pd.DataFrame({"average":av,"std dev":std, "percent error":percent_err}, index = acs)
    return av_tRNA

def heatmap2d(arr: np.ndarray):
    fig,ax=plt.subplots()
    a=ax.imshow(arr, cmap='coolwarm')
    fig.colorbar(a,ax=ax)
    plt.show()
    return fig,ax



#Import tRNA abundances for all samples in each tissues
tRNA_brain = pd.read_excel('tRNAabNormal_cells.xlsx', sheet_name = 'Normal brain',index_col=0)
tRNA_Bcell = pd.read_excel('tRNAabNormal_cells.xlsx', sheet_name = 'Normal B-cells',index_col=0)
tRNA_prostate = pd.read_excel('tRNAabNormal_cells.xlsx', sheet_name = 'Normal prostate',index_col=0)
tRNA_bladder = pd.read_excel('tRNAabNormal_cells.xlsx', sheet_name = 'Normal bladder cells',index_col=0)
tRNA_colon = pd.read_excel('tRNAabNormal_cells.xlsx', sheet_name = 'colon Adjacent Normal Mucosa',index_col=0)

#----------------------------------------------------------------------------#
## Make heatmap of relative expressions of tRNA data (raw)
#----------------------------------------------------------------------------#

#Convert pandas DataFrame to numpy arrays
n_bladder=tRNA_bladder.to_numpy()
n_Bcell=tRNA_Bcell.to_numpy()
n_brain=tRNA_brain.to_numpy()
n_colon=tRNA_colon.to_numpy()
n_prostate=tRNA_prostate.to_numpy()

bladder_range=np.arange(0,np.shape(n_bladder)[1])
bladder_range=np.arange(0,np.shape(n_bladder)[1])
bladder_range=np.arange(0,np.shape(n_bladder)[1])
bladder_range=np.arange(0,np.shape(n_bladder)[1])
bladder_range=np.arange(0,np.shape(n_bladder)[1])

#Concatenate numpy version of arrays into single array, 
#conserve rows as corresponding to anticodons 
n_tRNA=np.concatenate((n_brain,n_Bcell,n_prostate,n_bladder,n_colon),axis=1) 
#For each tRNA anticodon, calculate average abundance in all tissues
n_av=np.mean(n_tRNA,axis=1)
#Normalise by average tRNA abundance for fold change
n_tRNA=np.log2(n_tRNA/n_av[:,None]) #each row of n_tRNA divided by corresponding element in n_av

heatfig,heatax=heatmap2d(n_tRNA)
heatax.set_title("tRNA level fold-change in each tissue relative to average(log2 scale)")
y_pos=np.arange(len(tRNA_bladder.index))
heatax.set_yticks(y_pos)
heatax.set_yticklabels(tRNA_brain.index)

#Get number of samples for each tissue type
brain_samples=np.shape(n_brain)[1]
Bcell_samples=np.shape(n_Bcell)[1]
prostate_samples=np.shape(n_prostate)[1]
bladder_samples=np.shape(n_bladder)[1]
colon_samples=np.shape(n_colon)[1]
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
n_tRNA=np.concatenate((n_brain,n_Bcell,n_prostate,n_bladder,n_colon),axis=1) 
n_av=np.mean(n_brain,axis=1) #average abundance in brain for each anticodon

n_tRNA=np.log2(n_tRNA/n_av[:,None]) #each row of n_tRNA divided by corresponding element in n_av

heatfig,heatax=heatmap2d(n_tRNA)
heatax.set_title("tRNA level fold-change in each tissue relative to average(log2 scale)")
y_pos=np.arange(len(tRNA_bladder.index))
heatax.set_yticks(y_pos)
heatax.set_yticklabels(tRNA_brain.index)

#Get number of samples for each tissue type
brain_samples=np.shape(n_brain)[1]
Bcell_samples=np.shape(n_Bcell)[1]
prostate_samples=np.shape(n_prostate)[1]
bladder_samples=np.shape(n_bladder)[1]
colon_samples=np.shape(n_colon)[1]
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

#---------------------------------------------------------------------------#
## Find averages and add new columns to dataframe
#---------------------------------------------------------------------------#
av_tRNA=averagetRNAAbundance(tRNA_brain)
tRNA_brain = pd.concat([av_tRNA,tRNA_brain], axis=1)
av_tRNA=averagetRNAAbundance(tRNA_Bcell)
tRNA_Bcell = pd.concat([av_tRNA,tRNA_Bcell], axis=1)
av_tRNA=averagetRNAAbundance(tRNA_prostate)
tRNA_prostate = pd.concat([av_tRNA,tRNA_prostate], axis=1)
av_tRNA=averagetRNAAbundance(tRNA_bladder)
tRNA_bladder = pd.concat([av_tRNA,tRNA_bladder], axis=1)
av_tRNA=averagetRNAAbundance(tRNA_colon)
tRNA_colon = pd.concat([av_tRNA,tRNA_colon], axis=1)

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
n_tRNA_av=np.concatenate((tRNA_brain['average'].to_numpy(),
                          tRNA_Bcell['average'].to_numpy(),
                          tRNA_prostate['average'].to_numpy(),
                          tRNA_bladder['average'].to_numpy(),
                          tRNA_colon['average'].to_numpy())).reshape(-1,31).transpose()

n_av=np.mean(n_tRNA_av,axis=1)
n_tRNA_av=np.log10(n_tRNA_av/n_av[:,None]) #each row of n_tRNA divided by corresponding element in n_av, take log base 2
heatfig_av,heatax_av=heatmap2d(n_tRNA_av)
heatax_av.set_title("tRNA level fold-change in each tissue relative to average(log10 scale)")
y_pos=np.arange(len(tRNA_bladder.index))
heatax_av.set_yticks(y_pos)
heatax_av.set_yticklabels(tRNA_brain.index)
heatax_av.set_xticks(np.arange(np.shape(n_tRNA_av)[1]))
heatax_av.set_xticklabels(['Brain', 'B cells', 'Prostate', 'Bladder', 'Colon'],rotation=90)
#plt.savefig('tRNA average Abundances in tissues.png')

av_fold_change=np.mean(n_tRNA_av,axis=0)
stats.f_oneway(n_tRNA_av[:,0],n_tRNA_av[:,1],n_tRNA_av[:,2],n_tRNA_av[:,3],n_tRNA_av[:,4])
'''
#Save average tRNA abundances to file
saveToFile(tRNA_brain,'tRNA brain av.xlsx')
saveToFile(tRNA_colon,'tRNA colon av.xlsx')
saveToFile(tRNA_bladder,'tRNA bladder av.xlsx')
saveToFile(tRNA_prostate,'tRNA prostate av.xlsx')
saveToFile(tRNA_Bcell,'tRNA B cells av.xlsx')
'''