# -*- coding: utf-8 -*-
"""
Codon bias - adaptation comparison
Compare calculated average bias in tissue-enriched genes to 
the relative adaptiveness values of a codon (calculated as tool for tACI) in that tissue
Methods used:
    -Visually compare with a heatmap
    -Calculate correlation coefficient between bias and adaptation 
    as they vary in tissues for each codon

Here use simple non-weighted means for codon bias in tissue groups 
(all genes count the same regardless of how many of that codon they have)

Created on Sun Jun 14 21:11:06 2020
@author: Marianne
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
from sklearn import preprocessing
from usefulfunc import (uploadFile, saveToFile, heatmap2d)
import tAI_obj as t

# =============================================================================
# Import data
# =============================================================================
raw_bias = uploadFile("codon_bias.xlsx", index_col=0)

## Load tRNA abundance data
tRNA_brain = uploadFile('tRNA brain av.xlsx', index_col=0)
tRNA_Bcell = uploadFile('tRNA B cells av.xlsx', index_col=0)
tRNA_prostate = uploadFile('tRNA prostate av.xlsx', index_col=0)
tRNA_bladder = uploadFile('tRNA bladder av.xlsx', index_col=0)
tRNA_colon = uploadFile('tRNA colon av.xlsx', index_col=0)

## Create tACI objects to get relative adaptiveness values of codons
# Reverse complement anticodons in index using indexes in brain data
acs=list(tRNA_brain.index)
codons=list() 
for ac in acs:      
  #convert anticodon to codon it recognises exactly
  codons.append(t.tAI.revcomp(ac)) 
  
tRNA_brain["codons"]  = codons
# Get average abundance only by anticodons in Series object
S_tRNA=tRNA_brain.set_index('codons')["average"]
taci_brain=t.tAI(S_tRNA,bacteria=False) #create a tAI object with tRNA abundances to calculate the weights

tRNA_bladder["codons"]  = codons
# Get average abundance only by anticodons in Series object
S_tRNA=tRNA_bladder.set_index('codons')["average"]
taci_bladder=t.tAI(S_tRNA,bacteria=False) #create a tAI object with tRNA abundances to calculate the weights

tRNA_prostate["codons"]  = codons
# Get average abundance only by anticodons in Series object
S_tRNA=tRNA_prostate.set_index('codons')["average"]
taci_prostate=t.tAI(S_tRNA,bacteria=False) #create a tAI object with tRNA abundances to calculate the weights

tRNA_colon["codons"]  = codons
# Get average abundance only by anticodons in Series object
S_tRNA=tRNA_colon.set_index('codons')["average"]
taci_colon=t.tAI(S_tRNA,bacteria=False) #create a tAI object with tRNA abundances to calculate the weights

##Recompile gene groups and only keep the ones for which sequence available
#EV1Genes=pd.read_excel('Table_EV1.xlsx', sheet_name = 'C. Genes') # Import sheet "C. Genes" from excel file 
#EV1Genes=EV1Genes.set_index('Gene ID')
#
## Look at tissue-enriched genes
#genes=pd.concat([raw_bias,EV1Genes[['Tissue enriched','Group enriched','Tissue enhanced','Classification']]],
#                join='inner',axis=1) #dataframe for results using bladder tRNA levels
#groups=genes.groupby(['Classification'])
#enriched = genes.groupby(['Classification']).get_group('Tissue enriched')
#enriched=enriched.drop(['Group enriched', 'Tissue enhanced'],axis=1)

scaled_bias = raw_bias.drop(['ATG','TAA','TAG','TGA'],axis='columns', inplace=False) \
                        .loc[enriched.index,:].values
#scaled_data=preprocessing.StandardScaler().fit_transform(scaled_bias) #scaled_bias is a numpy array 
#scaled_bias=pd.DataFrame(scaled_data, 
#                       index=raw_bias.loc[enriched.index,:].index, #only keep tissue-enriched genes
#                       columns=raw_bias.drop(['ATG','TAA','TAG','TGA'],axis='columns', inplace=False).columns)
scaled_bias.dropna(inplace=True)
brain_bias=scaled_bias[enriched['Tissue enriched']=='Brain']
bladder_bias=scaled_bias[enriched['Tissue enriched']=='Urinary bladder']
prostate_bias=scaled_bias[enriched['Tissue enriched']=='Prostate']
colon_bias=scaled_bias[enriched['Tissue enriched']=='Colon']

#av_tissue_bias=np.concatenate((brain_bias.mean().to_numpy().reshape(-1,1)
#                ,bladder_bias.mean().to_numpy().reshape(-1,1)
#                ,prostate_bias.mean().to_numpy().reshape(-1,1)
#                ,colon_bias.mean().to_numpy().reshape(-1,1))
#                ,axis=1)
#
#codon_adapt=np.concatenate((taci_brain.weights.to_numpy().reshape(-1,1)
#                ,taci_bladder.weights.to_numpy().reshape(-1,1)
#                ,taci_prostate.weights.to_numpy().reshape(-1,1)
#                ,taci_colon.weights.to_numpy().reshape(-1,1))
#                ,axis=1)

##need to change codon adaptation and bias to same indexing!
#Also normalise codon adaptation to same scale as bias for easier comparison

df=pd.DataFrame({'brain bias': brain_bias.mean(), 'bladder bias':bladder_bias.mean(),'prostate bias':prostate_bias.mean(),'colon bias':colon_bias.mean(),
              'brain adaptation':taci_brain.weights,'bladder adaptation':taci_bladder.weights,
              'prostate adaptation':taci_prostate.weights,
              'colon adaptation':taci_colon.weights}, 
               index=taci_brain.weights.index)
np_adapt=df.loc[:,'brain adaptation':].to_numpy()
#np_adapt = preprocessing.StandardScaler().fit_transform(np_adapt)
np_bias=df.loc[:,:'colon bias'].to_numpy()

heatfig,heatax=heatmap2d(np.concatenate((np_bias,np_adapt),axis=1))
heatfig,heatax=plt.subplots(1,2)
a=heatax[0].imshow(np_bias, cmap='coolwarm')
b=heatax[1].imshow(np_adapt, cmap='coolwarm')
heatfig.colorbar(a,ax=heatax[0])
heatfig.colorbar(b,ax=heatax[1])
heatax[0].set_title("Average codon bias of tissue-specific genes")
heatax[1].set_title("Relative codon adaptiveness")
heatfig.suptitle("Tissue order: Brain, bladder, prostate, colon" )

y_pos=np.arange(len(df.index))
heatax[0].set_xlabel("Tissue samples")
heatax[0].set_ylabel("Codon")
heatax[0].set_yticks(y_pos)
heatax[0].set_yticklabels(df.index)
heatax[1].set_xlabel("Tissue samples")
heatax[1].set_ylabel("Codon")
heatax[1].set_yticks(y_pos)
heatax[1].set_yticklabels(df.index)
#plt.show()

preprocessing.StandardScaler().fit_transform(codon_adapt)

# =============================================================================
# Try to look at correlation between codon bias and relative adaptiveness value?
# =============================================================================
##Remember we have 4 tissues here
#Start with correlation for each codon
r=pd.Series(0.0, index=df.index) #correlation coefficient for each row
xx=df.iloc[:,:4].copy()
yy=df.iloc[:,4:].copy()
xy=pd.DataFrame(0.0, index=df.index, 
                columns=['brain', 'bladder', 'prostate', 'colon'])

for codon in df.index:    
    mean_bias = df.loc[codon].iloc[:4].mean()
    mean_adapt = df.loc[codon].iloc[4:].mean()
    xx.loc[codon]=df.loc[codon].iloc[:4]-mean_bias
    yy.loc[codon]=df.loc[codon].iloc[4:]-mean_adapt
    xy.loc[codon]=xx.loc[codon]*yy.loc[codon]
    r[codon]=xy.loc[codon].sum()/np.sqrt((xx.loc[codon]**2).sum()*(yy.loc[codon]**2).sum())



        
        
