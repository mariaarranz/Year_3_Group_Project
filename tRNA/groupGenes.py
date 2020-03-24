# -*- coding: utf-8 -*-
"""
Year 3 Group Project: Codons

Separate genes according to their enrichment status

Input: Table_EV1.xlsx, sheet 'C.Genes' from Wang et al, 2019 
        (containing protein abundance data)
    
Created on Wed Mar  4 14:12:25 2020

@author: Marianne
"""
import pandas as pd
import numpy as np
import tAI_obj as t #import all objects and functions defined in that file
import scipy.stats as stats
from matplotlib import pyplot as plt

tACI_values = t.uploadFile('tACI_values.xlsx', index_col=0) #Reupload file with tACI values
#Could also be using file with sequences
#seq=uploadFile('Table_EV4_excel.xlsx',index_col=1)


#Recompile gene groups and only keep the ones for which sequence available
EV1Genes=pd.read_excel('Table_EV1.xlsx', sheet_name = 'C. Genes') # Import sheet "C. Genes" from excel file 
EV1Genes=EV1Genes.set_index('Gene ID')

genes=pd.concat([tACI_values,EV1Genes[['Tissue enriched','Group enriched','Tissue enhanced','Classification']]],
                join='inner',axis=1) #dataframe for results using bladder tRNA levels
groups=genes.groupby(['Classification'])
enriched = genes.groupby(['Classification']).get_group('Tissue enriched')
enriched=enriched.drop(['Group enriched', 'Tissue enhanced'],axis=1)
enriched.head()

groups=enriched.groupby(['Tissue enriched'])


#Number of tissue enriched genes in each tissue
x_pos = np.arange(len(enriched['Tissue enriched'].value_counts()))
fig,ax=plt.subplots()
ax.bar(x_pos,enriched['Tissue enriched'].value_counts())
ax.set_xticks(x_pos)
ax.set_xticklabels(enriched['Tissue enriched'].value_counts().index, rotation=90)
ax.set_ylabel("Number of genes on group")
#ax.set_title('Number of tissue enriched genes per tissue')


#Whisker plots for tACI of each tissue-enriched gene group in each tissue
axs=enriched.boxplot(by='Tissue enriched', rot=90,return_type='axes')
for a in axs:
    a.set_ylabel('tACI')


#Whisker plots of tACI of each gene group in each tissue
groups.boxplot(rot=90)
plt.title('Whisker plots of tACI of each gene group in each tissue')

tissue_selection=['Brain','Prostate','Colon','Urinary bladder']
new_gp = pd.concat( [ groups.get_group(group) for group in tissue_selection]).groupby('Tissue enriched')
fig, ax=plt.subplots(nrows=2,ncols=2)
new_gp.boxplot(ax=ax)
for row in ax:
    for col in row:
        col.set_ylabel('tACI')
#Observation: all of them present same pattern in which tissue has larger or smaller tACI than the others
#Same deviation as when plot tACI for all tissues on same graph

##Start comparing tACI given bladder tRNA pool
stats.f_oneway(enriched['bladder taci'][enriched['Tissue enriched'] == 'Urinary bladder'], 
             enriched['bladder taci'][enriched['Tissue enriched'] == 'Prostate'],
             enriched['bladder taci'][enriched['Tissue enriched'] == 'Brain'],
             enriched['bladder taci'][enriched['Tissue enriched'] == 'Colon'],
             enriched['bladder taci'][enriched['Tissue enriched'] == 'Thyroid'])

taci_inbrain_all=tACI_values['brain taci']
taci_inbrain_enriched=enriched['brain taci'][enriched['Tissue enriched'] == 'Brain']

stats.ttest_ind(enriched['brain taci'][enriched['Tissue enriched'] == 'Liver'], taci_inbrain_enriched)