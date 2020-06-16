# -*- coding: utf-8 -*-
"""
PCA_codon_groups.py
Try to identify if different groups of genes vary in their codon usage on a PCA graph

Created on Sat Jun  6 10:52:25 2020

@author: Marianne
"""
import pandas as pd
import numpy as np
import scipy.stats as stats
from sklearn.decomposition import PCA
from sklearn import preprocessing
import matplotlib.pyplot as plt
import seaborn as sns 
from mpl_toolkits.mplot3d import Axes3D
#from usefulfunc import (uploadFile, SumSqEr,SumSqTr, SumSqTot)

# =============================================================================
# Import data into Python
# =============================================================================

raw_bias = uploadFile("codon_bias.xlsx", index_col=0)
#raw_bias.drop(['ATG','TGG'],axis='columns', inplace=True) #drop amino acids encoded by a single codon
# Drop codons with no synonyms and STOP codons 
raw_bias.drop(['ATG','TGG','TAA','TAG','TGA'],axis='columns', inplace=True)

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


# =============================================================================
# PCA calculations
# =============================================================================
## Version 1: follow PCA analysis method found at https://towardsdatascience.com/pca-using-python-scikit-learn-e653f8989e60
#pca = PCA() #create PCA object
#features= raw_bias.loc[enriched.index,:].values
#features=preprocessing.StandardScaler().fit_transform(features)
#components =  pca.fit_transform(features)

## Follow PCA analysis method used by Ruben
pca = PCA() #create PCA object
scaled_data=preprocessing.scale(raw_bias.loc[enriched.index,:])
scaled_df=pd.DataFrame(scaled_data, 
                       index=raw_bias.loc[enriched.index,:].index, #only keep tissue-enriched genes
                       columns=raw_bias.columns)
scaled_df.dropna(inplace=True)
scaled_data=scaled_df.to_numpy()
pca.fit(scaled_data)
pca_data=pca.transform(scaled_data)

#########################
# Draw a scree plot and a PCA plot
#########################
 
#The following code constructs the Scree plot
per_var = np.round(pca.explained_variance_ratio_* 100, decimals=1)
labels = ['PC' + str(x) for x in range(1, len(per_var)+1)]



fig, ax= plt.subplots()
ax.bar(x=range(1,len(per_var)+1), height=per_var, tick_label=labels)
ax.set_ylabel('Percentage of Explained Variance')
ax.set_xlabel('Principal Component')
ax.tick_params(axis='x', labelrotation=90)
ax.set_title('Scree Plot')
#ax.savefig('PCExplainedVariance.png')


# =============================================================================
# Visualise 2D projection
# =============================================================================
#Only keep tissue-enriched genes
pca_df=pd.DataFrame(pca_data, index=scaled_df.index, columns=labels)
pca_df=pd.concat([pca_df, enriched[['Tissue enriched']]], axis = 1, join='inner')

pd.plotting.scatter_matrix(pca_df.loc[:, "PC1":"PC7"], diagonal="kde")
plt.suptitle("Matrix scatter plot of the first 6 principal components of codon bias for tissue-enriched genes",
             fontsize = 20 )
plt.tight_layout()

#ax = fig.add_subplot(1,1,1)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_xlabel('Principal Component 1', fontsize = 15)
ax.set_ylabel('Principal Component 2', fontsize = 15)
ax.set_zlabel('Principal Component 3', fontsize = 15)
ax.set_title('3 component PCA', fontsize = 20)
#tissue_selection=['Lymph node', 'Gallbladder','Esophagus']
#tissue_selection=['Brain','Prostate','Colon','Urinary bladder'] #tissues for which have tRNA data in addition to abundance
tissue_selection = ['Brain','Testis','Liver','Heart','Fallopian tube']
#tissue_selection = list(np.unique(enriched['Tissue enriched'].values))
#colors = ['r', 'b', 'g','y','m']
colors = ['r', 'b', 'g','y','m']
for tissue, color in zip(tissue_selection,colors):
#for tissue in tissue_selection:
    indicesToKeep = pca_df['Tissue enriched'] == tissue
    ax.scatter(pca_df.loc[indicesToKeep, 'PC1']
               , pca_df.loc[indicesToKeep, 'PC2']
               #, zs=pca_df.loc[indicesToKeep, 'PC3']
               , c = color
               , s = 10)
ax.legend(tissue_selection)
ax.grid()

ax.scatter(pca_df['PC1'][pca_df['Tissue enriched']=='Liver'],
           pca_df['PC2'][pca_df['Tissue enriched']=='Liver'],
           pca_df['PC3'][pca_df['Tissue enriched']=='Liver'])


fig2 = plt.figure()
ax2 = fig2.add_subplot(111, projection='3d')
ax2.set_zlabel('Principal Component 3', fontsize = 15)
ax2.set_xlabel('Principal Component 4', fontsize = 15)
ax2.set_ylabel('Principal Component 5', fontsize = 15)
ax2.scatter(pca_df['PC4'][pca_df['Tissue enriched']=='Liver'],
           pca_df['PC5'][pca_df['Tissue enriched']=='Liver'],
           pca_df['PC3'][pca_df['Tissue enriched']=='Liver'],
           label='Liver')


##Investigate codon usage of the blobs on the 3rd PC direction

posblob_bias = raw_bias.loc[pca_df[pca_df['PC3']>=0].index]
posblob_bias['PC3 sign']='positive'
negblob_bias = raw_bias.loc[pca_df[pca_df['PC3']<0].index]
negblob_bias['PC3 sign']='negative'
blob_bias=pd.concat([posblob_bias, negblob_bias], axis=0)

blob_data=pd.melt(blob_bias, id_vars = ['PC3 sign'], value_vars=list(raw_bias.columns),
        var_name='codon', value_name='codon bias')
fig,ax=plt.subplots()
sns.catplot(x="codon", y="codon bias", hue="PC3 sign", kind="box", data=blob_data, ax=ax)
ax.tick_params(axis='x', labelrotation=90)



# =============================================================================
# ANOVA on pronciple component 1
# =============================================================================
##Compare codon usage between some subsets of tissue-enriched groups

#Start with tissues for which also have tRNA data
tissue_selection=['Brain','Prostate','Colon','Urinary bladder'] 
tissue_data = [pca_df['PC1'][pca_df['Tissue enriched'] == tissue]  \
                for tissue in tissue_selection]
print("Tissues:", *[tissue+' ' for tissue in tissue_selection])
print(stats.f_oneway(*tissue_data))
#If remove Brain:
tissue_data = [pca_df['PC1'][pca_df['Tissue enriched'] == tissue]  \
                for tissue in tissue_selection if tissue !='Brain']
print("Tissues:", *[tissue+' ' for tissue in tissue_selection if tissue !='Brain'])
print(stats.f_oneway(*tissue_data))


##Subsets of tissue-enriched groups with most genes in them
tissue_selection = ['Brain', 'Heart','Fallopian tube','Liver','Testis']
tissue_data = [pca_df['PC1'][pca_df['Tissue enriched'] == tissue]  \
                for tissue in tissue_selection]
print("Tissues:", *[tissue+' ' for tissue in tissue_selection])
print(stats.f_oneway(*tissue_data))
#If remove Brain:
tissue_data = [pca_df['PC1'][pca_df['Tissue enriched'] == tissue]  \
                for tissue in tissue_selection if tissue !='Brain']
print("Tissues:", *[tissue+' ' for tissue in tissue_selection if tissue !='Brain'])
print(stats.f_oneway(*tissue_data))


##Change to numpy array to calculate other statistics
np_tissues=[pca_df['PC1'][pca_df['Tissue enriched'] == tissue].to_numpy()  \
            for tissue in tissue_selection]

# =============================================================================
# Find which tissues may be statistically different from the others
# =============================================================================

# Start with ANOVA on all tissues
tissue_selection=list(np.unique(enriched['Tissue enriched'].values))

tissue_data=[pca_df['PC1'][pca_df['Tissue enriched'] == tissue]  \
                for tissue in list(np.unique(enriched['Tissue enriched'].values))]
np_tissues=[dataser.to_numpy() for dataser in tissue_data] #numpy version of list
Fstat, pval = stats.f_oneway(*tissue_data)
print(stats.f_oneway(*tissue_data))

eta_sq=SumSqTr(np_tissues)/SumSqTot(np_tissues)

#Exclude brain from analysis
tissue_data = [pca_df['PC1'][pca_df['Tissue enriched'] == tissue]  \
                for tissue in tissue_selection if tissue !='Brain']
print("Tissues:", *[tissue+' ' for tissue in tissue_selection if tissue !='Brain'])
print(stats.f_oneway(*tissue_data))


# =============================================================================
# Brain seems different from the others, so look at its codon usage
# =============================================================================
##Create boxplot for brain's codon usage vs other tissues
brain_bias=raw_bias[pca_df['Tissue enriched']=='Brain'].copy()
others_bias=raw_bias[pca_df['Tissue enriched']!='Brain'].copy()
brain_bias['Tissue']='Brain'
others_bias['Tissue']='Non-brain'

tissues_bias=pd.concat([brain_bias, others_bias], axis=0)
tissues_data=pd.melt(tissues_bias, id_vars = ['Tissue'], value_vars=list(raw_bias.columns),
        var_name='codon', value_name='codon bias')
fig,ax=plt.subplots()
sns.catplot(x="codon", y="codon bias", hue="Tissue", kind="box", data=tissues_data, ax=ax)
ax.tick_params(axis='x', labelrotation=90)

##To compare profile, draw heart's codon usage vs non-brain nor heart genes
heart_bias=raw_bias[pca_df['Tissue enriched']=='Brain'].copy()
others_bias= raw_bias[not (pca_df['Tissue enriched']=='Brain' or pca_df['Tissue enriched']=='Heart')].copy()
heart_bias['Tissue']='Heart'
others_bias['Tissue']='Non-brain nor heart'

tissues_bias=pd.concat([brain_bias, others_bias], axis=0)
tissues_data=pd.melt(tissues_bias, id_vars = ['Tissue'], value_vars=list(raw_bias.columns),
        var_name='codon', value_name='codon bias')
fig,ax=plt.subplots()
sns.catplot(x="codon", y="codon bias", hue="Tissue", kind="box", data=tissues_data, ax=ax)
ax.tick_params(axis='x', labelrotation=90)
