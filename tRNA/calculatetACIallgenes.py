# -*- coding: utf-8 -*-
"""
"calculatetACIallgenes.py"
Use tAI calcualtion object defined in tAI_ojb.py 
to calcualte tACI given varying tRNA pool 
for all genes for which have PTR values

Libraries required for stats analysis section:
researchpy
scipy
statsmodels

Created on Tue Mar  3 18:53:39 2020
@author: Marianne Buat
Group project 2019-2020: codon usage and regulation of gene expression
"""
#import numpy as np
import pandas as pd
import numpy as np
import tAI_obj as t #import all objects and functions defined in that file
import scipy.stats as stats
from matplotlib import pyplot as plt
import usefulfunc as uf #file where generic useful functions were defined

#import file with sequence, use gene ID as index
seq=uf.uploadFile('Table_EV4_excel.xlsx',index_col=1) 

#Load tRNA abundance data
tRNA_brain = uf.uploadFile('tRNA brain av.xlsx', index_col=0)
tRNA_Bcell = uf.uploadFile('tRNA B cells av.xlsx', index_col=0)
tRNA_prostate = uf.uploadFile('tRNA prostate av.xlsx', index_col=0)
tRNA_bladder = uf.uploadFile('tRNA bladder av.xlsx', index_col=0)
tRNA_colon = uf.uploadFile('tRNA colon av.xlsx', index_col=0)



##If need to reupload file after closing Python
tACI_values = uf.uploadFile('tACI_values.xlsx', index_col=0) 
'''
tACI_values=pd.DataFrame() 

## Start with brain
# Reverse complement anticodons in index using indexes in brain data
acs=list(tRNA_brain.index)
codons=list() 
for ac in acs:      
    #convert anticodon to codon it recognises exactly
    codons.append(t.tAI.revcomp(ac)) 

# Get average abundance only by anticodons in Series object
S_tRNA = pd.Series(tRNA_brain["average"], index=codons)
taci_brain=t.tAI(S_tRNA,bacteria=False) #create a tAI object with tRNA abundances to calculate the weights
    
                                       #keep default values for wobble pairing constrains
for ID in seq.index:
    CDS=seq.loc[ID, 'CDS_Sequence']
    tACI_values.loc[ID,'brain taci']=taci_brain.calc(CDS) #add tACI value to dataframe                                           

## Repeat for other tissues
# B cells
# Reverse complement anticodons in index using indexes in B cell data
acs=list(tRNA_Bcell.index)
codons=list() 
for ac in acs:      
    #convert anticodon to codon it recognises exactly
    codons.append(t.tAI.revcomp(ac)) 

# Get average abundance only by anticodons in Series object
S_tRNA = pd.Series(tRNA_Bcell["average"], index=codons)
taci_Bcell=t.tAI(S_tRNA,bacteria=False) #create a tAI object with tRNA abundances to calculate the weights
    
                                       #keep default values for wobble pairing constrains
for ID in seq.index:
    CDS=seq.loc[ID, 'CDS_Sequence']
    tACI_values.loc[ID,'B cell taci']=taci_Bcell.calc(CDS) #add tACI value to dataframe                  

# Prostate
# Reverse complement anticodons in index using indexes in prostate data
acs=list(tRNA_prostate.index)
codons=list() 
for ac in acs:      
    #convert anticodon to codon it recognises exactly
    codons.append(t.tAI.revcomp(ac)) 

# Get average abundance only by anticodons in Series object
S_tRNA = pd.Series(tRNA_prostate["average"], index=codons)
taci_prostate=t.tAI(S_tRNA,bacteria=False) #create a tAI object with tRNA abundances to calculate the weights
    
                                       #keep default values for wobble pairing constrains
for ID in seq.index:
    CDS=seq.loc[ID, 'CDS_Sequence']
    tACI_values.loc[ID,'prostate taci']=taci_prostate.calc(CDS) #add tACI value to dataframe                                           

# Bladder
# Reverse complement anticodons in index using indexes in bladder data
acs=list(tRNA_bladder.index)
codons=list() 
for ac in acs:      
    #convert anticodon to codon it recognises exactly
    codons.append(t.tAI.revcomp(ac)) 

# Get average abundance only by anticodons in Series object
S_tRNA = pd.Series(tRNA_bladder["average"], index=codons)
taci_bladder=t.tAI(S_tRNA,bacteria=False) #create a tAI object with tRNA abundances to calculate the weights
    
                                       #keep default values for wobble pairing constrains
for ID in seq.index:
    CDS=seq.loc[ID, 'CDS_Sequence']
    tACI_values.loc[ID,'bladder taci']=taci_bladder.calc(CDS) #add tACI value to dataframe   
                         
# colon
# Reverse complement anticodons in index using indexes in colon data
acs=list(tRNA_colon.index)
codons=list() 
for ac in acs:      
    #convert anticodon to codon it recognises exactly
    codons.append(t.tAI.revcomp(ac)) 

# Get average abundance only by anticodons in Series object
S_tRNA = pd.Series(tRNA_colon["average"], index=codons)
taci_colon=t.tAI(S_tRNA,bacteria=False) #create a tAI object with tRNA abundances to calculate the weights
    
                                       #keep default values for wobble pairing constrains
for ID in seq.index:
    CDS=seq.loc[ID, 'CDS_Sequence']
    tACI_values.loc[ID,'colon taci']=taci_colon.calc(CDS) #add tACI value to dataframe   

#uf.saveToFile(tACI_values, "tACI_values.xlsx")
'''

## Plotting things

x_pos=np.arange(len(tACI_values))
f,tax=plt.subplots()
tax.scatter(x_pos,tACI_values['brain taci'])
tax.plot(x_pos,[tACI_values['brain taci'].mean()]*len(tACI_values),color='red')


means=[tACI_values['brain taci'].mean(),
tACI_values['bladder taci'].mean(),
tACI_values['colon taci'].mean(),
tACI_values['B cell taci'].mean(),
tACI_values['prostate taci'].mean()]
print("tAI Means:{}".format(means))

std_devs=[tACI_values['brain taci'].std(),
tACI_values['bladder taci'].std(),
tACI_values['colon taci'].std(),
tACI_values['B cell taci'].std(),
tACI_values['prostate taci'].std()]
print("tAI std dev:{}".format(std_devs))

fig, ax = plt.subplots(nrows=2, ncols=3)
val=0
for row in ax:
    for col in row:
        if val>=len(tACI_values.columns):
            tACI_values.boxplot(ax=col)
            col.set_ylabel('tACI')
            col.set_title('Whisker plot of tACI in different tissues')
            break
        col.hist(tACI_values.iloc[:,val],bins=50)
        col.set_ylabel('Occurence')
        col.set_xlabel('tACI value')
        col.set_title('({})'.format(val+1))
        val+=1

plt.suptitle("Order: Brain, B cells, Prostate, Bladder, Colon")
plt.show()

#ANOVA to see whether tACI in one tissue consistently different to that in other tissues
ANOVA_taci0=stats.f_oneway(tACI_values['brain taci'],tACI_values['colon taci'],tACI_values['bladder taci'],tACI_values['B cell taci'],tACI_values['prostate taci'])
eta_taci0=uf.Eta_Squared(ANOVA_taci0[0],k=5,N=len(tACI_values)*5)

#----------------------------------------------------------------------------#
## Repeat tACI calculations but with shifted reading frame
#----------------------------------------------------------------------------#
##If need to reupload file after closing Python
tACI_shift1 = uf.uploadFile('tACI_shift1.xlsx', index_col=0) 
tACI_shift2 = uf.uploadFile('tACI_shift2.xlsx', index_col=0) 

'''
tACI_shift1=pd.DataFrame() 
tACI_shift2=pd.DataFrame() 


## Start with brain
# Reverse complement anticodons in index using indexes in brain data
acs=list(tRNA_brain.index)
codons=list() 
for ac in acs:      
    #convert anticodon to codon it recognises exactly
    codons.append(t.tAI.revcomp(ac)) 

# Get average abundance only by anticodons in Series object
S_tRNA = pd.Series(tRNA_brain["average"], index=codons)
taci_brain=t.tAI(S_tRNA,bacteria=False) #create a tAI object with tRNA abundances to calculate the weights
    
                                       #keep default values for wobble pairing constrains
for ID in seq.index:
    CDS=seq.loc[ID, 'CDS_Sequence'][1:] #shift reading frame by 1 nucleotide
    tACI_shift1.loc[ID,'brain taci']=taci_brain.calc(CDS) #add tACI value to dataframe   
    CDS=seq.loc[ID, 'CDS_Sequence'][2:] #shift reading frame by 1 nucleotide
    tACI_shift2.loc[ID,'brain taci']=taci_brain.calc(CDS)                                         

## Repeat for other tissues
# B cells
# Reverse complement anticodons in index using indexes in B cell data
acs=list(tRNA_Bcell.index)
codons=list() 
for ac in acs:      
    #convert anticodon to codon it recognises exactly
    codons.append(t.tAI.revcomp(ac)) 

# Get average abundance only by anticodons in Series object
S_tRNA = pd.Series(tRNA_Bcell["average"], index=codons)
taci_Bcell=t.tAI(S_tRNA,bacteria=False) #create a tAI object with tRNA abundances to calculate the weights
    
                                       #keep default values for wobble pairing constrains
for ID in seq.index:
    CDS=seq.loc[ID, 'CDS_Sequence'][1:]
    tACI_shift1.loc[ID,'B cell taci']=taci_Bcell.calc(CDS) #add tACI value to dataframe                  
    CDS=seq.loc[ID, 'CDS_Sequence'][2:]
    tACI_shift2.loc[ID,'B cell taci']=taci_Bcell.calc(CDS) #add tACI value to dataframe                  

# Prostate
# Reverse complement anticodons in index using indexes in prostate data
acs=list(tRNA_prostate.index)
codons=list() 
for ac in acs:      
    #convert anticodon to codon it recognises exactly
    codons.append(t.tAI.revcomp(ac)) 

# Get average abundance only by anticodons in Series object
S_tRNA = pd.Series(tRNA_prostate["average"], index=codons)
taci_prostate=t.tAI(S_tRNA,bacteria=False) #create a tAI object with tRNA abundances to calculate the weights
    
                                       #keep default values for wobble pairing constrains
for ID in seq.index:
    CDS=seq.loc[ID, 'CDS_Sequence'][1:]
    tACI_shift1.loc[ID,'prostate taci']=taci_prostate.calc(CDS) #add tACI value to dataframe                         
    CDS=seq.loc[ID, 'CDS_Sequence'][2:]
    tACI_shift2.loc[ID,'prostate taci']=taci_prostate.calc(CDS) #add tACI value to dataframe                                                  

# Bladder
# Reverse complement anticodons in index using indexes in bladder data
acs=list(tRNA_bladder.index)
codons=list() 
for ac in acs:      
    #convert anticodon to codon it recognises exactly
    codons.append(t.tAI.revcomp(ac)) 

# Get average abundance only by anticodons in Series object
S_tRNA = pd.Series(tRNA_bladder["average"], index=codons)
taci_bladder=t.tAI(S_tRNA,bacteria=False) #create a tAI object with tRNA abundances to calculate the weights
    
                                       #keep default values for wobble pairing constrains
for ID in seq.index:
    CDS=seq.loc[ID, 'CDS_Sequence'][1:]
    tACI_shift1.loc[ID,'bladder taci']=taci_bladder.calc(CDS) #add tACI value to dataframe   
    CDS=seq.loc[ID, 'CDS_Sequence'][2:]
    tACI_shift2.loc[ID,'bladder taci']=taci_bladder.calc(CDS)             
# colon
# Reverse complement anticodons in index using indexes in colon data
acs=list(tRNA_colon.index)
codons=list() 
for ac in acs:      
    #convert anticodon to codon it recognises exactly
    codons.append(t.tAI.revcomp(ac)) 

# Get average abundance only by anticodons in Series object
S_tRNA = pd.Series(tRNA_colon["average"], index=codons)
taci_colon=t.tAI(S_tRNA,bacteria=False) #create a tAI object with tRNA abundances to calculate the weights
    
                                       #keep default values for wobble pairing constrains
for ID in seq.index:
    CDS=seq.loc[ID, 'CDS_Sequence'][1:]
    tACI_shift1.loc[ID,'colon taci']=taci_colon.calc(CDS) #add tACI value to dataframe   
    CDS=seq.loc[ID, 'CDS_Sequence'][2:]
    tACI_shift2.loc[ID,'colon taci']=taci_colon.calc(CDS)
    


uf.saveToFile(tACI_shift1, "tACI_shift1.xlsx")
uf.saveToFile(tACI_shift2, "tACI_shift2.xlsx")
'''
## Look at some properties of shifted tACI
#Shifted reading frame by 1 nucleotide
means1=[tACI_shift1['brain taci'].mean(),
tACI_shift1['bladder taci'].mean(),
tACI_shift1['colon taci'].mean(),
tACI_shift1['B cell taci'].mean(),
tACI_shift1['prostate taci'].mean()]
print("tAI Means (1nt shift):{}".format(means1))

std_devs1=[tACI_shift1['brain taci'].std(),
tACI_shift1['bladder taci'].std(),
tACI_shift1['colon taci'].std(),
tACI_shift1['B cell taci'].std(),
tACI_shift1['prostate taci'].std()]
print("tAI std dev (1nt shift):{}".format(std_devs1))

fig1, ax1 = plt.subplots(nrows=2, ncols=3)
val=0
for row in ax1:
    for col in row:
        if val>=len(tACI_shift1.columns):
            tACI_shift1.boxplot(ax=col)
            col.set_ylabel('tACI')
            col.set_title('Whisker plot of tACI in different tissues, shift frame by 1 nt')
            break
        col.hist(tACI_shift1.iloc[:,val],bins=50)
        col.set_ylabel('Occurence')
        col.set_xlabel('tACI value')
        col.set_title('({})'.format(val+1))
        val+=1
plt.suptitle("Order: Brain, B cells, Prostate, Bladder, Colon")

#Shifted reading frame by 2 nucleotide

means2=[tACI_shift2['brain taci'].mean(),
tACI_shift2['bladder taci'].mean(),
tACI_shift2['colon taci'].mean(),
tACI_shift2['B cell taci'].mean(),
tACI_shift2['prostate taci'].mean()]
print("tAI Means (2nt shift):{}".format(means2))

std_devs2=[tACI_shift2['brain taci'].std(),
tACI_shift2['bladder taci'].std(),
tACI_shift2['colon taci'].std(),
tACI_shift2['B cell taci'].std(),
tACI_shift2['prostate taci'].std()]
print("tAI std dev (2nt shift):{}".format(std_devs2))

fig2, ax2 = plt.subplots(nrows=2, ncols=3)
val=0
for row in ax2:
    for col in row:
        if val>=len(tACI_shift2.columns):
            tACI_shift1.boxplot(ax=col)
            col.set_ylabel('tACI')
            col.set_title('Whisker plot of tACI in different tissues, shift frame by 2 nt')
            break
        col.hist(tACI_shift2.iloc[:,val],bins=50)
        col.set_ylabel('Occurence')
        col.set_xlabel('tACI value')
        col.set_title('({})'.format(val+1))
        val+=1

plt.suptitle("Order: Brain, B cells, Prostate, Bladder, Colon")
plt.show()

ANOVA_taci1=stats.f_oneway(tACI_shift1['brain taci'],tACI_shift1['colon taci'],tACI_shift1['bladder taci'],tACI_shift1['B cell taci'],tACI_shift1['prostate taci'])
eta_taci1=uf.Eta_Squared(ANOVA_taci1[0],k=5,N=len(tACI_values)*5)

ANOVA_taci2=stats.f_oneway(tACI_shift2['brain taci'],tACI_shift2['colon taci'],tACI_shift2['bladder taci'],tACI_shift2['B cell taci'],tACI_shift2['prostate taci'])
eta_taci2=uf.Eta_Squared(ANOVA_taci2[0],k=5,N=len(tACI_values)*5)

def avrg(my_list, *args, **kwargs):
    try:
        s=0
        for i in my_list:
            s+=i
        s/=len(my_list)
        return s
    except TypeError:
        print("TypeError, returning object instead")
        return my_list
def standDev(my_list):
    try:
        squared=[m**2 for m in my_list]
        std=len(my_list)/(len(my_list)-1)*(avrg(squared)-avrg(my_list)**2)
        return std
    except TypeError:
        print("Could not carry out calculation, returning 0 instead")
        return 0
        