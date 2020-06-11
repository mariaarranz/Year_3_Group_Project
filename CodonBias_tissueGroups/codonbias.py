# -*- coding: utf-8 -*-
"""
Codon bias
Count codon occurences of and relative abundance to other synonymous codons
Compare with shifted sequences and randomly-generated sequences

Requires these files in working directory:
    - 'Table_EV4_excel.xlsx' where the CDS sequences are conserved
    (or pre-import into dataframe variable 'seq')
    - usefulfunc.py for some generic functions in use in the code
    (or pre-run the file for using functions directly/import as uf from other file)

Created on Wed Jun  3 15:19:05 2020

@author: Marianne
"""
import pandas as pd
import numpy as np
#import scipy.stats as stats
#import itertools
import re
from matplotlib import pyplot as plt
#import usefulfunc as uf #file where generic useful functions were defined
#Or import only the functions we need:
#from usefulfunc import (uploadFile, saveToFile, 
#                        revcomp,randomSeqGenerator) 

# =============================================================================
# Define functions for use throughout file
# =============================================================================
def countTriplets(seq, codon_list=None):
    #Count number of unique repeat of 3 letters in a string
    #Input seq is a string, output is a Pandas Series with index the triplets found 
    #and entries the counts
    codon_id, codon_count = np.unique(re.findall('...', seq),
                                          return_counts=True)
   
    codon_count = pd.Series(codon_count, index=codon_id, name='codons')
    if codon_list!=None and isinstance(codon_list,list): #counts for predetermined list of triplets 
        # get values for standard codons only by joining with codon list series
        #change codons that were not found in sequence to zero counts
        codon_count=pd.Series(codon_count, index=codon_list).replace(np.NaN,0)       
        
    return codon_count

def calc_codonBias(CDS,genetic_code, return_count=False):
    """Calculate codon bias of a coding sequence
    codons is the list of valid codons that are counted
    """
    codon_list=list(get_all_values(genetic_code))
    CDS = CDS.upper().replace('U', 'T') #change to standard dNTP format
    codon_count = countTriplets(CDS, codon_list)    
    
    codon_bias=pd.Series(index=codon_list)
    for aa in genetic_code.keys():
        codon_count[aa]=codon_count[genetic_code[aa]].sum()
        for c in genetic_code[aa]:
            codon_bias[c]=codon_count[c]/codon_count[aa]
    if return_count:
        return codon_count, codon_bias
    else:
        return codon_bias
    

def get_all_values(d):
    #function I found to unravel all values contained within a nested dict
    #need to take list() of return value to have as a list
    #taken from https://stackoverflow.com/questions/7002429/how-can-i-extract-all-values-from-a-dictionary-in-python
    if isinstance(d, dict):
        for v in d.values():
            yield from get_all_values(v)
    elif isinstance(d, list):
        for v in d:
            yield from get_all_values(v)
    else:
        yield d 
        

    
# =============================================================================
# Main part of code
# =============================================================================

## Import file with sequence, use gene ID as index
#seq=uploadFile('Table_EV4_excel.xlsx',index_col=1) 

#Create dictionary for synonymous codons corresponding to each amino acid
genetic_code = {'Arg':['AGG','AGA','CGT','CGC','CGA','CGG'],
                'Leu':['TTA','TTG','CTT','CTC','CTA','CTG'],
                'Ser':['AGC','AGT','TCT','TCC','TCA','TCG'],
                'Ala':['GCA','GCC','GCG','GCT'],
                'Gly':['GGA','GGC','GGG','GGT'],
                'Pro':['CCT','CCC','CCA','CCG'],
                'Thr':['ACT','ACC','ACA','ACG'],
                'Val':['GTT','GTC','GTA','GTG'],
                'Ile':['ATT','ATC','ATA'],
                'Asp':['GAC','GAT'],
                'Glu':['GAA','GAG'],                
                'Phe':['TTT','TTC'],                
                'Tyr':['TAT','TAC'],
                'Cys':['TGT','TGC'],           
                'His':['CAT','CAC'],
                'Gln':['CAA','CAG'],
                'Asn':['AAT','AAC'],
                'Lys':['AAA','AAG'],                
                'Trp':['TGG'],
                'Met':['ATG'],
                'Stop':['TAA','TAG','TGA']}
all_codons=list(get_all_values(genetic_code)) #list of all possible 64 codons

##If for whatever reason want to order synonymous codons by alphabetic order
## Note it acts inplace (a.sort() will replace a by a sorted list and return None)
gencode=genetic_code.copy() #so do not alter original object 
                            #(for compatibility with order in existing objects)
for i in gencode:
    gencode[i].sort()
    print(gencode[i])
all_codons_sorted = list(get_all_values(gencode))
    
# =============================================================================
# Calculate codon bias for in-frame CDSs
# =============================================================================
## Calculate codon usage for all genes
codon_tables=pd.DataFrame(0,index=seq.index, columns=all_codons)
for ID in seq.index:
    CDS=seq.loc[ID, 'CDS_Sequence']
    CDS = CDS.upper().replace('U', 'T') #change to standard dNTP format
    codon_count = countTriplets(CDS, all_codons)
    codon_tables.loc[ID, codon_count.index]=codon_count

##Codon bias for each group of synonymous codons, calculate for each gene
codon_bias=codon_tables.copy()
for aa in genetic_code.keys():
    codon_tables[aa]=codon_tables[genetic_code[aa]].sum(axis=1)
    for c in genetic_code[aa]:
        codon_bias[c]=codon_tables[c]/codon_tables[aa]
        
#Find average bias over all sequences
codon_summed = codon_tables.sum() 
av_bias=pd.Series(index=all_codons)
for aa in genetic_code.keys():
    for c in genetic_code[aa]:
        av_bias[c]=codon_summed[c]/codon_summed[aa]

        
#Standard deviation in codon bias
sqrd_diff = codon_bias.copy()
#take squared difference of codon biases from mean for each gene and codon
for col in codon_bias.columns: 
    sqrd_diff[col]=np.square(sqrd_diff[col]-av_bias[col])
std_bias = np.sqrt(sqrd_diff.sum()/len(sqrd_diff)) #standard dev for each codon
    
##Save codon bias of genes       
#saveToFile(codon_bias,'codon_bias.xlsx')

# =============================================================================
## Repeat with shifted reading frames 
# Start shift from 2nd codon to avoid over-representing codons 
# starting in "TG."
# =============================================================================

## Shift reading frame by 1 
codon_tables_shift1=pd.DataFrame(0,index=seq.index, columns=all_codons)
for ID in seq.index:
    CDS=seq.loc[ID, 'CDS_Sequence'][4:]
    CDS= CDS.upper().replace('U', 'T') #change to standard dNTP format
    codon_count = countTriplets(CDS, all_codons)
    codon_tables_shift1.loc[ID, codon_count.index]=codon_count

codon_bias_shift1=codon_tables_shift1.copy()
for aa in genetic_code.keys():
    codon_tables_shift1[aa]=codon_tables_shift1[genetic_code[aa]].sum(axis=1)
    for c in genetic_code[aa]:
        codon_bias_shift1[c]=codon_tables_shift1[c]/codon_tables_shift1[aa]

#Find average bias over all sequences
codon_summed1 = codon_tables_shift1.sum() 
av_bias_shift1=pd.Series(index=all_codons)
for aa in genetic_code.keys():
    for c in genetic_code[aa]:
        av_bias_shift1[c]=codon_summed1[c]/codon_summed1[aa]
#Standard deviation in codon bias
sqrd_diff = codon_bias_shift1.copy()
#take squared difference of codon biases from mean for each gene and codon
for col in codon_bias_shift1.columns: 
    sqrd_diff[col]=np.square(sqrd_diff[col]-av_bias[col])
std_bias_shift1 = np.sqrt(sqrd_diff.sum()/len(sqrd_diff)) #standard dev for each codon

#saveToFile(codon_bias_shift1,'codon_bias_shift1.xlsx')       

## Shift reading frame by 2 
codon_tables_shift2=pd.DataFrame(0,index=seq.index, columns=all_codons)
for ID in seq.index:
    CDS=seq.loc[ID, 'CDS_Sequence'][5:]
    CDS=CDS.upper().replace('U', 'T') #change to standard dNTP format
    codon_count = countTriplets(CDS, all_codons)
    codon_tables_shift2.loc[ID, codon_count.index]=codon_count

codon_bias_shift2=codon_tables_shift2.copy()
for aa in genetic_code.keys():
    codon_tables_shift2[aa]=codon_tables_shift2[genetic_code[aa]].sum(axis=1)
    for c in genetic_code[aa]:
        codon_bias_shift2[c]=codon_tables_shift2[c]/codon_tables_shift2[aa]

#Find average bias over all sequences
codon_summed2 = codon_tables_shift2.sum() 
av_bias_shift2=pd.Series(index=all_codons)
for aa in genetic_code.keys():
    for c in genetic_code[aa]:
        av_bias_shift2[c]=codon_summed2[c]/codon_summed2[aa]
#Standard deviation in codon bias
sqrd_diff = codon_bias_shift2.copy()
#take squared difference of codon biases from mean for each gene and codon
for col in codon_bias_shift2.columns: 
    sqrd_diff[col]=np.square(sqrd_diff[col]-av_bias[col])
std_bias_shift2 = np.sqrt(sqrd_diff.sum()/len(sqrd_diff)) #standard dev for each codon

#saveToFile(codon_bias_shift2,'codon_bias_shift2.xlsx')     



# =============================================================================
# Generate random sequences with equivalent nucleotide content to our genes
# (null models)
# =============================================================================

## Right now generating equivalent sequences for each gene is commented 
# because there is a bug that I still need to fix
'''
with warnings.catch_warnings():
    warnings.simplefilter("error", RuntimeWarning)

    rand_equivalent = pd.DataFrame(index=seq.index, columns=['rand_CDS']+all_codons+list(genetic_code.keys()))
    rand_bias=pd.DataFrame(index=seq.index, columns=all_codons)
    for ID in seq.index:
        CDS=seq.loc[ID, 'CDS_Sequence'].upper().replace('U', 'T') #change to standard dNTP format
        nt, nt_content = np.unique(re.findall('.', CDS),
                                   return_counts=True) #count each nucleotide in sequence
        nt_content = pd.Series(nt_content, index=nt)
        nt_content = nt_content/nt_content.sum() #convert to frequency
        randCDS=randomSeqGenerator(len(CDS), [nt_content['A'], nt_content['T'],nt_content['C'], nt_content['G']]) 
        cc, cb=calc_codonBias(randCDS,genetic_code, return_count=True)
        rand_equivalent.loc[ID, all_codons+list(genetic_code.keys())]=cc
        rand_bias.loc[ID]=cb

codon_summed = rand_equivalent[all_codons+list(genetic_code.keys())].sum()
av_bias_rand=pd.Series(index=all_codons)
for aa in genetic_code.keys():
    for c in genetic_code[aa]:
        av_bias_rand[c]=codon_summed[c]/codon_summed[aa]
#Standard deviation in codon bias
sqrd_diff = rand_bias.copy()
#take squared difference of codon biases from mean for each gene and codon
for col in rand_bias.columns: 
    sqrd_diff[col]=np.square(sqrd_diff[col]-av_bias_rand[col])
std_bias_rand = np.sqrt(sqrd_diff.sum()/len(sqrd_diff)) #standard dev for each codon

'''  

## Generate a single, long random sequence representative of 
# the overall nucleotide content of the genome

nt_content=pd.Series([0,0,0,0], index=['A','T','C','G'])
for codon in all_codons:
    for letter in list(codon):
        nt_content[letter]+=codon_summed[codon]
nt_freq=nt_content/nt_content.sum()            

srand = randomSeqGenerator(500000, nt_freq[['A','T','C','G']].to_list()) 
srand_bias=calc_codonBias(srand,genetic_code)

## Random sequence if nucleotides were equally present
sranduni = randomSeqGenerator(500000)
sranduni_bias=calc_codonBias(sranduni,genetic_code)

# =============================================================================
# Visualise average codon bias
# =============================================================================


fig,ax=plt.subplots(3,1)
av = [av_bias,av_bias_shift1,av_bias_shift2]
std =[std_bias,std_bias_shift1,std_bias_shift2]
for i in range(3):
    x_pos=np.arange(len(av[i]))
    ax[i].bar(x_pos, av[i],
           yerr=std[i],
           label='{} shift'.format(i))
    ax[i].set_ylabel('codon bias')
    ax[i].set_xticks(x_pos)
    ax[i].set_xticklabels(av[i].index,rotation=90)
    ax[i].yaxis.grid(True)
    ax[i].legend()
fig.suptitle('Codon bias')
    

fig, ax = plt.subplots()
for i in range(3):
    x_pos=np.arange(len(av[i]))
    ax.errorbar(x_pos, av[i], 
           yerr=std[i], fmt = 'x',
           label='{} shift'.format(i))
ax.set_ylabel('codon bias')
ax.set_xticks(x_pos)
ax.set_xticklabels(av[i].index,rotation=90)
ax.yaxis.grid(True)
ax.set_title('Codon bias')
ax.legend()



# =============================================================================
# Visually check changing reading frame completely alters codon usage
# =============================================================================

## With codon frequency
fig,ax = plt.subplots()
desc_freq=codon_tables.sum()[all_codons]/(codon_tables.sum()[all_codons].sum())
desc_freq.sort_values(ascending=False, inplace=True)

desc_freq_shift1=pd.Series(codon_tables_shift1.sum()[all_codons]/(codon_tables_shift1.sum()[all_codons].sum()),index=desc_freq.index)
desc_freq_shift2=pd.Series(codon_tables_shift2.sum()[all_codons]/(codon_tables_shift2.sum()[all_codons].sum()),index=desc_freq.index)

ax.plot(desc_freq_shift1,'b*')
ax.plot(desc_freq_shift2,'r*')
ax.plot(desc_freq,'k*')

ax.set_xticklabels(desc_freq.index,rotation=90, fontsize = 12)

        
# =============================================================================
# Visualise changing codon bias vs null models   
# =============================================================================
fig, ax=plt.subplots()
# Plot bias in uniform and occurence-ajusted model random sequences
xpos = np.arange(len(all_codons))
ax.plot(xpos, sranduni_bias[all_codons], 'k-', label = "uniform model")
ax.plot(xpos, srand_bias[all_codons], 'b--', label = "content-ajusted model")
# Plot average bias in actual and shifted CDSs
ax.plot(xpos, av_bias[all_codons],'rx',label = "average bias")
#ax.errorbar(x_pos, av_bias, 
#           yerr=std_bias, fmt = 'rx',label = "average bias with sd")
#ax.plot(xpos, av_bias_shift1[all_codons],'gx', label = "1 shift")
#ax.plot(xpos, av_bias_shift2[all_codons],'yx',label = "2 shift")

ax.set_xticklabels(all_codons,rotation=90, fontsize = 12)
ax.set_ylabel('codon bias')
ax.set_xticks(x_pos)
ax.yaxis.grid(True)
ax.set_title('Average codon bias',fontsize = 20)
ax.legend()

codon_bias.boxplot(rot=90)              


