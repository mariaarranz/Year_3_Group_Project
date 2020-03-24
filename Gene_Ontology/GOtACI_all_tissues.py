#Code to get tACI values for all 4 tissues with tRNA data

import pandas as pd
import matplotlib.pyplot as plt
import itertools
import xlrd
import numpy as np
from pandas.plotting import boxplot_frame_groupby
import matplotlib as mpl

#READ IN DATA

#read in just gene sheet of table EV1 excel file
EV1data= pd.read_excel('Table_EV1.xlsx', sheet_name='C. Genes', index_col=0) 

#read in tACI values for all tissues
tAIdata = pd.read_excel('tACI_values.xlsx', index_col=0) 

#merge both tables to get tACI values for each gene
merged = pd.concat([EV1data[[ 'Gene name', 'Tissue enriched']], tAIdata], axis=1, join='inner')

#group genes by enrichment
tissuegroup = merged.groupby('Tissue enriched')

#-----------------------------------------------------------------------------------------------------------
#COLON
#get tRNA data correpsonding to each colon enriched gene
colon=tissuegroup.get_group('Colon')
colontACI=colon[['Gene name', 'brain taci', 'colon taci', 'prostate taci', 'bladder taci']]
#colontACI.to_excel(r'/Users/cerysbarclay/Desktop/colontaci.xlsx') 
#print this out and search these genes in GO term mapper. Save as new excel file
colontACI=colontACI.set_index('Gene name')

#read in colon ontology GO slim data
colonGOdata= pd.read_excel(r'/Users/cerysbarclay/Desktop/colonGOdata.xlsx')
colonGOdata=colonGOdata.set_index('GOID')
colonGOdata['split names']=colonGOdata["ANNOTATED_GENES"].str.split(', ')

#loop through GO terms and for every gene in each term average all tACI values
for row in colonGOdata.index:
        names=colonGOdata.loc[row, 'split names']
        braintACIs=colontACI.loc[names, 'brain taci']
        colonGOdata.loc[row, 'mean brain taci']=braintACIs.mean()

        colontACIs=colontACI.loc[names, 'colon taci']
        colonGOdata.loc[row, 'mean colon taci']=colontACIs.mean()

        bladdertACIs=colontACI.loc[names, 'bladder taci']
        colonGOdata.loc[row, 'mean bladder taci']=bladdertACIs.mean()

        prostatetACIs=colontACI.loc[names, 'prostate taci']
        colonGOdata.loc[row, 'mean prostate taci']=prostatetACIs.mean()
   
#Create new dataframe to store all tACI data in
tRNAdisplay=colonGOdata[['TERM']]
tRNAdisplay=tRNAdisplay.sort_values('TERM')
tRNAdisplay['colon with brain tRNA']=colonGOdata[['mean brain taci']]
tRNAdisplay['colon with colon tRNA']=colonGOdata[['mean colon taci']]
tRNAdisplay['colon with bladder tRNA']=colonGOdata[['mean bladder taci']]
tRNAdisplay['colon with prostate tRNA']=colonGOdata[['mean prostate taci']]

#-----------------------------------------------------------------------------------------------------------
#BRAIN

#get tRNA data correpsonding to each brain enriched gene
brain =tissuegroup.get_group('Brain')
braintACI=brain[['Gene name', 'brain taci', 'colon taci', 'prostate taci', 'bladder taci']]
#colontACI.to_excel(r'/Users/cerysbarclay/Desktop/colontaci.xlsx') 
#print this out and search these genes in GO term mapper. Save as new excel file
braintACI=braintACI.set_index('Gene name')

#read in brain ontology GO slim data
brainGOdata= pd.read_excel(r'/Users/cerysbarclay/Desktop/GOdata.xlsx')
brainGOdata=brainGOdata.set_index('GOID')
brainGOdata['split names']=brainGOdata["ANNOTATED_GENES"].str.split(', ')

#loop through GO terms and for every gene in each term average all tACI values
for row in brainGOdata.index:
        names=brainGOdata.loc[row, 'split names']
        #may be issue here
        braintACIs=braintACI.loc[names, 'brain taci']
        brainGOdata.loc[row, 'mean brain taci']=braintACIs.mean()

        colontACIs=braintACI.loc[names, 'colon taci']
        brainGOdata.loc[row, 'mean colon taci']=colontACIs.mean()

        bladdertACIs=braintACI.loc[names, 'bladder taci']
        brainGOdata.loc[row, 'mean bladder taci']=bladdertACIs.mean()

        prostatetACIs=braintACI.loc[names, 'prostate taci']
        brainGOdata.loc[row, 'mean prostate taci']=prostatetACIs.mean()
   
tRNAdisplay['brain with brain tRNA']=brainGOdata[['mean brain taci']]
tRNAdisplay['brain with colon tRNA']=brainGOdata[['mean colon taci']]
tRNAdisplay['brain with bladder tRNA']=brainGOdata[['mean bladder taci']]
tRNAdisplay['brain with prostate tRNA']=brainGOdata[['mean prostate taci']]

#-----------------------------------------------------------------------------------------------------------
#BLADDER

#get tRNA data corresponding to each bladder enriched gene
bladder=tissuegroup.get_group('Urinary bladder')
bladdertACI=bladder[['Gene name', 'brain taci', 'colon taci', 'prostate taci', 'bladder taci']]
#bladdertACI.to_excel(r'/Users/cerysbarclay/Desktop/bladdertaci.xlsx') 
#print this out and search these genes in GO term mapper. Save as new excel file
bladdertACI=bladdertACI.set_index('Gene name')

#read in bladder ontology GO slim data
bladderGOdata= pd.read_excel(r'/Users/cerysbarclay/Desktop/bladderGOdata.xlsx')
bladderGOdata=bladderGOdata.set_index('GOID')
bladderGOdata['split names']=bladderGOdata["ANNOTATED_GENES"].str.split(', ')

#loop through GO terms and for every gene in each term average all tACI values
for row in bladderGOdata.index:
        names=bladderGOdata.loc[row, 'split names']
 
        braintACIs=bladdertACI.loc[names, 'brain taci']
        bladderGOdata.loc[row, 'mean brain taci']=braintACIs.mean()

        colontACIs=bladdertACI.loc[names, 'colon taci']
        bladderGOdata.loc[row, 'mean colon taci']=colontACIs.mean()

        bladdertACIs=bladdertACI.loc[names, 'bladder taci']
        bladderGOdata.loc[row, 'mean bladder taci']=bladdertACIs.mean()

        prostatetACIs=bladdertACI.loc[names, 'prostate taci']
        bladderGOdata.loc[row, 'mean prostate taci']=prostatetACIs.mean()
   
tRNAdisplay['bladder with brain tRNA']=bladderGOdata[['mean brain taci']]
tRNAdisplay['bladder with colon tRNA']=bladderGOdata[['mean colon taci']]
tRNAdisplay['bladder with bladder tRNA']=bladderGOdata[['mean bladder taci']]
tRNAdisplay['bladder with prostate tRNA']=bladderGOdata[['mean prostate taci']]

#-----------------------------------------------------------------------------------------------------------
#PROSTATE

#get tRNA data corresponding to each prostate enriched gene
prostate=tissuegroup.get_group('Prostate')
prostatetACI=prostate[['Gene name', 'brain taci', 'colon taci', 'prostate taci', 'bladder taci']]
#prostatetACI.to_excel(r'/Users/cerysbarclay/Desktop/prostatetaci.xlsx') 
#print this out and search these genes in GO term mapper. Save as new excel file
prostatetACI=prostatetACI.set_index('Gene name')

#read in prostate ontology GO slim data
prostateGOdata= pd.read_excel(r'/Users/cerysbarclay/Desktop/prostateGOdata.xlsx')
prostateGOdata=prostateGOdata.set_index('GOID')
prostateGOdata['split names']=prostateGOdata["ANNOTATED_GENES"].str.split(', ')

#loop through GO terms and for every gene in each term average all tACI values
for row in prostateGOdata.index:
        names=prostateGOdata.loc[row, 'split names']
 
        braintACIs=prostatetACI.loc[names, 'brain taci']
        prostateGOdata.loc[row, 'mean brain taci']=braintACIs.mean()

        colontACIs=prostatetACI.loc[names, 'colon taci']
        prostateGOdata.loc[row, 'mean colon taci']=colontACIs.mean()

        bladdertACIs=prostatetACI.loc[names, 'bladder taci']
        prostateGOdata.loc[row, 'mean bladder taci']=bladdertACIs.mean()

        prostatetACIs=prostatetACI.loc[names, 'prostate taci']
        prostateGOdata.loc[row, 'mean prostate taci']=prostatetACIs.mean()
   

tRNAdisplay['prostate with brain tRNA']=prostateGOdata[['mean brain taci']]
tRNAdisplay['prostate with colon tRNA']=prostateGOdata[['mean colon taci']]
tRNAdisplay['prostate with bladder tRNA']=prostateGOdata[['mean bladder taci']]
tRNAdisplay['prostate with prostate tRNA']=prostateGOdata[['mean prostate taci']]

#convert number of genes to percentages
tRNAdisplay['number of colon genes']=(colonGOdata[['NUM_LIST_ANNOTATIONS']]/196)*100
tRNAdisplay['number of brain genes']=(brainGOdata[['NUM_LIST_ANNOTATIONS']]/1192)*100
tRNAdisplay['number of bladder genes']=(bladderGOdata[['NUM_LIST_ANNOTATIONS']]/154)*100
tRNAdisplay['number of prostate genes']=(prostateGOdata[['NUM_LIST_ANNOTATIONS']]/266)*100

tRNAdisplay.to_excel(r'/Users/cerysbarclay/Desktop/alltissuestRNAdisplay3.xlsx') 
