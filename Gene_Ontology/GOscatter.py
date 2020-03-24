
import pandas as pd
import matplotlib.pyplot as plt
import itertools
import xlrd
import numpy as np
from pandas.plotting import boxplot_frame_groupby
import matplotlib as mpl
import time

#READ IN DATA
data=pd.read_excel('GOtACIdata.xlsx', sheet_name='plain table', index_col=0) 

#create new dataframe to store data for plotting
colon=data[['% of colon genes', '%  of brain genes']]

#Add data relating to colon enriched genes to dataframe
colon['colonwithbraintrna']=data[['colon with brain tRNA']]
colon['colonwithcolontrna']=data[['colon with colon tRNA']]
colon['colonwithbladdertrna']=data[['colon with bladder tRNA']]
colon['colonwithprostatetrna']=data[['colon with prostate tRNA']]

#add this line if we were to look at brain enriched gene set
#colon['brainwithbraintrna']=data[['brain with brain tRNA']]

#Loop through each term and get number of genes assigned and mean tACI value
colonx=[]
colony=[]
for row in colon.index:
        trna=colon.loc[row, 'colonwithcolontrna']
        genes=colon.loc[row, '% of colon genes']
        colonx.append(genes)
        colony.append(trna)
        
#plot mean tACI against the number of genes in a sctater graph
plt.subplot(221)
plt.scatter(colonx, colony, c='green')
plt.xlabel('Genes assigned to GO term (%)')
plt.ylabel('Mean tACI')
plt.title('Colon tRNA')

#Repeat for each tRNA pool, adding each new graph as a subplot
brainx=[]
brainy=[]
for row in colon.index:
        trna=colon.loc[row, 'colonwithbraintrna']
        genes=colon.loc[row, '% of colon genes']
        brainx.append(genes)
        brainy.append(trna)

#plt.figure(2)
plt.subplot(222)
plt.scatter(brainx, brainy, c='blue')
plt.xlabel('Genes assigned to GO term (%)')
plt.ylabel('Mean tACI')
plt.title('Brain tRNA')

bladderx=[]
bladdery=[]
for row in colon.index:
        trna=colon.loc[row, 'colonwithbladdertrna']
        genes=colon.loc[row, '% of colon genes']
        bladderx.append(genes)
        bladdery.append(trna)

#plt.figure(3)
plt.subplot(223)
plt.scatter(bladderx, bladdery, c='red')
plt.xlabel('Genes assigned to GO term (%)')
plt.ylabel('Mean tACI')
plt.title('Bladder tRNA')


prostatex=[]
prostatey=[]
for row in colon.index:
        trna=colon.loc[row, 'colonwithprostatetrna']
        genes=colon.loc[row, '% of colon genes']
        prostatex.append(genes)
        prostatey.append(trna)

plt.subplot(224)
plt.scatter(prostatex, prostatey, c='orange')
plt.xlabel('Genes assigned to GO term (%)')
plt.ylabel('Mean tACI')
plt.title('Prostate tRNA')

#Display figure
plt.show()
