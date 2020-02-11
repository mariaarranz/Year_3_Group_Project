#PROGRAM JUST TO GET BOXPLOT FOR THE LIVER- after completion will put into for loop for all tissues

#Need to explore how showfliers is actually removing outliers

import pandas as pd
import matplotlib.pyplot as plt
import itertools
import xlrd
import numpy as np
from pandas.plotting import boxplot_frame_groupby
import matplotlib as mpl

#STEP1: READ IN DATA
data = pd.read_excel(r'/Users/cerysbarclay/Desktop/EV1.xlsx', sheet_name='C. Genes') #read in just gene sheet of table EV1 excel file (stored on h drive)
df=data.iloc[:,1:32] #Select only columns with abundance data

#STEP 2: GROUP BY ENRICHED TISSUE
grouped = df.groupby('Tissue enriched') #Group all genes by tissue they are enriched in 

for name, group in grouped: #print to check it works
   print(name)
   print(group)

#STEP 3: CREATE EXCEL FILE DISPLAYING LIVER DATA (ALL GENES ENRICHED IN LIVER WITH CORRESPONDING ABUNDANCES FOR EVERY TISSUE)
liver=grouped.get_group('Liver') #get gene name and abundance data just for liver
print(liver) #print to check 

liver_abundance=liver.drop(['Tissue enriched'], axis=1) #drop enriched column as no longer needed for boxplot

liver_abundance.to_excel(r'/Users/cerysbarclay/Desktop/Liver_data5.xlsx', sheet_name="Liver") 

#STEP 4: FOR EACH TISSUE AS A GROUP, CREATE A BOXPLOT OF ABUNDANCE DATA
dataforboxplot= liver_abundance.drop(['Gene name'], axis=1) #drop gene name and index columns not needed for boxplot

ax = dataforboxplot.boxplot(figsize=(20, 8), showfliers=False, rot=90) #Get rid of outliers

ax.set_xlabel('Tissue')
ax.set_ylabel('Abundance')
ax.set_title('Tissue Specific Abundance for Genes Enriched in Liver') 
plt.suptitle('')  # Getting rid of pandas-generated boxplot title

plt.show() #displays boxplot

#boxplot = dataforboxplot.boxplot(column=['Adrenal gland', 'Appendix', 'Brain'])
