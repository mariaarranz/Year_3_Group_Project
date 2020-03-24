# -*- coding: utf-8 -*-
"""
Year 3 Group project

"sumtRNAabundances.py"
Add together abundances of tRNAs with same anticodon 
but from different tRNA genes (measured with different probes)

Using tRNA data from Gingold et al (2014)
Here the anticodons have already been extracted from probe description, 
being careful to check for exeptions.

For input file "OutputNormal_Tissues_Only_noSeCMeti.xlsx", the abundances for 
initiator Methionine and for selenocystein encoded by stop TCA have been removed
(iMet could have been kept but the file anticodon modified to differentiate it 
from elongator Met)

After ended analysis, removed probe description columns 
(could have done that before feeding file to input)
Created on Wed Feb 19 16:06:22 2020

@author: Marianne
"""
import pandas as pd
def uploadFile(input_file, index_col=None, skiprows=None):
    """
    Upload input file into dataframe given arguments
    Mainly use to avoid super big function with too many input parameters
    Inputs:
        input_file: [str]: Filename of an Excel file contining the tRNA data 
        with one of the columns with anticodon-related information
        (please have them in the same folder as the code, or indicate full path)
        index_col [int]: Column (0-indexed) to use as row labels
        skiprows [list-like, int]: Rows to skip at the beginning (0-indexed)
    """
    df=pd.read_excel(input_file,index_col=index_col,skiprows = skiprows)
    return df

def saveToFile(df_in,output_file):
    ##Save updated output to file
    writer = pd.ExcelWriter('ac'+output_file)
    df_in.to_excel(writer)
    writer.save()
    print('Output written successfully to Excel File.')
    return

def removeDuplicates(listofElements):
    '''Remove duplicate elements from list
    
    Input : list to be screened for duplicates
    
    Alternatively, use list(set(aclist))
    but it doesn't conserve original list order
    
    Code from "Python : How to Remove Duplicates from a List" 
    Posted by Varun, April 20, 2018 on thispointer.com
    Available at:
        https://thispointer.com/python-how-to-remove-duplicates-from-a-list/
        [Accessed 12/02/2020]
    '''
    # Create an empty list to store unique elements
    uniqueList = []
    
    # Iterate over the original list and for each element
    # add it to uniqueList, if its not already there.
    for elem in listofElements:
        if elem not in uniqueList:
            uniqueList.append(elem)
    
    # Return the list of unique elements        
    return uniqueList


#Manually correct output data in Excel
#Remove Selenocysteine encoded by stop codon and initiator Methionine tRNAs
df_corrected=uploadFile("OutputNormal_Tissues_Only_noSeCMeti.xlsx", index_col=0)
acs=list(df_corrected.index)
acs=removeDuplicates(acs)
s_list=[]
for ac in acs:
    #use trick to avoid errors when there is a single entry 
    #and indexing dataframe returns a Pandas Series
    s=pd.DataFrame()
    s=s.append(df_corrected.loc[ac]) 
    s=s.sum(axis = 0, skipna = True) #sum the tRNA abundance values 
                                       
    s_list.append(s)#add summed values for the anticodon to dataframe
                          # (only select columns with abundance, not descriptors!) 
df_summed=pd.DataFrame(s_list, columns=s_list[0].index, index=acs)
saveToFile(df_summed, "SummedtRNAab"+input_file)