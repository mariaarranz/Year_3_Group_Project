# -*- coding: utf-8 -*-
"""
Year 3 Group project
"usefulfunc.py"
Useful functions that could be used in any other file

Created on Sun Mar 22 11:34:42 2020
@author: Marianne
"""

import pandas as pd
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt

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
    writer = pd.ExcelWriter(output_file)
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

def revcomp(nucstr):
        """Reverse complement base pairs in a string        
        """
        nucstr = nucstr.upper().replace('U', 'T') #change to standard dNTP format
        #take complement of each dNTP base, change to lower case to avoid replacing back
        comp = nucstr.replace('A','t').replace('T','A').replace('C','g').replace('G','C')
        comp=comp.upper()[::-1] #convert to upper case and revert sequence order
        return comp

def heatmap2d(arr: np.ndarray, colourscale='coolwarm'):
    #MAke a heatmap from a numpy array
    fig,ax=plt.subplots()
    a=ax.imshow(arr, cmap=colourscale)
    fig.colorbar(a,ax=ax)
    plt.show()
    return fig,ax


# Functions for ANOVA
def SumSqTr(treatment_list, N = None, k=None, nis=None):
    """
    Calculate the Sum of Squares for treatment for ANOVA
    treatment_list is a list of numpy arrays containing the data for each treatment,
    N the total number of samples, k the number of different treatments and 
    nis the number of samples per treatment
    If any of the last 3 parameters are not provided, they are automatically calculated from the data
    
    """
    if k==None: #number of different treatments
        k=len(treatment_list)
    if nis==None: #number of samples in each treatment
        nis=[]
        for t in treatment_list:
            nis.append(len(t))
    if N==None or sum(nis)!=N: #total number of samples
        N=sum(nis)
    
    tr_means=[np.mean(t) for t in treatment_list] #means of each treatment
    grand_mean=sum([nis[tr]*tr_means[tr] for tr in range(k)])/N #grand mean, 
            #equivalent to weighted average of treatment means 
            #or to average of all samples together regardless of treatment 
    SSTr=0
    for tr in range(0,k):
        SSTr += nis[tr]*(tr_means[tr]-grand_mean)**2
    return SSTr
        
    
def SumSqEr(treatment_list, N = None, k=None, nis=None):
    """Calculate the Sum of Squares due to errors for ANOVA
    treatment_list is a list of numpy arrays containing the data for each treatment,
    N the total number of samples, k the number of different treatments and 
    nis the number of samples per treatment
    If any of the last 3 parameters are not provided, they are automatically calculated from the data
    
    """
    if k==None:
        k=len(treatment_list)
    if nis==None:
        nis=[]
        for t in treatment_list:
            nis.append(len(t))
    if N==None or sum(nis)!=N:
        N=sum(nis)
    tr_means=[np.mean(t) for t in treatment_list] #means of each treatment
    SSE=0
    for tr in range(0,k):
        for sample in treatment_list[tr]:
            SSE += (sample-tr_means[tr])**2
    return SSE

def SumSqTot(treatment_list, N = None, k=None, nis=None):
    """
    Calculate the total sum of squares for ANOVA
    treatment_list is a list of numpy arrays containing the data for each treatment,
    N the total number of samples, k the number of different treatments and 
    nis the number of samples per treatment
    If any of the last 3 parameters are not provided, they are automatically calculated from the data
    
    """
    if nis==None:
        nis=[]
        for t in treatment_list:
            nis.append(len(t))
    if N==None or sum(nis)!=N:
        N=sum(nis)
    grand_mean=sum([sum(t) for t in treatment_list])/N
    return sum([sum([(sample-grand_mean)**2 for sample in t]) for t in treatment_list])
        
def Eta_Squared(F_stat,k,N):
    """
    Another way of calculating eta squared directly from the f statistic
    and the degrees of freedom (from the number of samples and treatments used)
    """
    return F_stat*(k-1)/(F_stat*(k-1)+N-k)