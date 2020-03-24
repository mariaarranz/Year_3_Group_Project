# -*- coding: utf-8 -*-
"""
Year 3 Group project

"draft_getAnticodonFromProbeName.py"
Draft and fragment testing file for extracting anticodons from a file

Retrieve the anticodon triplet from probe names in tRNA data files from Gingold et al (2014)
This is an example import function to be used and modified analytically for each file 
since each have slightly different formats and data

Input files: 
    (please have them in the same folder as the code, or indicate full path)
    By default 
    

Created on Sun Feb  9 19:03:18 2020

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

def extractAnticodon(df_in, anticodon_col, method,*methodArgs, **methodkwArgs):
    """
    Retrieve the anticodon triplets from identifier columns in the dataframe 
    containing the tRNA data. Returns a 
    (not adapted for reading from row indexes)
    
    Inputs:
        df_in [Pandas DataFrame]:datframe containing the tRNA data and a column 
        with string identifier from which the anticodon is to be extracted
        anticodon_col: expect the column label but also works with column number(0-based)
        method [function]: method to use to screen for the anticodon in a string
        *methodArgs and **methodkwArgs: variable number of arguments to be passed to anticodon extraction 
        method acstrExtractMethod 
        (respectively for positional and keyword arguments)
    Output: 
        (Modified) df_in: input dataframe with anticodons as row indexes
        acs [list]: list of anticodons strings found in each string of the anticodon
        column anticodon_col 
    
    Example:
        df_in loaded from "Normal_Tissues_Only.xlsx", 
        strings containing the anticodons are in the 1st column in the Excel file
        under the header "Anticodon" so pass "Anticodon" or 0 for anticodon_col (name prefered).
        According to the formatting of the string, need to pass a function 
        that will extract 3 letters from letter 3 to 6 (0-based) 
        (pass as lambda function / user-defined function with extra arguments if required)
    """
    #extract column with identifier strings
    #Account for user errors that may provide it as a number
    try:
        d=df_in[anticodon_col]
    except(KeyError): #anticodon column was not passed in expected format, so try with number
        d=df_in.iloc[:,anticodon_col]
        #print("caugth in time!")
        
    acs=d.apply(method,args=methodArgs,**methodkwArgs).tolist()
    df_in.index=acs
    return df_in,acs

def getAnticodonFromProbeName(input_file, output_file = None, anticodon_col=None,anticodon_pos=None, skiprow=None, zero_based = True):
    """
    Description:
    Retrieve the anticodon triplet from probe names in tRNA data files.
    It saves the result in an excel file that contains the same tRNA abundance
    information but with the anticodon as index column
    
    Input:
        input_file [str]: Filename of an Excel file contining the tRNA data 
        with one of the columns with anticodon-related information
        e.g. "Anticodon": AlaTGC-chr12.tRNA8 ("Normal_Tissues_Only.xlsx"),
        "tRNAgene": AlaTGC_chr12_tRNA8 ("Colon_cancer.xlsx")
        (please have them in the same folder as the code, or indicate full path)
        
        Optional parameters:
        output_file [str,None] (default None): 
            Name of the file to save the output data to.
            If output_file is None, then the output to the function is not saved 
            to an Excel file
        
        anticodon_col [string](default None): Indicates on which column of the excel 
        file the anticodon information to extract is
        
        anticodon_pos [non-negative int](default None):
            position of the first letter of the anticodon in a string (rest of the codon assumed to be the 2 following letters)
        zero_based [bool](default True): indicates whether the optional positions anticodon_col
        
            
    
    
            
    Output : a dataframe with the anticodons as index and the tRNA data as columns
    (names of the columns depend on input data)

    Details:
    """
    if zero_based==False:
        
        anticodon_pos-=1
        skiprow-=1
        
    df_in=pd.read_excel(input_file, skiprows = skiprow)
    
    if not anticodon_col==None:
        acstr=df_in.loc[anticodon_col]


    
def checkNon_negInt(my_input):
    """
    Throws a TypeError if my_input is not a non-negative integer      
    """
    if not isinstance(my_input,int):
        raise TypeError("Input is not an int")
    return

def acstrExtract1(idstr,pos):
    return(idstr[pos:pos+3])
    


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


"""
## Here testing space for command lines 
filename = "Normal_Tissues_Only.xlsx"

df_in=pd.read_excel(filename, skiprows = 1)
acs = df_in['Anticodon'].apply(extractAnticodon,args=[3]).tolist()


df_in.index=acs
print(df_in.head())

##Save updated output to file

writer = pd.ExcelWriter('ac'+filename)
df_in.to_excel(writer)
writer.save()
print('Output written successfully to Excel File.')


def testmultipleArgs(i,df_in,*a,**ka):
    print("test: {}".format(testmultipleArgs))
    print("1st positional argument: {}".format(i))
    print(a)
    print(ka)
    return df_in['Anticodon'].apply(acstrExtract1,args=a,**ka).tolist()
  
"""

## Integrate when using functions
input_file = "Normal_Tissues_Only.xlsx"
df_in=uploadFile(input_file,skiprows=1)

#Extract anticodons from info column and save output dataframe with ac as row indexes, 
#save list with anticodons
df_out,aclist = extractAnticodon(df_in, "Anticodon", acstrExtract1,pos=3)

noDup=removeDuplicates(aclist) #remove duplicates from anticodon list, 
                               #left with 32 entries instead of 155
                               #Could also have used list(set(aclist)), 
                               #but it doesn't conserve original anticodon order

saveToFile(df_out,"Output"+input_file)
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