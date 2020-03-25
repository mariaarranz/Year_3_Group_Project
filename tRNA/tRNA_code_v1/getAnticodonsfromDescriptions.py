# -*- coding: utf-8 -*-
"""
Year 3 Group project

"getAnticodonsfromDescriptions.py"
Retrieve the anticodon triplet from probe names in tRNA data files from Gingold et al (2014)

WARNING: may need to go back to output file and check for exceptions 
that may have been analysed incorrectly by the automated script, 
and correct them by hand

Functions included in this file:
    
 
Created on Tue Feb 11 20:52:20 2020

@author: Marianne
"""

import pandas as pd

# ------------------------ Function definitions ------------------------- #


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
    """Save updated output to file
    
    Inputs:
        df_in [pandas DataFrame]: DataFrame to be saved to file
        output_file [str]: name of file to save dataframe to 
    Outputs:
        Prints message to monitor
        Saves output file in folder 
    
    """
    writer = pd.ExcelWriter(output_file)
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
        
    #Extract anticodons using specified method 
    acs=d.apply(method,args=methodArgs,**methodkwArgs).tolist() #extract and convert to list
    df_in.index=acs #change dataframe index
    return df_in,acs

def acstrExtract1(idstr,pos):
    """Anticodon extraction method for "Normal_Tissues_Only.xlsx"
    
    Extract 3 letters after position pos (0-based) of anticodon description string
    Inputs:
        idstr: string to be scanned
        pos: index in string of 1st base of anticodon triplet (0-based)
    Outputs:
        Extracted anticodon 3-letter string 
    """ 
    return(idstr[pos:pos+3])


# -------------------- Integrate for file "Normal_Tissues_Only.xlsx" ----

## Integrate when using functions
input_file = "tRNAs/Normal_Tissues_Only.xlsx"
df_in=uploadFile(input_file,skiprows=1)

#Extract anticodons from info column and save output dataframe with ac as row indexes 
#Ditch list with anticodons if use '_' as 2nd var, save if use named variable
_,aalist = extractAnticodon(df_in, "Anticodon", acstrExtract1,pos=0) #also keep trace of amino acid encoded
df_out,aclist = extractAnticodon(df_in, "Anticodon", acstrExtract1,pos=3)

df_out['amino acid']=aalist
saveToFile(df_out,"tRNAs/tRNA probes with ac and aa.xlsx")

