#!/bin/python3
"""
Year 3 Group project
Retrieve gene ID in table EV1 from Wang et al (2019)
Simplest "dumb" way: 
    even though there maybe  multiple CDSs associated with gene,
    select (arbitrarily) the one from the 1st result

Choose which parts to comment in code depending on how advanced through gene table you are
 
Created on Tue Dec 24 2019 (Merry XMas!)

@author: MarianneBuat
Written with help and review from Martha CRUZ-LESBROS

"""

import pandas as pd
import requests, sys #to communicate with Ensembl servers
import json 


"""
Function definitions
"""
def ensemblGETrequest(ID, server = "https://rest.ensembl.org", partial_ext="", addi_ext = "" ):
    ext = partial_ext+ID+"?"+addi_ext
    while "Repeat the same request":
        try:
            r = requests.get(server+ext, headers={ "Content-Type" : "application/json"},timeout=100)
            r.raise_for_status()
        except (requests.exceptions.HTTPError,requests.exceptions.Timeout) as e: 
            print(e)
            return "error"
        else:
            #print("Successfully decoded")
            decoded=r.json()
            #print(repr(decoded)) #for testing 
            return decoded #return with result if successful

"""
Main code
"""
#Standard format of request to send to Ensembl, just changing ID each time
#Example: https://rest.ensembl.org/sequence/id/ENSG00000000003?species=human;type=cds;multiple_sequences=1
server = "https://rest.ensembl.org"
GETseq_ext="/sequence/id/"
GETseq_addi_ext = "species=human;type=cds;multiple_sequences=1"


#Import data from Excel sheets
#genes=pd.read_excel('Table_EV1.xlsx', sheet_name = 'C. Genes') # Import sheet "C. Genes" from excel file 
#outputDF = pd.DataFrame({'Gene_ID':genes.get('Gene ID'),'transcript_ID':pd.Series([]),'CDS':pd.Series([]), 'Liver_abundance':genes.get('Liver')})
outputDF=pd.read_excel('CDSFromGeneID.xlsx',index_col=0) #import pre-existing output file

# The ID 'ENST00000370378' at index 924 of sheet "C. Genes" is discontinued in the Ensembl database
#dec=ensemblGETrequest('ENST00000370378',server, GETseq_ext,GETseq_addi_ext)
#outputDF.loc[924,'Gene_ID']='ENSG00000189195'
#outputDF.loc[924,'CDS']=dec[0]['seq']
#outputDF.loc[924,'transcript_ID']=dec[0]['id']

indstop=8873 #last successful query or saved output
plante = []
#Loop until reach end
for ind in range(indstop,indstop+1000):
    ID = outputDF.loc[ind,'Gene_ID']
    decoded=ensemblGETrequest(ID,server, GETseq_ext,GETseq_addi_ext) #get sequences from server
    if decoded == "error":
        plante.append(ind)
        print("{}: error for query: {}".format(ind,ID))
    else:
        outputDF.loc[ind,'CDS']=decoded[0]['seq']
        outputDF.loc[ind,'transcript_ID']=decoded[0]['id']
    
        print("{}. Query: {}, ID: {}".format(ind,ID,decoded[0]['id']))
        #print("CDS: {}".format(decoded[0]['seq'])) #for testing
    indstop=ind


#Save to file
writer = pd.ExcelWriter('CDSFromGeneID.xlsx')
outputDF.to_excel(writer)
writer.save()
print('Output written successfully to Excel File.')

################
for ind in range(indstop,indstop+1000):
    ID = outputDF.loc[ind,'Gene_ID']
    decoded=ensemblGETrequest(ID,server, GETseq_ext,GETseq_addi_ext) #get sequences from server
    if decoded == "error":
        plante.append(ind)
        print("{}: error for query: {}".format(ind,ID))
    else:
        outputDF.loc[ind,'CDS']=decoded[0]['seq']
        outputDF.loc[ind,'transcript_ID']=decoded[0]['id']

        print("{}. Query: {}, ID: {}".format(ind,ID,decoded[0]['id']))
        #print("CDS: {}".format(decoded[0]['seq'])) #for testing
    indstop=ind


#Save to file
writer = pd.ExcelWriter('CDSFromGeneID.xlsx')
outputDF.to_excel(writer)
writer.save()
print('Output written successfully to Excel File.')


################
for ind in range(indstop,len(outputDF)):
    ID = outputDF.loc[ind,'Gene_ID']
    decoded=ensemblGETrequest(ID,server, GETseq_ext,GETseq_addi_ext) #get sequences from server
    if decoded == "error":
        plante.append(ind)
        print("{}: error for query: {}".format(ind,ID))
    else:
        outputDF.loc[ind,'CDS']=decoded[0]['seq']
        outputDF.loc[ind,'transcript_ID']=decoded[0]['id']

        print("{}. Query: {}, ID: {}".format(ind,ID,decoded[0]['id']))
        #print("CDS: {}".format(decoded[0]['seq'])) #for testing
    indstop=ind

#Save to file
writer = pd.ExcelWriter('CDSFromGeneID.xlsx')
outputDF.to_excel(writer)
writer.save()
print('Output written successfully to Excel File.')

###############
    
"""
#Save to file
writer = pd.ExcelWriter('CDSFromGeneIDfull.xlsx')
outputDF.to_excel(writer)
writer.save()
print('Full output written successfully to Excel File.')

print('liste des ind plantes')
[print(l) for l in plante]
"""
