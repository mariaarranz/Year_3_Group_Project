#!/bin/python3
"""
Year 3 Group project

"getCDSFromGeneID_timedout_20191231.py"

Retrieve gene ID in table EV1 from Wang et al (2019)
Simplest "dumb" way: 
    even though there maybe  multiple CDSs associated with gene,
    select (arbitrarily) the one from the 1st result
This part of the code is to retrieve sequence for IDs that timed out during first run.

Choose which parts to comment in code depending on how advanced through gene table you are
 
Created on Tue Dec 31 2019 (Happy new Year!)

@author: MarianneBuat

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
##Standard format of request to send to Ensembl, just changing ID each time
##Example: https://rest.ensembl.org/sequence/id/ENSG00000000003?species=human;type=cds;multiple_sequences=1
server = "https://rest.ensembl.org"
GETseq_ext="/sequence/id/"
GETseq_addi_ext = "species=human;type=cds;multiple_sequences=1"


##Import data from Excel sheets
outputDF=pd.read_excel('CDSFromGeneID_20191231_1228.xlsx',index_col=0) #import pre-existing output file

##Import text file with indexes of lines previous to timed out IDs
timedout1=pd.read_table("timeout_previous_lineout1.txt", sep="\n", header = None) 
timedout1=(timedout1[0]+1).tolist() # correct to index of for the timedout ID  
timedout2=pd.read_table("timeout_lineout2.txt", sep="\n", header = None) 
timedout2=timedout2[0].tolist()
timedout = timedout1+timedout2


while timedout!=[]: # continue looping as long as have requests that timed out
	plante = [] #to contain index of IDs that errored
	##Loop over all timed out IDs
	for ind in timedout:
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
		
	##Save updated output to file
	writer = pd.ExcelWriter('CDSFromGeneID_20191231_notimedouts.xlsx')
	outputDF.to_excel(writer)
	writer.save()
	print('Output written successfully to Excel File.')
	
	timedout=plante
	
print("No timed out requests")
 
"""
##Save to file
writer = pd.ExcelWriter('CDSFromGeneIDfull.xlsx')
outputDF.to_excel(writer)
writer.save()
print('Full output written successfully to Excel File.')

print('liste des ind plantes')
[print(l) for l in plante]
print(l)
"""
