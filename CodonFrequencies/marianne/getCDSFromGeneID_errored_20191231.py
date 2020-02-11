#!/bin/python3
"""
Year 3 Group project

"getCDSFromGeneID_errored_20191231.py"

Retrieve gene ID in table EV1 from Wang et al (2019)
Simplest "dumb" way: 
    even though there maybe  multiple CDSs associated with gene,
    select (arbitrarily) the one from the 1st result
This part of the code is to retrieve sequence for IDs that that caused an HTTP error during first run.

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
server = "https://rest.ensembl.org"
ov_trans_ext = "/overlap/translation/"
ov_trans_addi_ext="species=human"
GETseq_ext="/sequence/id/"
GETseq_addi_ext = "species=human;type=cds;multiple_sequences=1"

##Import data from Excel sheets
protGrps=pd.read_excel('Table_EV1.xlsx', sheet_name = 'B. Protein groups') # Import sheet "B. Protein Groups" from excel file 
key=[l.split(';') for l in protGrps['Gene ID (Ensembl)']]
key =[key[l][0] for l in range(len(key))]
protGrps=protGrps.set_index(pd.Series(key))
outputDF=pd.read_excel('CDSFromGeneID_20191231_notimedouts.xlsx',index_col=0) #import pre-existing output file

error1=pd.read_table("Error_previous_lineout1.txt", sep="\n", header = None) 
error1=(error1[0]+1).tolist() # correct to index of for the error ID  
error2=pd.read_table("Error_lineout2.txt", sep="\n", header = None) 
error2=error2[0].tolist()
error = error1+error2

trPlante=[]
seqPlante=[]
for ind in error:
    gID=outputDF.loc[ind,'Gene_ID']
    protIDstr = protGrps.loc[gID,'Protein ID (representatives)']
    protIDlist = protIDstr.split(';') #separate continuous ID string into list of IDs
    protID=protIDlist[0]
    decoded=ensemblGETrequest(protID,server,ov_trans_ext,ov_trans_addi_ext)
    if decoded == "error":
        trPlante.append(ind)
        print("{}: error for query: {}".format(ind,gID))
    else:
        transcriptID=decoded[0]['Parent'] #get transcript ID from any of the list elements in answer
        decoded=ensemblGETrequest(transcriptID,server, GETseq_ext,GETseq_addi_ext) #get sequences from server
        if decoded == "error":
            seqPlante.append(ind)
            print("{}: error for query: {}".format(ind,transcriptID))
        else:
            outputDF.loc[ind,'CDS']=decoded[0]['seq']
            outputDF.loc[ind,'transcript_ID']=decoded[0]['id']
            print("{}. Query: {}, ID: {}".format(ind,protID,decoded[0]['id']))

##Save updated output to file
writer = pd.ExcelWriter('CDSFromGeneID_20191231_noerrors.xlsx')
outputDF.to_excel(writer)
writer.save()
print('Output written successfully to Excel File.')



print('liste des ind plantes apres translation overlap:')
print(trPlante)
[print(l) for l in trPlante]
print('liste des ind plantes apres translation overlap:')
print(seqPlante)
[print(l) for l in seqPlante]


###############
##Be careful to decomment some parts if needed

removeIDs=[trPlante[0],trPlante[1],trPlante[4]]
trIDnew=[[trPlante[2],'ENSG00000283787','ENST00000640310'],[trPlante[3],'ENSG00000256806','ENST00000542475'],[13637,'ENSG00000284691','ENST00000498682'],[trPlante[-1],'ENSG00000171862','ENST00000371953']]
#outputDF.drop(removeIDs) #!!! Run only once!
plante = []
for t in trIDnew:
    outputDF.loc[t[0],'Gene_ID']=t[1]
    decoded=ensemblGETrequest(t[2],server, GETseq_ext,GETseq_addi_ext) #get sequences from server
    if decoded == "error":
            seqPlante.append(ind)
            print("{}: error for query: {}".format(ind,transcriptID))
    else:
        outputDF.loc[t[0],'CDS']=decoded[0]['seq']
        outputDF.loc[t[0],'transcript_ID']=decoded[0]['id']
        print("{}. Query: {}, ID: {}".format(t[0],t[2],decoded[0]['id']))

writer = pd.ExcelWriter('CDSFromGeneID_20191231_noerrors2.xlsx')
outputDF.to_excel(writer)
writer.save()
print('Output written successfully to Excel File.')


