# -*- coding: utf-8 -*-
"""
Created on Tue Dec 24 18:59:21 2019

@author: Marianne
"""

#import relevant libraries
import pandas as pd
import requests, sys #to communicate with Ensembl servers
import json 

"""
Function definitions
"""
def ensemblOverlapTranslation(protID, server = "https://rest.ensembl.org", partial_ext = "/overlap/translation/"):
    ext = partial_ext+protID+"?species=human"
 
    r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
    if not r.ok:
        r.raise_for_status()
        sys.exit()
    decoded=r.json()
    #print(repr(decoded)) #for testing 
    return decoded

def ensemblGETrequest(ID, server = "https://rest.ensembl.org", partial_ext="", addi_ext = "" ):
    ext = partial_ext+ID+"?"+addi_ext
    while True:    
        r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
        try:
            r.raise_for_status()
        except requests.exceptions.HTTPError as e: 
            print(e)
            print("I will ask you for input")
            n=float(input("Enter 1 to repeat request, anything else to exit program: "))
            if n==1:
                continue
            else:
                sys.exit()  
        #print("Successfully decoded")
        decoded=r.json()
        #print(repr(decoded)) #for testing 
        return decoded #return with result if successful
    
server = "https://rest.ensembl.org"
ov_trans_ext = "/overlap/translation/"
ov_trans_addi_ext="species=human"

seqGET_ext="/sequence/id/"
seqGETaddi_ext = "species=human;type=cds;multiple_sequences=1"
params = {'species':'human','type':'cds','multiple_sequences':1}

decID1=ensemblOverlapTranslation("ENSP00000229794")
ID1=decID1[0]['Parent']
decID2=ensemblOverlapTranslation("ENSP00000229795")
ID2=decID2[0]['Parent']

decprotID1=ensemblGETrequest(ID1,partial_ext=seqGET_ext,addi_ext=seqGETaddi_ext)
seq1=decprotID1[0]['seq']
decprotID2=ensemblGETrequest(ID2,partial_ext=seqGET_ext,addi_ext=seqGETaddi_ext)
seq2=decprotID2[0]['seq']

decgeneID=ensemblGETrequest(ID='ENSG00000112062',partial_ext=seqGET_ext,addi_ext=seqGETaddi_ext)

count=0
for ans in decgeneID:
    if ans['id']==ID1:
        print("Matching ID {}, index is {}".format(ID1,count) )
    elif ans['id']==ID2:
        print("Matching ID {}, index is {}".format(ID2,count) )    
    count+=1

decoded=ensemblGETrequest('ENSG00000069712',server, seqGET_ext,seqGETaddi_ext)
