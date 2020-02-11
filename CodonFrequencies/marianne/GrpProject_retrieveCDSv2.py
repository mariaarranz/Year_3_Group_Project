# -*- coding: utf-8 -*-
"""
Year 3 Project: Codons
Code to extract protein ID from table and 
retrieve the cds's from the Ensembl server.

Careful: loops through whole table so only use once on college computer!

Created on Mon Dec  9 17:54:22 2019

@author: Marianne Buat
"""
#import relevant libraries
import pandas as pd
import requests, sys #to communicate with Ensembl servers
import json 

###
def ensemblPostidRequest(id_List, t='cds'):
    #Uses the POST sequence/id of the Ensembl RESt API 
    #to Request multiple types of sequence by a stable identifier list.
    server = "https://rest.ensembl.org"
    ext = "/sequence/id"
    headers={ "Content-Type" : "application/json", "Accept" : "application/json","type":t}
    inputData = json.dumps({'ids':id_List,'type':[t]*len(id_List)}) #dictionary with ids and output type cds for all; 
    
    r = requests.post(server+ext, headers=headers, data=inputData) #request data from server        
    if not r.ok:
      r.raise_for_status()
      sys.exit()
    decoded =r.json() #read response from server 
    return decoded
###

df = pd.read_excel('Table_EV1.xlsx', sheet_name = 'C. Genes') # Import the excel file and call it xls_file

db_len = len(df)
#Create an output dataframe with only the parts we are interested in
geneIDs = df.get('Gene ID')
CDS=pd.Series([]) #Create new series to contain the sequences, will then be added to output dataframe
actualID=pd.Series([])

maxQueryLen = 49 #actually 50
server = "https://rest.ensembl.org"
ext = "/sequence/id"
headers={ "Content-Type" : "application/json", "Accept" : "application/json","type":'cds'}


#Loop through huge table to retrieve all sequences
ind = 0
while ind+maxQueryLen<db_len:
    gene_ids = geneIDs[ind:ind+maxQueryLen] #select 50 IDs from table
    id_list=gene_ids.values.tolist() #convert dataframe slice to list
     #convert to json string for compatibility (could also have used str())
    
    inputData = json.dumps({'ids':id_list,'type':['cds']*len(id_list)}) #dictionary with ids and output type cds for all; 
    
    r = requests.post(server+ext, headers=headers, data=inputData) #request data from server        
    if not r.ok:
      r.raise_for_status()
      sys.exit()
     
    decoded = r.json() #read response from server
    #decoded is a list of dictionaries
    for k in range(len(id_list)):
        CDS[ind+k] = decoded[k]['seq'] #store the retrieved sequences
        actualID[ind+k] = decoded[k]['id']
    ind+=maxQueryLen

#retrieve last ID's not included in loop 
gene_ids = geneIDs[ind:db_len]
id_list=gene_ids.values.tolist()
decoded = ensemblPostidRequest(id_list, t='cds')
for k in range(len(id_list)):
     CDS[ind+k] = decoded[k]['seq'] #store the retrieved sequences
     actualID[ind+k] = decoded[k]['id']
     
#Create an output dataframe for the results
outputDF = pd.DataFrame({'Gene ID':geneIDs,'answer ID':actualID,'CDS':CDS, 'Liver abundance':df['Liver']})
print(outputDF.head())                                                                

#Save to file
writer = pd.ExcelWriter('CDSQueryOutput.xlsx')
outputDF.to_excel(writer)
writer.save()
print('Output is written successfully to Excel File.')

