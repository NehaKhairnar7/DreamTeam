# -*- coding: utf-8 -*-
"""
Created on Wed Feb 12 17:44:36 2020

@author: sheri
"""

#Import packages
import csv 
import pandas as pd 
import requests 
import re 
import numpy as np
#-----------------------------------------------------------------------------------------------------------------#
#Read the kinase list csv file into a dataframe using pandas (pd)
df = pd.read_csv("clean_human_kinase.csv", index_col=0) 
#create an object that contains everything under the Uniprot Number
identifier = list(df.uniprot_number) 

#create a object containing the gene names
geneName=list(df.gene_name) 

#Need to produce 2 URLs for each protein, that contains the unique kinase 
#identifier number as well as the subcellular location from the database
#As there are sometimes 2 pages of results, need to make sure there is a URL 
#for page 1, and a URL for results on page 2
url1= 'https://www.ebi.ac.uk/QuickGO/services/annotation/search?includeFields=goName&geneProductId='
url2= '&aspect=cellular_component&limit=100&page=1'
url3= '&aspect=cellular_component&limit=100&page=2'


#Create empty list that will contain subcellular location information 
#from quickGO from page 1 of results
quickGODataList=[]

#Create empty list that will contain the subcellular location 
#information from quickGO
quickGODataList2=[]

#create an empty list for any protein names that may not be found using quickGO
errorList=[]

#for each kinase identifier, if no error is produced, merge the 2 urls for 
#page 1 results, separated by the unique kinase identifier name then merge the 
#2 urls for page 2 results
#append the information for each kinase to the lists
for i in identifier: 
    try: 
        url4=url1+i+url2 
        url5=url1+i+url3 
        quickGODataList.append(requests.get(url4).text) 
        quickGODataList2.append(requests.get(url5).text) 
    except:
        errorList.append(i)

#Make 2 empty lists, 1 for the list created from the results of searching the 
#UniprotDataList using regex1, the other empty list is for the results of 
#searching the results from regex1 using regex2
kinaseInfoList=[]
kinaseInfoList2=[]

#Make a regex that will location the cellular component information from the quickGO results
regex1=re.compile(r'"goName":"[A-Za-z]*\,*\s*\-*[a-z]*\,*\:*\-*\s*[A-Za-z]*\<*\-*\s*[A-Za-z]*\,*\s*[A-Za-z]*\-*\,*\s*[A-Za-z]*"')

#For the list from each page, earch the quickGODataList using regex1, append results to KinaseInfoList
for value in quickGODataList: 
    kinaseInfoList.append(regex1.findall(value))

for value in quickGODataList2: 
    kinaseInfoList2.append(regex1.findall(value))    

#create two lists that will contain the strings after unnecessary characters are removed
splitList=[]
splitList2=[]

#For each value in kinaseInfoList, remove extra characters that are not needed
for i in kinaseInfoList:
    i=str(i)
    splitList.append(i.replace('"','').replace('goName','').replace("[]", "").\
                     replace("]","").replace("[","").replace("u':","").\
                     replace("'","").replace(":",""))

#For each value in kinaseInfoList, remove extra characters that are not needed
for j in kinaseInfoList2:
    j=str(j)
    splitList2.append(j.replace('"','').replace('goName','').replace("[]", "")\
                      .replace("]","").replace("[","").replace("u':","").\
                      replace("'","").replace(":",""))

#Make a dictionary containing the gene name, Uniprot Kinase Number identifier 
#information and subcellular location information from each page of results
kinaseDict= {'Gene Name':geneName, 'Uniprot Number':identifier,
             'Subcellular Location1':splitList, 
             "Subcellular Location2":splitList2} 

#Use pandas to make a dataframe from kinaseDict
df=pd.DataFrame(kinaseDict) 

#Replace empty strings 'NaN' with 0
df = df.replace(np.nan, 0)
df['Subcellular Location']=df['Subcellular Location1'].astype(str)+','+ \
df['Subcellular Location2'].astype(str)

#Delete the columns for page 1 and 2
del df['Subcellular Location1']
del df['Subcellular Location2']

#Separate the list within cells of the dataframe so that each subcellular 
#location is a separate row
new_df=(df.set_index(['Gene Name','Uniprot Number'])
   .stack() 
   .str.split(',', expand=True) 
   .stack()
   .unstack(-2) 
   .reset_index(-1, drop=True)
   .reset_index()
)

#Remove any whitespace in the subcellular location column
new_df['Subcellular Location']=new_df['Subcellular Location'].str.strip()

#Tidy subcellular locations column by making all start of words capital letters
new_df['Subcellular Location']=new_df['Subcellular Location'].str.title()

#Drop rows where there are duplications in data
new_df2=new_df.drop_duplicates()

#Remove all rows where there are empty values in the Subcellular location column
final_df = new_df2[new_df2['Subcellular Location'] != '']

#Reindex dataframe
final_df.reset_index()

#Save the dataframe to a csv file
final_df.to_csv('Subcellular_location.csv')