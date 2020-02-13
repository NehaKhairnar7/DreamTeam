# -*- coding: utf-8 -*-
"""
Created on Fri Jan 31 18:37:16 2020

@author: sheri
"""
###User Data Analysis Script

#------------------------------------------------------------------------------------------------------------------

###Import packages into environment
from Database.kinase_functions import *
from sqlalchemy import create_engine, or_, and_
from sqlalchemy.orm import sessionmaker
from pprint import pprint
import csv 
import pandas as pd 
import re
import numpy as np
import math
from scipy.stats import norm
from bokeh.models import Span
from bokeh.resources import CDN
from bokeh.embed import file_html, components
from bokeh.plotting import figure, ColumnDataSource, show, output_file
from bokeh.models import HoverTool, WheelZoomTool, PanTool, BoxZoomTool, ResetTool, TapTool, SaveTool

#------------------------------------------------------------------------------------------------------------------

### Function that reads user data, filters the dataframe, matches the substrate
### with kinases from the database,and then carries out data analysis. 
def data_analysis(filename, p_val, CV, Sub):
    '''
    1. Read user data file, make sure there are only 7 columns and if CV 
       columns are missing, add two columns for CV, then rename all columns
    2. Clean and filter user data
    3. Take -log10 of p value and log2 of fold change and pass into 2 
       new columns
    4. Add column for kinases from database based on substrate 
    5. Filter by CV threshold set by user
    6. Carry out KSEA calculations
    :param filename: user data table
    :param p_val: user p value threshold
    :param CV: user coefficient of variance threshold
    :param Sub: user minimum number of substrates for kinase enrichment plot 
    :return final_substrate: Filtered dataframe of the input file.
    :return df_final3: Dataframe of kinase-substrate relationships and the 
     corresponding data from the input file.
    :return calculations_df: A dataframe of the data from kinase substrate 
     enrichment analysis calculations
    '''
    

    #Read user data file 
    df_input_original = pd.read_csv("instance/Data_Upload/"+ filename, \
                                    sep='\t')
    
    #Take first 7 columns from the user data file
    input_original_subset = df_input_original.iloc[:, 0:7]
    
    #Find number of columns in user data file
    col_number =  input_original_subset.shape[1]
    
    #If number of columns is 5, append CV columns as 0
    if col_number == 5:
        input_original_subset["control_cv"] = 0
        input_original_subset["condition_cv"] = 0  
    
    #Rename columns to standardised names
    df_cols=["Substrate", "control_mean", "inhibitor_mean", "fold_change", \
             "p_value", "ctrlCV", "treatCV" ]
    
    #Append renamed columns
    input_original_subset.columns=df_cols
    
    #Make the substrate column a string
    input_original_subset.iloc[:, 0] = \
    input_original_subset.iloc[:, 0].astype(str)

    #Drop all rows with 'None' in substrate for phosphosite
    input_original_subset = \
    (input_original_subset[~input_original_subset['Substrate'].\
                           str.contains("None")])

    #Make regex that will search for methionines in substrate column
    methionine_reg = r"\([M]\d+\)" 

    #Remove rows with methionine in the substrate column
    input_original_subset=(input_original_subset\
                           [~input_original_subset.Substrate.str.\
                            contains(methionine_reg)]) 
    
    #For CV columns with missing CV values, add a CV of 3
    input_original_subset['ctrlCV'] = input_original_subset['ctrlCV'].\
    replace(np.nan, 3)
    input_original_subset["treatCV"]=input_original_subset["treatCV"].\
    replace(np.nan, 3)
    
    #Make columns 2 to 7 float instead of string
    input_original_subset.iloc[:, 1:7] = input_original_subset.iloc[:, 1:7].\
    astype(float)

    #Need to separate the phosphosite from the substrate in the first column 
    #into 2 separate columns
    input_original_subset[['Substrate','Phosphosite']] = \
    input_original_subset.Substrate.str.split('\(|\)', expand=True).\
    iloc[:,[0,1]]

    #Drop rows where all columns are empty
    input_original_subset=input_original_subset.dropna(axis=1, how="all")
   
    #Take -log10 of the corrected p-value.
    uncorrected_p_values=input_original_subset.iloc[ :,4].astype(np.float64)
    log10_corrected_pvalue = (-np.log10(uncorrected_p_values))

    #Append -log10(P-values) to a new column in data frame.
    input_original_subset["-Log10 Corrected P-Value"]=log10_corrected_pvalue
    NegLog10Kinase=input_original_subset
    
    #Replace inf and -inf fold changes with 0
    NegLog10Kinase.loc[:, "fold_change"].replace([np.inf, -np.inf], 0,
                      inplace=True)
    
    #Take log2 of fold change column and append to new column
    log2FC=np.log2(NegLog10Kinase.iloc[:, 3])
    NegLog10Kinase["Log2 Fold Change"]=log2FC
    log2FCKinase = NegLog10Kinase
    
    #Replace Log2 fold change of 0 that is infinite with 0
    log2FCKinase.loc[:, "Log2 Fold Change"].replace([np.inf, -np.inf], 0,
                    inplace=True)
    final_substrate=log2FCKinase
   
    # Replace nan with 0.
    log2FCKinase.loc[:, "Log2 Fold Change"] =log2FCKinase.\
    loc[:,"Log2 Fold Change"].fillna(0)
    
    #Take each substrate and phosphosite from each row and append to a list
    Final_substrate=log2FCKinase
    Sub_phosp_list=[]
    for i, j, k in zip(log2FCKinase['Substrate'], 
                       log2FCKinase['Phosphosite'],
                       range(len(log2FCKinase))):
        Sub_phosp_list.append([])
        Sub_phosp_list[k].append(i)
        Sub_phosp_list[k].append(j)
    
    #Extract kinases from database using substrates and phosphosites from 
    #user input
    #Searches substrate in user input to see if it is in the database
    #if substrate is in the database, then append to list
    #Otherwise, append an empty list
    KinaseList=[]
    substrates={x: None for x in get_all_substrates_complete()}
    for i in Sub_phosp_list:
        Sub = i[0]
        Pho = i[1]
        if Sub in substrates:
            KinaseList.append(get_kinase_substrate_phosphosite(Sub, Pho))
        else:
            KinaseList.append([])
    
    #Make a new column for kinase        
    input_original_subset['Kinase']=KinaseList
    
    #As function to query kinases returns a dictionary of substrate, 
    #phosphosite and kinases, need to remove substrates and phosphosites
    #First, make each key of dictionary a column
    df_final = pd.concat([input_original_subset, input_original_subset\
                          ['Kinase'].apply(pd.Series)], axis = 1).\
        drop('Kinase', axis = 1)
    
    #Drop substrate and phosphosite columns
    df_final1=df_final.drop(['substrate', 'phosphosite'], axis=1)
    
    #Drop columns with empty rows
    df_final2=df_final1.dropna()
    
    #As for each substrate there is a sometimes a list of kinases, need to 
    #separate kinases into separate rows
    df_final2=df_final2.explode('kinase')
    
    #Drop rows where there is an empty row in kinase column
    df_final3=df_final2.dropna(subset = ["kinase"])
    
    #Convert CV to float
    CV=float(CV)
    
    #Filter rows by user inputted CV threshold
    df_final3= df_final3[(df_final3['ctrlCV'] <= CV) & 
                         (df_final3['treatCV'] <= CV)]  
    
    #Carry out KSEA calculations
    #mS = the mean fold change of a substrate set of a kinase
    #mP = the mean fold change of the entire data set
    #Delta = the standard deviation of the fold change of the dataset
    #m = the number of substrates in a substrate set for a kinase
    mS = df_final3.groupby('kinase')['Log2 Fold Change'].mean()
    mP = df_final3['Log2 Fold Change'].mean()
    delta=df_final3['Log2 Fold Change'].std()
    
    #for m, group phosphosites by kinase and find the length of each 
    #phosphosite set and append this length to new list
    m=[]
    Kinase_phosphosite=df_final3.groupby('kinase')['Phosphosite']
    for key, item in Kinase_phosphosite:
        m.append(len(item))
        
    #Calculate Z Scores
    Z_Scores=[]    
    for i, j in zip(mS, m):
        Z_Scores.append(((i-mP)*math.sqrt(j))/delta)
        
    #Calculate p values from z scores
    p_means=[]
    for i in Z_Scores:
        p_means.append(norm.sf(abs(i)))
    
    #Make each component of KSEA a different key and value in dictionary
    calculations_dict={'mS': mS, 'mP':mP, 'm':m, 'Delta':delta,
                       'Z_Scores':Z_Scores,"P_value":p_means}
    
    #Make calculations dictionary a dataframe
    calculations_df=pd.DataFrame(calculations_dict)
    
    #Reset dataframe index 
    calculations_df=calculations_df.reset_index(level=['kinase'])
    
    #Drop the kinase column of the final_substrate dataframe as it is empty 
    #and not needed for the sustrate volcano plot
    final_substrate=final_substrate.drop(['Kinase'], axis=1)
    
    #Return output
    return (calculations_df, final_substrate ,df_final3) 

#Make volcano plot of substrate data before matching with kinase
def VolcanoPlot_Sub(final_substrate, p_val, FC, CV):
    ''' 
    1. change user parameters to floats
    2. Add a new column for point colours, if point is above the -log(p value) 
       threshold  and below the -fold change threshold, colour red, if point is
       above -log(pvalue)threshold but above fold change threshold colour 
       green, otherwise keep grey .
    3. Add another new column for regulation status, based on same 
       thresholds as above.
    4. Use bokeh to make volcano plot of data
    5. Convert to html
    :param final_substrate:filtered dataframe of user data input
    :param p_val: user p value threshold
    :param FC: user fold change threshold 
    :param CV: user coefficient of variance threshold
    :return html: returns html of plot for embedding
    '''
    #Make user defined thresholds floats
    FC=float(FC)
    FC_N = -(float(FC))
    PV=-np.log10(float(p_val))
    
    #Add column 'colour' for colour points will be, specify conditions for 
    #each point colour
    final_substrate.loc[(final_substrate['Log2 Fold Change'] > FC) & 
                        (final_substrate['-Log10 Corrected P-Value'] > PV),
                        'colour' ] = "Green"  # upregulated
    
    final_substrate.loc[(final_substrate['Log2 Fold Change'] < FC_N) & 
                        (final_substrate['-Log10 Corrected P-Value'] > PV),
                        'colour' ] = "Red"   # downregulated
   
    final_substrate['colour'].fillna('grey', inplace=True)

    #add coloumn 'regulation' for regulation status, 'more abundant' will have 
    #same thresholds as 'green', 'less abundant' have same thresholds as red, 
    #and no change the same thresholds as grey.
    final_substrate.loc[(final_substrate['Log2 Fold Change'] > FC) & 
                        (final_substrate['-Log10 Corrected P-Value'] > PV),
                        'regulation' ] = "More abundant in treatment"  
    
    final_substrate.loc[(final_substrate['Log2 Fold Change'] < FC_N) & 
                        (final_substrate['-Log10 Corrected P-Value'] > PV), 
                        'regulation' ] = "Less abundant in treatment"   
    
    final_substrate['regulation'].fillna('No change', inplace=True)
    
    #set category for plot to substrates
    category = 'Substrate'
    
    #make category items each unique substrate
    category_items = final_substrate[category].unique()
    title="Summary of Proteins in Sample"

    #Make the source that final_sustrate dataframe
    source = ColumnDataSource(final_substrate)

    #Make hovertool, state what data from the dataframe will show when a 
    #point is hovered over
    hover = HoverTool(tooltips=[('Substrate', '@Substrate'),
                                ('Phosphosite', '@Phosphosite'),
                                ('Fold_change', '@{Log2 Fold Change}'),
                                ('p_value', '@{-Log10 Corrected P-Value}')])
    
    #Specify what other tools will be implemented in the plot
    tools = [hover, WheelZoomTool(), PanTool(), BoxZoomTool(), ResetTool(),
             SaveTool()]
    
    #make plot
    p = figure(tools=tools,title=title, plot_width=700,plot_height=500,
               toolbar_location='right',toolbar_sticky=False)
    
    #Add the datapoints using the colour from the dataframe and add legend w
    #ith regulation status
    p.circle(x = 'Log2 Fold Change', y = '-Log10 Corrected P-Value',
             source=source,size=10,color='colour', legend='regulation')
    
    #Make the threshold lines using user set parameters
    p_sig = Span(location=PV,dimension='width', line_color='black', 
                 line_dash='dashed', line_width=3)
    
    fold_sig_over = Span(location=FC,dimension='height', line_color='black',
                         line_dash='dashed', line_width=3)
    
    fold_sig_under=Span(location=FC_N,dimension='height', line_color='black',
                        line_dash='dashed', line_width=3)
    
    #Add the threshold lines 
    p.add_layout(p_sig)   
    p.add_layout(fold_sig_over)   
    p.add_layout(fold_sig_under)   
    
    #Add axis labels
    p.xaxis.axis_label = "Log2 Fold Change"
    p.yaxis.axis_label = "-Log10 Corrected P-Value"
            
    #output plot as html
    html=file_html(p, CDN, "Volcano Plot of Substrates" )
    
    #return plot as html
    return html

#Make kinase-substrate volcano plot
def VolcanoPlot(df_final3, p_val, FC, CV):
    ''' 
    1. change user parameters to floats
    2. Add a new column for point colours, if point is above the -log(p value) 
       threshold  and below the -fold change threshold, colour red, if point is
       above -log(pvalue)threshold but above fold change threshold colour 
       green, otherwise keep grey .
    3. Add another new column for regulation status, based on same 
       thresholds as above.
    4. Use bokeh to make volcano plot of data
    5. Convert to html
    :param df_final3: dataframe containing kinase-substrate data
    :param p_val: user p value threshold
    :param FC: user fold change threshold 
    :param CV: user coefficient of variance threshold
    :return html: returns html of plot for embedding
    '''
    
    #Make user defined thresholds floats
    FC=float(FC)
    FC_N=-(float(FC))
    PV=-np.log10(float(p_val))

    #Add column 'colour' for colour points will be, specify conditions 
    #for each point colour
    df_final3.loc[(df_final3['Log2 Fold Change'] > FC) & 
                  (df_final3['-Log10 Corrected P-Value'] > PV), 'colour' ] \
                  = "Green"
    df_final3.loc[(df_final3['Log2 Fold Change'] < FC_N)\
                  & (df_final3['-Log10 Corrected P-Value'] > PV), 'colour' ] \
                  = "Red"   # downregulated
    df_final3['colour'].fillna('grey', inplace=True)

    #add coloumn 'regulation' for regulation status, 'upregulated' will have 
    #same thresholds as 'green', 'downregulated' have same thresholds as red, 
    #and no change the same thresholds as grey.
    df_final3.loc[(df_final3['Log2 Fold Change'] > FC) \
                  & (df_final3['-Log10 Corrected P-Value'] > PV), 
                  'regulation' ] = "Upregulated"  
    df_final3.loc[(df_final3['Log2 Fold Change'] < FC_N) \
                  & (df_final3['-Log10 Corrected P-Value'] > PV),
                  'regulation' ] = "Downregulated"  
    df_final3['regulation'].fillna('No Change', inplace=True)
    
    #set category as substrates
    category = 'kinase'
    
    #set category items as each individual substrate
    category_items =df_final3[category].unique()
    
    #make title
    title="Kinase Activity Summary"
    
    #set df_final3 dataframe as data source
    source = ColumnDataSource(df_final3)
    
    #set hovertool data that will be shown when data point is hovered over
    hover = HoverTool(tooltips=[('Kinase','@kinase'),
                                ('Substrate', '@Substrate'),
                                ('Phosphosite', '@Phosphosite'),
                                ('Fold_change', '@{Log2 Fold Change}'),
                                ('p_value', '@{-Log10 Corrected P-Value}')])
            
    #set hover tools to be used in plot
    tools = [hover, WheelZoomTool(), PanTool(), BoxZoomTool(), ResetTool(), 
             SaveTool()]
    
    #make plot
    p = figure(tools=tools,title=title,plot_width=700,plot_height=500,
               toolbar_location='right', toolbar_sticky=False)
   
    #add points with colours corresponding to regulation status and add 
    #legend with regulation status
    p.circle(x = 'Log2 Fold Change', y = '-Log10 Corrected P-Value',
             source=source,size=10,color='colour', legend= 'regulation')
    
    #make threshold lines for different user parameters
    p_sig = Span(location=PV,dimension='width', line_color='black', 
                 line_dash='dashed', line_width=3)
    
    fold_sig_over=Span(location=FC,dimension='height', line_color='black',
                       line_dash='dashed', line_width=3)
    
    fold_sig_under=Span(location=FC_N,dimension='height', line_color='black',
                        line_dash='dashed', line_width=3)
    
    #add threshold lines
    p.add_layout(p_sig)   
    p.add_layout(fold_sig_over)   
    p.add_layout(fold_sig_under)  
        
    #Add axis labels
    p.xaxis.axis_label = "Log2 Fold Change"
    p.yaxis.axis_label = "-Log10 Corrected P-Value"
    
    #make html of plot
    html=file_html(p, CDN, "Volcano Plot of Filtered Kinases" )
    
    #return plot as html
    return html

def EnrichmentPlot(calculations_df, p_val, FC, CV, Sub):
    '''
    1. Make filtred dataframe from user defined minimum substrate number 
    2. Sort dataframe by z scores
    3. Use bokeh to create enrichment plot from calculations dataframe
    4. output plot as html
    :param calculations_df: dataframe containing data obtained during KSEA 
     calculations for each kinase
    :param p_val: user p value threshold
    :param FC: user defined fold change threshold
    :param CV: user coefficient of variance threshold
    :param Sub: user minimum number of substrates for kinase enrichment plot 
    :return html: html of enrichment plot
    '''
    
    #filter calculations_df by minimum substrate number as specified by user
    reduc_calculations_df=calculations_df[calculations_df['m']>= float(Sub)]
    
    #sort dataframe by z-scores
    reduc_calculations_df=reduc_calculations_df.sort_values(by='Z_Scores')

    #Add new columns for bar colour
    #orange if P value is below p value threshold, black if above p value threshold
    reduc_calculations_df.loc[(reduc_calculations_df['P_value'] < \
                               float(p_val)), 'colour'] = "Orange"  
    reduc_calculations_df.loc[(reduc_calculations_df['P_value'] > \
                               float(p_val)), 'colour' ] = "Black"
    
    #obtain kinases from dataframe and save as object
    kinase=reduc_calculations_df['kinase']
    
    #obtain z scores and save as object
    z_score=reduc_calculations_df['Z_Scores']
    
    #set source as the filtered dataframe
    source = ColumnDataSource(reduc_calculations_df)
    
    #define what hovertool will show when points are hovered over in plot
    hover = HoverTool(tooltips=[('Z-Score','@Z_Scores'),
                                ('Number of Substrates', '@m'),
                                ('P-value', '@P_value'),
                                ('Kinase', '@kinase')])

    #define what tools will be implemented in plot
    tools = [hover, WheelZoomTool(), PanTool(), BoxZoomTool(), ResetTool(), 
             SaveTool()]
    
    #make plot
    p = figure(tools=tools, y_range=kinase, 
               x_range=((z_score.min()-5), (z_score.max()+5)), 
               plot_width=600, plot_height=800, toolbar_location=None,
           title="Kinase Substrate Enrichment")
    
    #draw the bars, with the length of the bars corresponding the z scores and 
    #colour coresponding the the colour
    p.hbar(y="kinase", left=0, right='Z_Scores', height=0.3, color= 'colour', 
           source=source)
    
    #make axis label s
    p.xaxis.axis_label = "Z-Score)"
    p.yaxis.axis_label = "Kinase"
    
    #remove outlines
    p.outline_line_color = None
    p.ygrid.grid_line_color = None
    
    #convert plot to html
    html=file_html(p, CDN, "Ratio of Enrichment" )
    
    #return html of plot
    return html

def df2_html(calculations_df):
    '''
    1.convert the calculations_df dataframe to html
    :param calculations_df: the calculations dataframe created using the 
     KSEA calculations data
    :return df_final_html: the html for the calculations dataframe
    '''
    
    #convert the caluclations_df dataframe to html
    df_final_html=calculations_df.to_html()
    
    #return the html of dataframe
    return df_final_html