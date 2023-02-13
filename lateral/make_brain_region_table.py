#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 12 10:13:27 2022

@author: simonyj
"""

#%% Import tables used to make brain plots 

import pandas as pd

effect_sizes_path = "/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/parcels/lateral/ANOVA+contrast_effect_sizes.xlsx"
contrast_fdr_table = "/mnt/projects/VIA11/FREESURFER/Stats/Plots/lateral/lateral_combined_contrasts.xlsx" #this file is made from freesurfer_analysis_lateral.R and then freesurfer_analysis_lateral_plot.R

df1 = pd.read_excel(effect_sizes_path) 
df2 = pd.read_excel(contrast_fdr_table) 

df2['sex'] = df2['sex'].map({0:'female', 1:'male', 2:'both'})


new_df = pd.merge(df1, df2,  how='inner', left_on=["model_yvar","measure","hemisphere","sex","globvar_in_model"], right_on = ["Model_yvar","measure","hemisphere","sex","global_var_in_model"])

new_df = new_df[["model_yvar","measure","hemisphere","sex","global_var_in_model",'DOF_ttest',"tratio_BP-K","tratio_SZ-K","pval_BP-K","pval_SZ-K","fdr_pvals_BP-K","fdr_pvals_SZ-K","CohensD_BP-K","CohensD_SZ-K"]]
#df.year.astype(Int64)
#new_df['global_var_in_model'] = new_df['global_var_in_model'].map({0:'no', 1:'yes'})
new_df['global_var_in_model'] = new_df['global_var_in_model'].map({0:'without', 1:'with'})



#%% Use and reorder imported tables to make final tables for appendix

measures = ["volume", "area", "thickness"]
sexes = ["female", "male"]
glob = ["with","without"]
contrasts = ["BP-K","SZ-K"]


#extract the data for each table 
m = "volume"        #"volume", "area", "thickness"
s = "female"        #"female", "male"
g = "with"          #"with", "without"
c = "BP-K"   #"BP-K","SZ-K"

table_text = f"{s} {c} brain {m} difference from control from models {g} global covariate"

total_df = pd.DataFrame()

for m in measures:
    for s in sexes:
        for g in glob:
            for c in contrasts:
                
                excel_name = f"{s}_{c}_brain_{m}_models_{g}_glob.xlsx"
                
                
                df = new_df.loc[ (new_df["measure"].isin([m])) & (new_df["global_var_in_model"].isin([g])) & (new_df["sex"].isin([s]))]
                df = df.drop(df.filter(regex=c).columns, axis=1)
                df = df.drop(df.filter(regex=m).columns, axis=1)
                df = df.drop(df.filter(regex=s).columns, axis=1)
                df = df.drop(df.filter(regex=g).columns, axis=1)
                
                df["model_yvar"] = df["model_yvar"].str.replace(f"_{m}","")
                
                final_df = df 
                
                final_df["model_yvar"] = final_df["model_yvar"].str.replace("lh_","")
                final_df["model_yvar"] = final_df["model_yvar"].str.replace("rh_","")
                
                #reorder 
                names = list(final_df['model_yvar'])[:34] #names = list(set(list(final_df['model_yvar'])))
                hemi_order = ['lh','rh']
                final_df = final_df.astype({'model_yvar': pd.CategoricalDtype(names, ordered=True),
                                            'hemisphere': pd.CategoricalDtype(hemi_order, ordered=True)})
                final_df = final_df.sort_values(['model_yvar', 'hemisphere'])
                
                
                #final touches 
                final_df = final_df.drop(['measure','sex','global_var_in_model'],axis=1)
                final_df = final_df.round(3)
                final_df.columns = final_df.columns.str.rsplit('_').str.get(0)  # strip suffix at the right end only.
                
                #final_df.to_excel(f"/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/parcels/lateral/combined_tables/{excel_name}",index=False) 
                
                # add some extra information to final_df before appending it to total_df 
                final_df['measure'] = m
                final_df['sex'] = s
                final_df['glob'] = g
                final_df['contrast'] = c
                
                total_df = total_df.append(final_df, ignore_index=True)

#total_df = pd.DataFrame(total_list)
total_df.to_excel("/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/parcels/lateral/combined_tables/all_tables_combined.xlsx",index=False) 


#%% Add info & clean raw latex table code 
"""
df_lh = df.loc[ df["hemisphere"].isin(["lh"]) ]
df_rh = df.loc[ df["hemisphere"].isin(["rh"]) ]

df_lh["model_yvar"] = df_lh["model_yvar"].str.replace("lh_","")
df_rh["model_yvar"] = df_rh["model_yvar"].str.replace("rh_","")


keep_same = {"model_yvar","measure","hemisphere","sex","global_var_in_model"}
df_lh.columns = ['{}{}'.format(c, '' if c in keep_same else '_lh') for c in df_lh.columns]
df_rh.columns = ['{}{}'.format(c, '' if c in keep_same else '_rh') for c in df_rh.columns]


final_df = pd.merge(df_lh, df_rh,  how='left', left_on=["model_yvar","measure","sex","global_var_in_model"],right_on=["model_yvar","measure","sex","global_var_in_model"])

final_df = final_df.drop(['measure','sex','global_var_in_model','hemisphere_x','hemisphere_y'],axis=1)
final_df = final_df.round(3)
final_df.columns = final_df.columns.str.rsplit('_').str.get(0)  # strip suffix at the right end only.


latex_code_raw = final_df.to_latex(index=False)

sep = '\n'

latex_code = sep.join(['\\begin{table}[h]',latex_code_raw,'\caption{',table_text,'}','\end{table}'])
latex_code = latex_code.replace('\\toprule','&\multicolumn{5}{c|}{Left Hemisphere}& \multicolumn{5}{|c}{Right Hemisphere} \\\ \hline')
latex_code = latex_code.replace('\\midrule','\\hline')
latex_code = latex_code.replace('model','region')
latex_code = latex_code.replace('{lr','{l|r')
latex_code = latex_code.replace('rrrrrrrr','rrrrr|rrr')

latex_code = latex_code.replace('DOF','df')
latex_code = latex_code.replace('tratio','t')
latex_code = latex_code.replace('pval','p-value')
latex_code = latex_code.replace('fdr','FDR')

with open('/mnt/projects/VIA11/FREESURFER/Stats/Plots/lateral/lateral_stats_combined_latex_table.txt', 'w') as f:
    f.write(latex_code)"""
    
#replace  \toprule
# &\multicolumn{4}{c|}{Left Hemisphere}& \multicolumn{4}{|c}{Right Hemisphere} &\\ \hline

#with open('/mnt/projects/VIA11/FREESURFER/Stats/Plots/lateral/lateral_stats_combined_latex_table.txt', 'r') as f:
#    latex_code = f.read()

#latex_code = latex_code.replace("NaN","")