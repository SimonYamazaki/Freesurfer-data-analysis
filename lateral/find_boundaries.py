#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 29 11:57:09 2022

@author: simonyj
"""

import numpy as np
import matplotlib.pyplot as plt

from skimage import measure
import skimage.io
from skimage.color import rgb2gray
import skimage.segmentation

import glob 

import collections

import json

#%%

all_borders = collections.defaultdict(dict)

regions_manual = ["bankssts","caudal anterior cingulate","caudal middle frontal","cuneus","entorhinal",  
                   "fusiform","inferior parietal","inferior temporal","isthmus cingulate","lateral occipital",     
                   "lateral orbitofrontal","lingual","medial orbitofrontal","middle temporal","parahippocampal",      
                   "paracentral","pars opercularis","pars orbitalis","pars triangularis","pericalcarine", "postcentral","posterior cingulate","precentral","precuneus","rostral anterior cingulate",
                   "rostral middle frontal","superior frontal","superior parietal","superior temporal","supramarginal"  ,         
                   "frontal pole","temporal pole","transverse temporal","insula"]


show_plot = True

for region in regions_manual:    
    
    imgs = glob.glob(f"/mnt/projects/VIA11/FREESURFER/Stats/Plots/lateral/individual_regions2/{region}_*.png")
    imgs.sort()
    
    for im in imgs:
        border_type = None
        if all(x in im for x in ["lh_inside"]):#["_lh","_inside"]):
            border_type = "lh_inside"
        elif all(x in im for x in ["lh_outside"]):#["_lh","_outside"]):
            border_type = "lh_outside"
        elif all(x in im for x in ["rh_outside"]):#["_rh","_outside"]):
            border_type = "rh_outside"
        elif all(x in im for x in ["rh_inside"]):#["_rh","_inside"]):
            border_type = "rh_inside"
            
        img = skimage.io.imread(im)
        r = rgb2gray(img)
        
        #remove colorbar 
        r[0:250,500:700] = 1
        
        # Find contours at a constant value of 0.8
        contours = measure.find_contours(r,0.7)
        
        if len(contours)>0:
            all_borders[region][border_type] = contours[0]
        else:
            all_borders[region][border_type] = []
        
        # Display the image and plot all contours found
        if show_plot:
            fig, ax = plt.subplots()
            ax.imshow(r, cmap=plt.cm.gray)
            
            ax.axis('image')
            ax.set_xticks([])
            ax.set_yticks([])
            plt.show()
            

#%%
def convert(x):
    if hasattr(x, "tolist"):  # numpy arrays have this
        return {"$array": x.tolist()}  # Make a tagged object
    raise TypeError(x)


def deconvert(x):
    if len(x) == 1:  # Might be a tagged object...
        key, value = next(iter(x.items()))  # Grab the tag and value
        if key == "$array":  # If the tag is correct,
            return np.array(value)  # cast back to array
    return x

json_data = json.dumps(all_borders, default=convert)


with open('/mnt/projects/VIA11/FREESURFER/Stats/Plots/lateral/all_borders.json', 'w') as json_file:
    json.dump(json_data, json_file)
    
with open('/mnt/projects/VIA11/FREESURFER/Stats/Plots/lateral/all_borders.json') as json_file:
    json_data2 = json.load(json_file)
    
data2 = json.loads(json_data2, object_hook=deconvert)


#%%

# Display the image and plot all contours found
contrast = "SZ"
measure = "volume"
glob_var = "without_global_var"

img_name = "cohensd"
txt_name = "significant_regions"

img_path = glob.glob(f"/mnt/projects/VIA11/FREESURFER/Stats/Plots/lateral/{img_name}_{contrast}*_{measure}_{glob_var}_*.png")
sig_regions_path = glob.glob(f"/mnt/projects/VIA11/FREESURFER/Stats/Plots/lateral/{txt_name}_{contrast}*_{measure}_{glob_var}_*.txt")
img_path.sort()
sig_regions_path.sort()

fig = plt.figure(figsize=(60, 30))

# setting values to rows and column variables
rows = 2
columns = 4

#plot_order = [2,6,1,5,3,7,4,8]
plot_order = [2,1,3,4,6,5,7,8]


for n,(im, sig) in enumerate(zip(img_path, sig_regions_path)):
    fig.add_subplot(rows, columns, plot_order[n])
    
    r = skimage.io.imread(im)
    #r = rgb2gray(img)
    
    with open(sig) as f:
        sig_lines = f.read().splitlines()
    
    #if n not in (6,7):
        #remove colorbar 
    #r[0:240,550:700,:] = [255,255,255, 0.9]
        
    #fig, ax = plt.subplots()
    plt.imshow(r, cmap=plt.cm.gray)
    
    #for r, region_b in enumerate(all_borders):
    for region in regions_manual:
         
        if all(x in im for x in ["lh_inside"]):#["_lh","_inside"]):
            border_type = "lh_inside"
        elif all(x in im for x in ["lh_outside"]):#["_lh","_outside"]):
            border_type = "lh_outside"
        elif all(x in im for x in ["rh_outside"]):#["_rh","_outside"]):
            border_type = "rh_outside"
        elif all(x in im for x in ["rh_inside"]):#["_rh","_inside"]):
            border_type = "rh_inside"
            
        border = data2[region][border_type]
        
        if len(sig_lines[0])>0:
            sig_regions = [s.split('_')[1] for s in sig_lines]
        else:
            sig_regions = []
            
        if len(border)>0:
            if region.replace(" ","") in sig_regions:
                plt.plot(border[:, 1], border[:, 0], color='k', linewidth=1.5) #"#8C000F"
            else:
                plt.plot(border[:, 1], border[:, 0], color='k', linewidth=0.1) #0.3
        
            
    plt.axis('image')
    plt.axis('off')


plt.subplots_adjust(wspace=0, hspace=0)

plt.gcf().text(0.49, 0.53, "Hemisphere", fontsize=40) 
plt.gcf().text(0.31, 0.58, "Left", fontsize=40) 
plt.gcf().text(0.7, 0.58, "Right", fontsize=40) 

plt.gcf().text(0.49, 0.15, "Hemisphere", fontsize=40) 
plt.gcf().text(0.31, 0.2, "Left", fontsize=40) 
plt.gcf().text(0.7, 0.2, "Right", fontsize=40) 

plt.gcf().text(0.50, 0.85, "Female", fontsize=50) 
plt.gcf().text(0.50, 0.45, "Male", fontsize=50) 

plt.gcf().text(0.1, 0.9, f"Effect size of {contrast} brain {measure} difference from control from models {glob_var}", fontsize=50) 


plt.savefig('/mnt/projects/VIA11/FREESURFER/Stats/Plots/lateral/test.png',format='png') 


#%%

import pandas as pd

effect_sizes_path = "/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/parcels/lateral/ANOVA+contrast_effect_sizes.xlsx"
contrast_fdr_table = "/mnt/projects/VIA11/FREESURFER/Stats/Plots/lateral/lateral_combined_contrasts.xlsx"

df1 = pd.read_excel(effect_sizes_path) 
df2 = pd.read_excel(contrast_fdr_table) 

df2['sex'] = df2['sex'].map({0:'female', 1:'male', 2:'both'})


new_df = pd.merge(df1, df2,  how='inner', left_on=["model_yvar","measure","hemisphere","sex","globvar_in_model"], right_on = ["Model_yvar","measure","hemisphere","sex","global_var_in_model"])

new_df = new_df[["model_yvar","measure","hemisphere","sex","global_var_in_model","tratio_BP-K","tratio_SZ-K","pval_BP-K","pval_SZ-K","fdr_pvals_BP-K","fdr_pvals_SZ-K","CohensD_BP-K","CohensD_SZ-K"]]
#df.year.astype(Int64)
new_df['global_var_in_model'] = new_df['global_var_in_model'].map({0:'no', 1:'yes'})


#%%
#extract the data for each table 
m = "volume"
s = "female"
g = "no"

df = new_df.loc[ (new_df["measure"].isin([m])) & (new_df["global_var_in_model"].isin([g])) & (new_df["sex"].isin([s]))]

df_lh = df.loc[ df["hemisphere"].isin(["lh"]) ]
df_rh = df.loc[ df["hemisphere"].isin(["rh"]) ]

df_lh["model_yvar"] = df_lh["model_yvar"].str.replace("lh_","")
df_rh["model_yvar"] = df_rh["model_yvar"].str.replace("rh_","")


keep_same = {"model_yvar","measure","hemisphere","sex","global_var_in_model"}
df_lh.columns = ['{}{}'.format(c, '' if c in keep_same else '_lh') for c in df_lh.columns]
df_rh.columns = ['{}{}'.format(c, '' if c in keep_same else '_rh') for c in df_rh.columns]

final_df = pd.merge(df_lh, df_rh,  how='left', left_on=["model_yvar","measure","sex","global_var_in_model"],right_on=["model_yvar","measure","sex","global_var_in_model"])
#final_df = df_lh.merge(df_rh, on=["model_yvar","measure","hemisphere","sex","global_var_in_model"])


#%%

latex_code = final_df.to_latex(index=False)

with open('/mnt/projects/VIA11/FREESURFER/Stats/Plots/lateral/lateral_stats_combined_latex_table.txt', 'w') as f:
    f.write(latex_code)
    
with open('/mnt/projects/VIA11/FREESURFER/Stats/Plots/lateral/lateral_stats_combined_latex_table.txt', 'r') as f:
    latex_code = f.read()

#latex_code = latex_code.replace("NaN","")



"""
### with segmentation

contours = skimage.segmentation.find_boundaries(r,connectivity=1)
#bimg = r*np.logical_not(contours)
bimg = np.logical_not(contours)
fig, ax = plt.subplots()
ax.imshow(bimg, cmap=plt.cm.gray)

ax.axis('image')
ax.set_xticks([])
ax.set_yticks([])
plt.show()
"""