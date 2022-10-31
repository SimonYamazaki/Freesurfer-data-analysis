#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 12 10:11:59 2022

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

with open('/mnt/projects/VIA11/FREESURFER/Stats/Plots/lateral/all_borders.json') as json_file:
    json_data2 = json.load(json_file)
    
data2 = json.loads(json_data2, object_hook=deconvert)


contrast = "SZ"
measure = "thickness"
glob_var = "without_global_var"

c = ["SZ", "BP"]
m = ["volume", "area","thickness"]
g = ["without_global_var", "with_global_var"]

for contrast in c:
    for measure in m:
        for glob_var in g:
            
            print(contrast)
            print(measure)
            
            img_name = "cohensd"
            txt_name = "significant_regions"
            img_path = glob.glob(f"/mnt/projects/VIA11/FREESURFER/Stats/Plots/lateral/intermediate_plots/{img_name}_{contrast}*_{measure}_{glob_var}_*.png")
            sig_regions_path = glob.glob(f"/mnt/projects/VIA11/FREESURFER/Stats/Plots/lateral/intermediate_plots/{txt_name}_{contrast}*_{measure}_{glob_var}_*.txt")
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
                            plt.plot(border[:, 1], border[:, 0], color='k', linewidth=5) #"#8C000F"
                        else:
                            plt.plot(border[:, 1], border[:, 0], color='k', linewidth=0.3) #0.3
                      
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
            
            plt.savefig(f"/mnt/projects/VIA11/FREESURFER/Stats/Plots/lateral/final_plots/brain_{contrast}_{measure}_{glob_var}.png",format='png') 





#%%
from PIL import Image

img_paths = glob.glob("/mnt/projects/VIA11/FREESURFER/Stats/Plots/lateral/final_plots/*.png")
pdf_path = "/mnt/projects/VIA11/FREESURFER/Stats/Plots/lateral/final_plots/all_brain_maps2.pdf"

images = [Image.open(f).convert('RGB') for f in img_paths]

images[0].save(
    pdf_path, "PDF" ,resolution=100.0, save_all=True, append_images=images[1:]
)
