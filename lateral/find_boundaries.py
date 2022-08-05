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

all_borders = []

regions_manual = ["bankssts", "caudal anterior cingulate","caudal middle frontal","cuneus","entorhinal",  
                   "fusiform","inferior parietal","inferior temporal","isthmus cingulate","lateral occipital",     
                   "lateral orbitofrontal","lingual","medial orbitofrontal","middle temporal","parahippocampal",      
                   "paracentral","pars opercularis","pars orbitalis","pars triangularis","pericalcarine",  
                   "postcentral","posterior cingulate","precentral","precuneus","rostral anterior cingulate",
                   "rostral middle frontal","superior frontal","superior parietal","superior temporal","supramarginal"  ,         
                   "frontal pole","temporal pole","transverse temporal","insula"]

Dict = {}

for region in regions_manual:
    #region = regions_manual[0]
    
    region_borders = []
    
    imgs = glob.glob(f"/mnt/projects/VIA11/FREESURFER/Stats/Plots/lateral/individual_regions/{region}_*.png")
    imgs.sort()
    
    for im in imgs:
        img = skimage.io.imread(im)
        r = rgb2gray(img)
        
        #remove colorbar 
        r[0:250,500:700] = 1
        
        # Find contours at a constant value of 0.8
        contours = measure.find_contours(r,0.7)
        
        if len(contours)>0:
            #Dict[region][1]
            region_borders.append(contours[0])
            contours = [contours[0]]
            
            # Display the image and plot all contours found
            fig, ax = plt.subplots()
            ax.imshow(r, cmap=plt.cm.gray)
            
            for contour in contours:
                ax.plot(contour[:, 1], contour[:, 0], linewidth=2)
            
            ax.axis('image')
            ax.set_xticks([])
            ax.set_yticks([])
            plt.show()
        else:
            region_borders.append([])
            
            # Display the image and plot all contours found
            fig, ax = plt.subplots()
            ax.imshow(r, cmap=plt.cm.gray)
    
            ax.axis('image')
            ax.set_xticks([])
            ax.set_yticks([])
            plt.show()
    
    all_borders.append(region_borders)

#%%

# Display the image and plot all contours found
img_name = "cohensd_brain"
txt_name = "significant_regions"

img_path = glob.glob(f"/mnt/projects/VIA11/FREESURFER/Stats/Plots/lateral/{img_name}_*.png")
sig_regions_path = glob.glob(f"/mnt/projects/VIA11/FREESURFER/Stats/Plots/lateral/{txt_name}_*.txt")
img_path.sort()
sig_regions_path.sort()

fig = plt.figure(figsize=(20, 10))

# setting values to rows and column variables
rows = 2
columns = 4

plot_order = [2,6,1,5,3,7,4,8]

for n,(img, sig) in enumerate(zip(img_path, sig_regions_path)):
    fig.add_subplot(rows, columns, plot_order[n])
    
    r = skimage.io.imread(img)
    #r = rgb2gray(img)
    
    with open(sig) as f:
        sig_lines = f.read().splitlines()
    
    #fig, ax = plt.subplots()
    plt.imshow(r, cmap=plt.cm.gray)
    
    for r, region_b in enumerate(all_borders):
        
        if all(x in img for x in ["_lh_","_inside_"]): #"_lh_" and "_inside_" in img:
            border_type = 0
        elif all(x in img for x in ["_lh_","_outside_"]):
            border_type = 1
        elif all(x in img for x in ["_rh_","_outside_"]):
            border_type = 2
        elif all(x in img for x in ["_rh_","_inside_"]):
            border_type = 3
            
        border = region_b[border_type]
        
        if len(sig_lines[0])>0:
            sig_regions = [s.split('_')[1] for s in sig_lines]
        else:
            sig_regions = []
            
        if len(border)>0:
            if regions_manual[r].replace(" ","") in sig_regions:
                plt.plot(border[:, 1], border[:, 0], color='k', linewidth=1.5) #"#8C000F"
            else:
                plt.plot(border[:, 1], border[:, 0], color='w', linewidth=0.3) #0.3
            
            
    plt.axis('image')
    #plt.set_xticks([])
    #plt.set_yticks([])
    plt.axis('off')
    #plt.savefig(img,format='png')
    #plt.show()


plt.gcf().text(0.48, 0.53, "Hemisphere", fontsize=14) 
plt.gcf().text(0.3, 0.58, "Left", fontsize=14) 
plt.gcf().text(0.7, 0.58, "Right", fontsize=14) 

plt.gcf().text(0.48, 0.1, "Hemisphere", fontsize=14) 
plt.gcf().text(0.3, 0.15, "Left", fontsize=14) 
plt.gcf().text(0.7, 0.15, "Right", fontsize=14) 

plt.gcf().text(0.48, 0.85, "Female", fontsize=24) 
plt.gcf().text(0.48, 0.43, "Male", fontsize=24) 
plt.gcf().text(0.1, 0.9, "Effect size of BP brain volume difference from control from models without global variable", fontsize=24) 


plt.savefig('/mnt/projects/VIA11/FREESURFER/Stats/Plots/lateral/test.png',format='png') # BP_vol_cohensd_without_glob.png


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