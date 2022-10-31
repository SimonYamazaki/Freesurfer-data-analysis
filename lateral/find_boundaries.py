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