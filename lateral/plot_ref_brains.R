#load packages
library(data.table)
library(ggplot2)
library(ggdist)
library(gridExtra)
library(tidyr)
library(lsmeans)
library(grid)
library(ggnewscale)
library(writexl)
library(car)
library(NCmisc)
library(lsr)
library(readxl)

#set working directory for plots to be saved
setwd("/mnt/projects/VIA11/FREESURFER/Stats/Plots/lateral")


regions_manual = c("bankssts", "caudal anterior cingulate","caudal middle frontal","cuneus","entorhinal",  
                   "fusiform","inferior parietal","inferior temporal","isthmus cingulate","lateral occipital",     
                   "lateral orbitofrontal","lingual","medial orbitofrontal","middle temporal","parahippocampal",      
                   "paracentral","pars opercularis","pars orbitalis","pars triangularis","pericalcarine",  
                   "postcentral","posterior cingulate","precentral","precuneus","rostral anterior cingulate",
                   "rostral middle frontal","superior frontal","superior parietal","superior temporal","supramarginal"  ,         
                   "frontal pole","temporal pole","transverse temporal","insula")

h = c("lh","rh")
hemi = c("left","right")
camera = list()#c("left lateral","right lateral")
camera[[1]] = list(camera = list(eye = list(x = 2, y = 0, z = 0)))
camera[[2]] = list(camera = list(eye = list(x = -2, y = 0, z = 0)))


region = regions_manual[34]
region_vals = (region == regions_manual)*1

for (i in seq(1,2)){ #hemisphere
  #add cohensd to a dataframe to plot 
  someData = tibble(region = regions_manual,p = region_vals,)
  
  for (k in seq(1,2)){
    pls[[h[i]]][[side[k]]] <- ggseg3d(.data = someData,
                                     atlas = dk_3d,
                                     colour = "p", text = "p",
                                     surface = "inflated",
                                     hemisphere = c(hemi[i]),
                                     palette = c("#ffffff"=0,"#000000"=1),
                                     na.colour = "white",
                                     options.legend=list(tickfont=list(size=25),title=list(text="Cohens D"))) %>%
      layout(scene = camera[[k]]) %>% 
      remove_axes()
      #hide_colorbar()
  } #end side
} #end hemisphere

#plots 2 save
plots2save = c(paste("individual_regions/",region,"_plot1_ref.png",sep = ""), paste("individual_regions/",region,"_plot2_ref.png",sep = ""),
               paste("individual_regions/",region,"_plot3_ref.png",sep = ""), paste("individual_regions/",region,"_plot4_ref.png",sep = ""))

orca(pls[["lh"]][["right"]],plots2save[1])
orca(pls[["lh"]][["left"]],plots2save[2])
orca(pls[["rh"]][["right"]],plots2save[3])
orca(pls[["rh"]][["left"]],plots2save[4])


