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
library(htmlwidgets)
library(filesstrings)

#set working directory for plots to be saved
setwd("/mnt/projects/VIA11/FREESURFER/Stats/Plots/lateral")

#c("bankssts", "caudal anterior cingulate","caudal middle frontal","cuneus","entorhinal",  
#  "fusiform","inferior parietal","inferior temporal","isthmus cingulate","lateral occipital",     
#  "lateral orbitofrontal","lingual","medial orbitofrontal","middle temporal","parahippocampal",      
#  "paracentral","pars opercularis","pars orbitalis",
  
regions_manual = c("pars triangularis","pericalcarine","postcentral","posterior cingulate","precentral","precuneus","rostral anterior cingulate",
                   "rostral middle frontal","superior frontal","superior parietal","superior temporal","supramarginal"  ,         
                   "frontal pole","temporal pole","transverse temporal","insula")

h = c("lh","rh")
hemi = c("left","right")
side = c("inside","outside","outside","inside")
camera = list()#c("left lateral","right lateral")
camera[[1]] = list(camera = list(eye = list(x = 2, y = 0, z = 0)))
camera[[2]] = list(camera = list(eye = list(x = -2, y = 0, z = 0)))

Sys.setenv("plotly_username"="simonyj")
Sys.setenv("plotly_api_key"="IrXExPAPh7emxCFHdOBe")

Sys.setenv("plotly_username"="s174165")
Sys.setenv("plotly_api_key"="uldxFDzocuKNM5D4M4SB")


for (r in seq(1,length(regions_manual))){
  
region = regions_manual[r]
region_vals = (region == regions_manual)*1

print(region)

  for (i in seq(1,2)){ #hemisphere
    #add cohensd to a dataframe to plot 
    someData = tibble(region = regions_manual,p = region_vals,)
    
    for (k in seq(1,2)){
      if ( h[i] == "rh" ){
        kk=k+2
      }else {
        kk=k
      }
      
      pls[[regions_manual[r]]][[h[i]]][[side[kk]]] <- ggseg3d(.data = someData,
                                       atlas = dk_3d,
                                       colour = "p", text = "p",
                                       surface = "inflated",
                                       hemisphere = c(hemi[i]),
                                       palette = c("#ffffff"=0,"#000000"=1),
                                       na.colour = "white",
                                       options.legend=list(tickfont=list(size=20),title=list(text="Cohens D"))) %>%
        layout(scene = camera[[k]]) %>% 
        remove_axes() %>%
        hide_colorbar()
      
      #op <- options()
      #options(viewer = NULL)
      #out_file = gsub(" ","_",paste("'",paste(regions_manual[r],"ref_plot",h[i],side[kk],sep="_"),"'",sep=""))
      #pls[[regions_manual[r]]][[h[i]]][[side[kk]]] %>% htmlwidgets::onRender(
      #  paste("function(el, x) {var gd = document.getElementById(el.id);Plotly.downloadImage(gd, {format: 'png', width: 1000, height: 1000, filename: ",out_file,"});}",sep=""))
      #options(viewer = op$viewer)
      #Sys.sleep(8)
      #home_dir = "/mrhome/simonyj/Downloads"
      #save_folder = "/mnt/projects/VIA11/FREESURFER/Stats/Plots/lateral/individual_regions2"
      #file.copy(from = paste(home_dir,"/",gsub("'","",out_file),".png",sep=""), to = save_folder)#paste(save_folder,"/",gsub("'","",out_file),".png",sep=""))
      #file.remove(paste(home_dir,"/",gsub("'","",out_file),".png",sep=""))
      
      plotly_IMAGE(pls[[regions_manual[r]]][[h[i]]][[side[kk]]], width = 1000, height = 1000, format = "png", scale = 2,
                   out_file = paste("/mnt/projects/VIA11/FREESURFER/Stats/Plots/lateral/individual_regions2/",paste(regions_manual[r],"ref_plot",h[i],side[kk],sep="_"),".png",sep=""))
      
      #plots2save[[regions_manual[r]]][[h[i]]][[side[kk]]] = paste("individual_regions2/",paste(regions_manual[r],"ref_plot",h[i],side[kk],sep="_"),".png",sep="")
      
    } #end side
  } #end hemisphere
} #end region



#plots 2 save
#plots2save = c(paste("individual_regions2/",region,"_ref_plot_lh_outside.png",sep = ""), paste("individual_regions2/",region,"_ref_plot_lh_inside.png",sep = ""),
#               paste("individual_regions2/",region,"_ref_plot_rh_outside.png",sep = ""), paste("individual_regions2/",region,"_ref_plot_rh_inside.png",sep = ""))

#orca(pls[["lh"]][["outside"]],plots2save[1])
#orca(pls[["lh"]][["inside"]],plots2save[2])
#orca(pls[["rh"]][["outside"]],plots2save[3])
#orca(pls[["rh"]][["inside"]],plots2save[4])



#library(RSelenium)
#pls[["lh"]][["outside"]] %>%
#  export(file = "filename.svg",
#         selenium = RSelenium::rsDriver(browser = "firefox"))

#library(htmlwidgets)
# Save viewer settings (e.g. RStudio viewer pane)
#op <- options()
# Set viewer to web browser
#options(viewer = NULL)
# Use web browser to save image
#pls[["lh"]][["outside"]] %>% htmlwidgets::onRender(
#  "function(el, x) {
#  var gd = document.getElementById(el.id); 
#  Plotly.downloadImage(gd, {format: 'png', width: 1000, height: 1000, filename: 'plot'});
#  }"
#)
# Restore viewer to old setting (e.g. RStudio)
#options(viewer = op$viewer)



