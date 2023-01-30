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
library(ggseg3d)
library(plotly)


Sys.setenv("plotly_username" = "simonyj")
Sys.setenv("plotly_api_key" = "Fjx0xgvDeX23Z2DBu49P")

#save paths:
contrast_with_glob ="/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/parcels/lateral/lateral_Parcel_model_contrast_with_glob.xlsx"
contrast_without_glob = "/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/parcels/lateral/lateral_Parcel_model_contrast_without_glob.xlsx"


effect_sizes_path = "/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/parcels/lateral/ANOVA+contrast_effect_sizes.xlsx"

plot_save_dir = "/mnt/projects/VIA11/FREESURFER/Stats/Plots/lateral"



#####
# - Plots LSmeans based contrasts which compares BP and SZ to control

#save folder for plots:
save_folder = "/mnt/projects/VIA11/FREESURFER/Stats/Plots/lateral"

#prefix for the plot
LSmeans_prefix = "LSmean_difference_"


#set working directory for plots to be saved
setwd("/mnt/projects/VIA11/FREESURFER/Stats/Plots/lateral")




#load lsmeans and contrasts

#save paths:
contrast_path_with_glob ="/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/parcels/lateral/lateral_Parcel_model_contrast_with_glob.xlsx"
contrast_path_without_glob = "/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/parcels/lateral/lateral_Parcel_model_contrast_without_glob.xlsx"

contrast_table_with_glob <-  read_excel(contrast_path_with_glob)
contrast_table_without_glob <-  read_excel(contrast_path_without_glob)

contrast_table = rbind(contrast_table_with_glob, contrast_table_without_glob)

DFc = contrast_table



#LSmeans contrast plot

pivot_cols = c("Contrast_BP-K","Contrast_SZ-K")
DFp = DFc %>%
  pivot_longer(cols = pivot_cols,
               names_to = "diff_group",
               values_to = "diff_value",
               values_drop_na = FALSE)
DFp = as.data.frame(DFp)
DFp$diff_group[DFp$diff_group == "Contrast_BP-K"] = "BP"
DFp$diff_group[DFp$diff_group == "Contrast_SZ-K"] = "SZ"


pivot_cols_LCL = c("LCL_Contrast_BP-K","LCL_Contrast_SZ-K")
DF_LCL = DFc %>%
  pivot_longer(cols = pivot_cols_LCL,
               names_to = "diff_group_LCL",
               values_to = "diff_LCL",
               values_drop_na = FALSE)

pivot_cols_UCL = c("UCL_Contrast_BP-K","UCL_Contrast_SZ-K")
DF_UCL = DFc %>%
  pivot_longer(cols = pivot_cols_UCL,
               names_to = "diff_group_UCL",
               values_to = "diff_UCL",
               values_drop_na = FALSE)

DF_CL = data.frame(DF_LCL$diff_group_LCL, DF_LCL$diff_LCL, DF_UCL$diff_UCL)
names(DF_CL) = c("diff_group_CL","diff_LCL","diff_UCL")
DFp = cbind(DFp,DF_CL)

DFp$diff_value_P = 100*DFp$diff_value / DFp$K_LSmean
DFp$diff_LCL_P = 100*DFp$diff_LCL / DFp$K_LSmean
DFp$diff_UCL_P = 100*DFp$diff_UCL / DFp$K_LSmean


#make a new coloumn in the dataframe with a non zero value when the model
#contrast is significant. The value given is the y-axis position of a star being 
#drawn indicating significant contrast in the plot below
#DFp$BP_diff_sig = NA
#DFp$BP_diff_sig[DFp$BP_diff_pv < 0.05] = -4
#DFp$SZ_diff_sig = NA
#DFp$SZ_diff_sig[DFp$SZ_diff_pv < 0.05] = -5


#define relevant variables for naming
h = c("lh","rh")
m = c("area","thickness","volume")
groups = c("BP","SZ")
m = c("area","thickness","volume")
glob_names = c("Models WITHOUT global var","Models WITH global var")
glob = c("without_global_var","with_global_var")

sx = c("female","male")

sp=list()

for (j in seq(1,3)){ 
  for (g in seq(1,2)){
    for (i in seq(1,2)){  #hemisphere 
      for (s in seq(1,2)){  
        
        dfs = DFp[DFp$measure == m[j] & DFp$global_var_in_model == (g-1) & DFp$hemisphere == h[i] & DFp$sex == s-1,]
        
        sp[[2*s-2+i]]=ggplot(dfs, aes_string(group="diff_group",color="diff_group", x="Model_yvar", y="diff_value_P")) +
          labs(x="region",y="Difference from control [%]") +
          geom_line() + 
          ggtitle(paste(sx[s],h[i])) + 
          geom_hline(yintercept = 0) + 
          geom_errorbar(aes_string(width=0.3,group="diff_group_CL",color="diff_group_CL",ymin="diff_LCL_P", ymax="diff_UCL_P" )) +
          scale_color_manual("Group", values=c("#0072B2","#0072B2","#009E73", "#009E73"))+
          
          geom_point(aes_string(group="diff_group",color="diff_group",x="Model_yvar", y="diff_value_P")) +
          #geom_point(shape=8, aes_string(group="diff_group",x="model_yvar", y="BP_diff_sig"), color="#0072B2") +
          #geom_point(shape=8, aes_string(group="diff_group",x="model_yvar", y="SZ_diff_sig"), color="#009E73") +
          
          theme(axis.text.x = element_text(color = "black", size = 8, angle = 90)) + 
          ylim(-10, 10) +
          annotate("text", x = 8, y = 4,label = "* = Significant uncorrected contrast",family = "", fontface = 3, size=3)+
          
          {if (m[j]=="thickness") ylim(-5, 5)}
        
        #coord_flip()
      } 
    }
    
    top_title = paste("Lateral LSmean difference from control:",m[j],glob_names[g])
    ps=grid.arrange(grobs=sp, top=textGrob(top_title,gp=gpar(fontsize=20)))
    ggsave(paste(LSmeans_prefix,m[j],glob[g],".png",sep=""),ps,width = 15,height = 10)
    
  } #g
} #j







### plot the brains 
#remotes::install_github("LCBC-UiO/ggseg")


#### ggseg3d attempt
library(stringr)
library(ggplot2)
library(plotly)
library(ggseg3d)
#library(ggsegExtra) #remotes::install_github("ggseg/ggsegExtra")
#test <- make_ggseg3d_2_ggseg()

#col_names = names(datab)
#find names of all the regions based on coloumn names that include "lh_" 
#regions = col_names[grepl("^lh_", col_names)]
#regions <- unique(unlist(strsplit(regions,"_")))
#remove some specific coloumns
#col2rm = c(m,"WhiteSurfArea","lh","MeanThickness")
#regions = regions[!(regions %in% col2rm)]

regions_manual = c("bankssts", "caudal anterior cingulate","caudal middle frontal","cuneus","entorhinal",  
                   "fusiform","inferior parietal","inferior temporal","isthmus cingulate","lateral occipital",     
                   "lateral orbitofrontal","lingual","medial orbitofrontal","middle temporal","parahippocampal",      
                   "paracentral","pars opercularis","pars orbitalis","pars triangularis","pericalcarine",  
                   "postcentral","posterior cingulate","precentral","precuneus","rostral anterior cingulate",
                   "rostral middle frontal","superior frontal","superior parietal","superior temporal","supramarginal"  ,         
                   "frontal pole","temporal pole","transverse temporal","insula")


#load effect sizes
h = c("lh","rh")
m = c("area","thickness","volume")
sex = c("female","male","both") #3 sexes are defined, 0=female, 1=male, 2=both
glob = c("without_global_var","with_global_var")
hemi = c("left","right")
contrast = c("BP-K", "SZ-K") #"SZ-K"#

side = c("inside","outside","outside","inside")
camera = list()
camera[[1]] = list(camera = list(eye = list(x = 2, y = 0, z = 0)))
camera[[2]] = list(camera = list(eye = list(x = -2, y = 0, z = 0)))

plots2save = list()
pls = list()


#load and prepare contrast table
contrast_table_with_glob <-  read_excel(contrast_with_glob)
contrast_table_without_glob <-  read_excel(contrast_without_glob)
contrast_table = rbind(contrast_table_with_glob, contrast_table_without_glob)

#load effect size table 
eff_table <-  read_excel(effect_sizes_path)

#find biggest and smallest effect sizes
all_cohensd = as.numeric(rbind(eff_table$`CohensD_BP-K`, eff_table$`CohensD_SZ-K`))
max_cohens_D = max(all_cohensd)
min_cohens_D = min(all_cohensd)


#j = 3 #measure
#gg = 1 #global variable 
#contrast = "SZ-K"

for (cc in seq(1,2)){
  contrast_table[paste("fdr_pvals_",contrast[cc],sep="")] = NA
  
for (j in seq(1,3)){
for (gg in seq(1,2)){
for (ss in seq(1,2)){
  
#FDR correct p-values 
contrast_pval_coloumn = paste("pval_",contrast[cc],sep="")
raw_pvals <- contrast_table[c(contrast_table$measure == m[j] & contrast_table$sex == ss-1 & contrast_table$global_var_in_model == gg-1),][[contrast_pval_coloumn]]
fdr_pvals = p.adjust(raw_pvals, method = "fdr", n = length(raw_pvals))
contrast_table[paste("fdr_pvals_",contrast[cc],sep="")][c(contrast_table$measure == m[j] & contrast_table$sex == ss-1 & contrast_table$global_var_in_model == gg-1),] = fdr_pvals

for (i in seq(1,2)){ #hemisphere
  
  table_sex <- eff_table[c(eff_table$hemisphere == h[i] & eff_table$measure == m[j] & eff_table$sex == sex[ss] & eff_table$globvar_in_model == gg-1),]
  
  contrast_coloumn = paste("CohensD_",contrast[cc],sep="")
  cohensd = as.numeric(table_sex[[contrast_coloumn]])
  
  corrected_contrast_pval = contrast_table[paste("fdr_pvals_",contrast[cc],sep="")][c(contrast_table$measure == m[j] & contrast_table$sex == ss-1 & contrast_table$global_var_in_model == gg-1 & contrast_table$hemisphere == h[i]),]

  if ( any(as.numeric(unlist(corrected_contrast_pval)) < 0.05) ){
    significant_regions = contrast_table[c(contrast_table$measure == m[j] & contrast_table$sex == ss-1 & contrast_table$global_var_in_model == gg-1 & contrast_table$hemisphere == h[i]),][as.numeric(unlist(corrected_contrast_pval)) < 0.05,]$Model_yvar
  }  else {
    significant_regions = ""
  }
  #significant_text = paste(significant_regions, "<br>",collapse = '')
  #significant_text = paste("Significant regions:<br>",significant_text,collapse = '')
  
  #add cohensd to a dataframe to plot 
  someData = tibble(
    region = regions_manual,
    p = cohensd,
  )
  
  #plot cohens d
  for (k in seq(1,2)){
    if ( h[i] == "rh" ){
      kk=k+2
    }else {
      kk=k
    }

    pls[[h[i]]][[m[j]]][[sex[ss]]][[glob[gg]]][[side[kk]]] <- ggseg3d(.data = someData,
                                                                     atlas = dk_3d,
                                                                     colour = "p", text = "p",
                                                                     surface = "inflated",
                                                                     hemisphere = c(hemi[i]),
                                                                     palette = c("#0000ff"=min_cohens_D, "#ffffff"=0,"#ff0000"=max_cohens_D),
                                                                     options.legend=list(tickfont=list(size=25),tickvals = c(-0.6,-0.4,-0.2,0, 0.2,0.4, 0.6),title=list(text="Cohens D"))) %>%
      #pan_camera(camera[k]) %>%
      layout(scene = camera[[k]]) %>% 
      remove_axes() %>% 
      hide_colorbar()

    #plots2save[[h[i]]][[m[j]]][[sex[ss]]][[glob[gg]]][[side[kk]]] = paste(paste("cohensd","brain",h[i],side[kk],sex[ss],sep="_"),".png",sep="")# "cohensd_brain_lh_inside_male.png"
    
    img_path = paste("/mnt/projects/VIA11/FREESURFER/Stats/Plots/lateral/intermediate_plots/",paste("cohensd",contrast[cc],sex[ss],m[j],glob[gg],h[i],side[kk],sep="_"),".png",sep="")
    #plotly_IMAGE(pls[[h[i]]][[m[j]]][[sex[ss]]][[glob[gg]]][[side[kk]]], width = 1000, height = 1000, format = "png", scale = 2,out_file = img_path)
    
    txt2save = paste("intermediate_plots/significant","regions",contrast[cc],sex[ss],m[j],glob[gg],h[i],side[kk],".txt",sep="_")#"significant_regions_lh_inside_male.txt"
    #fileConn<-file(txt2save)
    #writeLines(significant_regions, fileConn)
    #close(fileConn)  
    
  } #end side
} #end hemisphere
} #sex

} #c 
} #g
} #j

#get high resolution colorbar
brain_colorbar <- ggseg3d(.data = someData,
                          atlas = dk_3d,
                          colour = "p", text = "p",
                          surface = "inflated",
                          hemisphere = c(hemi[i]),
                          palette = c("#0000ff"=min_cohens_D, "#ffffff"=0,"#ff0000"=max_cohens_D),
                          options.legend=list(tickfont=list(size=25),tickvals = c(-0.6,-0.4,-0.2,0, 0.2,0.4, 0.6),title=list(text="Cohens D"))) %>%
  layout(scene = camera[[k]]) %>% 
  remove_axes()

img_path = "/mnt/projects/VIA11/FREESURFER/Stats/Plots/lateral/final_plots/brain_colorbar.png"
plotly_IMAGE(brain_colorbar, width = 1000, height = 1000, format = "png", scale = 2,out_file = img_path)


contrast_fdr_table = "/mnt/projects/VIA11/FREESURFER/Stats/Plots/lateral/lateral_combined_contrasts.xlsx"
write_xlsx(contrast_table,contrast_fdr_table)


#https://chart-studio.plotly.com/settings/api#/
#https://plotly.com/r/reference/layout/coloraxis/

#these are in the same order as plot2save
#orca(pls[["lh"]][[m[j]]][["female"]][[glob[gg]]][["inside"]], plots2save[["lh"]][[m[j]]][["female"]][[glob[gg]]][["inside"]])
#orca(pls[["lh"]][[m[j]]][["female"]][[glob[gg]]][["outside"]], plots2save[["lh"]][[m[j]]][["female"]][[glob[gg]]][["outside"]])
#orca(pls[["rh"]][[m[j]]][["female"]][[glob[gg]]][["inside"]], plots2save[["rh"]][[m[j]]][["female"]][[glob[gg]]][["inside"]])
#orca(pls[["rh"]][[m[j]]][["female"]][[glob[gg]]][["outside"]], plots2save[["rh"]][[m[j]]][["female"]][[glob[gg]]][["outside"]])

#orca(pls[["lh"]][[m[j]]][["male"]][[glob[gg]]][["inside"]], plots2save[["lh"]][[m[j]]][["male"]][[glob[gg]]][["inside"]])
#orca(pls[["lh"]][[m[j]]][["male"]][[glob[gg]]][["outside"]], plots2save[["lh"]][[m[j]]][["male"]][[glob[gg]]][["outside"]])
#orca(pls[["rh"]][[m[j]]][["male"]][[glob[gg]]][["inside"]], plots2save[["rh"]][[m[j]]][["male"]][[glob[gg]]][["inside"]])
#orca(pls[["rh"]][[m[j]]][["male"]][[glob[gg]]][["outside"]], plots2save[["rh"]][[m[j]]][["male"]][[glob[gg]]][["outside"]])
#https://github.com/plotly/orca#installation

#plots= list()
#for (kk in seq(1,length(unlist(plots2save)))){
#  plots[[kk]] <- readPNG(paste(plot_save_dir,"/",str_split(unlist(plots2save)[kk],".png")[[1]][1],"_1.png",sep=""))
#} 


#top_title = paste(m[j],contrast,"contrast for models",glob[gg])

#ONLY RUN THIS WHEN ALL THE ORCA HAVE FINISHED
#png(paste(m[j],"_","brains_combined_zoom","_",contrast,"_",glob[gg],"_",".png",sep=""),width = 2000, height = 1200,)
#print( 
#  grid.arrange(rasterGrob(plots[[1]]),
#               rasterGrob(plots[[2]]),
#               rasterGrob(plots[[5]]),
#               rasterGrob(plots[[6]]), 
#               rasterGrob(plots[[3]]),
#               rasterGrob(plots[[4]]),
#               rasterGrob(plots[[7]]),
#               rasterGrob(plots[[8]]),
#               ncol=4,
               #bottom = "Something something",
#               top=textGrob(top_title,gp=gpar(fontsize=20)))
#)

#grid.text("Female", x = unit(0.5, "npc"), 
#          y = unit(.9, "npc"),gp = gpar(fontsize=30))
#grid.text("Male", x = unit(0.5, "npc"), 
#          y = unit(.4, "npc"),gp = gpar(fontsize=30))#, fontfamily="Times New Roman"))


#grid.text("Left                                                                                                                                                                                         Right", x = unit(0.5, "npc"), 
#          y = unit(.1, "npc"),gp = gpar(fontsize=20))
#grid.text("Hemisphere", x = unit(0.5, "npc"), 
#          y = unit(.07, "npc"),gp = gpar(fontsize=20))#, fontfamily="Times New Roman"))


#grid.text("Left                                                                                                                                                                                         Right", x = unit(0.5, "npc"), 
#          y = unit(.6, "npc"),gp = gpar(fontsize=20))
#grid.text("Hemisphere", x = unit(0.5, "npc"), 
#          y = unit(.57, "npc"),gp = gpar(fontsize=20))#, fontfamily="Times New Roman"))

#dev.off()




#### ggseg attempts 


#library(ggseg)
#library(patchwork)
#library(ggplot2)

#require(devtools)
#install_version("ggseg", version = "1.6.3", repos = "http://cran.us.r-project.org")

#remotes::install_github("ggseg/ggseg")

#plot(dk)
#dev.off()
#dev.set(5)
#ggseg()

#x11(type='Xlib')
#set_Polypath(FALSE)
#ggseg(usePolypath=FALSE)


#atlas_data = tibble(
#  subject = rep(c("bert"), 10),
#  label = rep(c("lh_b","lh_c","lh_e","lh_f","lh_i"), 2),
#  ThickAvg = sample(seq(0,.5,.001), 10),
#)

#model_yvars = col_names_list[[paste("rh","volume",sep = "_")]] 

#someData = tibble(
#  subject = seq(1,34,1),
#  label = model_yvars,
#  volume = sample(seq(0,.5,.001), 34),
#)



# Figure 4A
#fig4A <- someData %>%
#  # plot only subject bert
#  ggseg3d(mapping = aes(fill = volume)) +
#  labs(title = "from read_freesurfer_table", fill = "volume\n(mm∧3)") +
#  scale_fill_gradient(low = "firebrick", high = "goldenrod") +
#  guides(fill = guide_colourbar(barheight = 3))

# Figure 4B
#fig4B <- atlas_data %>%
#  # plot only subject bert
#  #filter(subject == "bert") %>%
#  ggseg(mapping = aes (fill = ThickAvg)) +
#  labs(title = "from read_atlas_files",
#       fill = "Thickness\n (mean)") +
#  scale_fill_gradient(low = "forestgreen", high = "goldenrod") +
#  guides(fill = guide_colourbar(barheight= 3))

#(fig4A / fig4B) +
#plot_annotation(tag_levels = "a")

#png(file="/mnt/projects/VIA11/FREESURFER/Stats/Plots/lateral/brains.png")
#fig4A +
#  plot_annotation()
#dev.off()




#ggseg3d(.data = someData,
#        atlas = dk_3d,
#        colour = "p", text = "p",
#        surface = "inflated",
#        hemisphere = c("left"),
#        palette = c("#ff0000"=min_cohens_D, "#ffffff"=0,"#0000ff"=max_cohens_D),
#        options.legend=list(title=list(text="Cohens D"))) %>%
#  pan_camera("left lateral") %>%
#  remove_axes() %>%
#  add_annotations(text = paste(m[j],contrast,"for",h[i],"of",sex[ss],glob[gg]),
#                  y = 1,
#                  yanchor = 'top',
#                  legendtitle = TRUE, showarrow = FALSE,
#                  font = list(family = 'sans serif',size = 20))



#sanity check female SZ contrast fdr corrections

#FDR correct p-values 
#ss=2
#i=1
#j = 3 #measure
#gg = 1 
#contrast = "SZ-K"

contrast_table_with_glob <-  read_excel(contrast_with_glob)
contrast_table_without_glob <-  read_excel(contrast_without_glob)
contrast_table = rbind(contrast_table_with_glob, contrast_table_without_glob)
ctable_both_hemi_sex <- contrast_table[c(contrast_table$measure == m[j] & contrast_table$sex == ss-1 & contrast_table$global_var_in_model == gg-1),]
contrast_pval_coloumn = paste("pval_",contrast[cc],sep="")
raw_pvals = as.numeric(ctable_both_hemi_sex[[contrast_pval_coloumn]])
fdr_pvals = p.adjust(raw_pvals, method = "fdr", n = length(raw_pvals))
ctable_both_hemi_sex[paste("fdr_pvals_",contrast[cc],sep="")] = fdr_pvals
table_sex <- eff_table[c(eff_table$hemisphere == h[i] & eff_table$measure == m[j] & eff_table$sex == sex[ss] & eff_table$globvar_in_model == gg-1),]

corrected_contrast_pval = ctable_both_hemi_sex[paste("fdr_pvals_",contrast[cc],sep="")][c(ctable_both_hemi_sex$hemisphere == h[i]),]
corrected_contrast_pval
#sex 0
#rh_caudalmiddlefrontal_area
#lh_lateralorbitofrontal_volume
#lh_postcentral_volume
#lh_precentral_volume
#lh_superiorfrontal_volume

#rh_caudalmiddlefrontal_volume
#rh_lateralorbitofrontal_volume***
#rh_medialorbitofrontal_volume
#rh_precentral_volume
#rh_precuneus_volume
#rh_superiorfrontal_volume



