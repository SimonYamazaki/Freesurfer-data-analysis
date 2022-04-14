#author: Simon Yamazaki Jensen

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



#This script is intended for computation of statistics on lateral regional brain 
#measures in multiple groups and comparing to a control group. The script generates:


#####
# - ANOVA tables with models that include a group/sex interaction 
#   an extra row of a model without the interaction is included if it turned out insignificant
#   saved in an excel sheet with a global covariate and an excel sheet without

#save paths:
GS_ANOVA_with_glob = "/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/parcels/lateral/Parcel_GS_ANOVA_pvals_with_glob.xlsx"
GS_ANOVA_without_glob = "/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/parcels/lateral/Parcel_GS_ANOVA_pvals_without_glob.xlsx"

ANOVA_with_glob = "/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/parcels/lateral/Parcel_ANOVA_pvals_with_glob.xlsx"
ANOVA_without_glob = "/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/parcels/lateral/Parcel_ANOVA_pvals_without_glob.xlsx"


#####
# - An excel sheet with model relevant contrasts for each of the models in the ANOVA tables
#   saved in an excel sheet with a global covariate and an excel sheet without

#save paths:
contrast_with_glob ="/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/parcels/lateral/Parcel_model_contrast_with_glob.xlsx"
contrast_without_glob = "/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/parcels/lateral/Parcel_model_contrast_without_glob.xlsx"


#####
# - Plots LSmeans based contrasts which compares BP and SZ to control

#save folder for plots:
save_folder = "/mnt/projects/VIA11/FREESURFER/Stats/Plots/lateral"

#prefix for the plot
LSmeans_prefix = "LSmean_difference_"


######
#data path
data_path = "/mnt/projects/VIA11/FREESURFER/Stats/Data/VIA11_allkey_160621_FreeSurfer_pruned_20220126.csv"

#load data
data_csv <- read.table(data_path, header = TRUE, sep = ",", dec = ".")


#inspect the head of data and summary
head(data_csv)
summary(data_csv)

#filter the data with include variable
# - extract rows with 1 in Include_FS_studies coloumn
data_csv_filtered <- data_csv[c(data_csv$Include_FS_studies == 1),]
data_csv_filtered <- data_csv_filtered[!is.na(data_csv_filtered$Include_FS_studies),]

# - extract rows with 1 in Include_FS_studies_euler_outliers_excluded
data_csv_filtered <- data_csv[c(data_csv$Include_FS_studies_euler_outliers_excluded == 1),]
data_csv_filtered <- data_csv_filtered[!is.na(data_csv_filtered$Include_FS_studies_euler_outliers_excluded),]

#rename to a shorter name for convenience
datab = data_csv_filtered

#set working directory for plots to be saved
setwd("/mnt/projects/VIA11/FREESURFER/Stats/Plots/lateral")




######################################
#    run models on each region
######################################

#### extract model yvars ####

#loop that extract all the relevant coloumns based on the regions variable
# extract model yvars 

col_names = names(datab)
col_names_list = list()
model_names = list()
h = c("lh","rh")
m = c("area","thickness","volume")
for (i in seq(1,2)){
  for (j in seq(1,3)){
    name = paste(h[i],m[j],sep="_")
    col_names_list[[name]] = col_names
    col_names_list[[name]] = col_names_list[[name]][grepl(paste("^",h[i],sep = ""), col_names_list[[name]])]
    col_names_list[[name]] = col_names_list[[name]][grepl(paste(m[j],'$',sep = ""), col_names_list[[name]])]
    col_names_list[[name]] = col_names_list[[name]][!grepl(paste("Area_area",'$',sep = ""), col_names_list[[name]])] #remove WhiteSurfaceArea
    model_names = c(model_names, col_names_list[[name]])
  }
}

#make new variables with shorter and contained names
# - tell r which variables are factors
datab$sex = as.factor(datab$Sex_child)
datab$group = as.factor(datab$HighRiskStatus)
datab$site = as.factor(datab$MR_Site)
datab$diag = as.factor(datab$Axis.1_diag_v11)
datab$age = as.numeric(datab$MRI_age)

#get a separate data frame for each sex
data_sex1 = datab[c(datab$sex == 1),]
data_sex0 = datab[c(datab$sex == 0),]
dataf = list(data_sex0,data_sex1)


###### Make inference with models

#initialize some variables for referencing
glob = c("Without global var","With global var")
DF = data.frame()

#loop that 
for (i in seq(1,2)){
  print(h[i])
for (j in seq(1,3)){
  print(m[j])
  model_yvars = col_names_list[[paste(h[i],m[j],sep = "_")]] 
  
  for (k in seq(1,length(model_yvars))){
    for (g in seq(1,2)){
      
      if (glob[g] == "Without global var"){
        #define model formulation for models without the global variable 
        f = paste(model_yvars[k],"~","+","group*sex","+","age","+","site","+","TotalEulerNumber")
      }
      else {
        #define model formulation for models with the global variable 
        if (m[j] == "area"){
          glob_var = "total_area"
          f = paste(model_yvars[k],"~","+","group*sex","+","age","+","site","+","TotalEulerNumber","+",glob_var)
        }
        else if (m[j] == "thickness"){
          glob_var = "mean_thickness"
          f = paste(model_yvars[k],"~","+","group*sex","+","age","+","site","+","TotalEulerNumber","+",glob_var)
        }
        else if (m[j] == "volume"){
          glob_var = "BrainTotalVol"
          f = paste(model_yvars[k],"~","+","group*sex","+","age","+","site","+","TotalEulerNumber","+",glob_var)
        }
      }
      
      #run the model
      model = lm(f,data=datab)
      xvars = attributes(Anova(model,type = "III"))$row.names
      
      #define the p-value of the group/sex interaction
      pv_group_sex = Anova(model,type = "III")$"Pr(>F)"[xvars=="group:sex"]
      
      
      #use the above defined model
      
      #if the group/sex interaction is insignificant 
      if (pv_group_sex > 0.05){
        #indicate the significance of the group/sex interation
        sign_GS = 0
        
        #save the model with group/sex interaction included
        model_gs = model
        
        #save the anova F-value for the group/sex interaction
        GS_F = Anova(model_gs,type = "III")$"F value"[xvars=="group:sex"]
        
        #save the model without the group/sex interaction
        model = update(model,~.-group:sex)
        
        #save both models in a list
        models = list(model_gs,model)
      }
      else{
        GS_F = Anova(model,type = "III")$"F value"[xvars=="group:sex"]
        sign_GS = 1
        models = list(model)
      }
      
      
      for (mi in seq(1,length(models))){
        mm = models[[mi]]
        group_F = Anova(mm,type = "III")$"F value"[xvars=="group"]
        sex_F = Anova(mm,type = "III")$"F value"[xvars=="sex"]
        age_F = Anova(mm,type = "III")$"F value"[xvars=="age"]
        site_F = Anova(mm,type = "III")$"F value"[xvars=="site"]
        EulerNumber_F = Anova(mm,type = "III")$"F value"[xvars=="TotalEulerNumber"]
        
        group_pv = Anova(mm,type = "III")$"Pr(>F)"[xvars=="group"]
        sex_pv = Anova(mm,type = "III")$"Pr(>F)"[xvars=="sex"]
        age_pv = Anova(mm,type = "III")$"Pr(>F)"[xvars=="age"]
        site_pv = Anova(mm,type = "III")$"Pr(>F)"[xvars=="site"]
        EulerNumber_pv = Anova(mm,type = "III")$"Pr(>F)"[xvars=="TotalEulerNumber"]
        
        if (glob[g] == "Without global var"){
          #no global measure model
          glob_var_F = NA
          glob_var_pv = NA
          glob_var = NA
        }       
        else{
          #global measure model
          glob_var_F = Anova(mm,type = "III")$"F value"[xvars==glob_var]
          glob_var_pv = Anova(mm,type = "III")$"Pr(>F)"[xvars==glob_var]
        } 
        
        
        if (mi == 1){
          rw = list(model_yvars[k], h[i],group_F, group_pv, 
                    sex_F, sex_pv, age_F, age_pv, 
                    site_F, site_pv, 
                    EulerNumber_F, EulerNumber_pv, 
                    glob_var_F, glob_var_pv, glob_var,
                    GS_F, pv_group_sex,
                    sign_GS, mi, g-1)
        }       
        else{
          rw = list(model_yvars[k], h[i], group_F, group_pv, 
                    sex_F, sex_pv, age_F, age_pv, 
                    site_F, site_pv, 
                    EulerNumber_F, EulerNumber_pv, 
                    glob_var_F, glob_var_pv, glob_var,
                    NA, NA,
                    sign_GS, mi-2, g-1)
        } 
        
        DF = rbindlist(list(DF, rw))
        
      } #for m
      
    } #g
    
  } #k
} #j
} #i

names(DF)<-c("model_yvar","hemisphere","group_Fval","group_pval",
             "sex_Fval","sex_pval","age_Fval","age_pval",
             "site_Fval","site_pval",
             "Eulernumber_Fval","Eulernumber_pval",
             "global_var_Fval","global_var_pval", "global_var_name",
             "Group_sex_Fval", "Group_sex_pval",
             "Significant_GS_interaction", "GS_in_model","global_var_in_model")


DF_xlsx_glob0 = DF[DF$global_var_in_model == 0, ]
DF_xlsx_glob1 = DF[DF$global_var_in_model == 1, ]

write_xlsx(DF_xlsx_glob1,GS_ANOVA_with_glob)
write_xlsx(DF_xlsx_glob0,GS_ANOVA_without_glob)




###### make inference with models
DF_xlsx = data.frame()
DFc = data.frame()
glob = c("without_global_var","with_global_var")

for (i in seq(1,2)){
  print(h[i])
for (j in seq(1,3)){
  print(m[j])
  model_yvars = col_names_list[[paste(h[i],m[j],sep = "_")]] 
  #model_yvars = unlist(col_names_list[[m[j]]] )
  for (k in seq(1,length(model_yvars))){
    for (s in seq(1,2)){
    for (g in seq(1,2)){
      
      if (glob[g] == "without_global_var"){
        #no global measure model
        f = paste(model_yvars[k],"~","group","+","age","+","site","+","TotalEulerNumber")
        model = lm(f,data=dataf[[s]])
        xvars = attributes(Anova(model,type = "III"))$row.names
        
      }
      else {
        #define model formulation for models with the global variable 
        if (m[j] == "area"){
          glob_var = "total_area"
          f2 = paste(model_yvars[k],"~","+","group","+","age","+","site","+","TotalEulerNumber","+",glob_var)
        }
        else if (m[j] == "thickness"){
          glob_var = "mean_thickness"
          f2 = paste(model_yvars[k],"~","+","group","+","age","+","site","+","TotalEulerNumber","+",glob_var)
        }
        else if (m[j] == "volume"){
          glob_var = "BrainTotalVol"
          f2 = paste(model_yvars[k],"~","+","group","+","age","+","site","+","TotalEulerNumber","+",glob_var)
        }
        #run model and define p-value for interaction
        model = lm(f2,data=datab)
        xvars = attributes(Anova(model,type = "III"))$row.names
      }
      
      
      #use the above defined model_ana 
      ls = lsmeans(model,pairwise~"group",adjust="none")
      c = ls$contrasts
      K_emm = summary(ls)$lsmeans$lsmean[2]
      
      #raw contrasts
      BP_diff = summary(c)$estimate[1]
      BP_diff_pv = summary(c)$"p.value"[1]
      BP_diff_tratio = summary(c)$"t.ratio"[1]
      BP_diff_LCL = confint(c)$lower.CL[1]
      BP_diff_UCL = confint(c)$upper.CL[1]
      
      
      SZ_diff = -summary(c)$estimate[3]
      SZ_diff_pv = summary(c)$"p.value"[3]
      SZ_diff_tratio = summary(c)$"t.ratio"[3]
      SZ_diff_LCL = -confint(c)$lower.CL[3]
      SZ_diff_UCL = -confint(c)$upper.CL[3]
      
      #for contrast xlsx
      rwc = list(model_yvars[k],m[j],
                      BP_diff,BP_diff_tratio,BP_diff_pv,
                      BP_diff_LCL, BP_diff_UCL,
                      SZ_diff,SZ_diff_tratio,SZ_diff_pv,
                      SZ_diff_LCL, SZ_diff_UCL,
                      K_emm,g-1,s-1,h[i])
      
      DFc = rbindlist(list(DFc, rwc))
      
      
      mm = model
      
      group_F = Anova(mm,type="III")$"F value"[xvars=="group"]
      age_F = Anova(mm,type = "III")$"F value"[xvars=="age"]
      site_F = Anova(mm,type = "III")$"F value"[xvars=="site"]
      EulerNumber_F = Anova(mm,type = "III")$"F value"[xvars=="TotalEulerNumber"]
      
      group_pv = Anova(mm,type = "III")$"Pr(>F)"[xvars=="group"]
      age_pv = Anova(mm,type = "III")$"Pr(>F)"[xvars=="age"]
      site_pv = Anova(mm,type = "III")$"Pr(>F)"[xvars=="site"]
      EulerNumber_pv = Anova(mm,type = "III")$"Pr(>F)"[xvars=="TotalEulerNumber"]
      
      if (glob[g] == "without_global_var"){
        #no global measure model
        glob_var_F = NA
        glob_var_pv = NA
      }       
      else{
        #global measure model
        glob_var_F = Anova(mm,type = "III")$"F value"[xvars==glob_var]
        glob_var_pv = Anova(mm,type = "III")$"Pr(>F)"[xvars==glob_var]
        
      } 
      
      #rows for ANOVA xlsx table, sex separated
      rw_xlsx = list(model_yvars[k],h[i], group_F, group_pv, 
                     age_F, age_pv, site_F, site_pv, 
                     EulerNumber_F, EulerNumber_pv, 
                     glob_var_F, glob_var_pv, glob_var,
                     g-1,s-1)
      
      DF_xlsx = rbindlist(list(DF_xlsx, rw_xlsx))
      
      
    } #g
    } #s 
  } #k
} #j
}#i


names(DFc) = c("Model_yvar","measure",
                    "Contrast_BP-K","tratio_BP-K","pval_BP-K",
                    "LCL_Contrast_BP-K", "UCL_Contrast_BP-K",
                    "Contrast_SZ-K","tratio_SZ-K","pval_SZ-K",
                    "LCL_Contrast_SZ-K", "UCL_Contrast_SZ-K",
                    "K_LSmean","global_var_in_model","sex","hemisphere")


names(DF_xlsx)<-c("Model_yvar","hemisphere","Group_Fval","Group_pval",
                  "Age_Fval","Age_pval","Site_Fval","Site_pval",
                  "EulerNumber_Fval","Eulernumber_pval",
                  "global_var_F","global_var_pv","global_var_name",
                  "global_var_in_model","sex")

DF_glob0 = DF_xlsx[DF_xlsx$global_var_in_model == 0, ]
DF_glob1 = DF_xlsx[DF_xlsx$global_var_in_model == 1, ]

write_xlsx(DF_glob1,ANOVA_with_glob)
write_xlsx(DF_glob0,ANOVA_without_glob)

DFc_glob0 = DFc[DFc$global_var_in_model == 0, ]
DFc_glob1 = DFc[DFc$global_var_in_model == 1, ]

write_xlsx(DFc_glob1,contrast_with_glob)
write_xlsx(DFc_glob0,contrast_without_glob)



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
groups = c("BP","SZ")
m = c("area","thickness","volume")
glob = c("Models WITHOUT global var","Models WITH global var")
sx = c("female","male")

sp=list()

for (j in seq(1,3)){
for (g in seq(1,2)){
for (i in seq(1,2)){  
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
} #g

top_title = paste("Lateral LSmean difference from control:",m[j],glob[g])
ps=grid.arrange(grobs=sp, top=textGrob(top_title,gp=gpar(fontsize=20)))
#ggsave(paste(LSmeans_prefix,m[j],".png",sep=""),ps,width = 15,height = 10)

} #j


