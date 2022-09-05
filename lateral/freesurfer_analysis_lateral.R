#author: Simon Yamazaki Jensen

#load packages
library(data.table)
library(ggplot2)
library(ggdist)
library(gridExtra)
library(tidyr)
library(lsmeans)
library(grid)
#library(ggnewscale)
library(writexl)
library(car)
library(NCmisc)
library(lsr)
library(readxl)



#This script is intended for computation of statistics on lateral regional brain 
#volumes in multiple groups and comparing to a control group on cortical structures
#The script generates:


#####
# - GS_ANOVA tables with models that include a group/sex interaction 
#   an extra row of a model without the interaction is included if the iteraction 
#   turned out insignificant
#   saved in an excel sheet with a global covariate and an excel sheet without

#save paths:
GS_ANOVA_with_glob = "/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/parcels/lateral/lateral_Parcel_GS_ANOVA_pvals_with_glob.xlsx"
GS_ANOVA_without_glob = "/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/parcels/lateral/lateral_Parcel_GS_ANOVA_pvals_without_glob.xlsx"

#####
# - ANOVA tables with models that are defined on sex separated data
#   one row for each model on each sex
#   saved in an excel sheet with a global covariate and an excel sheet without

#save paths
ANOVA_with_glob = "/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/parcels/lateral/lateral_Parcel_S_ANOVA_pvals_with_glob.xlsx"
ANOVA_without_glob = "/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/parcels/lateral/lateral_Parcel_S_ANOVA_pvals_without_glob.xlsx"


#####
# - An excel sheet with model relevant LSmean contrasts for each of the models
#   in the sex separated ANOVA tables
#   saved in an excel sheet with a global covariate and an excel sheet without

#save paths:
contrast_with_glob ="/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/parcels/lateral/lateral_Parcel_model_contrast_with_glob.xlsx"
contrast_without_glob = "/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/parcels/lateral/lateral_Parcel_model_contrast_without_glob.xlsx"


effect_sizes_path = "/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/parcels/lateral/ANOVA+contrast_effect_sizes.xlsx"

plot_save_dir = "/mnt/projects/VIA11/FREESURFER/Stats/Plots/lateral"



######
#data path
#data_path = "/mnt/projects/VIA11/FREESURFER/Stats/Data/VIA11_allkey_160621_FreeSurfer_pruned_20220126.csv"
data_path = "/mnt/projects/VIA11/FREESURFER/Stats/Data/VIA11_allkey_160621_FreeSurfer_pruned_20220509.csv"

#load data
data_csv <- read.table(data_path, header = TRUE, sep = ",", dec = ".")


#inspect the head of data and summary
head(data_csv)
summary(data_csv)

#filter the data with include variable
# - extract rows with 1 in Include_FS_studies coloumn
#data_csv_filtered <- data_csv[c(data_csv$Include_FS_studies == 1),]
#data_csv_filtered <- data_csv_filtered[!is.na(data_csv_filtered$Include_FS_studies),]

# - extract rows with 1 in Include_FS_studies_euler_outliers_excluded
#data_csv_filtered <- data_csv_filtered[c(data_csv_filtered$Include_FS_studies_euler_outliers_excluded == 1),]
#data_csv_filtered <- data_csv_filtered[!is.na(data_csv_filtered$Include_FS_studies_euler_outliers_excluded),]

# - extract rows with 1 in Include_FS_studies_euler_outliers_sibpairs_out
data_csv_filtered <- data_csv[c(data_csv$Include_FS_studies_euler_outliers_sibpairs_out == 1),]
data_csv_filtered <- data_csv_filtered[!is.na(data_csv_filtered$Include_FS_studies_euler_outliers_sibpairs_out),]


#rename to a shorter name for convenience
datab = data_csv_filtered

#set working directory for plots to be saved
setwd("/mnt/projects/VIA11/FREESURFER/Stats/Plots/lateral")


general_cov = c("age","site","TotalEulerNumber")
general_cov_formula = paste(general_cov,collapse = "+")
contrast_signs = c(1,1, -1) # for K/BP/SZ 



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
    col_names_list[[name]] = col_names_list[[name]][!grepl(paste("MeanThickness_thickness",'$',sep = ""), col_names_list[[name]])] #remove WhiteSurfaceArea
    model_names = c(model_names, col_names_list[[name]])
  }
}

#make new variables with shorter and contained names
# - tell r which variables are factors
datab$sex = as.factor(datab$Sex_child)
datab$group = as.factor(datab$HighRiskStatus_v11)
datab$site = as.factor(datab$MRI_site_v11)
datab$diag = as.factor(datab$ksads_any_diag_excl_elim_lft_v11) #Axis.1_diag_v11) 
datab$age = as.numeric(datab$MRI_age)

#get a separate data frame for each sex
data_sex1 = datab[c(datab$sex == 1),]
data_sex0 = datab[c(datab$sex == 0),]




###### Make inference with models
# generate ANOVA table with GS interaction 

#initialize some variables for referencing
glob = c("without_global_var","with_global_var")
DF = data.frame()
model_list = list()


#loop that 
for (j in seq(1,3)){
  for (i in seq(1,2)){
  print(h[i])
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
      model_list[[glob[g]]][[model_yvars[k]]][[h[i]]] = lm(f,data=datab)
      model = model_list[[glob[g]]][[model_yvars[k]]][[h[i]]]
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
      
      model_list[[glob[g]]][[model_yvars[k]]][[h[i]]] = model
      
      for (mi in seq(1,length(models))){
        mm = models[[mi]]
        subs = length(mm$residuals)
        
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
                    sign_GS, mi, g-1,subs)
        }       
        else{
          rw = list(model_yvars[k], h[i], group_F, group_pv, 
                    sex_F, sex_pv, age_F, age_pv, 
                    site_F, site_pv, 
                    EulerNumber_F, EulerNumber_pv, 
                    glob_var_F, glob_var_pv, glob_var,
                    NA, NA,
                    sign_GS, mi-2, g-1,subs)
        } 
        
        DF = rbindlist(list(DF, rw))
        
      } #for m
    } #g
    } #i
    
  } #k
} #j


names(DF)<-c("model_yvar","hemisphere","group_Fval","group_pval",
             "sex_Fval","sex_pval","age_Fval","age_pval",
             "site_Fval","site_pval",
             "Eulernumber_Fval","Eulernumber_pval",
             "global_var_Fval","global_var_pval", "global_var_name",
             "Group_sex_Fval", "Group_sex_pval",
             "Significant_GS_interaction", "GS_in_model","global_var_in_model","n_subs")


DF_xlsx_glob0 = DF[DF$global_var_in_model == 0, ]
DF_xlsx_glob1 = DF[DF$global_var_in_model == 1, ]

write_xlsx(DF_xlsx_glob1,GS_ANOVA_with_glob)
write_xlsx(DF_xlsx_glob0,GS_ANOVA_without_glob)




###### make inference with models

DF_xlsx = data.frame()
DFc = data.frame()
glob = c("without_global_var","with_global_var")
sex = c("female","male","both") #3 sexes are defined, 0=female, 1=male, 2=both
dataf = list(data_sex0,data_sex1,datab)

models = list()


for (j in seq(1,3)){
  for (i in seq(1,2)){
  print(h[i])
  print(m[j])
  model_yvars = col_names_list[[paste(h[i],m[j],sep = "_")]] 
  #model_yvars = unlist(col_names_list[[m[j]]] )
  
  for (k in seq(1,length(model_yvars))){
    for (s in seq(1,3)){
    for (g in seq(1,2)){
      
      if (glob[g] == "without_global_var"){
        #no global measure model
        f = paste(model_yvars[k],"~","group","+","age","+","site","+","TotalEulerNumber")
        glob_var = NA 
      }
      else {
        #define model formulation for models with the global variable 
        if (m[j] == "area"){
          glob_var = "total_area"
          f = paste(model_yvars[k],"~","+","group","+","age","+","site","+","TotalEulerNumber","+",glob_var)
        }
        else if (m[j] == "thickness"){
          glob_var = "mean_thickness"
          f = paste(model_yvars[k],"~","+","group","+","age","+","site","+","TotalEulerNumber","+",glob_var)
        }
        else if (m[j] == "volume"){
          glob_var = "BrainTotalVol"
          f = paste(model_yvars[k],"~","+","group","+","age","+","site","+","TotalEulerNumber","+",glob_var)
        }
      }
      
      if (sex[s] == "both"){
        f = paste(f,"+sex",sep = "")
      }
      
      models[[sex[s]]][[glob[g]]][[model_yvars[k]]][[h[i]]] = lm(f,data=dataf[[s]])
      model = models[[sex[s]]][[glob[g]]][[model_yvars[k]]][[h[i]]]
      xvars = attributes(Anova(model,type = "III"))$row.names
      subs = length(model$residuals)
      
      #use the above defined model_ana 
      ls = lsmeans(model,pairwise~"group",adjust="none")
      c = ls$contrasts
      diff_df = attributes(c)$dfargs$df
      
      BP_emm = summary(ls)$lsmeans$lsmean[1]
      K_emm = summary(ls)$lsmeans$lsmean[2]
      SZ_emm = summary(ls)$lsmeans$lsmean[3]
      
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
      rwc = list(model_yvars[k],m[j],BP_emm,K_emm,SZ_emm,
                      BP_diff,BP_diff_tratio,BP_diff_pv,
                      BP_diff_LCL, BP_diff_UCL,
                      SZ_diff,SZ_diff_tratio,SZ_diff_pv,
                      SZ_diff_LCL, SZ_diff_UCL,
                      diff_df,g-1,s-1,h[i])
      
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
                     g-1,s-1,subs)
      
      DF_xlsx = rbindlist(list(DF_xlsx, rw_xlsx))
      
      
    } #g
    }#i
    } #s 
  } #k
} #j


names(DFc) = c("Model_yvar","measure","BP_LSmean","K_LSmean","SZ_LSmean",
                    "Contrast_BP-K","tratio_BP-K","pval_BP-K",
                    "LCL_Contrast_BP-K", "UCL_Contrast_BP-K",
                    "Contrast_SZ-K","tratio_SZ-K","pval_SZ-K",
                    "LCL_Contrast_SZ-K", "UCL_Contrast_SZ-K",
                    "DOF_ttest","global_var_in_model","sex","hemisphere")

names(DF_xlsx)<-c("Model_yvar","hemisphere","Group_Fval","Group_pval",
                  "Age_Fval","Age_pval","Site_Fval","Site_pval",
                  "EulerNumber_Fval","Eulernumber_pval",
                  "global_var_F","global_var_pv","global_var_name",
                  "global_var_in_model","sex","n_subs")

DF_glob0 = DF_xlsx[DF_xlsx$global_var_in_model == 0, ]
DF_glob1 = DF_xlsx[DF_xlsx$global_var_in_model == 1, ]

write_xlsx(DF_glob1,ANOVA_with_glob)
write_xlsx(DF_glob0,ANOVA_without_glob)

DFc_glob0 = DFc[DFc$global_var_in_model == 0, ]
DFc_glob1 = DFc[DFc$global_var_in_model == 1, ]

write_xlsx(DFc_glob1,contrast_with_glob)
write_xlsx(DFc_glob0,contrast_without_glob)






## Effect sizes


glob = c("without_global_var","with_global_var")
sex = c("female","male","both") #3 sexes are defined, 0=female, 1=male, 2=both
glob_vars = c("BrainTotalVol","mean_thickness","total_area")


DF = data.frame()

#all non-changing variable names in models
all_var_names = c("group","sex",general_cov, glob_vars, "group:sex")


for (j in seq(1,3)){
  for (i in seq(1,2)){
    print(h[i])
    print(m[j])
    model_yvars = col_names_list[[paste(h[i],m[j],sep = "_")]]
for (k in seq(1,length(model_yvars))){
  for (g in seq(1,length(glob))){
    for (s in seq(1,length(sex))){
      
      if (sex[s] == "both"){
        model_eff = model_list[[glob[g]]][[model_yvars[k]]][[h[i]]]
        GS_in_model = ("group:sex" %in% attributes(Anova(model_eff,type = "III"))$row.names)*1
        model_type = "group + sex"
        if (GS_in_model == 1){
          model_type = "group*sex"
        }
      }
      else {
        model_type = NA
        model_eff = models[[sex[s]]][[glob[g]]][[model_yvars[k]]][[h[i]]]
      }
      
      GS_in_model = ("group:sex" %in% attributes(Anova(model_eff,type = "III"))$row.names)*1
      
      
      #get model terms 
      xvars = attributes(etaSquared(model_eff, type = 3, anova = T))$dimnames[[1]]
      
      #get etaSquared outputs
      yvars = attributes(etaSquared(model_eff, type = 3, anova = T))$dimnames[[2]]
      
      eta_vals = as.numeric(etaSquared(model_eff, type = 3, anova = T)[,yvars="eta.sq.part"])
      
      all_col_names = all_var_names
      
      var_idx = all_col_names %in% xvars
      var_idx2 = xvars %in% all_col_names
      row = rep(NA,length(all_col_names))
      xt = xvars[var_idx2]
      at = all_col_names[var_idx]
      et = eta_vals[var_idx2]
      row[var_idx] = et[sapply(at, function(x) { grep(paste("^",x,"$", sep=""), xt)})]
      
      
      #cohens D for contrasts
      emm_eff = emmeans(model_eff,specs="group")
      effs = eff_size(emm_eff, sigma=sigma(model_eff), edf=df.residual(model_eff))
      
      contrast_names = attributes(effs)$levels$contrast
      c_eff_vals = rep(NA,length(contrast_names))
      
      for (cc in seq(1,length(contrast_names))){
        c_eff_vals[cc] = contrast_signs[cc]*summary(effs)$effect.size[cc]
      }
      
      #rows for effect size xlsx table
      DF = rbindlist(list(DF, as.list(c(model_yvars[k],m[j],h[i],row,
                                        c_eff_vals,
                                        g-1,GS_in_model,sex[s],model_type))))
      
      
    } #end s
  } #end g
} #end k
} #end i
} #end j
paES_col_names = paste(all_var_names, "_par_eta_sq",sep="")


#flip the order of contrast names are defined if the sign is flipped
contrast_names = attributes(c)$levels$contrast

for (cn in seq(1,length(contrast_names))){
  if (contrast_signs[cn]<0){
    contrast_names[cn]=paste(strsplit(contrast_names[cn], split = " - ")[[1]][2],strsplit(contrast_names[cn], split = " - ")[[1]][1],sep = " - ")
  }    
} 

c_cols = expand.grid(c("CohensD"), contrast_names)
c_cols = paste(c_cols$Var1, c_cols$Var2,sep = "_")
c_cols = gsub(" ", "", c_cols)


names(DF) = c("model_yvar","measure","hemisphere",paES_col_names,
              c_cols,
              "globvar_in_model","GS_in_model","sex","model_type")

write_xlsx(DF,effect_sizes_path)









###   sanity checks 
model2 = lm(lh_inferiorparietal_volume ~ group*sex + age + TotalEulerNumber + site + BrainTotalVol, data=datab)
#model_IPvol = update(model_IPvol,~.-group:sex)
Anova(model2,type="III")

model0 = lm(lh_inferiorparietal_volume ~ group + age + TotalEulerNumber + site + BrainTotalVol, data=data_sex0)
model1 = lm(lh_inferiorparietal_volume ~ group + age + TotalEulerNumber + site, data=data_sex1)
Anova(model0,type="III")
Anova(model1,type="III")

lsmeans(model0,pairwise~"group",adjust="none")





