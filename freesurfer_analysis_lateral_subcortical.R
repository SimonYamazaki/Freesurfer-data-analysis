
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
#volumes in multiple groups and comparing to a control group on subcortical structures
#The script generates:


#####
# - GS_ANOVA tables with models that include a group/sex interaction 
#   an extra row of a model without the interaction is included if the iteraction 
#   turned out insignificant
#   saved in an excel sheet with a global covariate and an excel sheet without

#save paths:
GS_ANOVA_with_glob = "/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/parcels/lateral/subParcel_GS_ANOVA_pvals_with_glob.xlsx"
GS_ANOVA_without_glob = "/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/parcels/lateral/subParcel_GS_ANOVA_pvals_without_glob.xlsx"

#####
# - ANOVA tables with models that are defined on sex separated data
#   one row for each model on each sex
#   saved in an excel sheet with a global covariate and an excel sheet without

#save paths
ANOVA_with_glob = "/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/parcels/lateral/subParcel_ANOVA_pvals_with_glob.xlsx"
ANOVA_without_glob = "/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/parcels/lateral/subParcel_ANOVA_pvals_without_glob.xlsx"


#####
# - An excel sheet with model relevant LSmean contrasts for each of the models
#   in the sex separated ANOVA tables
#   saved in an excel sheet with a global covariate and an excel sheet without

#save paths:
contrast_with_glob ="/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/parcels/lateral/subParcel_model_contrast_with_glob.xlsx"
contrast_without_glob = "/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/parcels/lateral/subParcel_model_contrast_without_glob.xlsx"



######
#data path
data_path = "/mnt/projects/VIA11/FREESURFER/Stats/Allthevolumeswithz.csv"

#load data
data_csv <- read.table(data_path, header = TRUE, sep = ",", dec = ".")
datab = data_csv


######################################
#    run models on each region
######################################

#### extract model vars ####

#extract model response variables, i.e. LHS of the model formula 
col_names = names(data_csv)
model_yvars = col_names[grepl(paste("_volume",'$',sep = ""),col_names)]

#remove some specific coloumns from the list above
global_var = "TotalVolume_mm_volume"
model_yvars = model_yvars[!(model_yvars %in% global_var)] 


#make new variables with shorter and contained names
# - tell r which variables are factors
datab$sex = as.factor(datab$Sex_child)
datab$group = as.factor(datab$HighRiskStatus)
datab$site = as.factor(datab$MR_Site)
datab$age = as.numeric(datab$MRI_age)

#get a separate data frame for each sex
data_sex1 = datab[c(datab$sex == 1),]
data_sex0 = datab[c(datab$sex == 0),]

#add them to a list for easy referencing later 
dataf = list(data_sex0,data_sex1)


###### Make inference with models
# generate GS_ANOVA tables 

#initialize some variables for referencing
glob = c("Without global var","With global var")
DF = data.frame()


for (k in seq(1,length(model_yvars))){
  for (g in seq(1,2)){
    
    if (glob[g] == "Without global var"){
      #define model formulation for models without the global variable 
      f = paste(model_yvars[k],"~","+","group*sex","+","age","+","site","+","TotalEulerNumber")
    }
    else {
      glob_var = global_var
      #define model formulation for models with the global variable 
      f = paste(model_yvars[k],"~","group*sex","+","age","+","site","+","TotalEulerNumber","+",glob_var)
    }
    
    #run the model
    model = lm(f,data=datab)
    xvars = attributes(Anova(model,type = "III"))$row.names
    
    #define the p-value of the group/sex interaction
    GS_F = Anova(model,type = "III")$"F value"[xvars=="group:sex"]
    GS_pv = Anova(model,type = "III")$"Pr(>F)"[xvars=="group:sex"]
    
    
    #if the group/sex interaction is insignificant 
    if (GS_pv > 0.05){
      #indicate the significance of the group/sex interation
      sign_GS = 0
      
      #save the model with group/sex interaction included
      model_gs = model
      
      #save the model without the group/sex interaction
      model = update(model,~.-group:sex)
      
      #save both models in a list
      models = list(model_gs,model)
    }
    else{
      #if the group/sex interaction is significant only add the model with the 
      #interaction to the list "models"
      sign_GS = 1
      models = list(model)
    }
    
    
    for (mi in seq(1,length(models))){
      #the list "models" has either 1 or 2 elements. 
      #it has 2 elements if the group/sex interaction was insignificant, where the
      #second element is the model without the interaction 
      mm = models[[mi]]
      
      #define statistics to be included for each model
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
        #when there is no global covariate in the mode, the below elements of the 
        #model row should be empty
        glob_var_F = NA
        glob_var_pv = NA
        glob_var = NA
      }       
      else{
        #global measure model
        #the variables should be defined if the model includes the global covariate
        glob_var_F = Anova(mm,type = "III")$"F value"[xvars==glob_var]
        glob_var_pv = Anova(mm,type = "III")$"Pr(>F)"[xvars==glob_var]
      } 
      
      
      if (mi == 1){ 
        #the first model in the list "models" is always the one that includes 
        #a group/sex interaction
        
        #define what goes into each row in the excel of models
        rw = list(model_yvars[k], group_F, group_pv, 
                  sex_F, sex_pv, age_F, age_pv, 
                  site_F, site_pv, 
                  EulerNumber_F, EulerNumber_pv, 
                  glob_var_F, glob_var_pv, 
                  GS_F, GS_pv, glob_var,
                  sign_GS, mi, g-1)
      }       
      else{
        #if this part of the if statement is reached then it is a model without 
        #a group/sex interaction
        rw = list(model_yvars[k], group_F, group_pv, 
                  sex_F, sex_pv, age_F, age_pv, 
                  site_F, site_pv, 
                  EulerNumber_F, EulerNumber_pv, 
                  glob_var_F, glob_var_pv,
                  NA, NA, glob_var,
                  sign_GS, mi-2, g-1)
      } 
      
      #append the row for the current model to the dataframe 
      DF = rbindlist(list(DF, rw))
      
    } #for m
    
  } #g
  
} #k

#define coloumn names, i.e. what each element in the model rows represent 
names(DF)<-c("model_yvar","group_Fval","group_pval",
             "sex_Fval","sex_pval","age_Fval","age_pval",
             "site_Fval","site_pval",
             "Eulernumber_Fval","Eulernumber_pval",
             "global_var_Fval","global_var_pval", 
             "Group_sex_Fval", "Group_sex_pval", "global_var_name",
             "Significant_GS_interaction", "GS_in_model","global_var_in_model")

#split the dataframe of model statistics into one excel with the global variable
#and one without
DF_xlsx_glob0 = DF[DF$global_var_in_model == 0, ]
DF_xlsx_glob1 = DF[DF$global_var_in_model == 1, ]

#write the dataframes to excel files
write_xlsx(DF_xlsx_glob1,GS_ANOVA_with_glob)
write_xlsx(DF_xlsx_glob0,GS_ANOVA_without_glob)




###### make inference with models
# generate ANOVA tables for sex separated models
# generate model_contrast tables for sex separated models

#initialize dataframes and the global variable indicator
DF_xlsx = data.frame()
DFc = data.frame()
glob = c("without_global_var","with_global_var")

#loop that generates the ANOVA table dataframe and model contrast dataframe
for (k in seq(1,length(model_yvars))){
  for (s in seq(1,2)){
    for (g in seq(1,2)){
      
      if (glob[g] == "without_global_var"){
        #no global measure model
        f = paste(model_yvars[k],"~","group","+","age","+","site","+","TotalEulerNumber")
      }
      else {
        glob_var = global_var
        #define model formulation for models with the global variable 
        f = paste(model_yvars[k],"~","+","group","+","age","+","site","+","TotalEulerNumber","+",glob_var)
      }
      
      #run model and define p-value for interaction
      model = lm(f,data=dataf[[s]])
      xvars = attributes(Anova(model,type = "III"))$row.names
      
      #use the above defined model_ana 
      ls = lsmeans(model,pairwise~"group",adjust="none")
      c = ls$contrasts
      BP_emm = summary(ls)$lsmeans$lsmean[1]
      K_emm = summary(ls)$lsmeans$lsmean[2]
      SZ_emm = summary(ls)$lsmeans$lsmean[3]
      
      #raw contrasts
      BP_diff = summary(c)$estimate[1]
      BP_diff_pv = summary(c)$"p.value"[1]
      BP_diff_tratio = summary(c)$"t.ratio"[1]
      BP_diff_LCL = confint(c)$lower.CL[1]
      BP_diff_UCL = confint(c)$upper.CL[1]
      
      SZ_diff = -summary(c)$estimate[3] #note the negative signs because the contrast is displayed as K - SZ
      SZ_diff_pv = summary(c)$"p.value"[3]
      SZ_diff_tratio = summary(c)$"t.ratio"[3]
      SZ_diff_LCL = -confint(c)$lower.CL[3]
      SZ_diff_UCL = -confint(c)$upper.CL[3]
      
      #rows to be added to the model contrast dataframe
      rwc = list(model_yvars[k],BP_emm,K_emm,SZ_emm,
                 BP_diff,BP_diff_tratio,BP_diff_pv,
                 BP_diff_LCL, BP_diff_UCL,
                 SZ_diff,SZ_diff_tratio,SZ_diff_pv,
                 SZ_diff_LCL, SZ_diff_UCL,
                 g-1,s-1)
      
      DFc = rbindlist(list(DFc, rwc))
      
      
      #extract statistics to go into ANOVA table dataframe
      group_F = Anova(model,type="III")$"F value"[xvars=="group"]
      age_F = Anova(model,type = "III")$"F value"[xvars=="age"]
      site_F = Anova(model,type = "III")$"F value"[xvars=="site"]
      EulerNumber_F = Anova(model,type = "III")$"F value"[xvars=="TotalEulerNumber"]
      
      group_pv = Anova(model,type = "III")$"Pr(>F)"[xvars=="group"]
      age_pv = Anova(model,type = "III")$"Pr(>F)"[xvars=="age"]
      site_pv = Anova(model,type = "III")$"Pr(>F)"[xvars=="site"]
      EulerNumber_pv = Anova(model,type = "III")$"Pr(>F)"[xvars=="TotalEulerNumber"]
      
      if (glob[g] == "without_global_var"){
        #no global measure model
        glob_var_F = NA
        glob_var_pv = NA
        glob_var = NA
      }       
      else{
        #global measure model
        glob_var_F = Anova(model,type = "III")$"F value"[xvars==glob_var]
        glob_var_pv = Anova(model,type = "III")$"Pr(>F)"[xvars==glob_var]
      } 
      
      #rows for ANOVA xlsx table, sex separated
      rw_xlsx = list(model_yvars[k], group_F, group_pv, 
                     age_F, age_pv, site_F, site_pv, 
                     EulerNumber_F, EulerNumber_pv, 
                     glob_var_F, glob_var_pv, glob_var,
                     g-1,s-1)
      
      #add the above row to the dataframe 
      DF_xlsx = rbindlist(list(DF_xlsx, rw_xlsx))
      
      
    } #g
  } #s 
} #k


#define coloumn names
names(DFc) = c("Model_yvar","BP_LSmean","K_LSmean","SZ_LSmean",
               "Contrast_BP-K","tratio_BP-K","pval_BP-K",
               "LCL_Contrast_BP-K", "UCL_Contrast_BP-K",
               "Contrast_SZ-K","tratio_SZ-K","pval_SZ-K",
               "LCL_Contrast_SZ-K", "UCL_Contrast_SZ-K",
               "global_var_in_model","sex")

names(DF_xlsx)<-c("Model_yvar","Group_Fval","Group_pval",
                  "Age_Fval","Age_pval","Site_Fval","Site_pval",
                  "EulerNumber_Fval","Eulernumber_pval",
                  "global_var_F","global_var_pv","global_var_name",
                  "global_var_in_model","sex")

#write a separate excel file for models with the global covariate and without
DF_glob0 = DF_xlsx[DF_xlsx$global_var_in_model == 0, ]
DF_glob1 = DF_xlsx[DF_xlsx$global_var_in_model == 1, ]

write_xlsx(DF_glob1,ANOVA_with_glob)
write_xlsx(DF_glob0,ANOVA_without_glob)

DFc_glob0 = DFc[DFc$global_var_in_model == 0, ]
DFc_glob1 = DFc[DFc$global_var_in_model == 1, ]

write_xlsx(DFc_glob1,contrast_with_glob)
write_xlsx(DFc_glob0,contrast_without_glob)




## Sanity checks some of the models run and displayed in tables 

model_hp = lm(Hypothalamus_whole_volume ~ group*sex + age + TotalEulerNumber + site, data=datab)
Anova(model_hp,type="III")

model_hpg = lm(Hypothalamus_whole_volume ~ group*sex + age + TotalEulerNumber + site + TotalVolume_mm_volume, data=datab)
Anova(model_hpg,type="III")


model0 = lm(Hypothalamus_whole_volume ~ group + age + site + TotalEulerNumber + TotalVolume_mm_volume, data=data_sex0)
model1 = lm(Hypothalamus_whole_volume ~ group + age + site + TotalEulerNumber + TotalVolume_mm_volume, data=data_sex1)
Anova(model0,type="III")
Anova(model1,type="III")

lsmeans(model0,pairwise~"group",adjust="none")


model0 = lm(Hypothalamus_whole_volume ~ group + age + site + TotalEulerNumber, data=data_sex0)
model1 = lm(Hypothalamus_whole_volume ~ group + age + site + TotalEulerNumber, data=data_sex1)
Anova(model0,type="III")
Anova(model1,type="III")

lsmeans(model0,pairwise~"group",adjust="none")
