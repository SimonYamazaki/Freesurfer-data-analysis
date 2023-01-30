#author: Simon Yamazaki Jensen

#load packages
library(data.table)
library(ggplot2)
library(gridExtra)
library(lsmeans)
library(grid)
library(writexl)
library(car)
library(BayesFactor)
library(stringr)
library(lsr)

#funchir::stale_package_check('/mrhome/simonyj/Freesurfer-data-analysis/global_measures/global_measures_generalized.R')

#This script is intended for computation of statistics on global brain measures
#in multiple groups and comparing to a control group. The script generates:


#####
# - ANOVA tables with models that include a group/sex interaction 
#   an extra row of a model without the interaction is included if it turned out insignificant
#   saved in an excel sheet with results from models with covariate and an excel sheet without

#save paths:
GS_ANOVA_with_cov = "globvar_GS_ANOVA_pvals_with_glob.xlsx"
GS_ANOVA_without_cov = "globvar_GS_ANOVA_pvals_without_glob.xlsx"

ANOVA_with_cov = "globvar_S_ANOVA_pvals_with_glob.xlsx"
ANOVA_without_cov = "globvar_S_ANOVA_pvals_without_glob.xlsx"

#####
# - An excel sheet with model relevant contrasts for each of the models in the ANOVA tables
#   saved in an excel sheet with a covariate and an excel sheet without

#save paths:
contrast_with_cov = "globvar_S_Model_contrasts_with_glob.xlsx"
contrast_without_cov = "globvar_S_Model_contrasts_without_glob.xlsx"

#####
# - An excel sheet with model relevant effect sizes for each of the global models
#   saved in an excel sheet with a global covariate and an excel sheet without

effect_sizes_path = "ANOVA+contrast_effect_sizes.xlsx"

#file name postfix
ppwp_postfix = "_group_diff_pvalues_ICV"



#what data should the models be run for 
run_group = "K_BP_SZ" #choose from "K_BP_SZ", "BP_axis1", "SZ_axis1"
siblings = FALSE #if siblings should be included 
bf_iterations = 100000 #number of posterior samples used to compute bayesfactor


#data path
#data_path = "/mnt/projects/VIA11/FREESURFER/Stats/Data/VIA11_allkey_160621_FreeSurfer_pruned_20220126.csv"
data_path = "/mnt/projects/VIA11/FREESURFER/Stats/Data/VIA11_allkey_160621_FreeSurfer_pruned_20220509.csv"

#load data
data_csv <- read.table(data_path, header = TRUE, sep = ",", dec = ".")


#inspect the head of data and summary
head(data_csv)
summary(data_csv)


## Preprocess 

#filter the data with include variable
# - extract rows with 1 in Include_FS_studies coloumn
data_csv_filtered <- data_csv[c(data_csv$Include_FS_studies == 1),]
data_csv_filtered <- data_csv_filtered[!is.na(data_csv_filtered$Include_FS_studies),]


#preprocess and configure according to what data should be used
if (siblings){
  # - extract rows with 1 in Include_FS_studies_euler_outliers_excluded
  data_csv_filtered <- data_csv[c(data_csv$Include_FS_studies_euler_outliers_identified == 1),]
  data_csv_filtered <- data_csv_filtered[!is.na(data_csv_filtered$Include_FS_studies_euler_outliers_identified),]
  
  # what folder should the tables be saved in
  save_folder = paste("/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/global_variables/",run_group,"_with_siblings/",sep="")
  dir.create(file.path(dirname(save_folder), basename(save_folder)))
  
  #save folder for plots:
  plot_folder = paste("/mnt/projects/VIA11/FREESURFER/Stats/Plots/global_measures/",run_group,"_with_siblings/",sep="")
  dir.create(file.path(dirname(plot_folder), basename(plot_folder)))
  
  
}else {
  # - extract rows with 1 in Include_FS_studies_euler_outliers_sibpairs_out
  data_csv_filtered <- data_csv_filtered[c(data_csv_filtered$Include_FS_studies_euler_outliers_sibpairs_out == 1),]
  data_csv_filtered <- data_csv_filtered[!is.na(data_csv_filtered$Include_FS_studies_euler_outliers_sibpairs_out),]
  
  # what folder should the tables be saved in
  save_folder = paste("/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/global_variables/",run_group,"/",sep="")
  dir.create(file.path(dirname(save_folder), basename(save_folder)))
  
  #save folder for plots:
  plot_folder = paste("/mnt/projects/VIA11/FREESURFER/Stats/Plots/global_measures/",run_group,"/",sep="")
  dir.create(file.path(dirname(plot_folder), basename(plot_folder)))
  
}


#make new variables with shorter and contained names
# - tell r which variables are factors
datab = data_csv_filtered
datab$sex = as.factor(datab$Sex_child)
datab$site = as.factor(datab$MRI_site_v11)
datab$diag = as.factor(datab$ksads_any_diag_excl_elim_lft_v11) #Axis.1_diag_v11) 
datab$age = as.numeric(datab$MRI_age)


# Configure what groups should be used
# Define contrast sign
# check below what the contrast signs should be

if (run_group == "BP_axis1"){
  datab$group = datab$HRS_BP_K_axis1
  datab$group[datab$group == 3] = "K"
  datab$group[datab$group == 2] = "BP+"
  datab$group[datab$group == 1] = "BP-"
  datab$group = as.factor(datab$group)
  contrast_signs = c(1,1, 1) # for K/BP-/BP+
  }
if (run_group == "SZ_axis1"){
  datab$group = datab$HRS_SZ_K_axis1
  datab$group[datab$group == 3] = "K"
  datab$group[datab$group == 2] = "SZ+"
  datab$group[datab$group == 1] = "SZ-"
  datab$group = as.factor(datab$group)
  contrast_signs = c(-1,-1, 1) # for K/SZ-/SZ+ 
  }
if (run_group == "K_BP_SZ"){
  datab$group = as.factor(datab$HighRiskStatus_v11)
  contrast_signs = c(1,1, -1) # for K/BP/SZ 
}

#datab$group = as.factor(datab$HighRiskStatus_axis1_v11)

# - check how contrasts are made and how contrast signs should be
levels(datab$group)
m1 = lm(BrainTotalVol~group,data=datab)
lsmeans(m1,pairwise~"group", adjust="none")


#add a string encoding of sex for later convenience 
datab$sexs = datab$Sex_child
datab$sexs[datab$sexs == 0] = "female"
datab$sexs[datab$sexs == 1] = "male"

#set another working directory so that files are saved in this folder
setwd(save_folder)

#specify dependent variable (global variable) column names, i.e. left hand side of statistical models
y_vars = c("BrainTotalVol", "CortexVol", "total_area", "mean_thickness", "eICV_samseg")

#covariates that are included in all models 
#encoded as a formula string
general_cov = c("age","site","TotalEulerNumber")
general_cov_formula = paste(general_cov,collapse = "+")




##############################
# whole brain measure models #
##############################


### run models with eICV_samseg as covariate and without
# generate ANOVA tables with group/sex interaction
# each row is a modelin saved excel

#which global variables to model and compute statistics for
model_vars = y_vars

GS_pvals = list()
ps = list()
emm = list()
min_emm = list()
max_emm = list()
model = list()

DF = data.frame()

#loop that defines models 
for (i in seq(1,length(model_vars))){
  if (model_vars[i] == "eICV_samseg"){    #only run models without_eICV for the model with eICV as response variable
    glob = c("without_eICV")
  }
  else{
    glob = c("without_eICV","with_eICV") #run both for all other variables
  }
  
  for (g in seq(1,length(glob))){
    
    #define model formulas
    f_without_gs = paste(model_vars[i],"~","+group+sex","+",general_cov_formula,sep="")
    f = paste(f,"+group:sex",sep="")
    
    #add global variable 
    if (glob[g] == "with_eICV"){
      f = paste(f,"+eICV_samseg",sep="")
      f_without_gs = paste(f,"+eICV_samseg",sep="")
    }
    
    #run the actual model
    model[[glob[g]]][[model_vars[i]]] = lm(f,data=datab)  
    
    #save the anova F and p-value for the group/sex interaction
    GS_pv = Anova(model[[glob[g]]][[model_vars[i]]],type = "III")$"Pr(>F)"[xvars=="group:sex"]
    GS_F = Anova(model[[glob[g]]][[model_vars[i]]],type = "III")$"F value"[xvars=="group:sex"]
    
    #extract the variables that are modelling the global variable
    xvars = attributes(Anova(model[[glob[g]]][[model_vars[i]]],type = "III"))$row.names
    
    
    
    #if the group/sex interaction is insignificant 
    if (Anova(model[[glob[g]]][[model_vars[i]]],type = "III")$"Pr(>F)"[xvars=="group:sex"] > 0.05){
      #indicate the significance of the group/sex interation
      sign_GS = 0
      
      #save the model with group/sex interaction included
      model_gs = model[[glob[g]]][[model_vars[i]]]
      
      #save the model without the group/sex interaction
      model[[glob[g]]][[model_vars[i]]] = update(model[[glob[g]]][[model_vars[i]]],~.-group:sex)
      
      #save both models in a list
      models = list(model_gs,model[[glob[g]]][[model_vars[i]]])
      model_formula = list(f,f_without_gs)
      model_type = list("GS","No_interactions")
      
    }
    else{    #if the group/sex interaction is significant 
      
      #indicate the significance of the group/sex interation
      sign_GS = 1
      
      #save the model with the group/sex interaction
      models = list(model[[glob[g]]][[model_vars[i]]])
      model_formula = list(f)
      model_type = list("GS")
    }
    
    # loop that makes rows in the ANOVA table excel 
    # makes two rows if the group/sex interaction is insignificant 
    # one row with the interaction included and one row with out
    for (m in seq(1,length(models))){
      mm = models[[m]]
      subs = length(mm$residuals)
      mf = model_formula[[m]]
      
      #extract relevant statistics for each variable in the model 
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
      
      

      ############### Bayes Factor analysis 
      # for each variable in model
      
      #group
      bf1 = lmBF(formula = eval(parse(text=mf)), data=datab)
      bf2 = lmBF(formula = eval(parse(text=gsub("\\+group", "",mf))), data=datab)
      bf_interaction = recompute(bf1 / bf2, iterations = bf_iterations)
      bf_res = extractBF(bf_interaction, logbf = FALSE, onlybf = FALSE)
      bf_g = as.numeric(bf_res[1])
      bf_g_error = as.numeric(bf_res[2])
      
      # age
      bf1 = lmBF(formula = eval(parse(text=mf)), data=datab)
      bf2 = lmBF(formula = eval(parse(text=gsub("\\+age", "",mf))), data=datab)
      bf_interaction = recompute(bf1 / bf2, iterations = bf_iterations)
      bf_res = extractBF(bf_interaction, logbf = FALSE, onlybf = FALSE)
      bf_age = as.numeric(bf_res[1])
      bf_age_error = as.numeric(bf_res[2])
      
      # Eulernumer
      bf1 = lmBF(formula = eval(parse(text=mf)), data=datab)
      bf2 = lmBF(formula = eval(parse(text=gsub("\\+TotalEulerNumber", "",mf))), data=datab)
      bf_interaction = recompute(bf1 / bf2, iterations = bf_iterations)
      bf_res = extractBF(bf_interaction, logbf = FALSE, onlybf = FALSE)
      bf_euler = as.numeric(bf_res[1])
      bf_euler_error = as.numeric(bf_res[2])
      
      #Site
      bf1 = lmBF(formula = eval(parse(text=mf)), data=datab)
      bf2 = lmBF(formula = eval(parse(text=gsub("\\+site", "",mf))), data=datab)
      bf_interaction = recompute(bf1 / bf2, iterations = bf_iterations)
      bf_res = extractBF(bf_interaction, logbf = FALSE, onlybf = FALSE)
      bf_site = as.numeric(bf_res[1])
      bf_site_error = as.numeric(bf_res[2])
      
      #sex
      bf1 = lmBF(formula = eval(parse(text=mf)), data=datab)
      bf2 = lmBF(formula = eval(parse(text=gsub("\\+sex", "",mf))), data=datab)
      bf_interaction = recompute(bf1 / bf2, iterations = bf_iterations)
      bf_res = extractBF(bf_interaction, logbf = FALSE, onlybf = FALSE)
      bf_sex = as.numeric(bf_res[1])
      bf_sex_error = as.numeric(bf_res[2])
      
      #Group/sex BF
      bf1 = lmBF(formula = eval(parse(text=f)), data=datab)
      bf2 = lmBF(formula = eval(parse(text=f_without_gs)), data=datab)
      bf_interaction = recompute(bf1 / bf2, iterations = bf_iterations)
      bf_res = extractBF(bf_interaction, logbf = FALSE, onlybf = FALSE)
      bf_gs = as.numeric(bf_res[1])
      bf_gs_error = as.numeric(bf_res[2])
      
      ############# BF above 
      
      
      
      #if the model includes ICV, then also define the anova statistics of ICV
      if (glob[g] == "with_eICV"){
        ICV_F = Anova(mm,type = "III")$"F value"[xvars=="eICV_samseg"]
        ICV_pv = Anova(mm,type = "III")$"Pr(>F)"[xvars=="eICV_samseg"]
        
        #glob
        bf1 = lmBF(formula = eval(parse(text=mf)), data=datab)
        bf2 = lmBF(formula = eval(parse(text=gsub("\\+eICV_samseg", "",mf))), data=datab)
        bf_interaction = recompute(bf1 / bf2, iterations = bf_iterations)
        bf_res = extractBF(bf_interaction, logbf = FALSE, onlybf = FALSE)
        bf_glob = as.numeric(bf_res[1])
        bf_glob_error = as.numeric(bf_res[2])
      }
      else{ #else, dont define them if the model is without ICV
        ICV_F = NA
        ICV_pv = NA
        bf_glob = NA
        bf_glob_error = NA
      }       
      
      #if the model is with group/sex interaction 
      if (model_type[[m]]=="GS"){
        GS_in = 1
      }
      else if (model_type[[m]]=="No_interactions"){ 
        GS_F = NA
        GS_pv = NA
        GS_in = 0
        
        bf_gs = NA
        bf_gs_error = NA
      }
      
      #append the row to the dataframe 
      rw = c(model_vars[i],
             group_F, group_pv, bf_g, bf_g_error,
             sex_F, sex_pv, bf_sex, bf_sex_error,
             age_F, age_pv, bf_age, bf_age_error,
             site_F, site_pv, bf_site, bf_site_error,
             EulerNumber_F, EulerNumber_pv, bf_euler, bf_euler_error,
             ICV_F, ICV_pv, bf_glob, bf_glob_error,
             GS_F, GS_pv, bf_gs, bf_gs_error,
             sign_GS, GS_in, 
             g-1, model_type[[m]], subs)
      
      DF = rbindlist(list(DF, as.list(rw)))
      
    } # for m 
  } #g
} #end i

#define the coloum names of the dataframe which is converted to an excel 
col_names = c("Model_yvar",
              "Group_Fval","Group_pval","Group_bf", "Group_bf_error",
              "Sex_Fval","Sex_pval", "Sex_bf", "Sex_bf_error",
              "Age_Fval","Age_pval", "Age_bf", "Age_bf_error",
              "Site_Fval","Site_pval", "Site_bf", "Site_bf_error",
              "EulerNumber_Fval","Eulernumber_pval", "Eulernumber_bf", "Eulernumber_bf_error",
              "ICV_Fval","ICV_pval", "ICV_bf", "ICV_bf_error",
              "Group_sex_Fval","Group_sex_pval", "Group_sex_bf", "Group_sex_bf_error",
              "Significant_GS_interaction","GS_in_model",
              "ICV_in_model","model_type","n_subjects")

names(DF)<-col_names


#split up the dataframe in models with ICV and without
DF_xlsx_ICV0 = DF[DF$ICV_in_model == 0, ]
DF_xlsx_ICV1 = DF[DF$ICV_in_model == 1, ]

#save the two dataframes in the a specified save path 
write_xlsx(DF_xlsx_ICV1,paste(save_folder,GS_ANOVA_with_cov,sep = ""))
write_xlsx(DF_xlsx_ICV0,paste(save_folder,GS_ANOVA_without_cov,sep = ""))




#### CONTRAST XLSX with male / female split
# and built models[sex][glob][yvar]

#define some relevant naming 
model_vars = y_vars
sex = c("female","male","both") #3 sexes are defined, 0=female, 1=male, 2=both

#split up the data into male and female 
data_sex1 = datab[c(datab$sex == 1),]
data_sex0 = datab[c(datab$sex == 0),]
dataf = list(data_sex0, data_sex1, datab)

models = list()

DFc = data.frame()
DF_xlsx = data.frame()


# loop that runs the model on each global variable, each sex with and without ICV
for (k in seq(1,length(model_vars))){
  print(model_vars[k])
  
  if (model_vars[k] == "eICV_samseg"){ 
    glob = c("without_eICV")
  }
  else { 
    glob = c("without_eICV","with_eICV")
  }
  
  for (ss in seq(1,3)){ #use the 3 defined sexes to know what data to use
    
    # define an temporary sex indicator to index
    print(sex[ss])
    
    #extract the dataframe containing the relevant sex data
    df_sex = dataf[[ss]]
    
    
    for (gg in seq(1,length(glob))){
      
      if (glob[gg] == "without_eICV"){ #run models without ICV
        #define the models 
        if (sex[ss] == "both"){ #define a specific model for the case of both sex
          f = paste(model_vars[k],"~","+group","+","sex","+",general_cov_formula,sep="")
        }
        else { #another model if it is sex-separated
          f = paste(model_vars[k],"~","+group","+",general_cov_formula,sep="")
        }
        #run the model
        models[[sex[ss]]][[glob[gg]]][[model_vars[k]]] = lm(f,data=df_sex)
        print(glob[gg])
        
        #save the current model into a temporary variable to make code less cluttered
        model_ana = models[[sex[ss]]][[glob[gg]]][[model_vars[k]]]
        xvars = attributes(Anova(model_ana,type = "III"))$row.names
      }
      
      else { #run models with ICV
        #define models
        if (sex[ss] == "both"){
          f = paste(model_vars[k],"~","+group","+","sex","+",general_cov_formula,"+","eICV_samseg",sep="")
        }
        else {
          f = paste(model_vars[k],"~","+group","+",general_cov_formula,"+","eICV_samseg",sep="")
        }        
        models[[sex[ss]]][[glob[gg]]][[model_vars[k]]] = lm(f,data=df_sex)
        model_ana = models[[sex[ss]]][[glob[gg]]][[model_vars[k]]]
        xvars = attributes(Anova(model_ana,type = "III"))$row.names
      }
      
      
      #use the above defined model_ana 
      #compute lsmeans pairwise for each group with no correction on p-values
      ls = lsmeans(model_ana,pairwise~"group", adjust="none")
      
      #compute contrast statistics
      c = ls$contrasts
      
      group_levels = attributes(ls$lsmeans)$levels$group
      stopifnot(length(contrast_signs) == length(group_levels))
      
      contrast_names = attributes(c)$levels$contrast
      
      rwc_xlsx = list(model_vars[k])
      
      for (gl in seq(1,length(group_levels))){
        emm = summary(ls)$lsmeans$lsmean[gl]
        rwc_xlsx = append(rwc_xlsx,list(emm))
      }
      
      for (cc in seq(1,length(contrast_names))){
        diff = contrast_signs[cc]*summary(c)$estimate[cc]
        diff_pv = summary(c)$"p.value"[cc]
        diff_tratio = summary(c)$"t.ratio"[cc]
        diff_LCL = contrast_signs[cc]*confint(c)$lower.CL[cc]
        diff_UCL = contrast_signs[cc]*confint(c)$upper.CL[cc]
        
        cgroup1 = str_split(contrast_names[cc]," - ")[[1]][1]
        cgroup2 = str_split(contrast_names[cc]," - ")[[1]][2]
        df_cgroups_sex = df_sex[df_sex$group==cgroup1 | df_sex$group==cgroup2,]
        
        bf1 = lmBF(formula = eval(parse(text=f)), data=df_cgroups_sex)
        bf2 = lmBF(formula = eval(parse(text=gsub("\\+group", "",f))), data=df_cgroups_sex)
        bf = recompute(bf1 / bf2, iterations = bf_iterations)
        bf_res = extractBF(bf, logbf = FALSE, onlybf = FALSE)
        bf_contrast = as.numeric(bf_res[1])
        bf_contrast_error = as.numeric(bf_res[2])
        
        rwc_xlsx = append(rwc_xlsx,list(diff, diff_tratio, diff_pv,
                                        bf_contrast,bf_contrast_error,
                                        diff_LCL,diff_UCL))

      }
      diff_df = attributes(c)$dfargs$df
      rwc_xlsx = append(rwc_xlsx,list(diff_df,gg-1,ss-1))
      DFc = rbindlist(list(DFc, rwc_xlsx))
      
      if (sex[ss] != "both"){
        mm=model_ana
        
        #ANOVA sex separated
        group_F = Anova(mm,type="III")$"F value"[xvars=="group"]
        age_F = Anova(mm,type = "III")$"F value"[xvars=="age"]
        site_F = Anova(mm,type = "III")$"F value"[xvars=="site"]
        EulerNumber_F = Anova(mm,type = "III")$"F value"[xvars=="TotalEulerNumber"]
        
        group_pv = Anova(mm,type = "III")$"Pr(>F)"[xvars=="group"]
        age_pv = Anova(mm,type = "III")$"Pr(>F)"[xvars=="age"]
        site_pv = Anova(mm,type = "III")$"Pr(>F)"[xvars=="site"]
        EulerNumber_pv = Anova(mm,type = "III")$"Pr(>F)"[xvars=="TotalEulerNumber"]
        
        #group BF
        bf1 = lmBF(formula = eval(parse(text=f)), data=df_sex)
        bf2 = lmBF(formula = eval(parse(text=gsub("\\+group", "",f))), data=df_sex)
        bf = recompute(bf1 / bf2, iterations = bf_iterations)
        bf_res = extractBF(bf, logbf = FALSE, onlybf = FALSE)
        bf_group = as.numeric(bf_res[1])
        bf_group_error = as.numeric(bf_res[2])
        
        #age BF
        bf1 = lmBF(formula = eval(parse(text=f)), data=df_sex)
        bf2 = lmBF(formula = eval(parse(text=gsub("\\+age", "",f))), data=df_sex)
        bf = recompute(bf1 / bf2, iterations = bf_iterations)
        bf_res = extractBF(bf, logbf = FALSE, onlybf = FALSE)
        bf_age = as.numeric(bf_res[1])
        bf_age_error = as.numeric(bf_res[2])
        
        #site BF
        bf1 = lmBF(formula = eval(parse(text=f)), data=df_sex)
        bf2 = lmBF(formula = eval(parse(text=gsub("\\+site", "",f))), data=df_sex)
        bf = recompute(bf1 / bf2, iterations = bf_iterations)
        bf_res = extractBF(bf, logbf = FALSE, onlybf = FALSE)
        bf_site = as.numeric(bf_res[1])
        bf_site_error = as.numeric(bf_res[2])
        
        #euler BF
        bf1 = lmBF(formula = eval(parse(text=f)), data=df_sex)
        bf2 = lmBF(formula = eval(parse(text=gsub("\\+TotalEulerNumber", "",f))), data=df_sex)
        bf = recompute(bf1 / bf2, iterations = bf_iterations)
        bf_res = extractBF(bf, logbf = FALSE, onlybf = FALSE)
        bf_euler = as.numeric(bf_res[1])
        bf_euler_error = as.numeric(bf_res[2])
        
        if (glob[gg] == "without_eICV"){
          #no global measure model
          glob_var_F = NA
          glob_var_pv = NA
          glob_var_bf = NA
          glob_var_bf_error = NA
        }       
        else{
          #global measure model
          glob_var_F = Anova(mm,type = "III")$"F value"[xvars=="eICV_samseg"]
          glob_var_pv = Anova(mm,type = "III")$"Pr(>F)"[xvars=="eICV_samseg"]
          
          bf1 = lmBF(formula = eval(parse(text=f)), data=df_sex)
          bf2 = lmBF(formula = eval(parse(text=gsub("\\+eICV_samseg", "",f))), data=df_sex)
          bf = recompute(bf1 / bf2, iterations = bf_iterations)
          bf_res = extractBF(bf, logbf = FALSE, onlybf = FALSE)
          bf_glob = as.numeric(bf_res[1])
          bf_glob_error = as.numeric(bf_res[2])
        } 
        
        
        #rows for sex separated ANOVA xlsx table
        rw_xlsx = list(model_vars[k], 
                       group_F, group_pv, bf_group, bf_group_error,
                       age_F, age_pv, bf_age, bf_age_error,
                       site_F, site_pv, bf_site, bf_site_error,
                       EulerNumber_F, EulerNumber_pv, bf_euler, bf_euler_error,
                       glob_var_F, glob_var_pv, bf_glob, bf_glob_error, "eICV_samseg",
                       gg-1,ss-1)
        
        DF_xlsx = rbindlist(list(DF_xlsx, rw_xlsx))
        
      } #fi both sex
      
    } #g 
  } #s
} #k

names(DF_xlsx)<-c("Model_yvar",
                  "Group_Fval","Group_pval","Group_BF","Group_BF_error",
                  "Age_Fval","Age_pval", "Age_BF","Age_BF_error",
                  "Site_Fval","Site_pval", "Site_BF","Site_BF_error",
                  "EulerNumber_Fval","Eulernumber_pval", "Eulernumber_BF","Eulernumber_BF_error",
                  "global_var_F","global_var_pv","global_var_BF","global_var_BF_error", "global_var_name",
                  "global_var_in_model","sex")

#flip the order of contrast names are defined if the sign is flipped
contrast_names = attributes(c)$levels$contrast

for (cn in seq(1,length(contrast_names))){
  if (contrast_signs[cn]<0){
    contrast_names[cn]=paste(strsplit(contrast_names[cn], split = " - ")[[1]][2],strsplit(contrast_names[cn], split = " - ")[[1]][1],sep = " - ")
  }    
}  
   
c_cols = expand.grid(c("contrast","t-ratio","pval","BF","BF_error","LCL","UCL"), contrast_names)
c_cols = paste(c_cols$Var1, c_cols$Var2,sep = "_")
c_cols = gsub(" ", "", c_cols)

group_lsmeans = paste("LSmean_",group_levels,sep = "")

#same excel saving procedure as above 
col_names = c("Model_yvar", group_lsmeans, c_cols,
              "ttest_DOF","global_var_in_model","sex")
names(DFc)<-col_names



DF_xlsx_glob0 = DFc[DFc$global_var_in_model == 0, ]
DF_xlsx_glob1 = DFc[DFc$global_var_in_model == 1, ]

write_xlsx(DF_xlsx_glob1,paste(save_folder,contrast_with_cov,sep = ""))
write_xlsx(DF_xlsx_glob0,paste(save_folder,contrast_without_cov,sep = ""))


DF_xlsx_glob0 = DF_xlsx[DF_xlsx$global_var_in_model == 0, ]
DF_xlsx_glob1 = DF_xlsx[DF_xlsx$global_var_in_model == 1, ]

write_xlsx(DF_xlsx_glob1,paste(save_folder,ANOVA_with_cov,sep = ""))
write_xlsx(DF_xlsx_glob0,paste(save_folder,ANOVA_without_cov,sep = ""))


### plot ppwp (contrast p-value plot)
# and the LSmeans 

for (i in seq(1,length(model_vars))){
  ps = list()
  emm = list()
  min_emm = list()
  max_emm = list()
  print(model_vars[i])
  
  #define what models LSmeans should be plotted for 
  if (model_vars[i] == "eICV_samseg"){
    glob = c("without_eICV")
    sexplot = c("female","male")
  }
  else if (model_vars[i] == "mean_thickness"){
    glob = c("without_eICV","with_eICV")
    sexplot = c("both")
  }
  else{
    glob = c("without_eICV","with_eICV")
    sexplot = c("female","male")
  }
  
  # loop to ppwp
  for (g in seq(1,length(glob))){
    for (s in seq(1,length(sexplot))){
      
      emm[[s]] = emmeans(models[[sexplot[s]]][[glob[g]]][[model_vars[i]]],specs = "group")
      plot_idx = s+length(glob)*g-length(glob)
      ps[[plot_idx]] = pwpp(emm[[s]],adjust="none",sort = FALSE) +
        labs(x="Uncorrected P-value") +
        ggtitle(paste(glob[g])) +
        geom_vline(xintercept = 0.05,linetype="dashed") +
        coord_flip()
      
      min_emm[[plot_idx]] = min(summary(emm[[s]])$emmean)
      max_emm[[plot_idx]] = max(summary(emm[[s]])$emmean)
    } #end s
  } #end g
  
  # plot LSmeans
  for (g in seq(1,length(glob))){
    for (s in seq(1,length(sexplot))){
      emm[[s]] = emmeans(models[[sexplot[s]]][[glob[g]]][[model_vars[i]]],specs = "group")
      plot_idx = s+length(glob)*g-length(glob)
      ps[[plot_idx +4]] = plot(emm[[s]]) +
        aes(color=group) +
        #facet_grid(cols =vars(sex)) +
        ggtitle(paste(sex[s])) +
        scale_x_continuous(limits = c(min(unlist(min_emm)), max(unlist(max_emm)))) +
        scale_y_discrete(limits=levels(datab$group)) +
        coord_flip()
    } #end g
  } #end s
  ps<-ps[!sapply(ps,is.null)]
  top_title = paste(model_vars[i]," models")
  ga=grid.arrange(grobs=ps, ncol=min(length(glob),length(sexplot))*2,top=textGrob(top_title,gp=gpar(fontsize=20)))
  
  #save the plot here
  ggsave(paste(plot_folder,model_vars[i],ppwp_postfix,".png",sep=""),ga,width = 10,height = 10)
} #end i








#### Effect size ####


#model_eff = lm(BrainTotalVol ~ group*sex + age + site + TotalEulerNumber, data=datab)
model_eff = lm(BrainTotalVol ~ group + age + site + TotalEulerNumber, data=dataf[[1]])
emm_eff = emmeans(model_eff,specs="group")

#cohens d
eff_size(emm_eff, sigma=sigma(model_eff), edf=df.residual(model_eff))

#etaSquared(model_eff, type = 3, anova = T)[xvars=="group",yvars="eta.sq.part"]



model_vars = c("BrainTotalVol", "CortexVol", "total_area", "mean_thickness")

glob = c("without_eICV","with_eICV")
sex = c("female","male","both") #3 sexes are defined, 0=female, 1=male, 2=both


DF = data.frame()

#all non-changing variable names in models
all_var_names = c("group","sex",general_cov, "eICV_samseg", "group:sex")


for (k in seq(1,length(model_vars))){
  for (g in seq(1,length(glob))){
    for (s in seq(1,length(sex))){
      
      if (sex[s] == "both"){
        model_eff = model[[glob[g]]][[model_vars[k]]]
        GS_in_model = ("group:sex" %in% attributes(Anova(model_eff,type = "III"))$row.names)*1
        model_type = "group + sex"
        if (GS_in_model == 1){
          model_type = "group*sex"
        }
      }
      else {
        model_type = NA
        model_eff = models[[sex[s]]][[glob[g]]][[model_vars[k]]]
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
      DF = rbindlist(list(DF, as.list(c(model_vars[k],row,
                                        c_eff_vals,
                                        g-1,GS_in_model,sex[s],model_type))))
      
      
    } #end s
  } #end g
} #end i

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


names(DF) = c("model_yvar",paES_col_names,
              c_cols,
              "globvar_in_model","GS_in_model","sex","model_type")

write_xlsx(DF,effect_sizes_path)




###############################################################################

#####                           Testing area                              ##### 


#             EVERYTHING BELOW IS ONLY FOR VALIDATION OF METHODS              #

###############################################################################


model_bvol_glob = lm(BrainTotalVol ~ group*sex + age + TotalEulerNumber + site, data=datab) #eICV_samseg
Anova(model_bvol_glob,type="III")
lsmeans(model_bvol_glob,pairwise ~ group)

model_bvol_glob2 = lm(BrainTotalVol ~ group + age + TotalEulerNumber + site, data=data_sex0) #eICV_samseg
Anova(model_bvol_glob2,type="III")
attributes(lsmeans(model_bvol_glob2,pairwise ~ group)$contrasts)$dfargs$df


model_bvol_glob2 = lm(BrainTotalVol ~ group + age + TotalEulerNumber + site, data=data_sex1) #eICV_samseg
Anova(model_bvol_glob2,type="III")

#model_bvol_glob2 = lm(BrainTotalVol ~ group + sex + age + TotalEulerNumber + site, data=datab) #eICV_samseg
#anova(model_bvol_glob, model_bvol_glob2)
#model_bvol_glob = update(model_bvol_glob,~.-group:sex)
#etaSquared(model_bvol_glob, type = 3, anova = T)

emm_eff = emmeans(model_bvol_glob,specs="group")

#cohens d
eff_size(emm_eff, sigma=sigma(model_bvol_glob), edf=df.residual(model_bvol_glob))



#### Bayesian factor analysis 

###find ud af om bayes factor er log - det er ikke log 
#priors ? modified separate g-priors on categorical, regular g-priors with all g's same on continues 
#what is the prior odds?
#få styr på proportional error 
#få styr på multiple comparisons tests - side 19 i bayestestR
#change prior scales???


## Each effect in anova model 

#Group/sex
bf1 = lmBF(BrainTotalVol ~ group + sex + group:sex + age + TotalEulerNumber + site, data=datab)
bf2 = lmBF(BrainTotalVol ~ group + sex +             age + TotalEulerNumber + site, data=datab)
bf_interaction = bf1 / bf2
newbf = recompute(bf_interaction, iterations = 10000)
newbf
bf = extractBF(newbf, logbf = FALSE, onlybf = FALSE)
#8.075819 ±0.53%
as.numeric(bf[1])
as.numeric(bf[2])

#age
bf1 = lmBF(mean_thickness ~ group + age + TotalEulerNumber + site, data=data_sex0)
bf2 = lmBF(mean_thickness ~ group + TotalEulerNumber + site, data=data_sex0)
bf1 / bf2

newbf = recompute(bf1 / bf2, iterations = 500000)
newbf

classic = lm(mean_thickness ~ group + age + TotalEulerNumber + site, data=data_sex0)
Anova(classic, type="III")


#Euler number
bf1 = lmBF(BrainTotalVol ~ group + age + TotalEulerNumber + site, data=data_sex0)
bf2 = lmBF(BrainTotalVol ~ group + age + site, data=data_sex0)
bf1 / bf2

#Site
bf1 = lmBF(CortexVol ~ group + age + TotalEulerNumber + site, data=data_sex0)
bf2 = lmBF(CortexVol ~ group + age + TotalEulerNumber, data=data_sex0)
bf1 / bf2
newbf = recompute(bf1 / bf2, iterations = 50000)
newbf

classic = lm(CortexVol ~ group + age + TotalEulerNumber + site, data=data_sex0)
Anova(classic, type="III")

#Sex #THIS ANALYSIS IS INVALID
bf1 = lmBF(BrainTotalVol ~ group + sex + age + TotalEulerNumber + site, data=datab)
bf2 = lmBF(BrainTotalVol ~ group       + age + TotalEulerNumber + site, data=datab)
bf1 / bf2


#sex 1
bf1 = lmBF(BrainTotalVol ~ group + age + TotalEulerNumber + site, data=data_sex1)
bf2 = lmBF(BrainTotalVol ~         age + TotalEulerNumber + site, data=data_sex1)
bf_group = bf1 / bf2
bf_group

#sex 0
bf1 = lmBF(BrainTotalVol ~ group + age + TotalEulerNumber + site, data=data_sex0)
bf2 = lmBF(BrainTotalVol ~         age + TotalEulerNumber + site, data=data_sex0)
bf_group = bf1 / bf2
bf_group



bf1 = lmBF(CortexVol ~ group + age + TotalEulerNumber + site, data=data_sex1)
bf2 = lmBF(CortexVol ~         age + TotalEulerNumber + site, data=data_sex1)
bf1/bf2

#Compute group contrasts for a certain model

#get bayes factor for the classical test
data_BP_K_sex0 = data_sex0[!data_sex0$group=="SZ",]
bf1 = lmBF(CortexVol ~ group + age + TotalEulerNumber + site, data=data_BP_K_sex0)
bf2 = lmBF(CortexVol ~         age + TotalEulerNumber + site, data=data_BP_K_sex0)
bf1/bf2
recompute(bf1/bf2,iterations=1000000)


#this is the equivalent classical tests
bf1 = lm(CortexVol ~ group + age + TotalEulerNumber + site, data=data_BP_K_sex0)
bf2 = lm(CortexVol ~         age + TotalEulerNumber + site, data=data_BP_K_sex0)
anova(bf2,bf1)
lsmeans(bf1,pairwise ~ group,adjust="none")


data_SZ_K_sex1 = data_sex1[!data_sex1$group=="BP",]
bf1 = lmBF(CortexVol ~ group + age + TotalEulerNumber + site, data=data_SZ_K_sex1)
bf2 = lmBF(CortexVol ~         age + TotalEulerNumber + site, data=data_SZ_K_sex1)
bf1/bf2


#this is the equivalent classical tests
bf1 = lm(CortexVol ~ group + age + TotalEulerNumber + site, data=data_SZ_K_sex1)
bf2 = lm(CortexVol ~         age + TotalEulerNumber + site, data=data_SZ_K_sex1)
bf3 = lm(CortexVol ~ 1, data=data_SZ_K_sex1)
anova(bf2,bf1)
anova(bf3,bf1)
lsmeans(bf1,pairwise ~ group,adjust="none")




# Following tutorial from: http://bayesfactor.blogspot.com/2015/01/multiple-comparisons-with-bayesfactor-2.html
#data_BP_K = datab[!datab$group=="SZ",]
#data_BP_K_sex0 = data_sex0[!data_sex0$group=="SZ",]
#data_SZ_K_sex0 = data_sex0[!data_sex0$group=="BP",]

#bf1 = lm(BrainTotalVol ~ group + age + TotalEulerNumber + site, data=data_BP_K)
#lsmeans(bf1,pairwise ~ group, adjust="none")

#bf1 = lmBF(formula=BrainTotalVol ~ group + age + TotalEulerNumber + site, data=data_BP_K_sex0)
#samples = posterior(bf1, iterations = 10000)
##mean(samples[, "group-BP"]) - mean(samples[, "group-K"])
#consistent = ( samples[, "group-BP"] > samples[, "group-K"] )
#posterior_prob = (sum(consistent) / 10000)
#bf_restriction_against_full = posterior_prob / (1 / 2) # bf comparing the restriction against full model (in which all 3 means are unequal)
#bf_restriction_against_full

#the “full” model (that all three means are unequal, group1 != group2 != group3)
#the “null” model (all means are equal, group1=group2=group3)

#bf_full_against_null = as.vector(bf1)
#bf_restriction_against_null = bf_restriction_against_full * bf_full_against_null
#bf_restriction_against_null


# SZ - K contrast
#bf1 = lm(BrainTotalVol ~ group + age + TotalEulerNumber + site, data=data_SZ_K_sex0)
#lsmeans(bf1,pairwise ~ group, adjust="none")

#bf1 = lmBF(formula=BrainTotalVol ~ group + age + TotalEulerNumber + site, data=data_SZ_K_sex0)
#samples = posterior(bf1, iterations = 10000)
#consistent = ( samples[, "group-SZ"] > samples[, "group-K"] )
#posterior_prob = (sum(consistent) / 10000)
#bf_restriction_against_full = posterior_prob / (1 / 2) # bf comparing the restriction against full model (in which all 3 means are unequal)
#bf_restriction_against_full



## Use residuals and ttestBF 
#data_BP_sex0 = data_sex0[data_sex0$group=="BP",]
#BP_model = lm(BrainTotalVol ~ 1 + age + TotalEulerNumber + site, data=data_BP_sex0)
#res_BP = BP_model$residuals
#BP_level = res_BP + as.numeric(BP_model$coefficients)[1]

#data_K_sex0 = data_sex0[data_sex0$group=="K",]
#K_model = lm(BrainTotalVol ~  1 + age + TotalEulerNumber + site, data=data_K_sex0)
#res_K = K_model$residuals
#K_level = res_K + as.numeric(K_model$coefficients)[1]

#ttestBF(x=K_level, y=BP_level)
#t.test(x=K_level, y=BP_level, var.eq=TRUE)



## BayestestR emm BF
#library(rstanarm)
#library(bayestestR)
#stan_model <- stan_glm(BrainTotalVol ~ group + age + TotalEulerNumber + site, data=datab,
#                        prior = normal(0,1),                          # as per Rouder et al., 2012
#                        prior_intercept = student_t(3,0,10),          # weakly informative
#                        prior_aux = exponential(.1),
#                        chains = 4, iter = 10000)                 # weakly informative
#c_color_main <- pairs(emmeans(stan_model, ~ group))
#describe_posterior(c_color_main,
#                   estimate = "median", dispersion = TRUE,
#                   ci = .9, ci_method = "hdi",
#                   test = c("bayesfactor"),
#                   bf_prior = stan_model)



## Bayesian LSmeans
#source("https://gist.githubusercontent.com/mattansb/f383af8c0dfe88922fb5c75e8572d03e/raw/9616f4a2c959716a1e82bb4a20fd60a8fe8d0f3e/BFBayesFactor_2_emmeans.R")
#options(contrasts=c('contr.sum', 'contr.poly'))
#library(magrittr)
#library(purrr)
#fit_BF <- anovaBF(BrainTotalVol ~ group + site, data = datab)
#emmeans(fit_BF[4],pairwise~group)


## ttestBF with predictions, however WITHOUT lsmeans but regular means
#data_BP_K = datab[!datab$group=="SZ",]
#predict(model_bvol_glob, newdata = data.frame(group = c("K","BP","SZ"), TotalEulerNumber = mean(datab$TotalEulerNumber), age = mean(datab$age) ))
#model_bvol_glob = lm(BrainTotalVol ~ group + age + TotalEulerNumber + site, data=data_sex0) #eICV_samseg
#predict(model_bvol_glob )
#BP_vals = predict(model_bvol_glob )[data_sexdata.frame(group="BP",age=0,TotalEulerNumber=0,site=factor(1))]
#SZ_vals = predict(model_bvol_glob )[data_sex0$group=="SZ"]
#K_vals = predict(model_bvol_glob )[data_sex0$group=="K"]
#mean(BP_vals)
#mean(SZ_vals)
#mean(K_vals)
#BP_raw_mean_c = mean(BP_vals) - mean(K_vals)
#SZ_raw_mean_c = mean(SZ_vals) - mean(K_vals)
#ls=lsmeans(model_bvol_glob,pairwise ~ group,adjust="none")
#ls
#diff_percent = (BP_raw_mean_c  - summary(ls)$contrasts$estimate[1])*100 / summary(ls)$contrasts$estimate[1]
#diff_percent
#ttestBF(x=as.numeric(BP_vals), y=as.numeric(K_vals), iterations=1000)


## estimate the standard errors on the lsmean
#sbp = sd(data_sex0$BrainTotalVol[data_sex0$group=="BP"]) #/ sqrt(length(data_sex0$group))
#ssz = sd(data_sex0$BrainTotalVol[data_sex0$group=="SZ"])
#sk = sd(data_sex0$BrainTotalVol[data_sex0$group=="K"])
#s = c(sbp, ssz, sk)
#n = c(sum(data_sex0$group=="BP"), sum(data_sex0$group=="SZ"), sum(data_sex0$group=="K"))
#s.p <- sqrt( sum((n - 1) * s^2) / sum(n - 1) ) 
#s.p / sqrt(n)


## post hoc t-test with t-statistic
#model_bvol_glob = lm(BrainTotalVol ~ group + sex + age + TotalEulerNumber + site, data=datab) #eICV_samseg
#model_bvol_glob = lm(BrainTotalVol ~ group + age + TotalEulerNumber + site, data=data_sex0) #eICV_samseg
#Anova(model_bvol_glob,type="III")
#ls=lsmeans(model_bvol_glob,pairwise ~ group,adjust="none")
# the standard error of the contrast test
#sqrt( sigma(model_bvol_glob)**2 * (1/sum(data_sex0$group=="BP") + 1/sum(data_sex0$group=="K")) )
# how lsmean does the test
#dof = length(data_sex0$group) - 3 - 1 - 1 - 1
#2*pt(q=2.555192, df=dof, lower.tail=FALSE)
#t_ratios = summary(ls)$contrasts$"t.ratio"
#dof = model_bvol_glob$df.residual
#looks like the degrees of freedom is n1+n2-2
# or its n1-1
#result <- ttest.tstat(t = t_ratios[1], n1 = dof+1)#n1 = sum(datab$group=="BP"), n2 = sum(datab$group=="K"))
#exp(result[['bf']])
#result <- ttest.tstat(t = t_ratios[1], n1 = sum(datab$group=="BP"), n2 = sum(datab$group=="K"))
#exp(result[['bf']])
#result[['properror']]
#result <- ttest.tstat(t = t_ratios[2], n1 = dof+1) #n1 = sum(datab$group=="BP"), n2 = sum(datab$group=="SZ"))
#exp(result[['bf']])
#result[['properror']]
#result <- ttest.tstat(t = t_ratios[3], n1 = dof+1) #n1 = sum(datab$group=="SZ"), n2 = sum(datab$group=="K"))
#exp(result[['bf']])
#result[['properror']]


## Something something 
#library(emmeans)
#model1 = lm(BrainTotalVol ~ group + age + TotalEulerNumber + site, data=data_sex0)
#emm = emmeans(model1,specs="group")
#equivalence_test(contrast(emm))
#library(brms)
#model <- brms::brm(BrainTotalVol ~ group + age + TotalEulerNumber + site, data=data_sex0)
#equivalence_test(model)
