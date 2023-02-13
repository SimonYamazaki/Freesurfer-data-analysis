

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
library(rlist)

#This script is intended for computation of statistics on behavioural measures
#based on GLMs on data from all groups (to potentially find an overall group effect)

#models on the form below are run:
#Global brain measure ~ behavioral measure + group*sex + age + site + eulernumber 

# models with one of the following behavioural measures as covariate is produced: 
#("CBCL_ext_cg_v11","CBCL_int_cg_v11","CBCL_totsc_cg_v11","CGASx_v11")
#with the following global brain measures as response variables: 
#("BrainTotalVol", "CortexVol", "total_area", "mean_thickness", "eICV_samseg")

#the following covariates are included in all models by default
#("age","site","TotalEulerNumber")


#returns: 
#relevant statistics from GLMs run for each combination of response variable and behavioural measures and saved into an excel sheets 
#files generated are: 

#GS models (with group-sex interaction)
# - excel sheet with variable level GLM effects with global covariate (eICV_samseg)
# - excel sheet with variable level GLM effects without global covariate (eICV_samseg)

#S models
# - excel sheet with variable level GLM effects run separately on data from each sex with global covariate (eICV_samseg)
# - excel sheet with variable level GLM effects run separately on data from each sex without global covariate (eICV_samseg)

# - excel sheet with group contrasts for GLMs (both GS and S models from above) without global covariate (eICV_samseg)
# - excel sheet with group contrasts for GLMs (both GS and S models from above) with global covariate (eICV_samseg)
# - excel sheet with variable level GLM effect sizes

#For file names and directories saved into, refer to the variables save_folder, GS_B_model_ANOVA_with_glob, GS_B_model_ANOVA_without_glob, B_model_ANOVA_with_glob, B_model_ANOVA_without_glob, B_model_contrast_with_glob, B_model_contrast_without_glob, B_effect_sizes_path
#data is loaded from the path specified in data_path

#includes sheets for whole group (with G/S interaction)



#save folder for tables and plots:
save_folder = "/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/global_variables/behav/B/"
#save_folder_plot = "/mnt/projects/VIA11/FREESURFER/Stats/Plots/global_measures/behav/B/"

#create the folders if they dont exist
dir.create(file.path(dirname(save_folder), basename(save_folder)))
#dir.create(file.path(dirname(save_folder_plot), basename(save_folder_plot)))


#File names
GS_B_model_ANOVA_with_glob = "globvar_GS_ANOVA_pvals_with_glob.xlsx"
GS_B_model_ANOVA_without_glob = "globvar_GS_ANOVA_pvals_without_glob.xlsx"

#B model sheets for each sex, a row for each sex
B_model_ANOVA_with_glob = "globvar_S_ANOVA_pvals_with_glob.xlsx"
B_model_ANOVA_without_glob = "globvar_S_ANOVA_pvals_without_glob.xlsx"

#contrast of the GS and S models above
B_model_contrast_with_glob = "globvar_S_Model_contrasts_with_glob.xlsx"
B_model_contrast_without_glob = "globvar_S_Model_contrasts_without_glob.xlsx"

B_effect_sizes_path = "ANOVA+contrast_effect_sizes.xlsx"


bf_iterations = 10000


#data path
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

# - extract rows with 1 in Include_FS_studies_euler_outliers_sibpairs_out
data_csv_filtered <- data_csv_filtered[c(data_csv_filtered$Include_FS_studies_euler_outliers_sibpairs_out == 1),]
data_csv_filtered <- data_csv_filtered[!is.na(data_csv_filtered$Include_FS_studies_euler_outliers_sibpairs_out),]


#make new variables with shorter and contained names
# - tell r which variables are factors
datab = data_csv_filtered
datab$sex = as.factor(datab$Sex_child)
datab$site = as.factor(datab$MRI_site_v11)
datab$diag = as.factor(datab$ksads_any_diag_excl_elim_lft_v11)
datab$age = as.numeric(datab$MRI_age_v11)

#what groups should be used
datab$group = as.factor(datab$HighRiskStatus_v11)


#add a string encoding of sex for later convenience 
datab$sexs = datab$Sex_child
datab$sexs[datab$sexs == 0] = "female"
datab$sexs[datab$sexs == 1] = "male"




###### for loop that generates tables of global brain measure as y-variable 
#with each of the behavioural variables as covariates


#specify dependent variable (global variable) column names, i.e. left hand side of statistical models
model_vars = c("BrainTotalVol", "CortexVol", "total_area", "mean_thickness", "eICV_samseg")
behav_vars = c("CBCL_ext_cg_v11","CBCL_int_cg_v11","CBCL_totsc_cg_v11","CGASx_v11")

#covariates that are included in all models 
#encoded as a formula string
general_cov = c("age","site","TotalEulerNumber")
general_cov_formula = paste(general_cov,collapse = " + ")


contrast_signs = c(1,1, -1) # for K/BP/SZ 



´



## B models whole group


model = list()

DF = data.frame()

#loop that defines models 
for (i in seq(1,length(model_vars))){
  if (model_vars[i] == "eICV_samseg"){    #only run models without_eICV for the model on eICV_samseg
    glob = c("without_eICV")
  }
  else{
    glob = c("without_eICV","with_eICV") #run both for all other variables
  }
  print(model_vars[i])
  
  for (g in seq(1,length(glob))){
    
    #model formulas WITHOUT behav
    if (glob[g] == "with_eICV"){
      #f = paste(model_vars[i],"~","+","group*sex","+",general_cov_formula,"+","eICV_samseg") #model with ICV
      f = paste(model_vars[i],"~","+group+sex+group:sex","+",general_cov_formula,"+","eICV_samseg",sep="") #model with ICV
      f_without_gs = paste(model_vars[i],"~","+","group+sex","+",general_cov_formula,"+","eICV_samseg",sep="") #model with ICV
    }
    else{
      #f = paste(model_vars[i],"~","+","group*sex","+",general_cov_formula) #model without ICV
      f = paste(model_vars[i],"~","+group+sex+group:sex","+",general_cov_formula,sep="") #model without ICV
      f_without_gs = paste(model_vars[i],"~","+","group+sex","+",general_cov_formula,sep="") #model with ICV
    }
    
    #loop that adds and analyze models WITH behav
    for (b in seq(1,length(behav_vars))){
      print(behav_vars[b])
      
      beh = behav_vars[b]
      ff = paste(f,"+",beh,sep="")
      
      #run the actual model
      model_behav = lm(ff,data=datab)
      
      xvars = attributes(Anova(model_behav,type = "III"))$row.names
      
      GS_F = Anova(model_behav,type = "III")$"F value"[xvars=="group:sex"]
      GS_pv = Anova(model_behav,type = "III")$"Pr(>F)"[xvars=="group:sex"]
      
      
      #if the group/sex interaction is insignificant
      #update the current model_behav if insignificant, and save the old one
      if (Anova(model_behav,type = "III")$"Pr(>F)"[xvars=="group:sex"] > 0.05){
        #indicate the significance of the group/sex interaction
        sign_GS = 0
        
        model_gs = model_behav
        
        model_behav = update(model_behav,~.-group:sex)
        
        models = list(model_gs, model_behav)
        model_formula = list(f,f_without_gs)
        model_type = list("GS","No_interactions")
        
      }
      else{    #if the group/sex interaction is significant - just keep the model
        #indicate the significance of the group/sex interaction
        sign_GS = 1
        models = list(model_behav)
        model_formula = list(f)
        model_type = list("GS")
      }
      
      
      xvars = attributes(Anova(model_behav,type = "III"))$row.names
      

      # loop that makes rows in the ANOVA table excel 
      # makes two rows if the group/sex interaction is insignificant 
      # one row with the interaction included and one row without
      # makes 3 rows if the group/behav is insignificant 
      
      for (m in seq(1,length(models))){
        mm = models[[m]]
        mf = model_formula[[m]]
        
        xvars = attributes(Anova(mm,type = "III"))$row.names
        
        #extract relevant statistics for each variable in the model 
        behav_F = Anova(mm,type = "III")$"F value"[xvars==beh]
        group_F = Anova(mm,type = "III")$"F value"[xvars=="group"]
        sex_F = Anova(mm,type = "III")$"F value"[xvars=="sex"]
        age_F = Anova(mm,type = "III")$"F value"[xvars=="age"]
        site_F = Anova(mm,type = "III")$"F value"[xvars=="site"]
        EulerNumber_F = Anova(mm,type = "III")$"F value"[xvars=="TotalEulerNumber"]
        
        behav_pv = Anova(mm,type = "III")$"Pr(>F)"[xvars==beh]
        group_pv = Anova(mm,type = "III")$"Pr(>F)"[xvars=="group"]
        sex_pv = Anova(mm,type = "III")$"Pr(>F)"[xvars=="sex"]
        age_pv = Anova(mm,type = "III")$"Pr(>F)"[xvars=="age"]
        site_pv = Anova(mm,type = "III")$"Pr(>F)"[xvars=="site"]
        EulerNumber_pv = Anova(mm,type = "III")$"Pr(>F)"[xvars=="TotalEulerNumber"]
        
      
        ############### Bayes Factor analysis 
        
        #behav
        bf1 = lmBF(formula = eval(parse(text=mf)), data=datab)
        bf2 = lmBF(formula = eval(parse(text=gsub(paste("\\+",beh,sep=""), "",mf))), data=datab)
        bf_interaction = recompute(bf1 / bf2, iterations = bf_iterations)
        bf_res = extractBF(bf_interaction, logbf = FALSE, onlybf = FALSE)
        bf_behav = as.numeric(bf_res[1])
        bf_behav_error = as.numeric(bf_res[2])
        
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
        
        subs = length(model_behav$residuals)
        
        #the first index in m is always the model with a group/sex and group/behav interaction
        if (model_type[[m]]=="GS"){
          GS_in = 1
          
          #Group/sex BF
          bf1 = lmBF(formula = eval(parse(text=f)), data=datab)
          bf2 = lmBF(formula = eval(parse(text=f_without_gs)), data=datab)
          bf_interaction = recompute(bf1 / bf2, iterations = bf_iterations)
          bf_res = extractBF(bf_interaction, logbf = FALSE, onlybf = FALSE)
          bf_gs = as.numeric(bf_res[1])
          bf_gs_error = as.numeric(bf_res[2])
          
        }
        else if (model_type[[m]]=="No_interactions"){ #if m is 3 then the model is does not include group/sex and group/behav interaction
          GS_F = NA
          GS_pv = NA
          GS_in = 0
          
          bf_gs = NA
          bf_gs_error = NA
        }
        
        model[[glob[g]]][[behav_vars[b]]][[model_vars[i]]] = model_behav
        
        #append the row to the dataframe 
        #DF = rbindlist(list(DF, rw))
        DF = rbindlist(list(DF, as.list(c(model_vars[i], beh,
                                          group_F, group_pv, bf_g, bf_g_error,
                                          sex_F, sex_pv, bf_sex, bf_sex_error,
                                          age_F, age_pv, bf_age, bf_age_error,
                                          site_F, site_pv, bf_site, bf_site_error,
                                          EulerNumber_F, EulerNumber_pv, bf_euler, bf_euler_error,
                                          ICV_F, ICV_pv, bf_glob, bf_glob_error,
                                          behav_F, behav_pv, bf_behav, bf_behav_error,
                                          
                                          GS_F, GS_pv, bf_gs, bf_gs_error,
                                          sign_GS,GS_in, 
                                          
                                          g-1,model_type[[m]],subs))))
        
      } # for m 
      
    } # b
  } #g
} #end i


#define the coloum names of the dataframe which is converted to an excel 
col_names = c("Model_yvar", "behav_var",
              "Group_Fval","Group_pval", "Group_bf", "Group_bf_error",
              "Sex_Fval","Sex_pval", "Sex_bf", "Sex_bf_error",
              "Age_Fval","Age_pval", "Age_bf", "Age_bf_error",
              "Site_Fval","Site_pval", "Site_bf", "Site_bf_error",
              "EulerNumber_Fval","Eulernumber_pval", "Eulernumber_bf", "Eulernumber_bf_error",
              "ICV_Fval","ICV_pval", "ICV_bf", "ICV_bf_error",
              "behav_Fval","behav_pval", "behav_bf", "behav_bf_error",
              "Group_sex_Fval","Group_sex_pval", "Group_sex_bf", "Group_sex_bf_error",
              "Significant_GS_interaction","GS_in_model",
              "ICV_in_model","model_type","n_subjects")

names(DF)<-col_names


DF_xlsx_glob0 = DF[DF$ICV_in_model == 0, ]
DF_xlsx_glob1 = DF[DF$ICV_in_model == 1, ]

write_xlsx(DF_xlsx_glob1,paste(save_folder,GS_B_model_ANOVA_with_glob,sep=""))
write_xlsx(DF_xlsx_glob0,paste(save_folder,GS_B_model_ANOVA_without_glob,sep=""))




#B models for each sex

#### CONTRAST XLSX with male / female split
# and built models[sex][glob][yvar]

#define some relevant naming 
glob = c("without_eICV","with_eICV")
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
  
  
  #loop that adds and analyze models WITH behav
  for (b in seq(1,length(behav_vars))){
    print(behav_vars[b])
    
    beh = behav_vars[b]
    
    for (ss in seq(1,3)){ #use the 3 defined sexes to know what data to use
      
      # define an temporary sex indicator to index
      print(sex[ss])
      
      #extract the dataframe containing the relevant sex data
      df_sex = dataf[[ss]]
      
      
      for (gg in seq(1,2)){
        f = paste(model_vars[k],"~","+group","+",beh,"+",general_cov_formula,sep="")
        
        if (glob[gg] == "with_eICV"){ #run models without ICV
          f = paste(f,"+eICV_samseg",sep="")
        }
        
        if (sex[ss] == "both"){ #define a specific model for the case of both sex
          f = paste(f,"+sex",sep="")
        }
        
        #run the model
        models[[sex[ss]]][[glob[gg]]][[behav_vars[b]]][[model_vars[k]]] = lm(f,data=df_sex)
        
        #save the current model into a temporary variable to make code less cluttered
        model_ana = models[[sex[ss]]][[glob[gg]]][[behav_vars[b]]][[model_vars[k]]]
        xvars = attributes(Anova(model_ana,type = "III"))$row.names
        
        
        
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
        
        
        
        #ANOVA sex separated
        mm=model_ana
          
        behav_F = Anova(mm,type = "III")$"F value"[xvars==beh]
        group_F = Anova(mm,type="III")$"F value"[xvars=="group"]
        age_F = Anova(mm,type = "III")$"F value"[xvars=="age"]
        site_F = Anova(mm,type = "III")$"F value"[xvars=="site"]
        EulerNumber_F = Anova(mm,type = "III")$"F value"[xvars=="TotalEulerNumber"]
        
        behav_pv = Anova(mm,type = "III")$"Pr(>F)"[xvars==beh]
        group_pv = Anova(mm,type = "III")$"Pr(>F)"[xvars=="group"]
        age_pv = Anova(mm,type = "III")$"Pr(>F)"[xvars=="age"]
        site_pv = Anova(mm,type = "III")$"Pr(>F)"[xvars=="site"]
        EulerNumber_pv = Anova(mm,type = "III")$"Pr(>F)"[xvars=="TotalEulerNumber"]
        
        ### Bayes Factor analysis 
        
        #behav BF
        bf1 = lmBF(formula = eval(parse(text=f)), data=df_sex)
        bf2 = lmBF(formula = eval(parse(text=gsub(paste("\\+",beh,sep=""), "",f))), data=df_sex)
        bf = recompute(bf1 / bf2, iterations = bf_iterations)
        bf_res = extractBF(bf, logbf = FALSE, onlybf = FALSE)
        bf_behav = as.numeric(bf_res[1])
        bf_behav_error = as.numeric(bf_res[2])
        
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
          bf_glob = NA
          bf_glob_error = NA
        }       
        else{
          #global measure model
          glob_var_F = Anova(mm,type = "III")$"F value"[xvars=="eICV_samseg"]
          glob_var_pv = Anova(mm,type = "III")$"Pr(>F)"[xvars=="eICV_samseg"]
          
          #glob BF
          bf1 = lmBF(formula = eval(parse(text=f)), data=df_sex)
          bf2 = lmBF(formula = eval(parse(text=gsub("\\+eICV_samseg", "",f))), data=df_sex)
          bf = recompute(bf1 / bf2, iterations = bf_iterations)
          bf_res = extractBF(bf, logbf = FALSE, onlybf = FALSE)
          bf_glob = as.numeric(bf_res[1])
          bf_glob_error = as.numeric(bf_res[2])
        } 
        
        if (sex[ss] == "both"){
          sex_F = Anova(mm,type = "III")$"F value"[xvars=="sex"]
          sex_pv = Anova(mm,type = "III")$"Pr(>F)"[xvars=="sex"]
          model_type = "group + sex"
          
          #sex BF
          bf1 = lmBF(formula = eval(parse(text=f)), data=df_sex)
          bf2 = lmBF(formula = eval(parse(text=gsub("\\+sex", "",f))), data=df_sex)
          bf = recompute(bf1 / bf2, iterations = bf_iterations)
          bf_res = extractBF(bf, logbf = FALSE, onlybf = FALSE)
          bf_sex = as.numeric(bf_res[1])
          bf_sex_error = as.numeric(bf_res[2])
        }
        else {
          sex_F = NA
          sex_pv = NA
          bf_sex = NA
          bf_sex_error = NA
          model_type = "group"
        }
        
        #contrast table
        rwc_xlsx = append(rwc_xlsx,list(gg-1,ss-1,model_type))
        DFc = rbindlist(list(DFc, rwc_xlsx))
        
        #rows for ANOVA xlsx table
        rw_xlsx = list(model_vars[k], beh, 
                       group_F, group_pv, bf_group, bf_group_error,
                       sex_F, sex_pv, bf_sex, bf_sex_error,
                       age_F, age_pv, bf_age, bf_age_error,
                       site_F, site_pv, bf_site, bf_site_error,
                       EulerNumber_F, EulerNumber_pv, bf_euler, bf_euler_error,
                       behav_F, behav_pv, bf_behav, bf_behav_error,
                       glob_var_F, glob_var_pv, bf_glob, bf_glob_error,
                       "eICV_samseg",
                       gg-1,ss-1,model_type)
        
        DF_xlsx = rbindlist(list(DF_xlsx, rw_xlsx))
        
      } # glob 
    } #sex 
  } #behav
} #model_vars

names(DF_xlsx)<-c("Model_yvar","Behav_var",
                  "Group_Fval","Group_pval","Group_bf", "Group_bf_error",
                  "Sex_Fval","Sex_pval", "Sex_bf", "Sex_bf_error",
                  "Age_Fval","Age_pval","Age_bf", "Age_bf_error",
                  "Site_Fval","Site_pval","Site_bf", "Site_bf_error",
                  "EulerNumber_Fval","Eulernumber_pval","Eulernumber_bf", "Eulernumber_bf_error",
                  "Behav_Fval","Behav_pval","Behav_bf", "Behav_bf_error",
                  "global_var_F","global_var_pv","global_var_bf", "global_var_bf_error",
                  "global_var_name",
                  "global_var_in_model","sex","model_type")


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
              "global_var_in_model","sex","Model_type")
names(DFc)<-col_names



DF_xlsx_glob0 = DFc[DFc$global_var_in_model == 0, ]
DF_xlsx_glob1 = DFc[DFc$global_var_in_model == 1, ]

write_xlsx(DF_xlsx_glob1,paste(save_folder,B_model_contrast_with_glob,sep = ""))
write_xlsx(DF_xlsx_glob0,paste(save_folder,B_model_contrast_without_glob,sep = ""))


DF_xlsx_glob0 = DF_xlsx[DF_xlsx$global_var_in_model == 0, ]
DF_xlsx_glob1 = DF_xlsx[DF_xlsx$global_var_in_model == 1, ]

write_xlsx(DF_xlsx_glob1,paste(save_folder,B_model_ANOVA_with_glob,sep = ""))
write_xlsx(DF_xlsx_glob0,paste(save_folder,B_model_ANOVA_without_glob,sep = ""))





#### Effect size ####


library(lsr)

#model_eff = lm(BrainTotalVol ~ group*sex + age + site + TotalEulerNumber, data=datab)
#model_eff = lm(BrainTotalVol ~ group + age + site + TotalEulerNumber, data=dataf[[1]])
#emm_eff = emmeans(model_eff,specs="group")

#cohens d
#eff_size(emm_eff, sigma=sigma(model_eff), edf=df.residual(model_eff))

#etaSquared(model_eff, type = 3, anova = T)[xvars=="group",yvars="eta.sq.part"]

#model_vars = c("BrainTotalVol", "CortexVol", "total_area", "mean_thickness")

#glob = c("without_eICV","with_eICV")
sex = c("female","male","both") #3 sexes are defined, 0=female, 1=male, 2=both


DF = data.frame()


#all non-changing variable names in models
all_var_names = c("group","sex",general_cov, "eICV_samseg", "group:sex")


for (k in seq(1,length(model_vars))){
  if (model_vars[k] == "eICV_samseg"){    #only run models without_eICV for the model on eICV_samseg
    glob = c("without_eICV")
  }
  else{
    glob = c("without_eICV","with_eICV") #run both for all other variables
  }
  for (b in seq(1,length(behav_vars))){
    for (g in seq(1,length(glob))){
      for (s in seq(1,length(sex))){
        
        all_col_names = c(all_var_names, behav_vars[b])
        
        if (sex[s] == "both"){
          model_eff = model[[glob[g]]][[behav_vars[b]]][[model_vars[k]]]
          GS_in_model = ("group:sex" %in% attributes(Anova(model_eff,type = "III"))$row.names)*1
          
          model_type = "group + sex"
          if (GS_in_model == 1){
            model_type = "group*sex"
          }
        }
        else {
          model_type = NA
          GS_in_model = 0 #no GS interaction is possible for sex separated models
          model_eff = models[[sex[s]]][[glob[g]]][[behav_vars[b]]][[model_vars[k]]]
        }
        
        
        #get model terms 
        xvars = attributes(etaSquared(model_eff, type = 3, anova = T))$dimnames[[1]]
        
        #get etaSquared outputs
        yvars = attributes(etaSquared(model_eff, type = 3, anova = T))$dimnames[[2]]
        
        eta_vals = as.numeric(etaSquared(model_eff, type = 3, anova = T)[,yvars="eta.sq.part"])
        
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
        DF = rbindlist(list(DF, as.list(c(model_vars[k], behav_vars[b],
                                          row,
                                          c_eff_vals,
                                          g-1,GS_in_model,sex[s],model_type))))
        
      } #end s
    } #end g
  }
} #end i


all_col_names = c(all_var_names, "behav")

paES_col_names = paste(all_col_names, "_par_eta_sq",sep="")


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


names(DF) = c("model_yvar","behav_var",paES_col_names,
              c_cols,
              "globvar_in_model","GS_in_model","sex","model_type")

write_xlsx(DF,paste(save_folder,B_effect_sizes_path,sep = ""))









# test some of the models that are to be run 
#with global covariate 
model_bvol_glob = lm(BrainTotalVol ~ group*sex + age + TotalEulerNumber + site + eICV_samseg, data=datab)
Anova(model_bvol_glob,type="III")


model_bvol_glob = lm(mean_thickness ~ group*sex + group*CBCL_ext_cg_v11 + age + TotalEulerNumber + site + eICV_samseg, data=datab)
Anova(model_bvol_glob,type="III")
model_bvol_glob = update(model_bvol_glob,~.-group:sex)
Anova(model_bvol_glob,type="III")
#model_bvol_glob = update(model_bvol_glob,~.-eval(parse(text=paste("group:","CBCL_ext",sep=""))))
#model_bvol_glob = update(model_bvol_glob,~.-group:CBCL_ext)
model_bvol_glob = update(model_bvol_glob, formula=drop.terms(model_bvol_glob$terms, grep( "group:CBCL_ext_cg_v11", attr( model_bvol_glob$terms, "term.labels") ), keep.response=TRUE) )
Anova(model_bvol_glob,type="III")



