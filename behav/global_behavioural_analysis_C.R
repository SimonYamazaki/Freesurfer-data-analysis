

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

#This script is intended for computation of statistics on global brain measures
#in multiple groups and comparing to a control group with a behavioral covariate


#save folder for tables and plots:
save_folder = "/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/global_variables/behav/C/"
save_folder_plot = "/mnt/projects/VIA11/FREESURFER/Stats/Plots/global_measures/behav/C/"
#create the folders if they dont exist
dir.create(file.path(dirname(save_folder), basename(save_folder)))
dir.create(file.path(dirname(save_folder_plot), basename(save_folder_plot)))



#C
#includes sheets for whole group (with G/S interaction) and with group/behav interaction
GS_C_model_ANOVA_with_glob= "globvar_GS_GB_ANOVA_pvals_with_glob.xlsx"
GS_C_model_ANOVA_without_glob= "globvar_GS_GB_ANOVA_pvals_without_glob.xlsx"

#C model sheets for each sex, a row for each sex
C_model_ANOVA_with_glob = "globvar_S_GB_ANOVA_pvals_with_glob.xlsx"
C_model_ANOVA_without_glob= "globvar_S_GB_ANOVA_pvals_without_glob.xlsx"

#contrast of the GS and S models above
C_model_contrast_with_glob = "globvar_S_Model_contrasts_with_glob.xlsx"
C_model_contrast_without_glob = "globvar_S_Model_contrasts_without_glob.xlsx"

C_effect_sizes_path = "ANOVA+contrast_effect_sizes.xlsx"



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







## C models whole group


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
      f = paste(model_vars[i],"~","+","group*sex","+",general_cov_formula,"+","eICV_samseg") #model with ICV
    }
    else{
      f = paste(model_vars[i],"~","+","group*sex","+",general_cov_formula) #model without ICV
    }
    
    #loop that adds and analyze models WITH behav
    for (b in seq(1,length(behav_vars))){
      print(behav_vars[b])
      
      beh = behav_vars[b]
      ff = paste(f,"+",beh,"*group")
      
      #run the actual model
      model_behav = lm(ff,data=datab)
      
      xvars = attributes(Anova(model_behav,type = "III"))$row.names
      
      GS_F = Anova(model_behav,type = "III")$"F value"[xvars=="group:sex"]
      GS_pv = Anova(model_behav,type = "III")$"Pr(>F)"[xvars=="group:sex"]
      
      GB_F = Anova(model_behav,type = "III")$"F value"[xvars==paste("group:",beh,sep="")]
      GB_pv = Anova(model_behav,type = "III")$"Pr(>F)"[xvars==paste("group:",beh,sep="")]
      
      
      #if the group/sex interaction is insignificant
      #update the current model_behav if insignificant, and save the old one
      if (Anova(model_behav,type = "III")$"Pr(>F)"[xvars=="group:sex"] > 0.05){
        #indicate the significance of the group/sex interaction
        sign_GS = 0
        
        model_gs = model_behav
        
        model_behav = update(model_behav,~.-group:sex)
        
        models = list(model_gs, model_behav)
        model_type = list("GS+GB","GB")
        
      }
      else{    #if the group/sex interaction is significant - just keep the model
        #indicate the significance of the group/sex interaction
        sign_GS = 1
        models = list(model_behav)
        model_type = list("GS+GB")
      }
      
      
      xvars = attributes(Anova(model_behav,type = "III"))$row.names
      
      #group/behav interaction
      #update the current model_behav if insignificant, and save the old one
      if (Anova(model_behav,type = "III")$"Pr(>F)"[xvars==paste("group:",beh,sep="")] > 0.05){

        sign_GB = 0
        
        model_behav = update(model_behav, formula=drop.terms(model_behav$terms, grep( paste("group:",beh,sep=""), attr( model_behav$terms, "term.labels") ), keep.response=TRUE) )

        models = list.append(models,model_behav)
        
        if (sign_GS == 1){        
          model_type = list.append(model_type, "GS")
        }
        else{          
          model_type = list.append(model_type, "No_interactions")
        }
      }
      else{    #if the group/sex interaction is significant 
        # we dont add any new model - last added model is the appropriate one 
        sign_GB = 1
      }
      xvars = attributes(Anova(model_behav,type = "III"))$row.names
      
      # loop that makes rows in the ANOVA table excel 
      # makes two rows if the group/sex interaction is insignificant 
      # one row with the interaction included and one row without
      # makes 3 rows if the group/behav is insignificant 
      
      for (m in seq(1,length(models))){
        mm = models[[m]]
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
        
        #if the model includes ICV, then also define the anova statistics of ICV
        if (glob[g] == "with_eICV"){
          ICV_F = Anova(mm,type = "III")$"F value"[xvars=="eICV_samseg"]
          ICV_pv = Anova(mm,type = "III")$"Pr(>F)"[xvars=="eICV_samseg"]
        }
        else{ #else, dont define them if the model is without ICV
          ICV_F = NA
          ICV_pv = NA
        }   
        
        subs = length(model_behav$residuals)
        
        #the first index in m is always the model with a group/sex and group/behav interaction
        if (model_type[[m]]=="GS+GB"){
          GB_in = 1
          GS_in = 1
        }
        else if (model_type[[m]]=="GB"){ #the model does not include group/sex interaction
          GB_F = Anova(mm,type = "III")$"F value"[xvars==paste("group:",beh,sep="")]
          GB_pv = Anova(mm,type = "III")$"Pr(>F)"[xvars==paste("group:",beh,sep="")]
          GS_F = NA
          GS_pv = NA
          GB_in = 1
          GS_in = 0
        }
        else if (model_type[[m]]=="GS"){ #the model does not include group/behav interaction
          GB_F = NA
          GB_pv = NA
          GS_F = Anova(mm,type = "III")$"F value"[xvars=="group:sex"]
          GS_pv = Anova(mm,type = "III")$"Pr(>F)"[xvars=="group:sex"]
          GB_in = 0
          GS_in = 1
        }
        else if (model_type[[m]]=="No_interactions"){ #if m is 3 then the model is does not include group/sex and group/behav interaction
          GB_F = NA
          GB_pv = NA
          GS_F = NA
          GS_pv = NA
          GB_in = 0
          GS_in = 0
        }
        
        model[[glob[g]]][[behav_vars[b]]][[model_vars[i]]] = model_behav
        
        #append the row to the dataframe 
        #DF = rbindlist(list(DF, rw))
        DF = rbindlist(list(DF, as.list(c(model_vars[i], beh,
                                          group_F, group_pv, sex_F, sex_pv,
                                          age_F, age_pv, site_F, site_pv, 
                                          EulerNumber_F, EulerNumber_pv, 
                                          ICV_F, ICV_pv, behav_F, behav_pv,
                                          
                                          GB_F, GB_pv, GS_F, GS_pv, 
                                          sign_GB,sign_GS,GB_in,GS_in, 
                                          
                                          g-1,model_type[[m]],subs))))
        
      } # for m 
      
    } # b
  } #g
} #end i


#define the coloum names of the dataframe which is converted to an excel 
col_names = c("Model_yvar", "behav_var",
              "Group_Fval","Group_pval","Sex_Fval","Sex_pval",
              "Age_Fval","Age_pval","Site_Fval","Site_pval",
              "EulerNumber_Fval","Eulernumber_pval",
              "ICV_Fval","ICV_pval","behav_Fval","behav_pval",
              "group_behav_Fval","group_behav_pval","Group_sex_Fval","Group_sex_pval",
              "Significant_GB_interaction","Significant_GS_interaction","GB_in_model","GS_in_model",
              "ICV_in_model","model_type","n_subjects")

names(DF)<-col_names


DF_xlsx_glob0 = DF[DF$ICV_in_model == 0, ]
DF_xlsx_glob1 = DF[DF$ICV_in_model == 1, ]

write_xlsx(DF_xlsx_glob1,paste(save_folder,GS_C_model_ANOVA_with_glob,sep=""))
write_xlsx(DF_xlsx_glob0,paste(save_folder,GS_C_model_ANOVA_without_glob,sep=""))




#C models for each sex

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
        f = paste(model_vars[k],"~","group","+",beh,"+",general_cov_formula)
        
        if (glob[gg] == "with_eICV"){ #run models without ICV
          f = paste(f,"+ eICV_samseg")
        }
        
        if (sex[ss] == "both"){ #define a specific model for the case of both sex
          f = paste(f,"+ sex")
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
          
          rwc_xlsx = append(rwc_xlsx,list(diff, diff_tratio, diff_pv,
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
        
        if (glob[gg] == "without_eICV"){
          #no global measure model
          glob_var_F = NA
          glob_var_pv = NA
        }       
        else{
          #global measure model
          glob_var_F = Anova(mm,type = "III")$"F value"[xvars=="eICV_samseg"]
          glob_var_pv = Anova(mm,type = "III")$"Pr(>F)"[xvars=="eICV_samseg"]
          
        } 
        
        if (sex[ss] == "both"){
          sex_F = Anova(mm,type = "III")$"F value"[xvars=="sex"]
          sex_pv = Anova(mm,type = "III")$"Pr(>F)"[xvars=="sex"]
          model_type = "group + sex"
        }
        else {
          sex_F = NA
          sex_pv = NA
          model_type = "group"
        }
        
        #contrast table
        rwc_xlsx = append(rwc_xlsx,list(gg-1,ss-1,model_type))
        DFc = rbindlist(list(DFc, rwc_xlsx))
        
        #rows for ANOVA xlsx table
        rw_xlsx = list(model_vars[k], beh, 
                       group_F, group_pv, sex_F, sex_pv,
                       age_F, age_pv, site_F, site_pv, 
                       EulerNumber_F, EulerNumber_pv, 
                       behav_F, behav_pv,
                       glob_var_F, glob_var_pv, "eICV_samseg",
                       gg-1,ss-1,model_type)
        
        DF_xlsx = rbindlist(list(DF_xlsx, rw_xlsx))
          
      } # glob 
    } #sex 
  } #behav
} #model_vars

names(DF_xlsx)<-c("Model_yvar","Behav_var",
                  "Group_Fval","Group_pval","Sex_Fval","Sex_pval",
                  "Age_Fval","Age_pval","Site_Fval","Site_pval",
                  "EulerNumber_Fval","Eulernumber_pval",
                  "Behav_Fval","Behav_pval",
                  "global_var_F","global_var_pv","global_var_name",
                  "global_var_in_model","sex","model_type")


#flip the order of contrast names are defined if the sign is flipped
contrast_names = attributes(c)$levels$contrast

for (cn in seq(1,length(contrast_names))){
  if (contrast_signs[cn]<0){
    contrast_names[cn]=paste(strsplit(contrast_names[cn], split = " - ")[[1]][2],strsplit(contrast_names[cn], split = " - ")[[1]][1],sep = " - ")
  }    
}  

c_cols = expand.grid(c("contrast","t-ratio","pval","LCL","UCL"), contrast_names)
c_cols = paste(c_cols$Var1, c_cols$Var2,sep = "_")
c_cols = gsub(" ", "", c_cols)

group_lsmeans = paste("LSmean_",group_levels)

#same excel saving procedure as above 
col_names = c("Model_yvar", group_lsmeans, c_cols,
              "global_var_in_model","sex","model_type")
names(DFc)<-col_names



DF_xlsx_glob0 = DFc[DFc$global_var_in_model == 0, ]
DF_xlsx_glob1 = DFc[DFc$global_var_in_model == 1, ]

write_xlsx(DF_xlsx_glob1,paste(save_folder,C_model_contrast_with_glob,sep = ""))
write_xlsx(DF_xlsx_glob0,paste(save_folder,C_model_contrast_without_glob,sep = ""))


DF_xlsx_glob0 = DF_xlsx[DF_xlsx$global_var_in_model == 0, ]
DF_xlsx_glob1 = DF_xlsx[DF_xlsx$global_var_in_model == 1, ]

write_xlsx(DF_xlsx_glob1,paste(save_folder,C_model_ANOVA_with_glob,sep = ""))
write_xlsx(DF_xlsx_glob0,paste(save_folder,C_model_ANOVA_without_glob,sep = ""))





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

#behav interaction names
behav_inter = expand.grid(c("group:"), behav_vars)
behav_inter = paste(behav_inter$Var1, behav_inter$Var2,sep = "")
behav_inter = gsub(" ", "", behav_inter)

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
        
        all_col_names = c(all_var_names, behav_vars[b], paste("group:",behav_vars[b],sep=""))
        
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
        #row[var_idx] = et[as.numeric(sapply(at, function(x) { grep(x, xt) }))]
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


all_col_names = c(all_var_names, "behav", "behav:group")

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

write_xlsx(DF,paste(save_folder,C_effect_sizes_path,sep = ""))









# test some of the models that are to be run 
#with global covariate 
model_bvol_glob = lm(BrainTotalVol ~ group*sex + age + TotalEulerNumber + site + eICV_samseg, data=datab)
Anova(model_bvol_glob,type="III")


model_bvol_glob = lm(BrainTotalVol ~ group*sex + group*CBCL_ext_cg_v11 + age + TotalEulerNumber + site + eICV_samseg, data=datab)
#model_bvol_glob = lm(BrainTotalVol ~ group + group*CBCL_ext_cg_v11 + age + TotalEulerNumber + site + eICV_samseg, data=dataf[[1]])

Anova(model_bvol_glob,type="III")
model_bvol_glob = update(model_bvol_glob,~.-group:sex)
Anova(model_bvol_glob,type="III")
#model_bvol_glob = update(model_bvol_glob,~.-eval(parse(text=paste("group:","CBCL_ext",sep=""))))
#model_bvol_glob = update(model_bvol_glob,~.-group:CBCL_ext)
model_bvol_glob = update(model_bvol_glob, formula=drop.terms(model_bvol_glob$terms, grep( "group:CBCL_ext_cg_v11", attr( model_bvol_glob$terms, "term.labels") ), keep.response=TRUE) )
Anova(model_bvol_glob,type="III")

etaSquared(model_bvol_glob, type = 3, anova = T)

