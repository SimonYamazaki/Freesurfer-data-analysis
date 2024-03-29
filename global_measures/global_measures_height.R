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


#This script is intended for computation of statistics on global brain measures
#measures in multiple groups based on GLM. Also includes post-hoc comparisons to the control group.

#the following covariates are included in all models by default
#("age","site","TotalEulerNumber","height")


#returns: 
#files generated are: 

#GS models (with group-sex interaction)
# - excel sheet with variable level GLM effects with global covariate (eICV_samseg)
# - excel sheet with variable level GLM effects without global covariate (eICV_samseg)

#S models
# - excel sheet with variable level GLM effects run separately on data from each sex with global covariate (eICV_samseg)
# - excel sheet with variable level GLM effects run separately on data from each sex without global covariate (eICV_samseg)

#post-hoc
# - excel sheet with post-hoc group contrasts for GLMs (both GS and S models from above)without global covariate (eICV_samseg)
# - excel sheet with post-hoc group contrasts for GLMs (both GS and S models from above) with global covariate (eICV_samseg)
# - excel sheet with GLM effect sizes

# - group difference from control lsmeans and p-value plots 

#For file names and directories saved into, refer to the variables GS_ANOVA_with_glob, GS_ANOVA_without_cov, ANOVA_with_cov, ANOVA_without_cov, contrast_with_cov and contrast_without_cov
#and the save_folder and plot_folder definitions in this script
#data is loaded from the path specified in data_path


# what folder should the tables be saved in
save_folder = "/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/global_variables/height/"

#save paths:
GS_ANOVA_with_cov = "globvar_GS_ANOVA_pvals_with_glob.xlsx"
GS_ANOVA_without_cov = "globvar_GS_ANOVA_pvals_without_glob.xlsx"

S_ANOVA_with_cov = "globvar_S_ANOVA_pvals_with_glob.xlsx"
S_ANOVA_without_cov = "globvar_S_ANOVA_pvals_without_glob.xlsx"

contrast_with_cov = "globvar_S_Model_contrasts_with_glob.xlsx"
contrast_without_cov = "globvar_S_Model_contrasts_without_glob.xlsx"

effect_sizes_path = "ANOVA+contrast_effect_sizes.xlsx"



#save folder for plots:
plot_folder = "/mnt/projects/VIA11/FREESURFER/Stats/Plots/global_measures/height/"

#file name postfix
ppwp_postfix = "_group_diff_pvalues_ICV"


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
datab$height = as.numeric(datab$MRI_height_v11)


#what groups should be used
datab$group = as.factor(datab$HighRiskStatus_v11)

#add a string encoding of sex for later convenience 
datab$sexs = datab$Sex_child
datab$sexs[datab$sexs == 0] = "female"
datab$sexs[datab$sexs == 1] = "male"

#set another working directory so that files are saved in this folder
setwd(save_folder)

#specify dependent variable (global variable) column names, i.e. left hand side of statistical models
y_vars = c("BrainTotalVol", "CortexVol", "total_area", "mean_thickness", "eICV_samseg")

#covariates that are included in all models 
general_cov = c("age","site","TotalEulerNumber","height")
general_cov_formula = paste(general_cov,collapse = " + ")

##############################
# whole brain measure models #
##############################


### run models with eICV_samseg as covariate and without
# generate ANOVA tables with group/sex interaction
# each row is a model 

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
  if (model_vars[i] == "eICV_samseg"){    #only run models without_eICV for the model on eICV_samseg
    glob = c("without_eICV")
  }
  else{
    glob = c("without_eICV","with_eICV") #run both for all other variables
  }
  
  for (g in seq(1,length(glob))){
    if (glob[g] == "with_eICV"){
      f = paste(model_vars[i],"~","+","group*sex","+",general_cov_formula,"+","eICV_samseg") #model with ICV
    }
    else{
      f = paste(model_vars[i],"~","+","group*sex","+",general_cov_formula) #model without ICV
    }
    #run the actual model
    model[[glob[g]]][[i]] = lm(f,data=datab)  
    
    #extract the variables that are modelling the global variable
    xvars = attributes(Anova(model[[glob[g]]][[i]],type = "III"))$row.names
    
    #if the group/sex interaction is insignificant 
    if (Anova(model[[glob[g]]][[i]],type = "III")$"Pr(>F)"[xvars=="group:sex"] > 0.05){
      #indicate the significance of the group/sex interation
      sign_GS = 0
      
      #save the model with group/sex interaction included
      model_gs = model[[glob[g]]][[i]]
      
      #save the anova F and p-value for the group/sex interaction
      GS_pvals = Anova(model_gs,type = "III")$"Pr(>F)"[xvars=="group:sex"]
      GS_F = Anova(model_gs,type = "III")$"F value"[xvars=="group:sex"]
      
      #save the model without the group/sex interaction
      model[[glob[g]]][[i]] = update(model[[glob[g]]][[i]],~.-group:sex)
      
      #save both models in a list
      models = list(model_gs,model[[glob[g]]][[i]])
      
    }
    else{    #if the group/sex interaction is significant 
      
      #indicate the significance of the group/sex interation
      sign_GS = 1
      
      #save the anova F and p-value for the group/sex interaction
      GS_pvals = Anova(model[[glob[g]]][[i]],type = "III")$"Pr(>F)"[xvars=="group:sex"]
      GS_F = Anova(model[[glob[g]]][[i]],type = "III")$"F value"[xvars=="group:sex"]
      
      #save the model with the group/sex interaction
      models = list(model[[glob[g]]][[i]])
    }
    
    # loop that makes rows in the ANOVA table excel 
    # makes two rows if the group/sex interaction is insignificant 
    # one row with the interaction included and one row with out
    for (m in seq(1,length(models))){
      mm = models[[m]]
      
      #extract relevant statistics for each variable in the model 
      group_F = Anova(mm,type = "III")$"F value"[xvars=="group"]
      sex_F = Anova(mm,type = "III")$"F value"[xvars=="sex"]
      age_F = Anova(mm,type = "III")$"F value"[xvars=="age"]
      site_F = Anova(mm,type = "III")$"F value"[xvars=="site"]
      height_F = Anova(mm,type = "III")$"F value"[xvars=="height"]
      EulerNumber_F = Anova(mm,type = "III")$"F value"[xvars=="TotalEulerNumber"]
      
      group_pv = Anova(mm,type = "III")$"Pr(>F)"[xvars=="group"]
      sex_pv = Anova(mm,type = "III")$"Pr(>F)"[xvars=="sex"]
      age_pv = Anova(mm,type = "III")$"Pr(>F)"[xvars=="age"]
      site_pv = Anova(mm,type = "III")$"Pr(>F)"[xvars=="site"]
      height_pv = Anova(mm,type = "III")$"Pr(>F)"[xvars=="height"]
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
      
      #the first index in m is always the with a group/sex interaction
      if (m==1){
        rw = list(model_vars[i], 
                  group_F, group_pv, sex_F, sex_pv,
                  age_F, age_pv, site_F, site_pv, 
                  EulerNumber_F, EulerNumber_pv,
                  height_F, height_pv,
                  ICV_F, ICV_pv, GS_F, GS_pvals,
                  sign_GS,m,g-1)
      }
      else{ #if m is 2 then the model is does not include the interaction
        rw = list(model_vars[i], 
                  group_F, group_pv, sex_F, sex_pv, 
                  age_F, age_pv, site_F, site_pv, 
                  EulerNumber_F, EulerNumber_pv, 
                  height_F, height_pv,
                  ICV_F, ICV_pv, NA,   NA,
                  sign_GS,m-2,g-1)
      }
      
      #append the row to the dataframe 
      DF = rbindlist(list(DF, rw))
      
    } # for m 
  } #g
} #end i

#define the coloum names of the dataframe which is converted to an excel 
col_names = c("Model_yvar",
              "Group_Fval","Group_pval","Sex_Fval","Sex_pval",
              "Age_Fval","Age_pval","Site_Fval","Site_pval",
              "EulerNumber_Fval","Eulernumber_pval",
              "Height_Fval", "Height_pval",
              "ICV_Fval","ICV_pval","Group_sex_Fval","Group_sex_pval",
              "Significant_GS_interaction","GS_in_model","eICV_in_model")
names(DF)<-col_names

#A potential way of reducing significant digits 
#df_with_ICV[,2:ncol(df_with_ICV)] <- signif(df_with_ICV[,2:ncol(df_with_ICV)],digits=3)
#df_without_ICV[,2:ncol(df_without_ICV)] <- signif(df_without_ICV[,2:ncol(df_without_ICV)],digits=3)

#split up the dataframe in models with ICV and without
DF_xlsx_ICV0 = DF[DF$eICV_in_model == 0, ]
DF_xlsx_ICV1 = DF[DF$eICV_in_model == 1, ]

#save the two dataframes in the a specified save path 
write_xlsx(DF_xlsx_ICV1,paste(save_folder,GS_ANOVA_with_cov,sep = ""))
write_xlsx(DF_xlsx_ICV0,paste(save_folder,GS_ANOVA_without_cov,sep = ""))




#### CONTRAST XLSX with male / female split
# and built models[sex][glob][yvar]

#define some relevant naming 
model_vars = y_vars
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
  for (ss in seq(1,3)){ #use the 3 defined sexes to know what data to use
    
    # define an temporary sex indicator to index
    print(sex[ss])
    
    #extract the dataframe containing the relevant sex data
    df_sex = dataf[[ss]]
    
    #only continue with both sex data when the global variable is mean_thickness 
    #if (model_vars[k] == "mean_thickness" & sex[ss] != "both"){
    #  next
    #}
    if (model_vars[k] != "mean_thickness" & sex[ss] == "both"){
      next
    }
    
    
    for (gg in seq(1,2)){
      
      if (glob[gg] == "without_eICV"){ #run models without ICV
        #define the models 
        if (sex[ss] == "both"){ #define a specific model for the case of both sex
          f = paste(model_vars[k],"~","group","+","sex","+",general_cov_formula)
        }
        else { #another model if it is sex-separated
          f = paste(model_vars[k],"~","group","+",general_cov_formula)
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
          f2 = paste(model_vars[k],"~","group","+","sex","+",general_cov_formula,"+","eICV_samseg")
        }
        else {
          f2 = paste(model_vars[k],"~","group","+",general_cov_formula,"+","eICV_samseg")
        }        
        models[[sex[ss]]][[glob[gg]]][[model_vars[k]]] = lm(f2,data=df_sex)
        model_ana = models[[sex[ss]]][[glob[gg]]][[model_vars[k]]]
        xvars = attributes(Anova(model_ana,type = "III"))$row.names
      }
      
      
      #use the above defined model_ana 
      #compute lsmeans pairwise for each group with no correction on p-values
      ls = lsmeans(model_ana,pairwise~"group", adjust="none")
      
      #compute contrast statistics
      c = ls$contrasts
      
      #extract the lsmean of the control group
      K_emm = summary(ls)$lsmeans$lsmean[2]
      
      #estimates and statistics on raw contrasts
      #for bipolar group
      BP_diff = summary(c)$estimate[1]
      BP_diff_pv = summary(c)$"p.value"[1]
      BP_diff_tratio = summary(c)$"t.ratio"[1]
      BP_diff_LCL = confint(c)$lower.CL[1]
      BP_diff_UCL = confint(c)$upper.CL[1]
      
      #for Skizo 
      SZ_diff = -summary(c)$estimate[3]
      SZ_diff_pv = summary(c)$"p.value"[3]
      SZ_diff_tratio = summary(c)$"t.ratio"[3]
      SZ_diff_LCL = -confint(c)$lower.CL[3]
      SZ_diff_UCL = -confint(c)$upper.CL[3]
      
      
      #compute percent difference from control
      BP_diff_P = 100*BP_diff /K_emm
      BP_diff_LCL_P = 100*BP_diff_LCL /K_emm
      BP_diff_UCL_P = 100*BP_diff_UCL /K_emm
      
      SZ_diff_P = 100*SZ_diff /K_emm
      SZ_diff_LCL_P = 100*SZ_diff_LCL /K_emm
      SZ_diff_UCL_P = 100*SZ_diff_UCL /K_emm
      
      
      #rows for contrast xlsx table
      rwc_xlsx = list(model_vars[k],K_emm,
                      BP_diff,BP_diff_tratio,BP_diff_pv,
                      SZ_diff,SZ_diff_tratio,SZ_diff_pv,
                      gg-1,ss-1)
      
      DFc = rbindlist(list(DFc, rwc_xlsx))
      
      
      if (sex[ss] != "both"){
        mm=model_ana
        
        #ANOVA sex separated
        group_F = Anova(mm,type="III")$"F value"[xvars=="group"]
        age_F = Anova(mm,type = "III")$"F value"[xvars=="age"]
        site_F = Anova(mm,type = "III")$"F value"[xvars=="site"]
        height_F = Anova(mm,type = "III")$"F value"[xvars=="height"]
        EulerNumber_F = Anova(mm,type = "III")$"F value"[xvars=="TotalEulerNumber"]
        
        group_pv = Anova(mm,type = "III")$"Pr(>F)"[xvars=="group"]
        age_pv = Anova(mm,type = "III")$"Pr(>F)"[xvars=="age"]
        site_pv = Anova(mm,type = "III")$"Pr(>F)"[xvars=="site"]
        height_pv = Anova(mm,type = "III")$"Pr(<F)"[xvars=="height"]
        EulerNumber_pv = Anova(mm,type = "III")$"Pr(>F)"[xvars=="TotalEulerNumber"]
        
        if (glob[g] == "without_eICV"){
          #no global measure model
          glob_var_F = NA
          glob_var_pv = NA
        }       
        else{
          #global measure model
          glob_var_F = Anova(mm,type = "III")$"F value"[xvars=="eICV_samseg"]
          glob_var_pv = Anova(mm,type = "III")$"Pr(>F)"[xvars=="eICV_samseg"]
          
        } 
        
        
        #rows for ANOVA xlsx table
        rw_xlsx = list(model_vars[k], group_F, group_pv, 
                       age_F, age_pv, site_F, site_pv, 
                       EulerNumber_F, EulerNumber_pv, 
                       height_F, height_pv,
                       glob_var_F, glob_var_pv, "eICV_samseg",
                       gg-1,ss-1)
        
        DF_xlsx = rbindlist(list(DF_xlsx, rw_xlsx))
        
      } #fi both sex
      
    } #g 
  } #s
} #k

names(DF_xlsx)<-c("Model_yvar","Group_Fval","Group_pval",
                  "Age_Fval","Age_pval","Site_Fval","Site_pval",
                  "EulerNumber_Fval","Eulernumber_pval",
                  "height_Fval","height_pval",
                  "global_var_F","global_var_pv","global_var_name",
                  "global_var_in_model","sex")

#same excel saving procedure as above 
col_names = c("Model_yvar", "K_LSmean",
              "Contrast_BP-K","tratio_BP-K","pval_BP-K",
              "Contrast_SZ-K","tratio_SZ-K","pval_SZ-K",
              "global_var_in_model","sex")
names(DFc)<-col_names

DF_xlsx_glob0 = DFc[DFc$global_var_in_model == 0, ]
DF_xlsx_glob1 = DFc[DFc$global_var_in_model == 1, ]

write_xlsx(DF_xlsx_glob1,paste(save_folder,contrast_with_cov,sep = ""))
write_xlsx(DF_xlsx_glob0,paste(save_folder,contrast_without_cov,sep = ""))


DF_xlsx_glob0 = DF_xlsx[DF_xlsx$global_var_in_model == 0, ]
DF_xlsx_glob1 = DF_xlsx[DF_xlsx$global_var_in_model == 1, ]

write_xlsx(DF_xlsx_glob1,paste(save_folder,S_ANOVA_with_cov,sep = ""))
write_xlsx(DF_xlsx_glob0,paste(save_folder,S_ANOVA_without_cov,sep = ""))


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
#Anova(model_eff,type = "III")

#model_eff = lm(BrainTotalVol ~ group + age + site + TotalEulerNumber, data=dataf[[1]])

#emm_eff = emmeans(model_eff,specs="group")

#cohens d
#eff_size(emm_eff, sigma=sigma(model_eff), edf=df.residual(model_eff))


library(lsr)
#etaSquared(model_eff, type = 3, anova = T)

#xvars = attributes(etaSquared(model_eff, type = 3, anova = T))$dimnames[[1]]
#yvars = attributes(etaSquared(model_eff, type = 3, anova = T))$dimnames[[2]]

#etaSquared(model_eff, type = 3, anova = T)[xvars=="group",yvars="eta.sq.part"]



model_vars = c("BrainTotalVol", "CortexVol", "total_area", "mean_thickness")

glob = c("without_eICV","with_eICV")
sex = c("female","male","both") #3 sexes are defined, 0=female, 1=male, 2=both


DF = data.frame()

all_var_names = c("group",general_cov,
                  "eICV_samseg","group:sex")

for (k in seq(1,length(model_vars))){
  for (g in seq(1,length(glob))){
    for (s in seq(1,length(sex))){
      
      if (sex[s] == "both"){
        model_eff = model[[glob[g]]][[k]]
      }
      else {
        model_eff = models[[sex[s]]][[glob[g]]][[model_vars[k]]]
      }
      
      xvars = attributes(etaSquared(model_eff, type = 3, anova = T))$dimnames[[1]]
      yvars = attributes(etaSquared(model_eff, type = 3, anova = T))$dimnames[[2]]
      
      eta_vals = as.numeric(etaSquared(model_eff, type = 3, anova = T)[,yvars="eta.sq.part"])
      
      var_idx = all_var_names %in% xvars
      
      row = rep(NA,length(all_col_names))
      row[var_idx] = eta_vals
      
      
      emm_eff = emmeans(model_eff,specs="group")
      effs = eff_size(emm_eff, sigma=sigma(model_eff), edf=df.residual(model_eff))
      effs_BP = summary(effs)$effect.size[1]
      effs_SZ = -summary(effs)$effect.size[3]
      
      GS_in_model = ("group:sex" %in% attributes(Anova(model_eff,type = "III"))$row.names)*1
      
      #rows for effect size xlsx table
      DF = rbindlist(list(DF, as.list(c(model_vars[k],row,
                                        effs_BP,effs_SZ,
                                        g-1,GS_in_model,sex[s]))))
      
      
    } #end s
  } #end g
} #end i

all_col_names = paste(all_var_names, "_par_eta_sq",sep="")

names(DF) = c("model_yvar",all_col_names,
              "BP-K_contrast_cohensD","SZ-K_contrast_cohensD",
              "globvar_in_model","GS_in_model","sex")

write_xlsx(DF,effect_sizes_path)

