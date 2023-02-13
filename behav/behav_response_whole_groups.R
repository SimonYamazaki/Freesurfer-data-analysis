
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
#based on GLMs on data from all groups (which may show an overall group effect)

# models with one of the following behavioural measures as response variables is produced: 
#("CBCL_ext_cg_v11","CBCL_int_cg_v11","CBCL_totsc_cg_v11","CGASx_v11")
#with the following global brain measures as covariates: 
#("BrainTotalVol", "CortexVol", "total_area", "mean_thickness", "eICV_samseg")

#the following covariates are included in all models by default
#("age","site","TotalEulerNumber")


#returns: 
#relevant statistics from GLMs run for each combination of response variable and global brain measures and saved into an excel sheets 
#files generated are: 
# - excel sheet with variable level GLM effects without global covariate (eICV_samseg)
# - excel sheet with variable level GLM effects with global covariate (eICV_samseg)
# - excel sheet group contrasts for GLMs without global covariate (eICV_samseg)
# - excel sheet group contrasts for GLMs with global covariate (eICV_samseg)

#For file names and directories saved into, refer to the variables save_folder, model_ANOVA_with_glob, model_ANOVA_without_glob, model_contrast_with_glob, model_contrast_without_glob

#data is loaded from the path specified in data_path



#save folder for tables:
save_folder = "/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/global_variables/behav/behav_response/"

model_ANOVA_with_glob = "behavvar_GS_ANOVA_pvals_with_glob.xlsx"
model_ANOVA_without_glob = "behavvar_GS_ANOVA_pvals_without_glob.xlsx"

model_contrast_with_glob = "behavvar_Model_contrasts_with_glob.xlsx"
model_contrast_without_glob = "behavvar_Model_contrasts_without_glob.xlsx"



#data path
data_path = "/mnt/projects/VIA11/FREESURFER/Stats/Data/VIA11_allkey_160621_FreeSurfer_pruned_20220509.csv"

#load data
data_csv <- read.table(data_path, header = TRUE, sep = ",", dec = ".")



#inspect the head of data and summary
head(data_csv)
summary(data_csv)

#filter the data with include variable
# - extract rows with 1 in Include_FS_studies coloumn
data_csv_filtered <- data_csv[c(data_csv$Include_FS_studies == 1),]
data_csv_filtered <- data_csv_filtered[!is.na(data_csv_filtered$Include_FS_studies),]


# - extract rows with 1 in Include_FS_studies_euler_outliers_sibpairs_out
data_csv_filtered <- data_csv[c(data_csv$Include_FS_studies_euler_outliers_sibpairs_out == 1),]
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


#### models run with behav as response variable 

model_vars = behav_vars
glob = c("without_eICV","with_eICV")

model = list()

DF = data.frame()

#loop that defines models 
for (i in seq(1,length(model_vars))){
  
  for (g in seq(1,length(glob))){
    
    if (glob[g] == "with_eICV"){
      f = paste(model_vars[i],"~","+","group*sex","+","age","+","site","+","TotalEulerNumber","+","eICV_samseg") #model with ICV
    }
    else{
      f = paste(model_vars[i],"~","+","group*sex","+","age","+","site","+","TotalEulerNumber") #model without ICV
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
      EulerNumber_F = Anova(mm,type = "III")$"F value"[xvars=="TotalEulerNumber"]
      
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
      
      #the first index in m is always the with a group/sex interaction
      if (m==1){
        rw = list(model_vars[i], 
                  group_F, group_pv, sex_F, sex_pv,
                  age_F, age_pv, site_F, site_pv, 
                  EulerNumber_F, EulerNumber_pv, 
                  ICV_F, ICV_pv, GS_F, GS_pvals,
                  sign_GS,m,g-1,
                  subs)
      }
      else{ #if m is 2 then the model is does not include the interaction
        rw = list(model_vars[i], 
                  group_F, group_pv, sex_F, sex_pv, 
                  age_F, age_pv, site_F, site_pv, 
                  EulerNumber_F, EulerNumber_pv, 
                  ICV_F, ICV_pv, NA,   NA,
                  sign_GS,m-2,g-1,
                  subs)
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
              "ICV_Fval","ICV_pval","Group_sex_Fval","Group_sex_pval",
              "Significant_GS_interaction","GS_in_model","ICV_in_model",
              "n_subs")
names(DF)<-col_names



DF_xlsx_glob0 = DF[DF$ICV_in_model == 0, ]
DF_xlsx_glob1 = DF[DF$ICV_in_model == 1, ]

write_xlsx(DF_xlsx_glob1,paste(save_folder,model_ANOVA_with_glob,sep=""))
write_xlsx(DF_xlsx_glob0,paste(save_folder,model_ANOVA_without_glob,sep=""))





### contrast for behav variables as y-variable 

model_yvars = behav_vars
DF = data.frame()

for (k in seq(1,length(model_yvars))){
  for (g in seq(1,2)){
    
    modelb = model[[glob[g]]][[k]]
    
    #use the above defined model_ana 
    ls = lsmeans(modelb,pairwise~"group",adjust="none")
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
    
    
    SZ_diff = -summary(c)$estimate[3]
    SZ_diff_pv = summary(c)$"p.value"[3]
    SZ_diff_tratio = summary(c)$"t.ratio"[3]
    SZ_diff_LCL = -confint(c)$lower.CL[3]
    SZ_diff_UCL = -confint(c)$upper.CL[3]
    
    
    #percent difference
    BP_diff_P = 100*BP_diff /K_emm
    BP_diff_LCL_P = 100*BP_diff_LCL /K_emm
    BP_diff_UCL_P = 100*BP_diff_UCL /K_emm
    
    SZ_diff_P = 100*SZ_diff /K_emm
    SZ_diff_LCL_P = 100*SZ_diff_LCL /K_emm
    SZ_diff_UCL_P = 100*SZ_diff_UCL /K_emm
    
    
    #for contrast xlsx
    rwc_xlsx = list(model_yvars[k],BP_emm,K_emm,SZ_emm,
                    BP_diff,BP_diff_tratio,BP_diff_pv,
                    SZ_diff,SZ_diff_tratio,SZ_diff_pv,
                    g-1)
    
    DF = rbindlist(list(DF, rwc_xlsx))
    
    
  } #g
  
} #k

names(DF) = c("Model_yvar","BP_LSmean","K_LSmean","SZ_LSmean",
              "Contrast_BP-K","tratio_BP-K","pval_BP-K",
              #"LCL_Contrast_BP-K", "UCL_Contrast_BP-K",
              "Contrast_SZ-K","tratio_SZ-K","pval_SZ-K",
              #"LCL_Contrast_SZ-K", "UCL_Contrast_SZ-K",
              "global_var_in_model")


DF_xlsx_glob0 = DF[DF$global_var_in_model == 0, ]
DF_xlsx_glob1 = DF[DF$global_var_in_model == 1, ]

write_xlsx(DF_xlsx_glob1,paste(save_folder,model_contrast_with_glob,sep=""))
write_xlsx(DF_xlsx_glob0,paste(save_folder,model_contrast_without_glob,sep=""))



# test some of the models that are to be run 

#with global covariate 
model_bvol_glob = lm(BrainTotalVol ~ group*sex + age + TotalEulerNumber + site + eICV_samseg, data=datab)
Anova(model_bvol_glob,type="III")

model_bvol_glob = lm(BrainTotalVol ~ group*sex + group*CBCL_ext + age + TotalEulerNumber + site + eICV_samseg, data=datab)
Anova(model_bvol_glob,type="III")
model_bvol_glob = update(model_bvol_glob,~.-group:sex)
Anova(model_bvol_glob,type="III")
#model_bvol_glob = update(model_bvol_glob,~.-eval(parse(text=paste("group:","CBCL_ext",sep=""))))
model_bvol_glob = update(model_bvol_glob,~.-group:CBCL_ext)
#model_bvol_glob = update(model_bvol_glob, formula=drop.terms(model_bvol_glob$terms, grep( "group:CBCL", attr( model_bvol_glob$terms, "term.labels") ), keep.response=TRUE) )
Anova(model_bvol_glob,type="III")


model_bvol_glob = lm(BrainTotalVol ~ group*sex + group*CBCL_ext + age + TotalEulerNumber + site + eICV_samseg, data=datab)
Anova(model_bvol_glob,type="III")
model_bvol_glob = update(model_bvol_glob,~.-group:sex)
Anova(model_bvol_glob,type="III")
#model_bvol_glob = update(model_bvol_glob,~.-eval(parse(text=paste("group:","CBCL_ext",sep=""))))
model_bvol_glob = update(model_bvol_glob,~.-group:CBCL_ext)
#model_bvol_glob = update(model_bvol_glob, formula=drop.terms(model_bvol_glob$terms, grep( "group:CBCL", attr( model_bvol_glob$terms, "term.labels") ), keep.response=TRUE) )
Anova(model_bvol_glob,type="III")


model_CBCL_ext = lm(CBCL_ext ~ group*sex + age + TotalEulerNumber + site, data=datab)
Anova(model_CBCL_ext,type="III")
model_CBCL_ext = update(model_CBCL_ext,~.-group:sex)
Anova(model_CBCL_ext,type="III")
lsmeans(model_CBCL_ext,pairwise~"group", adjust="none")

#model=aov(YIELD~VARIETY) #Build a model with the normal ANOVA command
#res=model$residuals #Create an object of the residuals of Y
#shapiro.test(res)