

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
#in multiple groups and comparing to a control group. The script generates:

#save paths:
GS_ANOVA_with_cov = "/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/global_variables/globvar_GS_ANOVA_pvals_with_glob.xlsx"
GS_ANOVA_without_cov = "/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/global_variables/globvar_GS_ANOVA_pvals_without_glob.xlsx"

ANOVA_with_cov = "/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/global_variables/globvar_ANOVA_pvals_with_glob.xlsx"
ANOVA_without_cov = "/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/global_variables/globvar_ANOVA_pvals_without_glob.xlsx"

#save paths:
contrast_with_cov = "/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/global_variables/globvar_S_Model_contrasts_with_glob.xlsx"
contrast_without_cov = "/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/global_variables/globvar_S_Model_contrasts_without_glob.xlsx"


#save folder for plots:
save_folder = "/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/global_variables/behav/"
save_folder_plot = "/mnt/projects/VIA11/FREESURFER/Stats/Plots/global_measures/behav/"


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


#make new variables with shorter and contained names
# - tell r which variables are factors
datab = data_csv_filtered
datab$sex = as.factor(datab$Sex_child)
datab$group = as.factor(datab$HighRiskStatus)
datab$site = as.factor(datab$MR_Site)
datab$diag = as.factor(datab$Axis.1_diag_v11)
datab$age = as.numeric(datab$MRI_age)


#add a string encoding of sex for later convenience 
datab$sexs = datab$Sex_child
datab$sexs[datab$sexs == 0] = "female"
datab$sexs[datab$sexs == 1] = "male"

#set another working directory so that files are saved in this folder
#setwd(save_folder)


#with global covariate 
model_bvol_glob = lm(BrainTotalVol ~ group*sex + age + TotalEulerNumber + site + eICV_samseg, data=datab)
Anova(model_bvol_glob,type="III")

model_bvol_glob = lm(BrainTotalVol ~ group*sex + group*CBCL_ext + age + TotalEulerNumber + site + eICV_samseg, data=datab)
Anova(model_bvol_glob,type="III")
model_bvol_glob = update(model_bvol_glob,~.-group:sex)
Anova(model_bvol_glob,type="III")
#model_bvol_glob = update(model_bvol_glob,~.-eval(parse(text=paste("group:","CBCL_ext",sep=""))))
#model_bvol_glob = update(model_bvol_glob,~.-group:CBCL_ext)
model_bvol_glob = update(model_bvol_glob, formula=drop.terms(model_bvol_glob$terms, grep( "group:CBCL", attr( model_bvol_glob$terms, "term.labels") ), keep.response=TRUE) )
Anova(model_bvol_glob,type="III")


model_CBCL_ext = lm(CBCL_ext ~ group*sex + age + TotalEulerNumber + site + eICV_samseg, data=datab)
Anova(model_CBCL_ext,type="III")
model_CBCL_ext = update(model_CBCL_ext,~.-group:sex)
Anova(model_CBCL_ext,type="III")





### check normal distribution of behavioural data in each group
#hist(datab$CBCL_ext[datab$group == "SZ"])
#hist(datab$CBCL_ext[datab$group == "BP"])
#hist(datab$CBCL_ext[datab$group == "K"])


#model=aov(YIELD~VARIETY) #Build a model with the normal ANOVA command
#res=model$residuals #Create an object of the residuals of Y
#shapiro.test(res)


###### the for loop


#specify dependent variable (global variable) column names, i.e. left hand side of statistical models
model_vars = c("BrainTotalVol", "CortexVol", "total_area", "mean_thickness", "eICV_samseg")
behav_vars = c("CBCL_ext","CBCL_int","CBCL_totsc","CGAS")


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
      
      if (glob[g] == "with_eICV"){
        f = paste(model_vars[i],"~","+","group*sex","+","age","+","site","+","TotalEulerNumber","+","eICV_samseg") #model with ICV
      }
      else{
        f = paste(model_vars[i],"~","+","group*sex","+","age","+","site","+","TotalEulerNumber") #model without ICV
      }
      
      for (b in seq(1,length(behav_vars))){
        print(behav_vars[b])
        
        beh = behav_vars[b]
        ff = paste(f,"+",beh,"*group")
        
        #run the actual model
        model_behav = lm(ff,data=datab)
        xvars = attributes(Anova(model_behav,type = "III"))$row.names
        
        #if the group/sex interaction is insignificant 
        if (Anova(model_behav,type = "III")$"Pr(>F)"[xvars=="group:sex"] > 0.05){
          #indicate the significance of the group/sex interation
          sign_GS = 0
          model_behav = update(model_behav,~.-group:sex)
        }
        else{    #if the group/sex interaction is significant 
          #indicate the significance of the group/sex interation
          sign_GS = 1
        }
        xvars = attributes(Anova(model_behav,type = "III"))$row.names
        
        if (Anova(model_behav,type = "III")$"Pr(>F)"[xvars==paste("group:",beh,sep="")] > 0.05){
          sign_GB = 0

          GB_F = Anova(model_behav,type = "III")$"F value"[xvars==paste("group:",beh,sep="")]
          GB_pv = Anova(model_behav,type = "III")$"Pr(>F)"[xvars==paste("group:",beh,sep="")]
          
          model_gb = model_behav
          model_behav = update(model_behav, formula=drop.terms(model_behav$terms, grep( paste("group:",beh,sep=""), attr( model_behav$terms, "term.labels") ), keep.response=TRUE) )
          
          models = list(model_gb,model_behav)
        }
        else{    #if the group/sex interaction is significant 
          sign_GB = 1

          GB_F = Anova(model_behav,type = "III")$"F value"[xvars==paste("group:",beh,sep="")]
          GB_pv = Anova(model_behav,type = "III")$"Pr(>F)"[xvars==paste("group:",beh,sep="")]
          
          models = list(model_behav)
        }
        xvars = attributes(Anova(model_behav,type = "III"))$row.names
        
      # loop that makes rows in the ANOVA table excel 
      # makes two rows if the group/sex interaction is insignificant 
      # one row with the interaction included and one row with out
      
      for (m in seq(1,length(models))){
        mm = models[[m]]
        
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
        
        #the first index in m is always the with a group/sex interaction
        if (m==1){
          rw = list(model_vars[i], beh,
                    group_F, group_pv, sex_F, sex_pv,
                    age_F, age_pv, site_F, site_pv, 
                    EulerNumber_F, EulerNumber_pv, 
                    ICV_F, ICV_pv, #GS_F, GS_pvals,
                    behav_F, behav_pv, GB_F, GB_pv,
                    sign_GB,sign_GS,1,m,g-1,
                    subs)
        }
        else{ #if m is 2 then the model is does not include the interaction
          rw = list(model_vars[i], beh,
                    group_F, group_pv, sex_F, sex_pv, 
                    age_F, age_pv, site_F, site_pv, 
                    EulerNumber_F, EulerNumber_pv, 
                    ICV_F, ICV_pv, #NA,   NA,
                    behav_F, behav_pv, NA, NA,
                    sign_GB,sign_GS,0,m-2,g-1,
                    subs)
        }
        
        #append the row to the dataframe 
        DF = rbindlist(list(DF, rw))
        
      } # for m 
        
        png(file=paste(save_folder_plot,model_vars[i],beh,"model_diagnostics.png",sep = "_"))
        model_diag = model_behav
        mixed_res = rstudent(model_diag)
        par(mar=c(1,1,1,1))
        par(mfrow=c(1,3))
        qqnorm(mixed_res,main = NULL)
        qqline(mixed_res,main = NULL)
        title("qq-plot of residuals")
        plot(mixed_res ~ fitted(model_diag),xlab="Fitted",ylab="Standardized residuals")
        title("std residuals vs fitted values")
        boxplot(mixed_res ~ datab$group[!is.na(datab[[beh]])])
        title("residual variance in groups")
        par(mfrow=c(1,1))
        dev.off()
    } # b
  } #g
} #end i

#define the coloum names of the dataframe which is converted to an excel 
col_names = c("Model_yvar", "behav_var",
              "Group_Fval","Group_pval","Sex_Fval","Sex_pval",
              "Age_Fval","Age_pval","Site_Fval","Site_pval",
              "EulerNumber_Fval","Eulernumber_pval",
              "ICV_Fval","ICV_pval",#"Group_sex_Fval","Group_sex_pval",
              "behav_Fval","behav_pval","group_behav_Fval","group_behav_pval",
              "Significant_GB_interaction","Significant_GS_interaction","GB_in_model","GS_in_model","ICV_in_model",
              "n_subjects")

names(DF)<-col_names


DF_xlsx_glob0 = DF[DF$ICV_in_model == 0, ]
DF_xlsx_glob1 = DF[DF$ICV_in_model == 1, ]

write_xlsx(DF_xlsx_glob1,paste(save_folder,"globvar_GS_GB_ANOVA_pvals_with_glob.xlsx",sep=""))
write_xlsx(DF_xlsx_glob0,paste(save_folder,"globvar_GS_GB_ANOVA_pvals_without_glob.xlsx",sep=""))




#### behav as response variable 

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
      
      #the first index in m is always the with a group/sex interaction
      if (m==1){
        rw = list(model_vars[i], 
                  group_F, group_pv, sex_F, sex_pv,
                  age_F, age_pv, site_F, site_pv, 
                  EulerNumber_F, EulerNumber_pv, 
                  ICV_F, ICV_pv, GS_F, GS_pvals,
                  sign_GS,m,g-1)
      }
      else{ #if m is 2 then the model is does not include the interaction
        rw = list(model_vars[i], 
                  group_F, group_pv, sex_F, sex_pv, 
                  age_F, age_pv, site_F, site_pv, 
                  EulerNumber_F, EulerNumber_pv, 
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
              "ICV_Fval","ICV_pval","Group_sex_Fval","Group_sex_pval",
              "Significant_GS_interaction","GS_in_model","ICV_in_model")
names(DF)<-col_names



DF_xlsx_glob0 = DF[DF$ICV_in_model == 0, ]
DF_xlsx_glob1 = DF[DF$ICV_in_model == 1, ]

write_xlsx(DF_xlsx_glob1,paste(save_folder,"behavvar_GS_ANOVA_pvals_with_glob.xlsx",sep=""))
write_xlsx(DF_xlsx_glob0,paste(save_folder,"behavvar_GS_ANOVA_pvals_without_glob.xlsx",sep=""))










