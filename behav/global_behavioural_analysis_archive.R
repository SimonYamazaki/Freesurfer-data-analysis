

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
#in multiple groups and comparing to a control group. The script generates:


#save folder for tables and plots:
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





###### for loop that generates tables of global brain measure as y-variable 
#with each of the behavioural variables as covariates


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
        #f = paste(model_vars[i],"~","+","age","+","site","+","TotalEulerNumber","+","eICV_samseg") #model with ICV
      }
      else{
        f = paste(model_vars[i],"~","+","group*sex","+","age","+","site","+","TotalEulerNumber") #model without ICV
        #f = paste(model_vars[i],"~","+","age","+","site","+","TotalEulerNumber") #model without ICV
      }
      
      for (b in seq(1,length(behav_vars))){
        print(behav_vars[b])
        
        beh = behav_vars[b]
        ff = paste(f,"+",beh,"*group")
        #ff = paste(f,"+",beh,"*sex")
        
        
        #run the actual model
        model_behav = lm(ff,data=datab)
        
        xvars = attributes(Anova(model_behav,type = "III"))$row.names
        
        GS_F = Anova(model_behav,type = "III")$"F value"[xvars=="group:sex"]
        GS_pv = Anova(model_behav,type = "III")$"Pr(>F)"[xvars=="group:sex"]

        GB_F = Anova(model_behav,type = "III")$"F value"[xvars==paste("group:",beh,sep="")]
        GB_pv = Anova(model_behav,type = "III")$"Pr(>F)"[xvars==paste("group:",beh,sep="")]

        # GS_F = Anova(model_behav,type = "III")$"F value"[xvars==paste(beh,":sex",sep="")]
        # GS_pv = Anova(model_behav,type = "III")$"Pr(>F)"[xvars==paste(beh,":sex",sep="")]

        
        #if the group/sex interaction is insignificant
        if (Anova(model_behav,type = "III")$"Pr(>F)"[xvars=="group:sex"] > 0.05){
          #indicate the significance of the group/sex interation
          sign_GS = 0

          model_gs = model_behav

          model_behav = update(model_behav,~.-group:sex)

          models = list(model_gs, model_behav)

        }
        else{    #if the group/sex interaction is significant
          #indicate the significance of the group/sex interation
          sign_GS = 1
          models = list(model_behav)
        }
        
        #sign_GS = 2
        #models = list(model_behav)
        
        xvars = attributes(Anova(model_behav,type = "III"))$row.names
        
        #group/behav interaction
        if (Anova(model_behav,type = "III")$"Pr(>F)"[xvars==paste("group:",beh,sep="")] > 0.05){
        #if (Anova(model_behav,type = "III")$"Pr(>F)"[xvars==paste(beh,":sex",sep="")] > 0.05){
            
          sign_GB = 0

          #model_gb = model_behav
          model_behav = update(model_behav, formula=drop.terms(model_behav$terms, grep( paste("group:",beh,sep=""), attr( model_behav$terms, "term.labels") ), keep.response=TRUE) )
          #model_behav = update(model_behav, formula=drop.terms(model_behav$terms, grep( paste(beh,":sex",sep=""), attr( model_behav$terms, "term.labels") ), keep.response=TRUE) )
          
          models = list.append(models,model_behav)
        }
        else{    #if the group/sex interaction is significant 
          sign_GB = 1

        }
        xvars = attributes(Anova(model_behav,type = "III"))$row.names
        
      # loop that makes rows in the ANOVA table excel 
      # makes two rows if the group/sex interaction is insignificant 
      # one row with the interaction included and one row with out
      
      for (m in seq(1,length(models))){
        mm = models[[m]]
        xvars = attributes(Anova(mm,type = "III"))$row.names

        #extract relevant statistics for each variable in the model 
        behav_F = Anova(mm,type = "III")$"F value"[xvars==beh]
        #group_F = Anova(mm,type = "III")$"F value"[xvars=="group"]
        sex_F = Anova(mm,type = "III")$"F value"[xvars=="sex"]
        age_F = Anova(mm,type = "III")$"F value"[xvars=="age"]
        site_F = Anova(mm,type = "III")$"F value"[xvars=="site"]
        EulerNumber_F = Anova(mm,type = "III")$"F value"[xvars=="TotalEulerNumber"]
        
        behav_pv = Anova(mm,type = "III")$"Pr(>F)"[xvars==beh]
        #group_pv = Anova(mm,type = "III")$"Pr(>F)"[xvars=="group"]
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
        
        #GB_F = NA
        #GB_pv = NA
        #group_F = NA
        #group_pv = NA
        
        #the first index in m is always the model with a group/sex and group/behav interaction
        if (m==1){
          rw = list(model_vars[i], beh,
                    group_F, group_pv, sex_F, sex_pv,
                    age_F, age_pv, site_F, site_pv, 
                    EulerNumber_F, EulerNumber_pv, 
                    ICV_F, ICV_pv, behav_F, behav_pv, 
                    GB_F, GB_pv, GS_F, GS_pv,
                    sign_GB,sign_GS,1,1,g-1,
                    subs)
        }
        else if (m==2){ #if m is 2 then the model is does not include group/sex interaction
          GB_F = Anova(mm,type = "III")$"F value"[xvars==paste("group:",beh,sep="")]
          GB_pv = Anova(mm,type = "III")$"Pr(>F)"[xvars==paste("group:",beh,sep="")]
          
          rw = list(model_vars[i], beh,
                    group_F, group_pv, sex_F, sex_pv, 
                    age_F, age_pv, site_F, site_pv, 
                    EulerNumber_F, EulerNumber_pv, 
                    ICV_F, ICV_pv, behav_F, behav_pv, 
                    GB_F, GB_pv, GS_F, GS_pv,#NA, NA, 
                    sign_GB,sign_GS,1,0,g-1,
                    subs)
        }
        else if (m==3){ #if m is 3 then the model is does not include group/sex and group/behav interaction
          rw = list(model_vars[i], beh,
                    group_F, group_pv, sex_F, sex_pv, 
                    age_F, age_pv, site_F, site_pv, 
                    EulerNumber_F, EulerNumber_pv, 
                    ICV_F, ICV_pv, behav_F, behav_pv, 
                    NA, NA, NA,NA,
                    sign_GB,sign_GS,0,0,g-1,
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
              "ICV_Fval","ICV_pval","behav_Fval","behav_pval",
              "group_behav_Fval","group_behav_pval","Group_sex_Fval","Group_sex_pval",
              "Significant_GB_interaction","Significant_GS_interaction","GB_in_model","GS_in_model","ICV_in_model",
              "n_subjects")

# col_names = c("Model_yvar", "behav_var",
#               "Group_Fval","Group_pval","Sex_Fval","Sex_pval",
#               "Age_Fval","Age_pval","Site_Fval","Site_pval",
#               "EulerNumber_Fval","Eulernumber_pval",
#               "ICV_Fval","ICV_pval","behav_Fval","behav_pval",
#               "group_behav_Fval","group_behav_pval","behav_sex_Fval","behav_sex_pval",
#               "Significant_GB_interaction","Significant_GS_interaction","GB_in_model","GS_in_model","ICV_in_model",
#               "n_subjects")

names(DF)<-col_names


DF_xlsx_glob0 = DF[DF$ICV_in_model == 0, ]
DF_xlsx_glob1 = DF[DF$ICV_in_model == 1, ]

#write_xlsx(DF_xlsx_glob1,paste(save_folder,"globvar_GS_GB_ANOVA_pvals_with_glob.xlsx",sep=""))
#write_xlsx(DF_xlsx_glob0,paste(save_folder,"globvar_GS_GB_ANOVA_pvals_without_glob.xlsx",sep=""))




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

write_xlsx(DF_xlsx_glob1,paste(save_folder,"behavvar_GS_ANOVA_pvals_with_glob.xlsx",sep=""))
write_xlsx(DF_xlsx_glob0,paste(save_folder,"behavvar_GS_ANOVA_pvals_without_glob.xlsx",sep=""))






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

write_xlsx(DF_xlsx_glob1,paste(save_folder,"behavvar_Model_contrasts_with_glob.xlsx",sep=""))
write_xlsx(DF_xlsx_glob0,paste(save_folder,"behavvar_Model_contrasts_without_glob.xlsx",sep=""))





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




