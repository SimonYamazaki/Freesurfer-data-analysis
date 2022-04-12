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


#This script is intended for computation of statistics on bilateral regional brain 
#measures in multiple groups and comparing to a control group. The script generates:


#####
# - ANOVA tables with models that include a group/sex interaction 
#   an extra row of a model without the interaction is included if it turned out insignificant
#   saved in an excel sheet with a global covariate and an excel sheet without

#save paths:
GS_ANOVA_with_glob = "/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/parcels/Parcel_GS_ANOVA_pvals_with_glob.xlsx"
GS_ANOVA_without_glob = "/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/parcels/Parcel_GS_ANOVA_pvals_without_glob.xlsx"


#####
# - An excel sheet with model relevant contrasts for each of the models in the ANOVA tables
#   saved in an excel sheet with a global covariate and an excel sheet without

#save paths:
contrast_with_glob ="/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/parcels/Parcel_thickness_model_contrast_with_glob.xlsx"
contrast_without_glob = "/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/parcels/Parcel_thickness_model_contrast_without_glob.xlsx"


#####
# - Plots LSmeans based contrasts which compares BP and SZ to control 
#   one plot which includes model contrasts from model with global covarate and without
#   no split into sex, however contrasts are based on LSmean sex

#prefix on violin plot
LSmeans_prefix = "LSmean_difference_thickness_gender_combined"


#save folder for plots:
save_folder = "/mnt/projects/VIA11/FREESURFER/Stats/Plots/global_measures"

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
setwd("/mnt/projects/VIA11/FREESURFER/Stats/Plots/bilateral")




######################################
#    run models on each region
######################################

#### extract model yvars and make the bilateral coloumns ####

#compute statistics on these measures
m = c("area","thickness","volume")
col_names = names(datab)

#find names of all the regions based on coloumn names that include "lh_" 
regions = col_names[grepl("^lh_", col_names)]
regions <- unique(unlist(strsplit(regions,"_")))

#remove some specific coloumns
col2rm = c(m,"WhiteSurfArea","lh","MeanThickness")
regions = regions[!(regions %in% col2rm)] 
col_names_list = list()

#loop that extract all the relevant coloumns based on the regions variable
#and computes new bilateral coloumns, based on the mean value of the hemispheres
for (j in seq(1,length(m))){
  col_names_list[[m[j]]] = list()
  
  for (i in seq(1,length(regions))){
    mm = paste(regions[i],m[j],sep = "_")
    bmm = paste("bi_",mm,sep = "")
    col_names_list[[m[j]]] = c(col_names_list[[m[j]]],bmm)
    cols = col_names[grepl(paste("",mm,sep = "_"), col_names)]
    datab[bmm] = rowMeans(datab[cols])
  }
}

#make new variables with shorter and contained names
# - tell r which variables are factors
datab$sex = as.factor(datab$Sex_child)
datab$group = as.factor(datab$HighRiskStatus)
datab$site = as.factor(datab$MR_Site)
datab$diag = as.factor(datab$Axis.1_diag_v11)
datab$age = as.numeric(datab$MRI_age)



###### Make inference with models

#initialize some variables
glob = c("Without global var","With global var")
DF = data.frame()

#loop that 
for (j in seq(1,3)){
  print(m[j])
  model_yvars = unlist(col_names_list[[m[j]]] )
  for (k in seq(1,length(model_yvars))){
    for (g in seq(1,2)){
      
      if (glob[g] == "Without global var"){
        #define model formulation for models without the global variable 
        f = paste(model_yvars[k],"~","+","group*sex","+","age","+","site","+","TotalEulerNumber")
        
        #run the model
        model = lm(f,data=datab)
        xvars = attributes(Anova(model,type = "III"))$row.names
        
        #define the p-value of the group/sex interaction
        pv_group_sex = Anova(model,type = "III")$"Pr(>F)"[xvars=="group:sex"]
        
      }
      else {
        #define model formulation for models with the global variable 
        if (m[j] == "area"){
          glob_var = "total_area"
          f2 = paste(model_yvars[k],"~","+","group*sex","+","age","+","site","+","TotalEulerNumber","+",glob_var)
        }
        else if (m[j] == "thickness"){
          glob_var = "mean_thickness"
          f2 = paste(model_yvars[k],"~","+","group*sex","+","age","+","site","+","TotalEulerNumber","+",glob_var)
        }
        else if (m[j] == "volume"){
          glob_var = "BrainTotalVol"
          f2 = paste(model_yvars[k],"~","+","group*sex","+","age","+","site","+","TotalEulerNumber","+",glob_var)
        }
        #run model and define p-value for interaction
        model = lm(f2,data=datab)
        xvars = attributes(Anova(model,type = "III"))$row.names
        pv_group_sex = Anova(model,type = "III")$"Pr(>F)"[xvars=="group:sex"]
      }
    

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
          rw = list(model_yvars[k], group_F, group_pv, 
                    sex_F, sex_pv, age_F, age_pv, 
                    site_F, site_pv, 
                    EulerNumber_F, EulerNumber_pv, 
                    glob_var_F, glob_var_pv, glob_var,
                    GS_F, pv_group_sex,
                    sign_GS, mi, g-1)
        }       
        else{
          rw = list(model_yvars[k], group_F, group_pv, 
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


names(DF)<-c("model_yvar","group_Fval","group_pval",
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



#compute LSmean contrast for each region thickness with lsmean sex
#make model contrast table to be saved as excel
# only for thickness!!!

###### make inference with models
DF_thickness = data.frame()
DF_thickness_xlsx = data.frame()


for (j in seq(2,2)){
  print(m[j])
  model_yvars = unlist(col_names_list[[m[j]]] )
  for (k in seq(1,length(model_yvars))){
    for (g in seq(1,2)){

      if (glob[g] == "Without global var"){
        #no global measure model
        f = paste(model_yvars[k],"~","group+sex","+","age","+","site","+","TotalEulerNumber")
        model = lm(f,data=datab)
        xvars = attributes(Anova(model,type = "III"))$row.names
        
      }
      else {
        #global measure model
        if (m[j] == "thickness"){
          f2 = paste(model_yvars[k],"~","group+sex","+","age","+","site","+","TotalEulerNumber","+","mean_thickness")
        }
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
      
      
      #percent difference
      BP_diff_P = 100*BP_diff /K_emm
      BP_diff_LCL_P = 100*BP_diff_LCL /K_emm
      BP_diff_UCL_P = 100*BP_diff_UCL /K_emm
      
      SZ_diff_P = 100*SZ_diff /K_emm
      SZ_diff_LCL_P = 100*SZ_diff_LCL /K_emm
      SZ_diff_UCL_P = 100*SZ_diff_UCL /K_emm
      
      
      #for plotting 
      rw = list(model_yvars[k], K_emm,
                BP_diff_P, BP_diff_pv,
                BP_diff_LCL_P, BP_diff_UCL_P, 
                SZ_diff_P, SZ_diff_pv,
                SZ_diff_LCL_P, SZ_diff_UCL_P,
                g-1)
      
      DF_thickness = rbindlist(list(DF_thickness, rw))
      
      
      #for contrast xlsx
      rwc_xlsx = list(model_yvars[k],
                      BP_diff,BP_diff_tratio,BP_diff_pv,
                      SZ_diff,SZ_diff_tratio,SZ_diff_pv,
                      K_emm,g-1)
      
      DF_thickness_xlsx = rbindlist(list(DF_thickness_xlsx, rwc_xlsx))
      
      
    } #g
    
  } #k
} #j

names(DF_thickness)<-c("model_yvar","K_emm",
                       "BP_diff_P", "BP_diff_pv",
                       "BP_diff_LCL_P", "BP_diff_UCL_P", 
                       "SZ_diff_P", "SZ_diff_pv",
                       "SZ_diff_LCL_P", "SZ_diff_UCL_P",
                       "global_var_in_model")

names(DF_thickness_xlsx) <- c("Model_yvar",
                            "Contrast_BP-K","tratio_BP-K","pval_BP-K",
                            "Contrast_SZ-K","tratio_SZ-K","pval_SZ-K",
                            "K_LSmean","global_var_in_model")

DF_xlsx_glob0 = DF_thickness_xlsx[DF_thickness_xlsx$global_var_in_model == 0, ]
DF_xlsx_glob1 = DF_thickness_xlsx[DF_thickness_xlsx$global_var_in_model == 1, ]

write_xlsx(DF_xlsx_glob1,contrast_with_glob)
write_xlsx(DF_xlsx_glob0,contrast_without_glob)



#thickness plot

pivot_cols = c("BP_diff_P","SZ_diff_P")
DFp = DF_thickness %>%
  pivot_longer(cols = pivot_cols,
               names_to = "diff_group",
               values_to = "diff_value_P",
               values_drop_na = FALSE)
DFp = as.data.frame(DFp)
DFp$diff_group[DFp$diff_group == "BP_diff_P"] = "BP"
DFp$diff_group[DFp$diff_group == "SZ_diff_P"] = "SZ"


pivot_cols_LCL = c("BP_diff_LCL_P","SZ_diff_LCL_P")
DF_LCL = DF_thickness %>%
  pivot_longer(cols = pivot_cols_LCL,
               names_to = "diff_group_LCL",
               values_to = "diff_LCL_P",
               values_drop_na = FALSE)

pivot_cols_UCL = c("BP_diff_UCL_P","SZ_diff_UCL_P")
DF_UCL = DF_thickness %>%
  pivot_longer(cols = pivot_cols_UCL,
               names_to = "diff_group_UCL",
               values_to = "diff_UCL_P",
               values_drop_na = FALSE)
DF_CL = data.frame(DF_LCL$diff_group_LCL, DF_LCL$diff_LCL_P, DF_UCL$diff_UCL_P)
names(DF_CL) = c("diff_group_CL","diff_LCL_P","diff_UCL_P")
DFp = cbind(DFp,DF_CL)


#make a new coloumn in the dataframe with a non zero value when the model
#contrast is significant. The value given is the y-axis position of a star being 
#drawn indicating significant contrast in the plot below
DFp$BP_diff_sig = NA
DFp$BP_diff_sig[DFp$BP_diff_pv < 0.05] = -4
DFp$SZ_diff_sig = NA
DFp$SZ_diff_sig[DFp$SZ_diff_pv < 0.05] = -5


#define relevant variables for naming
groups = c("BP","SZ")
m = c("area","thickness","volume")
glob = c("Models WITHOUT global var","Models WITH global var")
sx = c("female","male")

sp=list()

for (g in seq(1,2)){  
  
  dfs = DFp[c(DFp$global_var_in_model == (g-1)),]
  
  sp[[g]]=ggplot(dfs, aes_string(group="diff_group",color="diff_group", x="model_yvar", y="diff_value_P")) +
    labs(x="region",y="Difference from control [%]") +
    geom_line() + 
    ggtitle(paste(glob[g])) + 
    geom_hline(yintercept = 0) + 
    geom_errorbar(aes_string(width=0.3,group="diff_group_CL",color="diff_group_CL",ymin="diff_LCL_P", ymax="diff_UCL_P" )) +
    scale_color_manual("Group", values=c("#0072B2","#0072B2","#009E73", "#009E73"))+
    
    geom_point(aes_string(group="diff_group",color="diff_group",x="model_yvar", y="diff_value_P")) +
    geom_point(shape=8, aes_string(group="diff_group",x="model_yvar", y="BP_diff_sig"), color="#0072B2") +
    geom_point(shape=8, aes_string(group="diff_group",x="model_yvar", y="SZ_diff_sig"), color="#009E73") +
    
    theme(axis.text.x = element_text(color = "black", size = 8, angle = 90)) + 
    ylim(-10, 10) +
    {if (m[j]=="thickness") ylim(-5, 5)}
  #coord_flip()
  
}
top_title = paste("Bilateral LSmean difference from control from both sex: thickness")
ps=grid.arrange(grobs=sp, top=textGrob(top_title,gp=gpar(fontsize=20)))
ggsave(paste(LSmeans_prefix,".png",sep=""),ps,width = 10,height = 10)







##### Test individual region statistics, FOR SANITY CHECKS ##### 



model_IPvol = lm(bi_inferiorparietal_volume ~ group*sex + age + TotalEulerNumber + site, data=datab)
#model_IPvol = update(model_IPvol,~.-group:sex)
lsmeans(model_IPvol,pairwise~"group",by="sex",adjust="none")


model1 = lm(bi_superiorparietal_volume ~ group*sex + age + site + TotalEulerNumber + BrainTotalVol, data=datab)
model1 = lm(bi_superiorparietal_volume ~ group+sex + age + site + TotalEulerNumber + BrainTotalVol, data=datab)
Anova(model1,type="III")
model2 = lm(bi_superiorparietal_volume ~ group*sex + age + site + TotalEulerNumber, data=datab)
Anova(model2,type="III")


model1 = lm(bi_medialorbitofrontal_area ~ group*sex + age + TotalEulerNumber + total_area + site, data=datab)
model1 = lm(bi_medialorbitofrontal_area ~ group+sex + age + TotalEulerNumber + total_area + site, data=datab)
Anova(model1,type="III")

model1 = lm(bi_medialorbitofrontal_area ~ group*sex + age + TotalEulerNumber + site, data=datab)
model1 = lm(bi_medialorbitofrontal_area ~ group+sex + age + TotalEulerNumber + site, data=datab)
Anova(model1,type="III")

model1 = lm(bi_fusiform_area ~ group*sex + age + TotalEulerNumber + site, data=datab)
model1 = lm(bi_fusiform_area ~ group+sex + age + TotalEulerNumber + site, data=datab)
Anova(model1,type="III")


model1 = lm(bi_medialorbitofrontal_thickness ~ group*sex + age + TotalEulerNumber + mean_thickness + site, data=datab)
Anova(model1,type="III")
model1 = lm(bi_medialorbitofrontal_thickness ~ group+sex + age + TotalEulerNumber + mean_thickness + site, data=datab)
Anova(model1,type="III")
lsmeans(model1,pairwise~"group",adjust="none")

model1 = lm(bi_medialorbitofrontal_thickness ~ group*sex + age + TotalEulerNumber + site, data=datab)
Anova(model1,type="III")
model1 = lm(bi_medialorbitofrontal_thickness ~ group+sex + age + TotalEulerNumber + site, data=datab)
Anova(model1,type="III")
lsmeans(model1,pairwise~"group",adjust="none")


model1 = lm(bi_rostralanteriorcingulate_thickness ~ group*sex + age + site + TotalEulerNumber + mean_thickness, data=datab)
Anova(model1,type="III")
model2 = lm(bi_rostralanteriorcingulate_thickness ~ group + sex + age + site + TotalEulerNumber  + mean_thickness, data=datab)
Anova(model2,type="III")
lsmeans(model1,pairwise~"group",adjust="none")


