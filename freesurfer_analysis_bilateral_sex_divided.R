
setwd("/mnt/projects/VIA11/FREESURFER/Stats/Data")

#load packages
library(data.table)
library(ggplot2)
library(ggdist)
library(gridExtra)
library(tidyr)
library(lsmeans)
library(grid)
library(car)
library(writexl)


#load data
data_csv <- read.table("VIA11_allkey_160621_FreeSurfer_pruned_20220126.csv", header = TRUE, sep = ",", dec = ".")


#filter the data with include variable
data_csv_filtered <- data_csv[c(data_csv$Include_FS_studies == 1),]
data_csv_filtered <- data_csv_filtered[!is.na(data_csv_filtered$Include_FS_studies),]

data_csv_filtered <- data_csv[c(data_csv$Include_FS_studies_euler_outliers_excluded == 1),]
data_csv_filtered <- data_csv_filtered[!is.na(data_csv_filtered$Include_FS_studies_euler_outliers_excluded),]

#tell r which variables are factors
datab = data_csv_filtered


setwd("/mnt/projects/VIA11/FREESURFER/Stats/Plots/bilateral")


######################################
#    run models on each region
######################################

# extract model yvars 

#make the bilateral coloumns
m = c("area","thickness","volume")
col_names = names(datab)
regions = col_names[grepl("^lh_", col_names)]
regions <- unique(unlist(strsplit(regions,"_")))
col2rm = c(m,"WhiteSurfArea","lh","MeanThickness")
regions = regions[!(regions %in% col2rm)] 
col_names_list = list()

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

datab$sex = as.factor(datab$Sex_child)
datab$group = as.factor(datab$HighRiskStatus)
datab$site = as.factor(datab$MR_Site)
datab$diag = as.factor(datab$Axis.1_diag_v11)
datab$age = as.numeric(datab$MRI_age)


#get a separate data frame for each sex
data_sex1 = datab[c(datab$sex == 1),]
data_sex0 = datab[c(datab$sex == 0),]
dataf = list(data_sex0,data_sex1)



###### make inference with models
models = list()
models_glob = list()

DF = data.frame()
DF_xlsx = data.frame()
DFc_xlsx = data.frame()

for (j in seq(1,3)){
  print(m[j])
  model_yvars = unlist(col_names_list[[m[j]]] )
  for (k in seq(1,length(model_yvars))){
    for (s in seq(1,2)){
      for (g in seq(0,1)){
      
      if (g == 0){
        #no global measure model
        glob_var = NA
        f = paste(model_yvars[k],"~","+","group","+","age","+","site","+","TotalEulerNumber")
        models[[m[j]]][[model_yvars[k]]] = lm(f,data=dataf[[s]])
        
        model_ana = models[[m[j]]][[model_yvars[k]]]
        xvars = attributes(Anova(model_ana,type = "III"))$row.names
        
      }
      else {
        #global measure model
        if (m[j] == "area"){
          glob_var = "total_area"
          f2 = paste(model_yvars[k],"~","group","+","age","+","site","+","TotalEulerNumber","+",glob_var)
        }
        else if (m[j] == "thickness"){
          glob_var = "mean_thickness"
          f2 = paste(model_yvars[k],"~","group","+","age","+","site","+","TotalEulerNumber","+",glob_var)
        }
        else if (m[j] == "volume"){
          glob_var = "BrainTotalVol"
          f2 = paste(model_yvars[k],"~","group","+","age","+","site","+","TotalEulerNumber","+",glob_var)
        }
        models_glob[[m[j]]][[model_yvars[k]]] = lm(f2,data=dataf[[s]])
        model_ana = models_glob[[m[j]]][[model_yvars[k]]]
        xvars = attributes(Anova(model_ana,type = "III"))$row.names
      }
      
      
      #use the above defined model_ana 
      
      ls = lsmeans(model_ana,pairwise~"group", adjust="none")
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
      
      mm = model_ana
      
      group_F = Anova(mm,type="III")$"F value"[xvars=="group"]
      age_F = Anova(mm,type = "III")$"F value"[xvars=="age"]
      site_F = Anova(mm,type = "III")$"F value"[xvars=="site"]
      EulerNumber_F = Anova(mm,type = "III")$"F value"[xvars=="TotalEulerNumber"]
      
      group_pv = Anova(mm,type = "III")$"Pr(>F)"[xvars=="group"]
      age_pv = Anova(mm,type = "III")$"Pr(>F)"[xvars=="age"]
      site_pv = Anova(mm,type = "III")$"Pr(>F)"[xvars=="site"]
      EulerNumber_pv = Anova(mm,type = "III")$"Pr(>F)"[xvars=="TotalEulerNumber"]
      
      if (g == 0){
        #no global measure model
        glob_var_F = NA
        glob_var_pv = NA
      }       
      else{
        #global measure model
        glob_var_F = Anova(mm,type = "III")$"F value"[xvars==glob_var]
        glob_var_pv = Anova(mm,type = "III")$"Pr(>F)"[xvars==glob_var]
        
      } 
      
      
      #rows for ANOVA xlsx table
      rw_xlsx = list(model_yvars[k], group_F, group_pv, 
                     age_F, age_pv, site_F, site_pv, 
                     EulerNumber_F, EulerNumber_pv, 
                     glob_var_F, glob_var_pv, glob_var,
                     g,s-1)
      
      DF_xlsx = rbindlist(list(DF_xlsx, rw_xlsx))
      
      
      #rows for contrast xlsx table
      rwc_xlsx = list(model_yvars[k],
                  BP_diff,BP_diff_tratio,BP_diff_pv,
                  SZ_diff,SZ_diff_tratio,SZ_diff_pv,
                  K_emm,g,s-1)
      
      DFc_xlsx = rbindlist(list(DFc_xlsx, rwc_xlsx))
      
      
      
      #rows for plotting
      rw = list(model_yvars[k], K_emm,
                BP_diff, BP_diff_pv, 
                BP_diff_LCL, BP_diff_UCL, 
                SZ_diff, SZ_diff_pv,
                SZ_diff_LCL, SZ_diff_UCL, 
                BP_diff_P, BP_diff_LCL_P, BP_diff_UCL_P, 
                SZ_diff_P, SZ_diff_LCL_P, SZ_diff_UCL_P, 
                g,(s-1))
      
      DF = rbindlist(list(DF, rw))
      
    } #g 
    } #s
  } #k
} #j


names(DF)<-c("model_yvar","K_emm",
             "BP_diff","BP_diff_pv",
             "BP_diff_LCL","BP_diff_UCL",
             "SZ_diff","SZ_diff_pv",
             "SZ_diff_LCL","SZ_diff_UCL",
             "BP_diff_P", "BP_diff_LCL_P", "BP_diff_UCL_P", 
             "SZ_diff_P", "SZ_diff_LCL_P", "SZ_diff_UCL_P", 
             "global_var_in_model","sex")

names(DF_xlsx)<-c("Model_yvar","Group_Fval","Group_pval",
                  "Age_Fval","Age_pval","Site_Fval","Site_pval",
                  "EulerNumber_Fval","Eulernumber_pval",
                  "global_var_F","global_var_pv","global_var_name",
                  "global_var_in_model","sex")

DF_xlsx_glob0 = DF_xlsx[DF_xlsx$global_var_in_model == 0, ]
DF_xlsx_glob1 = DF_xlsx[DF_xlsx$global_var_in_model == 1, ]

#df_with_ICV[,2:ncol(df_with_ICV)] <- signif(df_with_ICV[,2:ncol(df_with_ICV)],digits=3)
#df_without_ICV[,2:ncol(df_without_ICV)] <- signif(df_without_ICV[,2:ncol(df_without_ICV)],digits=3)

write_xlsx(DF_xlsx_glob1,"/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/parcels/Parcel_ANOVA_pvals_with_glob.xlsx")
write_xlsx(DF_xlsx_glob0,"/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/parcels/Parcel_ANOVA_pvals_without_glob.xlsx")


names(DFc_xlsx) = c("Model_yvar",
               "Contrast_BP-K","tratio_BP-K","pval_BP-K",
               "Contrast_SZ-K","tratio_SZ-K","pval_SZ-K",
               "K_LSmean","global_var_in_model","sex")

DFc_xlsx_glob0 = DFc_xlsx[DFc_xlsx$global_var_in_model == 0, ]
DFc_xlsx_glob1 = DFc_xlsx[DFc_xlsx$global_var_in_model == 1, ]

write_xlsx(DFc_xlsx_glob1,"/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/parcels/Parcel_model_contrast_with_glob.xlsx")
write_xlsx(DFc_xlsx_glob0,"/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/parcels/Parcel_model_contrast_without_glob.xlsx")





#change data frame for plotting

pivot_cols = c("BP_diff_P","SZ_diff_P")
DFp = DF %>%
  pivot_longer(cols = pivot_cols,
               names_to = "diff_group",
               values_to = "diff_value_P",
               values_drop_na = FALSE)
DFp = as.data.frame(DFp)
DFp$diff_group[DFp$diff_group == "BP_diff_P"] = "BP"
DFp$diff_group[DFp$diff_group == "SZ_diff_P"] = "SZ"

pivot_cols_LCL = c("BP_diff_LCL_P","SZ_diff_LCL_P")
DF_LCL = DF %>%
  pivot_longer(cols = pivot_cols_LCL,
               names_to = "diff_group_LCL",
               values_to = "diff_LCL_P",
               values_drop_na = FALSE)

pivot_cols_UCL = c("BP_diff_UCL_P","SZ_diff_UCL_P")
DF_UCL = DF %>%
  pivot_longer(cols = pivot_cols_UCL,
               names_to = "diff_group_UCL",
               values_to = "diff_UCL_P",
               values_drop_na = FALSE)
DF_CL = data.frame(DF_LCL$diff_group_LCL, DF_LCL$diff_LCL_P, DF_UCL$diff_UCL_P)
names(DF_CL) = c("diff_group_CL","diff_LCL_P","diff_UCL_P")
DFp = cbind(DFp,DF_CL)


DFp$BP_diff_sig = NA
DFp$BP_diff_sig[DFp$BP_diff_pv < 0.05] = -15
DFp$BP_diff_sig[DFp$BP_diff_pv < 0.05 & grepl(paste("_thickness",sep = ""), DFp$model_yvar)] = -6

DFp$SZ_diff_sig = NA
DFp$SZ_diff_sig[DFp$SZ_diff_pv < 0.05] = -17
DFp$SZ_diff_sig[DFp$SZ_diff_pv < 0.05 & grepl(paste("_thickness",sep = ""), DFp$model_yvar)] = -5



###### mean difference plot
groups = c("BP","SZ")
#P = c("_P","")
Plab = c("Difference from control [%]","Orignial units")
m = c("area","thickness","volume")
glob = c("Models WITHOUT global var","Models WITH global var")
sx = c("female","male")

sp=list()
pp = 1

for (j in seq(1,length(m))){
  
  data_DF = DFp[grepl(paste("_",m[j],sep = ""), DFp$model_yvar),]
  
  for (g in seq(1,2)){  
    
    dfs = data_DF[c(DFp$global_var_in_model == (g-1)),]
    
    for (s in seq(1,length(sx))){
      
      dfss = dfs[c(dfs$sex == s-1),]
      
      sp[[2*(g)-2+s]]=ggplot(dfss, aes_string(group="diff_group",color="diff_group", x="model_yvar", y="diff_value_P")) +
        labs(x="region",y=Plab[pp]) +
        geom_line() + 
        ggtitle(paste(glob[g],":",sx[s])) + 
        geom_hline(yintercept = 0) + 
        geom_errorbar(aes(width=0.3,group=diff_group_CL, color=diff_group_CL, ymin=diff_LCL_P, ymax=diff_UCL_P )) +
        scale_color_manual("Group", values=c("#0072B2","#0072B2","#009E73", "#009E73"))+
        
        geom_point(aes_string(group="diff_group",color="diff_group",x="model_yvar", y="diff_value_P")) +
        geom_point(shape=8, aes_string(group="diff_group",x="model_yvar", y="BP_diff_sig"), color="#0072B2") +
        geom_point(shape=8, aes_string(group="diff_group",x="model_yvar", y="SZ_diff_sig"),color="#009E73") +
        theme(axis.text.x = element_text(color = "#993333", size = 8, angle = 90)) +
        ylim(-17, 17) +
        annotate("text", x = 6, y = 13,label = "* = Significant contrast",family = "", fontface = 3, size=3)+
        {if (m[j]=="thickness") ylim(-6, 6)}
      {if (sx[s]=="male") theme(axis.text.x = element_text(color = "blue", size = 8, angle = 90))}
    } #end for s
    
  }
  top_title = paste("Bilateral LSmean difference from control: ",m[j])
  ps=grid.arrange(grobs=sp, top=textGrob(top_title,gp=gpar(fontsize=20)))
  ggsave(paste("LSmean_difference_",m[j],".png",sep=""),ps,width = 15,height = 10)
  
}





### SANITY CHECKS OF MODELS BEING RUN ####


#with glob var
model_ba_sex1 = lm(bi_bankssts_area ~ group + age + site + TotalEulerNumber + total_area, data=data_sex1)
model_ba_sex0 = lm(bi_bankssts_area ~ group + age + site + TotalEulerNumber + total_area, data=data_sex0)
Anova(model_ba_sex0,type="III")
Anova(model_ba_sex1,type="III")
lsmeans(model_ba_sex0,pairwise~"group",adjust="none")
lsmeans(model_ba_sex1,pairwise~"group",adjust="none")

model_IPvol_sex1 = lm(bi_frontalpole_volume ~ group + age + site + TotalEulerNumber + BrainTotalVol, data=data_sex1)
model_IPvol_sex0 = lm(bi_frontalpole_volume ~ group + age + site + TotalEulerNumber + BrainTotalVol, data=data_sex0)
Anova(model_IPvol_sex1,type="III")
Anova(model_IPvol_sex0,type="III")
lsmeans(model_IPvol_sex0,pairwise~"group",adjust="none")
lsmeans(model_IPvol_sex1,pairwise~"group",adjust="none")


#without
model_ba_sex1 = lm(bi_bankssts_area ~ group + age + site + TotalEulerNumber, data=data_sex1)
model_ba_sex0 = lm(bi_bankssts_area ~ group + age + site + TotalEulerNumber, data=data_sex0)
Anova(model_ba_sex0,type="III")
Anova(model_ba_sex1,type="III")
lsmeans(model_ba_sex0,pairwise~"group",adjust="none")
lsmeans(model_ba_sex1,pairwise~"group",adjust="none")

model_IPvol_sex1 = lm(bi_frontalpole_volume ~ group + age + site + TotalEulerNumber, data=data_sex1)
model_IPvol_sex0 = lm(bi_frontalpole_volume ~ group + age + site + TotalEulerNumber, data=data_sex0)
Anova(model_IPvol_sex1,type="III")
Anova(model_IPvol_sex0,type="III")
lsmeans(model_IPvol_sex0,pairwise~"group",adjust="none")
lsmeans(model_IPvol_sex1,pairwise~"group",adjust="none")


model1 = lm(bi_parsopercularis_area ~ group + age + TotalEulerNumber + total_area + site, data=data_sex1)
model0 = lm(bi_parsopercularis_area ~ group + age + TotalEulerNumber + total_area + site, data=data_sex0)
lsmeans(model1,pairwise~"group",adjust="none")
lsmeans(model0,pairwise~"group",adjust="none")


ls1 = lsmeans(model_IPvol_sex1,pairwise~"group",adjust="none")
ls0 = lsmeans(model_IPvol_sex0,pairwise~"group",adjust="none")
ls1
ls0
confint(ls1)
confint(ls0)


modelm = lm(bi_insula_volume ~ group*sex + age + site + TotalEulerNumber + BrainTotalVol, data = datab)
xvars = attributes(Anova(modelm,type="II"))$row.names
Anova(modelm,type="II")$"F value"[xvars=="group"]
xvars = attributes(Anova(modelm,type="III"))$row.names
Anova(modelm,type="III")$"F value"[xvars=="group"]




