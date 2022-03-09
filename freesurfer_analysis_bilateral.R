
setwd("/mnt/projects/VIA11/FREESURFER/Stats/Data")

#load packages
library(data.table)
library(ggplot2)
library(ggdist)
library(gridExtra)
library(lmerTest)
library(lme4)
library(tidyr)
library(dplyr)
library(stringr)
library(see)
library(lsmeans)
library(grid)
library(Cairo)
library(grDevices)
library(ggpcp)


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



##### Test individual region values
model_IPvol = lm(bi_inferiorparietal_volume ~ group*sex + age + TotalEulerNumber + site, data=datab)
anova(model_IPvol)
#model_IPvol = update(model_IPvol,~.-group:sex)
lsmeans(model_IPvol,pairwise~"group",by="sex",adjust="none")


model_IPvol = lm(bi_inferiorparietal_volume ~ group*sex + age + TotalEulerNumber + site, data=datab)
anova(model_IPvol)

lsmeans(model_IPvol,pairwise~"group",by="sex",adjust="none")

data_sex1 = datab[c(datab$sex == 1),]
data_sex0 = datab[c(datab$sex == 0),]
model_IPvol_sex1 = lm(bi_inferiorparietal_volume ~ group + age + TotalEulerNumber + site, data=data_sex1)
model_IPvol_sex0 = lm(bi_inferiorparietal_volume ~ group + age + TotalEulerNumber + site, data=data_sex0)
lsmeans(model_IPvol_sex1,pairwise~"group",adjust="none")
lsmeans(model_IPvol_sex0,pairwise~"group",adjust="none")



model_SPvol = lm(bi_superiorparietal_volume ~ group + sex + age + TotalEulerNumber + site, data=datab)
anova(model_SPvol)
lsmeans(model_SPvol,pairwise~"group",by="sex",adjust="none")




model_MOFs = lm(bi_medialorbitofrontal_thickness ~ group + sex + age + TotalEulerNumber + eICV_samseg + site, data=datab)
anova(model_MOFs)
lsmeans(model_MOFs,pairwise~"group",by="sex",adjust="none")

model_MOFs_sex0 = lm(bi_medialorbitofrontal_thickness ~ group + age + TotalEulerNumber + eICV_samseg + site, data=data_sex0)
model_MOFs_sex1 = lm(bi_medialorbitofrontal_thickness ~ group + age + TotalEulerNumber + eICV_samseg + site, data=data_sex1)
lsmeans(model_MOFs_sex0,pairwise~"group",adjust="none")
lsmeans(model_MOFs_sex1,pairwise~"group",adjust="none")

model_MOF = lm(bi_medialorbitofrontal_thickness ~ group + age + TotalEulerNumber + eICV_samseg + site, data=datab)
anova(model_MOF)
lsmeans(model_MOF,pairwise~"group",adjust="none")


###### make inference with models
models = list()
models_glob = list()

DF = data.frame()
DFs = data.frame()
DF_glob = data.frame()
DFs_glob = data.frame()


for (j in seq(1,3)){
  print(m[j])
  model_yvars = unlist(col_names_list[[m[j]]] )
  for (k in seq(1,length(model_yvars))){
    for (mi in seq(1,2)){
      
      if (mi == 1){
        #no global measure model
        f = paste(model_yvars[k],"~","+","group*sex","+","age","+","site","+","TotalEulerNumber")
        models[[m[j]]][[model_yvars[k]]] = lm(f,data=datab)
        
        model_ana = models[[m[j]]][[model_yvars[k]]]
        xvars = attributes(anova(model_ana))$row.names
        pv_group_sex = anova(model_ana)$"Pr(>F)"[xvars=="group:sex"]
        
        if (anova(model_ana)$"Pr(>F)"[xvars=="group:sex"] > 0.05){
          model_ana = update(model_ana,~.-group:sex)
        } 
        models[[m[j]]][[model_yvars[k]]] = model_ana
        
      }
      else {
        #global measure model
        if (m[j] == "area"){
          f2 = paste(model_yvars[k],"~","+","group*sex","+","age","+","site","+","TotalEulerNumber","+","total_area")
        }
        else if (m[j] == "thickness"){
          f2 = paste(model_yvars[k],"~","+","group*sex","+","age","+","site","+","TotalEulerNumber","+","mean_thickness")
        }
        else if (m[j] == "volume"){
          f2 = paste(model_yvars[k],"~","+","group*sex","+","age","+","site","+","TotalEulerNumber","+","BrainTotalVol")
        }
        models_glob[[m[j]]][[model_yvars[k]]] = lm(f2,data=datab)
        model_ana = models_glob[[m[j]]][[model_yvars[k]]]
        xvars = attributes(anova(model_ana))$row.names
        pv_group_sex = anova(model_ana)$"Pr(>F)"[xvars=="group:sex"]
          
        if (anova(model_ana)$"Pr(>F)"[xvars=="group:sex"] > 0.05){
          model_ana = update(model_ana,~.-group:sex)
        } 
        models_glob[[m[j]]][[model_yvars[k]]] = model_ana
      }
      
      
      #use the above defined model_ana 
      
      ls = lsmeans(model_ana,pairwise~"group", by = "sex", adjust="none")
      c = ls$contrasts
      
      K_emm = summary(ls)$lsmeans$lsmean[summary(ls$lsmeans)$sex == 0][2]
      sex1_K_emm = summary(ls)$lsmeans$lsmean[summary(ls$lsmeans)$sex == 1][2]
      
      #raw contrasts
      BP_diff = summary(c)$estimate[summary(c)$sex == 0][1]
      BP_diff_LCL = confint(c)$lower.CL[confint(c)$sex == 0][1]
      BP_diff_UCL = confint(c)$upper.CL[confint(c)$sex == 0][1]
      
      SZ_diff = -summary(c)$estimate[summary(c)$sex == 0][3]
      SZ_diff_LCL = -confint(c)$lower.CL[confint(c)$sex == 0][3]
      SZ_diff_UCL = -confint(c)$upper.CL[confint(c)$sex == 0][3]
      
      sex1_BP_diff = summary(c)$estimate[summary(c)$sex == 1][1]
      sex1_BP_diff_LCL = confint(c)$lower.CL[confint(c)$sex == 1][1]
      sex1_BP_diff_UCL = confint(c)$upper.CL[confint(c)$sex == 1][1]
      
      sex1_SZ_diff = -summary(c)$estimate[summary(c)$sex == 1][3]
      sex1_SZ_diff_LCL = -confint(c)$lower.CL[confint(c)$sex == 1][3]
      sex1_SZ_diff_UCL = -confint(c)$upper.CL[confint(c)$sex == 1][3]
      
      
      #percent difference
      BP_diff_P = 100*summary(c)$estimate[summary(c)$sex == 0][1] /K_emm
      BP_diff_LCL_P = 100*confint(c)$lower.CL[confint(c)$sex == 0][1] /K_emm
      BP_diff_UCL_P = 100*confint(c)$upper.CL[confint(c)$sex == 0][1] /K_emm
      
      SZ_diff_P = -100*summary(c)$estimate[summary(c)$sex == 0][3] /K_emm
      SZ_diff_LCL_P = -100*confint(c)$lower.CL[confint(c)$sex == 0][3] /K_emm
      SZ_diff_UCL_P = -100*confint(c)$upper.CL[confint(c)$sex == 0][3] /K_emm
      
      sex1_BP_diff_P = 100*summary(c)$estimate[summary(c)$sex == 1][1] /sex1_K_emm
      sex1_BP_diff_LCL_P = 100*confint(c)$lower.CL[confint(c)$sex == 1][1] /sex1_K_emm
      sex1_BP_diff_UCL_P = 100*confint(c)$upper.CL[confint(c)$sex == 1][1] /sex1_K_emm
      
      sex1_SZ_diff_P = -100*summary(c)$estimate[summary(c)$sex == 1][3] /sex1_K_emm
      sex1_SZ_diff_LCL_P = -100*confint(c)$lower.CL[confint(c)$sex == 1][3] /sex1_K_emm
      sex1_SZ_diff_UCL_P = -100*confint(c)$upper.CL[confint(c)$sex == 1][3] /sex1_K_emm
      
      
      pv_group = anova(model_ana)$"Pr(>F)"[1]
      
      mm = model_ana
      
      group_F = anova(mm)$"F value"[xvars=="group"]
      sex_F = anova(mm)$"F value"[xvars=="sex"]
      age_F = anova(mm)$"F value"[xvars=="age"]
      site_F = anova(mm)$"F value"[xvars=="site"]
      EulerNumber_F = anova(mm)$"F value"[xvars=="TotalEulerNumber"]
      
      group_pv = anova(mm)$"Pr(>F)"[xvars=="group"]
      sex_pv = anova(mm)$"Pr(>F)"[xvars=="sex"]
      age_pv = anova(mm)$"Pr(>F)"[xvars=="age"]
      site_pv = anova(mm)$"Pr(>F)"[xvars=="site"]
      EulerNumber_pv = anova(mm)$"Pr(>F)"[xvars=="TotalEulerNumber"]
      
      
      #rows for xlsx table
      rw = list(model_yvars[k],group_pv,pv_group_sex, K_emm, sex1_K_emm,
                BP_diff, BP_diff_LCL, BP_diff_UCL, 
                SZ_diff, SZ_diff_LCL, SZ_diff_UCL, 
                sex1_BP_diff, sex1_BP_diff_LCL, sex1_BP_diff_UCL, 
                sex1_SZ_diff, sex1_SZ_diff_LCL, sex1_SZ_diff_UCL,
                BP_diff_P, BP_diff_LCL_P, BP_diff_UCL_P, 
                SZ_diff_P, SZ_diff_LCL_P, SZ_diff_UCL_P, 
                sex1_BP_diff_P, sex1_BP_diff_LCL_P, sex1_BP_diff_UCL_P, 
                sex1_SZ_diff_P, sex1_SZ_diff_LCL_P, sex1_SZ_diff_UCL_P,mi)
      
      
      #rows for data plotting data frame
      rws0 = list(model_yvars[k], pv_group_sex, K_emm, sex1_K_emm, 
                  BP_diff, BP_diff_LCL, BP_diff_UCL, 
                  SZ_diff, SZ_diff_LCL, SZ_diff_UCL,
                  BP_diff_P, BP_diff_LCL_P, BP_diff_UCL_P, 
                  SZ_diff_P, SZ_diff_LCL_P, SZ_diff_UCL_P,0)
      
      rws1 = list(model_yvars[k], pv_group_sex, K_emm, sex1_K_emm, 
                  sex1_BP_diff, sex1_BP_diff_LCL, sex1_BP_diff_UCL,
                  sex1_SZ_diff, sex1_SZ_diff_LCL, sex1_SZ_diff_UCL,
                  sex1_BP_diff_P, sex1_BP_diff_LCL_P, sex1_BP_diff_UCL_P,
                  sex1_SZ_diff_P, sex1_SZ_diff_LCL_P, sex1_SZ_diff_UCL_P,1)
      
      
      if (mi == 1){
        DF = rbindlist(list(DF, rw))
        DFs = rbindlist(list(DFs, rws0))
        DFs = rbindlist(list(DFs, rws1))
      }
      else{
        model_LRT_pv = anova(models_glob[[m[j]]][[model_yvars[k]]],models[[m[j]]][[model_yvars[k]]])$"Pr(>Chisq)"[2]
        rw = append(rw,model_LRT_pv)
        
        DF_glob = rbindlist(list(DF_glob, rw))
        DFs_glob = rbindlist(list(DFs_glob, rws0))
        DFs_glob = rbindlist(list(DFs_glob, rws1))
      }
      
    } #mi 
    
  } #k
} #j


DF = data.frame(DF)
DFs = data.frame(DFs)
DF_glob = data.frame(DF_glob)
DFs_glob = data.frame(DFs_glob)

names(DF)<-c("model_yvar","Group_p_value","Group_sex_p_value","K_emm", "sex1_K_emm",
             "BP_diff","BP_diff_LCL","BP_diff_UCL",
             "SZ_diff","SZ_diff_LCL","SZ_diff_UCL",
             "sex1_BP_diff","sex1_BP_diff_LCL","sex1_BP_diff_UCL",
             "sex1_SZ_diff","sex1_SZ_diff_LCL","sex1_SZ_diff_UCL",
             "BP_diff_P", "BP_diff_LCL_P", "BP_diff_UCL_P", 
             "SZ_diff_P", "SZ_diff_LCL_P", "SZ_diff_UCL_P", 
             "sex1_BP_diff_P", "sex1_BP_diff_LCL_P", "sex1_BP_diff_UCL_P", 
             "sex1_SZ_diff_P", "sex1_SZ_diff_LCL_P", "sex1_SZ_diff_UCL_P","model_type")
names(DFs)<-c("model_yvar","Group_sex_p_value","K_emm", "sex1_K_emm",
              "BP_diff","BP_diff_LCL","BP_diff_UCL",
              "SZ_diff","SZ_diff_LCL","SZ_diff_UCL",
              "BP_diff_P", "BP_diff_LCL_P", "BP_diff_UCL_P",
              "SZ_diff_P", "SZ_diff_LCL_P", "SZ_diff_UCL_P", "sex")
names(DF_glob)<-c("model_yvar","Group_p_value","Group_sex_p_value","K_emm", "sex1_K_emm",
                  "BP_diff","BP_diff_LCL","BP_diff_UCL",
                  "SZ_diff","SZ_diff_LCL","SZ_diff_UCL",
                  "sex1_BP_diff","sex1_BP_diff_LCL","sex1_BP_diff_UCL",
                  "sex1_SZ_diff","sex1_SZ_diff_LCL","sex1_SZ_diff_UCL",
                  "BP_diff_P", "BP_diff_LCL_P", "BP_diff_UCL_P", 
                  "SZ_diff_P", "SZ_diff_LCL_P", "SZ_diff_UCL_P", 
                  "sex1_BP_diff_P", "sex1_BP_diff_LCL_P", "sex1_BP_diff_UCL_P", 
                  "sex1_SZ_diff_P", "sex1_SZ_diff_LCL_P", "sex1_SZ_diff_UCL_P","model_type")
names(DFs_glob)<-c("model_yvar","Group_sex_p_value","K_emm", "sex1_K_emm",
                   "BP_diff","BP_diff_LCL","BP_diff_UCL",
                   "SZ_diff","SZ_diff_LCL","SZ_diff_UCL",
                   "BP_diff_P","BP_diff_LCL_P","BP_diff_UCL_P",
                   "SZ_diff_P","SZ_diff_LCL_P","SZ_diff_UCL_P","sex")

DFs$sex = as.factor(DFs$sex)
DFs_glob$sex = as.factor(DFs_glob$sex)




#change data frame for plotting

pivot_cols = c("BP_diff_P","SZ_diff_P")
DFp = DFs %>%
  pivot_longer(cols = pivot_cols,
               names_to = "diff_group",
               values_to = "diff_value_P",
               values_drop_na = FALSE)
DFp = as.data.frame(DFp)
DFp$diff_group[DFp$diff_group == "BP_diff_P"] = "BP"
DFp$diff_group[DFp$diff_group == "SZ_diff_P"] = "SZ"

DFp_glob = DFs_glob %>%
  pivot_longer(cols = pivot_cols,
               names_to = "diff_group",
               values_to = "diff_value_P",
               values_drop_na = FALSE)

DFp_glob = as.data.frame(DFp_glob)
DFp_glob$diff_group[DFp_glob$diff_group == "BP_diff_P"] = "BP"
DFp_glob$diff_group[DFp_glob$diff_group == "SZ_diff_P"] = "SZ"


#for non-glob
pivot_cols_LCL = c("BP_diff_LCL_P","SZ_diff_LCL_P")
DF_LCL = DFs %>%
  pivot_longer(cols = pivot_cols_LCL,
               names_to = "diff_group_LCL",
               values_to = "diff_LCL_P",
               values_drop_na = FALSE)

pivot_cols_UCL = c("BP_diff_UCL_P","SZ_diff_UCL_P")
DF_UCL = DFs %>%
  pivot_longer(cols = pivot_cols_UCL,
               names_to = "diff_group_UCL",
               values_to = "diff_UCL_P",
               values_drop_na = FALSE)
DF_CL = data.frame(DF_LCL$diff_group_LCL, DF_LCL$diff_LCL_P, DF_UCL$diff_UCL_P)
names(DF_CL) = c("diff_group_CL","diff_LCL_P","diff_UCL_P")
DFp = cbind(DFp,DF_CL)


#for glob
pivot_cols_LCL = c("BP_diff_LCL_P","SZ_diff_LCL_P")
DF_LCL = DFs_glob %>%
  pivot_longer(cols = pivot_cols_LCL,
               names_to = "diff_group_LCL",
               values_to = "diff_LCL_P",
               values_drop_na = FALSE)

pivot_cols_UCL = c("BP_diff_UCL_P","SZ_diff_UCL_P")
DF_UCL = DFs_glob %>%
  pivot_longer(cols = pivot_cols_UCL,
               names_to = "diff_group_UCL",
               values_to = "diff_UCL_P",
               values_drop_na = FALSE)
DF_CL = data.frame(DF_LCL$diff_group_LCL, DF_LCL$diff_LCL_P, DF_UCL$diff_UCL_P)
names(DF_CL) = c("diff_group_CL","diff_LCL_P","diff_UCL_P")
DFp_glob = cbind(DFp_glob,DF_CL)


DFp$Group_sex_sig = NA
DFp_glob$Group_sex_sig = NA

DFp$Group_sex_sig[DFp$Group_sex_p_value < 0.05] = -15
DFp_glob$Group_sex_sig[DFp_glob$Group_sex_p_value < 0.05] = -15



###### mean difference plot
groups = c("BP","SZ")
P = c("_P","")
Plab = c("Difference from control [%]","Orignial units")
m = c("area","thickness","volume")
glob = c("Models WITHOUT global var","Models WITH global var")
sx = c("female","male")

sp=list()

pp = 1
  
for (j in seq(1,length(m))){
  
  dfs = DFp[grepl(paste("_",m[j],sep = ""), DFp$model_yvar),]
  dfs_glob = DFp_glob[grepl(paste("_",m[j],sep = ""), DFp$model_yvar),]
  data_DF = list(dfs,dfs_glob)
  
  for (g in seq(1,2)){  
    
    dfs = data_DF[[g]]
    
    for (s in seq(1,length(sx))){
      
      dfss = dfs[c(dfs$sex == s-1),]
      
      sp[[2*g-2+s]]=ggplot(dfss, aes_string(group="diff_group",color="diff_group", x="model_yvar", y="diff_value_P")) +
        labs(x="region",y=Plab[pp]) +
        geom_line() + 
        ggtitle(paste(glob[g],":",sx[s])) + 
        geom_hline(yintercept = 0) + 
        geom_errorbar(aes(width=0.3,group=diff_group_CL, color=diff_group_CL, ymin=diff_LCL_P, ymax=diff_UCL_P )) +
        scale_color_manual("Group", values=c("#0072B2","#0072B2","#009E73", "#009E73"))+
        
        geom_point(aes_string(group="diff_group",color="diff_group",x="model_yvar", y="diff_value_P")) +
        geom_point(shape=8, aes_string(group="diff_group",x="model_yvar", y="Group_sex_sig")) +
        theme(axis.text.x = element_text(color = "#993333", size = 8, angle = 90)) +
        ylim(-15, 15) +
        annotate("text", x = 6, y = 13,label = "* = Significant G/S interaction",family = "", fontface = 3, size=3)+
        {if (m[j]=="thickness") ylim(-3, 3)}
        {if (sx[s]=="male") theme(axis.text.x = element_text(color = "blue", size = 8, angle = 90))}
    } #end for s

  }
  top_title = paste("Bilateral LSmean difference from control: ",m[j])
  ps=grid.arrange(grobs=sp, top=textGrob(top_title,gp=gpar(fontsize=20)))
  ggsave(paste("LSmean_difference_",m[j],".png",sep=""),ps,width = 15,height = 10)
  
}





#mean difference for thickness without sex separation
# only for thickness!!!

###### make inference with models
models = list()
models_glob = list()

DF_thickness = data.frame()
DF_glob_thickness = data.frame()

for (j in seq(2,2)){
  print(m[j])
  model_yvars = unlist(col_names_list[[m[j]]] )
  for (k in seq(1,length(model_yvars))){
    for (mi in seq(1,2)){
      
      if (mi == 1){
        #no global measure model
        f = paste(model_yvars[k],"~","+","group*sex","+","age","+","site","+","TotalEulerNumber")
        models[[m[j]]][[model_yvars[k]]] = lm(f,data=datab)
        
        model_ana = models[[m[j]]][[model_yvars[k]]]
        xvars = attributes(anova(model_ana))$row.names
        pv_group_sex = anova(model_ana)$"Pr(>F)"[xvars=="group:sex"]
        
        if (anova(model_ana)$"Pr(>F)"[xvars=="group:sex"] > 0.05){
          model_ana = update(model_ana,~.-group:sex)
        } 
        models[[m[j]]][[model_yvars[k]]] = model_ana
        
      }
      else {
        #global measure model
        if (m[j] == "thickness"){
          f2 = paste(model_yvars[k],"~","+","group*sex","+","age","+","site","+","TotalEulerNumber","+","mean_thickness")
        }

        models_glob[[m[j]]][[model_yvars[k]]] = lm(f2,data=datab)
        model_ana = models_glob[[m[j]]][[model_yvars[k]]]
        xvars = attributes(anova(model_ana))$row.names
        pv_group_sex = anova(model_ana)$"Pr(>F)"[xvars=="group:sex"]
        
        if (anova(model_ana)$"Pr(>F)"[xvars=="group:sex"] > 0.05){
          model_ana = update(model_ana,~.-group:sex)
        } 
        models_glob[[m[j]]][[model_yvars[k]]] = model_ana
      }
      
      
      #use the above defined model_ana 
      ls = lsmeans(model_ana,pairwise~"group",adjust="none")
      c = ls$contrasts
      K_emm = summary(ls)$lsmeans$lsmean[2]

      #percent difference
      BP_diff_P = 100*summary(c)$estimate[1] /K_emm
      BP_diff_LCL_P = 100*confint(c)$lower.CL[1] /K_emm
      BP_diff_UCL_P = 100*confint(c)$upper.CL[1] /K_emm
      
      SZ_diff_P = -100*summary(c)$estimate[3] /K_emm
      SZ_diff_LCL_P = -100*confint(c)$lower.CL[3] /K_emm
      SZ_diff_UCL_P = -100*confint(c)$upper.CL[3] /K_emm
      
      pv_group = anova(model_ana)$"Pr(>F)"[1]
      
      rw = list(model_yvars[k],pv_group,pv_group_sex, K_emm,
                BP_diff_P, BP_diff_LCL_P, BP_diff_UCL_P, 
                SZ_diff_P, SZ_diff_LCL_P, SZ_diff_UCL_P,mi)
      
      if (mi == 1){
        DF_thickness = rbindlist(list(DF_thickness, rw))
      }
      else{
        model_LRT_pv = anova(models_glob[[m[j]]][[model_yvars[k]]],models[[m[j]]][[model_yvars[k]]])$"Pr(>Chisq)"[2]
        rw = append(rw,model_LRT_pv)
        
        DF_glob_thickness = rbindlist(list(DF_glob_thickness, rw))
      }
      
    } #mi 
    
  } #k
} #j

names(DF_thickness)<-c("model_yvar","Group_p_value","Group_sex_p_value","K_emm",
             "BP_diff_P", "BP_diff_LCL_P", "BP_diff_UCL_P", 
             "SZ_diff_P", "SZ_diff_LCL_P", "SZ_diff_UCL_P","model_type")
names(DF_glob_thickness)<-c("model_yvar","Group_p_value","Group_sex_p_value","K_emm",
                  "BP_diff_P", "BP_diff_LCL_P", "BP_diff_UCL_P", 
                  "SZ_diff_P", "SZ_diff_LCL_P", "SZ_diff_UCL_P", "model_type")


pivot_cols = c("BP_diff_P","SZ_diff_P")
DFp = DF_thickness %>%
  pivot_longer(cols = pivot_cols,
               names_to = "diff_group",
               values_to = "diff_value_P",
               values_drop_na = FALSE)
DFp = as.data.frame(DFp)
DFp$diff_group[DFp$diff_group == "BP_diff_P"] = "BP"
DFp$diff_group[DFp$diff_group == "SZ_diff_P"] = "SZ"

DFp_glob = DF_glob_thickness %>%
  pivot_longer(cols = pivot_cols,
               names_to = "diff_group",
               values_to = "diff_value_P",
               values_drop_na = FALSE)

DFp_glob = as.data.frame(DFp_glob)
DFp_glob$diff_group[DFp_glob$diff_group == "BP_diff_P"] = "BP"
DFp_glob$diff_group[DFp_glob$diff_group == "SZ_diff_P"] = "SZ"




sp=list()

pp = 1

data_DF = list(DFp,DFp_glob)

for (g in seq(1,2)){  
  
  dfs = data_DF[[g]]
    
    sp[[g]]=ggplot(dfs, aes_string(group="diff_group",color="diff_group", x="model_yvar", y="diff_value_P")) +
      labs(x="region",y=Plab[pp]) +
      geom_line() + 
      ggtitle(paste(glob[g])) + 
      geom_hline(yintercept = 0) + 
      #geom_errorbar(aes_string(width=0.3,group="diff_group_CL",color="diff_group_CL",ymin="diff_LCL_P", ymax="diff_UCL_P" )) +
      geom_point(aes_string(group="diff_group",color="diff_group",x="model_yvar", y="diff_value_P")) +
      theme(axis.text.x = element_text(color = "black", size = 8, angle = 90)) + 
      ylim(-10, 10) +
      {if (m[j]=="thickness") ylim(-3, 3)}
      #coord_flip()
  
}
top_title = paste("Bilateral LSmean difference from control: thickness")
ps=grid.arrange(grobs=sp, top=textGrob(top_title,gp=gpar(fontsize=20)))
ggsave(paste("LSmean_difference_thickness_gender_combined",".png",sep=""),ps,width = 10,height = 10)








################ MAYBE USE LATER ##############


#plot variance in each region for each measure

mm = "area"
mc = paste(col_names_list[[mm]])
dft = datab %>%
  pivot_longer(cols = mc,
               names_to = "NAME_region_measure_h",
               values_to = "VALUE_region_measure_h",
               values_drop_na = FALSE)
dft = as.data.frame(dft)

dft$NORM_region_measure_h = dft$VALUE_region_measure_h

for (i in seq(1,length(col_names_list[[mm]]))){
  avg = mean(dft$VALUE_region_measure_h[ dft$NAME_region_measure_h  == col_names_list[[mm]][i] ])
  dft$NORM_region_measure_h = ( dft$VALUE_region_measure_h[ dft$NAME_region_measure_h  == col_names_list[[mm]][i] ] - avg ) / avg
}

p_var = list()

#plot data
p_var[[1]]=with(dft,
                ggplot() +
                  aes_string(x = "NAME_region_measure_h", color = "group", group = "group", y = "VALUE_region_measure_h") +
                  geom_jitter(width = 0.4, size=0.1) + 
                  stat_summary(fun = mean, geom = "point",size=2) +
                  stat_summary(fun = mean, geom = "point",size=2,pch=21,colour="black") +
                  scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) + 
                  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
                  labs(y = paste("Region",m), x = "brain region") +
                  ggtitle(paste(mm))
)


p_var[[2]]= with(dft,
                 ggplot() +
                   aes(x = factor(NAME_region_measure_h), y = 100*NORM_region_measure_h, color = group) +
                   geom_violin(position = "identity",alpha=0.3) +
                   geom_jitter(width = 0.4, size=0.1) + 
                   geom_hline(yintercept = 0) + 
                   stat_summary(fun = mean, geom = "point",size=3,aes(colour = group)) +
                   #stat_summary(fun = mean, geom = "point",size=3,pch=21) +
                   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
                   labs(y = "Percent difference from mean [%]", x = "brain region") +
                   ggtitle(paste(mm))
)
ps=grid.arrange(grobs=p_var)
#ggsave(paste("original_unit_all_region_variance_",mm,".png",sep=""),ps,width = 14,height = 10)

ps=grid.arrange(p_var[[2]])
#ggsave(paste("all_region_variance_",mm,".png",sep=""),ps,width = 14,height = 6)


p_var2 = list()
half_names = c("bankssts","caudalanteriorcingulate","caudalmiddlefrontal","cuneus","entorhinal",
               "fusiform","inferiorparietal","inferiortemporal","isthmuscingulate","lateraloccipital",
               "lateralorbitofrontal","lingual","medialorbitofrontal","middletemporal","parahippocampal",
               "paracentral","parsopercularis")

half_names2 = c("parsorbitalis","parstriangularis","pericalcarine","postcentral","posteriorcingulate",
                "precentral","precuneus","rostralanteriorcingulate","rostralmiddlefrontal",
                "superiorfrontal","superiorparietal","superiortemporal","supramarginal","frontalpole",
                "temporalpole","transversetemporal","insula")

dft2 = dft[!grepl(half_names[1], dft$NAME_region_measure_h,fixed=TRUE),]
dft3 = dft[!grepl(half_names2[1], dft$NAME_region_measure_h,fixed=TRUE),]

for (i in seq(2,length(half_names))){
  dft2 = dft2[!grepl(half_names[i], dft2$NAME_region_measure_h,fixed=TRUE),]
  dft3 = dft3[!grepl(half_names2[i], dft3$NAME_region_measure_h,fixed=TRUE),]
}


#group sex interaction
DF_list = list()
DF_list[[1]] = DF
DF_list[[2]] = DF_glob

pgs = list()

for (g in seq(1,2)){
  df = DF_list[[g]]
  for (i in seq(1,3)){
      #df1 = df[grepl(paste("^",h[j],sep = ""), df$model_yvar),]
      df1 = df[grepl(paste("_",m[i],sep = ""), df$model_yvar),]
      
      pgs[[i*2+g-2]] =ggplot(df1, aes(as.factor( model_yvar ), as.numeric(Group_sex_p_value))) + 
        geom_point() + 
        geom_hline(yintercept = 0.05) + 
        labs(y = "G/S anova p-value", x = "brain region") +
        ggtitle(m[i]) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  }
}
top_title = "Model WITHOUT global var                 Model WITH global var"
ps=grid.arrange(grobs=pgs,top=textGrob(top_title,gp=gpar(fontsize=20)))
#ggsave(paste(glob[g],"_group:sex_interaction_pvalue",".png",sep=""),ps,width = 10,height = 12)

#sanity check above results with manual models 
model_test1 = lm(bi_parahippocampal_area ~ group*sex + age + site + TotalEulerNumber, data=datab)
DF$Group_sex_p_value[grepl(paste("bi_parahippocampal_area"), DF$model_yvar)]
anova(model_test1)

model_test2 = lm(bi_parahippocampal_area ~ group*sex + age + site + TotalEulerNumber + total_area, data=datab)
DF_glob$Group_sex_p_value[grepl(paste("bi_parahippocampal_area"), DF_glob$model_yvar)]
anova(model_test2)



###### mean difference plot
groups = c("BP","SZ")
P = c("_P","")
Plab = c("Difference from control [%]","Orignial units")
h = c("lh","rh")
m = c("area","thickness","volume")
glob = c("Model WITHOUT global var","Model WITH global var")

datafs = list(DFs,DFs_glob)
sp=list()

pp = 1

for (g in seq(1,2)){
  
    for (j in seq(1,3)){
      
      for (k in seq(1,2)){
        dfs = datafs[[k]]
        
        #df1 = dfs[grepl(paste("^",h[i],sep = ""), dfs$model_yvar),]
        df1 = dfs[grepl(paste("_",m[j],sep = ""), dfs$model_yvar),]
        
        nam = paste(groups[g],"_diff",P[pp],sep = "")
        df1[[nam]] = as.numeric(df1[[nam]])
        nam = paste(groups[g],"_diff_LCL",P[pp],sep = "")
        df1[[nam]] = as.numeric(df1[[nam]])
        nam = paste(groups[g],"_diff_UCL",P[pp],sep = "")
        df1[[nam]] = as.numeric(df1[[nam]])
        df1[["model_yvar"]] = as.factor(df1[["model_yvar"]])
        
        
        sp[[k]]=ggplot(df1, aes_string(group="sex",color="sex", x="model_yvar", y=paste(groups[g],"_diff",P[pp],sep = ""))) +
          labs(x="region",y=Plab[pp]) +
          see::geom_violinhalf(trim=FALSE) + #aes(group = sex)
          geom_line() + 
          ggtitle(paste(glob[k])) + 
          geom_hline(yintercept = 0) + 
          geom_errorbar(aes_string(width=0.3,group="sex",color="sex",ymin=paste(groups[g],"_diff_LCL",P[pp],sep = ""), ymax=paste(groups[g],"_diff_UCL",P[pp],sep = ""))) +
          geom_point(aes_string(group="sex",color="sex",x="model_yvar", y=paste(groups[g],"_diff",P[pp],sep = ""))) +
          theme(axis.text.x = element_text(face = "bold", color = "#993333", 
                                           size = 12, angle = 0),
                axis.text.y = element_text(face = "bold", color = "blue", 
                                           size = 10, angle = 0)) +
          ylim(-20, 20) +
          {if (k==1) theme(legend.position = "none")} +
          {if (k==2) theme(axis.text.y=element_blank())} +
          {if (k==2) xlab(NULL)} +
          coord_flip()
        
      }
      top_title = paste(m[j],"for group",groups[g])
      lay <- rbind(c(1,1,1,1,1,1,1,2,2,2,2),
                   c(1,1,1,1,1,1,1,2,2,2,2),
                   c(1,1,1,1,1,1,1,2,2,2,2))
      ps=grid.arrange(grobs=sp, layout_matrix=lay, top=textGrob(top_title,gp=gpar(fontsize=20)))
      #ggsave(paste(groups[g],"_",m[j],"_mean_difference",".png",sep=""),ps,width = 10,height = 12)
    }
}



#sanity check
model_test3 = lmer(lh_parahippocampal_volume ~ group + sex + age + BrainTotalVol + (1| site), data=datab)
lsmeans(model_test3,pairwise~"group")
model_bvol2 = lmer(mean_thickness ~ -1 + group + sex + age + (1| site), data=datab)
anova(model_bvol2)
lsmeans(model_bvol2,pairwise~"group",by="sex")
plot(lsmeans(model_bvol2,pairwise~"group",by="sex"),comparison=TRUE)







#############  MANOVA   ##############
#y_vars = col_names_list[[paste("lh","area",sep = "_")]] 

# Global measures # 
y_vars = c("BrainTotalVol", "CortexVol", "total_area", "mean_thickness")
y_vars_idx = (names(datab) %in% y_vars == TRUE)
my_vars <- data.matrix(datab[,c(y_vars_idx)])

mmodel = manova(my_vars ~ group*sex + age + site + TotalEulerNumber, data = datab)
anova(mmodel)
mmodel = manova(my_vars ~ group + sex + age + site + TotalEulerNumber, data = datab)
anova(mmodel)
summary.aov(mmodel)


y_vars = c("CortexVol", "total_area", "mean_thickness")
y_vars_idx = (names(datab) %in% y_vars == TRUE)
my_vars <- data.matrix(datab[,c(y_vars_idx)])

mmodel = manova(my_vars ~ group*sex + age + site + TotalEulerNumber, data = datab)
anova(mmodel)
mmodel = manova(my_vars ~ group + sex + age + site + TotalEulerNumber, data = datab)
anova(mmodel)
summary.aov(mmodel)



PH = lsmeans(mmodel,pairwise~"group", by = "sex")
plot(PH,comparison=TRUE,xlab="Multivariate lsmean")



# 34 regions for each measure
y_vars = col_names_list[[paste("lh","area",sep = "_")]] 
y_vars_idx = (names(datab) %in% y_vars == TRUE)
my_vars <- data.matrix(datab[,c(y_vars_idx)])

mmodel = manova(my_vars ~ group*sex + age + site + eICV_samseg, data = datab)
anova(mmodel)
mmodel = manova(my_vars ~ group + sex + age + site + eICV_samseg, data = datab)
anova(mmodel)
summary.aov(mmodel)




#test individual models 
model_vars = c("BrainTotalVol", "CortexVol", "total_area", "eICV_samseg")
model = list()
for (i in seq(1,length(model_vars))){
  f = paste(model_vars[i],"~","-1","+","group*sex","+","group:age","+","(1| site)")
  model[[i]] = lmer(f,data=datab)
}
model[[i+1]] = lmer(mean_thickness ~ -1 + group + sex + group:age + (1| site),data=datab)
model_vars = c(model_vars,"mean_thickness")

i = 1
PH = lsmeans(model[[i]],pairwise~"group", by = "sex")
plot(PH,comparison=TRUE,xlab=paste(model_vars[i],"lsmean"))


#sanity check
model_check = lm(total_area ~ group*sex + age + site, data = datab)
anova(model_check)

