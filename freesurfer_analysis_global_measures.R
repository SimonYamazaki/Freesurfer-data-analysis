
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
library(ggnewscale)
library(forcats)
library(writexl)
library(car)

#load data
data_csv <- read.table("VIA11_allkey_160621_FreeSurfer_pruned_20220126.csv", header = TRUE, sep = ",", dec = ".")

#inspect the head of data and summary
head(data_csv)
summary(data_csv)

#filter the data with include variable
data_csv_filtered <- data_csv[c(data_csv$Include_FS_studies == 1),]
data_csv_filtered <- data_csv_filtered[!is.na(data_csv_filtered$Include_FS_studies),]

data_csv_filtered <- data_csv[c(data_csv$Include_FS_studies_euler_outliers_excluded == 1),]
data_csv_filtered <- data_csv_filtered[!is.na(data_csv_filtered$Include_FS_studies_euler_outliers_excluded),]


#tell r which variables are factors
datab = data_csv_filtered
datab$sex = as.factor(datab$Sex_child)
datab$group = as.factor(datab$HighRiskStatus)
datab$site = as.factor(datab$MR_Site)
datab$diag = as.factor(datab$Axis.1_diag_v11)
datab$age = as.numeric(datab$MRI_age)


setwd("/mnt/projects/VIA11/FREESURFER/Stats/Plots/global_measures")

y_vars = c("BrainTotalVol", "CortexVol", "total_area", "mean_thickness", "eICV_samseg")


#plot variance in each region for each measure
dft = datab %>%
  pivot_longer(cols = y_vars,
               names_to = "NAME_measure",
               values_to = "VALUE_measure",
               values_drop_na = FALSE)
dft = as.data.frame(dft)

dft$NORM_measure = dft$VALUE_measure
dft$NORM_measure_K = NA #dftm$VALUE_measure
dft$NORM_measure_K_mt = NA #dftm$VALUE_measure


for (i in seq(1,length(y_vars))){
    avg = mean(dft$VALUE_measure[ dft$NAME_measure  == y_vars[i]])
    avg_K = mean(dft$VALUE_measure[ dft$NAME_measure  == y_vars[i] & dft$group  == "K"])
    
    dft$NORM_measure[ dft$NAME_measure  == y_vars[i] ] = ( dft$VALUE_measure[ dft$NAME_measure  == y_vars[i] ] - avg) / avg
    
    dft$NORM_measure_K[ dft$NAME_measure  == y_vars[i] ] = ( dft$VALUE_measure[ dft$NAME_measure  == y_vars[i] ] - avg_K ) / avg_K

}
  
p_var = list()

#original units
p_var[[1]]=with(dft,
                ggplot() +
                  aes_string(x = "NAME_measure", color = "group", group = "group", y = "VALUE_measure") +
                  geom_jitter(width = 0.3, size=0.1) + 
                  stat_summary(fun = mean, geom = "point",size=2) +
                  stat_summary(fun = mean, geom = "point",size=2,pch=21,colour="black") +
                  scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) + 
                  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
                  labs(y = paste("Region"), x = "brain region") +
                  ggtitle("Global measures")
)
grid.arrange(p_var[[1]])

#raw means
p_var[[2]]= with(dft,
                 ggplot() +
                   aes(x = factor(NAME_measure), y = 100*NORM_measure, color = group) +
                   geom_violin(position = "identity",alpha=0.3) +
                   geom_jitter(width = 0.3, size=0.1) + 
                   geom_hline(yintercept = 0) + 
                   #geom_boxplot(aes(colour = group),width = 0.1) +
                   stat_summary(fun = mean, geom = "point",size=3,aes(colour = group)) +
                   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
                   labs(y = "Percent difference from mean [%]", x = "brain region") +
                   ggtitle("Global measures")
)
grid.arrange(p_var[[2]])


## raw contrast means
p_var[[3]]= with(dft,
                 ggplot() +
                   aes(x = factor(NAME_measure), y = 100*NORM_measure_K, color = group) +
                   geom_violin(position = "identity",alpha=0.3) +
                   geom_jitter(width = 0.3, size=0.1) + 
                   geom_hline(yintercept = 0) + 
                   
                   ggnewscale::new_scale_colour() +
                   stat_summary(mapping = aes(x = factor(NAME_measure),y = 100*NORM_measure_K, group=group, colour=group),fun = "mean", geom = "point", size=2) +
                    
                   scale_colour_manual("Mean contrast", values = c("red", "green", "blue")) + 
                   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
                   labs(y = "Percent difference from control mean [%]", x = "brain region") +
                   ggtitle("Raw global diff from control")
)

ps=grid.arrange(p_var[[3]])




##############################
# whole brain measure models #
##############################


#WITH global variable covariate 
data_BP_K = datab[!datab$group=="SZ",]
model_bvol_glob = lm(BrainTotalVol ~ group + sex + age + TotalEulerNumber + eICV_samseg + site, data=data_BP_K)
anova(model_bvol_glob)
model_bvol_glob = update(model_bvol_glob,~.-group:sex)
anova(model_bvol_glob)
lsmeans(model_bvol_glob,pairwise~"group",adjust="none")



datab$eICV_samseg_demean = datab$eICV_samseg - mean(datab$eICV_samseg)

datab$ICV_ratio = datab$BrainTotalVol / datab$eICV_samseg

model_bvol_glob = lm(ICV_ratio ~ group*sex + age + site, data=datab)
anova(model_bvol_glob)

model_bvol_glob = lm(BrainTotalVol ~ -1 + group*sex + age + TotalEulerNumber + eICV_samseg + site, data=datab)
model_bvol_glob = lm(BrainTotalVol ~ group*sex + age + TotalEulerNumber + site, data=datab)
#model_bvol_glob = lm(BrainTotalVol ~ group*sex + age + eICV_samseg + site, data=datab)
#model_bvol_glob = lm(BrainTotalVol ~ group*sex + age + TotalEulerNumber + site, data=datab)
model_bvol_glob = update(model_bvol_glob,~.-group:sex)

anova(model_bvol_glob)
Anova(model_bvol_glob,type="II")
Anova(model_bvol_glob,type="III")

drop1(model_bvol_glob, test="F")


lsmeans(model_bvol_glob,pairwise~"group",by="sex",adjust="none")


anova(model_bvol_glob)
summary(model_bvol_glob)
model_bvol_glob = update(model_bvol_glob,~.-group:sex)

model_bvol_glob = lm(BrainTotalVol ~ group*sex + age + TotalEulerNumber + eICV_samseg + site, data=datab)
anova(model_bvol_glob)

model_bvol_glob = lm(BrainTotalVol ~ group + sex + age + TotalEulerNumber + site, data=datab)
anova(model_bvol_glob)
summary(model_bvol_glob)
lsmeans(model_bvol_glob,pairwise~"group",adjust="none")



micv = as.numeric(mean(datab$eICV_samseg))
licv = as.numeric(quantile(datab$eICV_samseg)[2])
uicv = as.numeric(quantile(datab$eICV_samseg)[4])

micv = as.numeric(mean(datab$TotalEulerNumber))
licv = as.numeric(quantile(datab$TotalEulerNumber)[2])
uicv = as.numeric(quantile(datab$TotalEulerNumber)[4])

ls = lsmeans(model_bvol_glob,"group",by="TotalEulerNumber",at = list(TotalEulerNumber = c(licv,micv,uicv)),adjust="none")
plot(ls)
#lsmeans(model_bvol_glob,pairwise~"group",by="eICV_samseg",at = list(eICV_samseg = c(licv,micv,uicv)),adjust="none")
#lsmeans(model_bvol_glob,pairwise~"group",adjust="none")

emm = emmeans(model_bvol_glob,specs = "group",adjust="none")
pwpp(emm)


model_cvol_glob = lm(CortexVol ~ group*sex + age + eICV_samseg + site, data=datab)
anova(model_cvol_glob)
model_cvol_glob = update(model_cvol_glob,~.-group:sex)
anova(model_cvol_glob)

model_area_glob = lm(total_area ~ group*sex + age + eICV_samseg + site, data=datab)
anova(model_area_glob)
model_area_glob = update(model_area_glob,~.-group:sex)
anova(model_area_glob)

model_mt_glob = lm(mean_thickness ~ group*sex + age + eICV_samseg + site, data=datab)
anova(model_mt_glob)
model_mt_glob = update(model_mt_glob,~.-group:sex)
anova(model_mt_glob)
lsmeans(model_mt_glob,pairwise~"group")


#WITHOUT global variable models 
model_bvol = lm(BrainTotalVol ~ group*sex + age + site, data=datab)
anova(model_bvol)
summary(model_bvol)

model_cvol = lm(CortexVol ~ group*sex + age + site, data=datab)
anova(model_cvol)


model_area = lm(total_area ~ group*sex + age + site + TotalEulerNumber, data=datab)
anova(model_area)
ls = lsmeans(model_area,pairwise~"group", by = "sex",adjust="none")
confint(ls)


model_mt = lm(mean_thickness ~ group*sex  + age + site  + TotalEulerNumber, data=datab)
anova(model_mt)
model_mt = update(model_mt,~.-group:sex)
anova(model_mt)
lsmeans(model_mt,pairwise~"group", by = "sex",adjust="none")


model_icv = lm(eICV_samseg ~ group*sex + age + site, data=datab)
anova(model_icv)



### run models with eICV_samseg as covariate and without

model_vars = c("BrainTotalVol", "CortexVol", "total_area","mean_thickness", "eICV_samseg")
glob = c("with_eICV", "without_eICV")
GS_pvals = list()
ps = list()
emm = list()
min_emm = list()
max_emm = list()
model = list()

df_with_ICV = data.frame()
df_without_ICV = data.frame()


for (i in seq(1,length(model_vars))){
  if (model_vars[i] == "eICV_samseg"){
    glob = c("without_eICV")
  }
  else{
    glob = c("with_eICV", "without_eICV")
  }

  for (g in seq(1,length(glob))){
    if (glob[g] == "with_eICV"){
      f = paste(model_vars[i],"~","+","group*sex","+","age","+","site","+","TotalEulerNumber","+","eICV_samseg")
    }
    else{
      f = paste(model_vars[i],"~","+","group*sex","+","age","+","site","+","TotalEulerNumber")
    }
    model[[glob[g]]][[i]] = lm(f,data=datab)
    xvars = attributes(Anova(model[[glob[g]]][[i]],type = "III"))$row.names
      
    if (Anova(model[[glob[g]]][[i]],type = "III")$"Pr(>F)"[xvars=="group:sex"] > 0.05){
      model_gs = model[[glob[g]]][[i]]
      GS_pvals[[glob[g]]][[model_vars[i]]] = Anova(model_gs,type = "III")$"Pr(>F)"[xvars=="group:sex"]
      GS_F = Anova(model_gs,type = "III")$"F value"[xvars=="group:sex"]
      
      model[[glob[g]]][[i]] = update(model[[glob[g]]][[i]],~.-group:sex)
      sign_GS = 0
      models = list(model_gs,model[[glob[g]]][[i]])
    }
    else{
      GS_pvals[[glob[g]]][[model_vars[i]]] = Anova(model[[glob[g]]][[i]],type = "III")$"Pr(>F)"[xvars=="group:sex"]
      GS_F = Anova(model[[glob[g]]][[i]],type = "III")$"F value"[xvars=="group:sex"]
      sign_GS = 1
      models = list(model[[glob[g]]][[i]])
    }
    
    for (m in seq(1,length(models))){
      mm = models[[m]]
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
      
      if (glob[g] == "with_eICV"){
        ICV_F = Anova(mm,type = "III")$"F value"[xvars=="eICV_samseg"]
        ICV_pv = Anova(mm,type = "III")$"Pr(>F)"[xvars=="eICV_samseg"]
        
        if (m==1){
          rw = list(model_vars[i], group_F, group_pv, sex_F, sex_pv, age_F, age_pv, site_F, site_pv, EulerNumber_F, EulerNumber_pv, ICV_F, ICV_pv, GS_F, GS_pvals[[glob[g]]][[model_vars[i]]],sign_GS,m)
        }
        else{
          rw = list(model_vars[i], group_F, group_pv, sex_F, sex_pv, age_F, age_pv, site_F, site_pv, EulerNumber_F, EulerNumber_pv, ICV_F, ICV_pv, NA,   NA,                                  sign_GS,m-2)
        }
        df_with_ICV = rbindlist(list(df_with_ICV, rw))
      }
      else{
        if (m==1){
          rw = list(model_vars[i], group_F, group_pv, sex_F, sex_pv, age_F, age_pv, site_F, site_pv, EulerNumber_F, EulerNumber_pv, GS_F, GS_pvals[[glob[g]]][[model_vars[i]]],sign_GS,m)
        }
        else{
          rw = list(model_vars[i], group_F, group_pv, sex_F, sex_pv, age_F, age_pv, site_F, site_pv, EulerNumber_F, EulerNumber_pv, NA,NA,sign_GS,m-2)
        }
        df_without_ICV = rbindlist(list(df_without_ICV, rw))
      }
      
    }
    
    emm[[g]] = emmeans(model[[glob[g]]][[i]],specs = "group",by="sex")
    ps[[g]] = pwpp(emm[[g]],adjust="none",sort = FALSE) +
              #aes(y = factor("group", level = c("BP","K","SZ"))) +
              labs(x="Uncorrected P-value") +
              ggtitle(paste(glob[g])) +
              geom_vline(xintercept = 0.05,linetype="dashed") +
              coord_flip()
    
    min_emm[[g]] = min(summary(emm[[g]])$emmean)
    max_emm[[g]] = max(summary(emm[[g]])$emmean)
  } #g
  
  for (g in seq(1,2)){
  ps[[g+2]] = plot(emm[[g]]) +
              aes(color=group) +
              facet_grid(cols =vars(sex)) +
              scale_x_continuous(limits = c(min(unlist(min_emm)), max(unlist(max_emm)))) +
              scale_y_discrete(limits=c("BP","K","SZ")) +
              coord_flip()
  }
  
  top_title = paste(model_vars[i]," models")
  ga=grid.arrange(grobs=ps, ncol=2,top=textGrob(top_title,gp=gpar(fontsize=20)))
  ggsave(paste(model_vars[i],"_group_diff_pvalues_ICV",".png",sep=""),ga,width = 10,height = 10)
}

col_names = c("Model_yvar","Group_Fval","Group_pval","Sex_Fval","Sex_pval","Age_Fval","Age_pval","Site_Fval","Site_pval","EulerNumber_Fval","Eulernumber_pval","Group_sex_Fval","Group_sex_pval","Significant_GS_interaction","model_with_GS")
col_names_ICV = c("Model_yvar","Group_Fval","Group_pval","Sex_Fval","Sex_pval","Age_Fval","Age_pval","Site_Fval","Site_pval","EulerNumber_Fval","Eulernumber_pval","ICV_Fval","ICV_pval","Group_sex_Fval","Group_sex_pval","Significant_GS_interaction","model_with_GS")

names(df_with_ICV)<-col_names_ICV
names(df_without_ICV)<-col_names

df_with_ICV[,2:ncol(df_with_ICV)] <- signif(df_with_ICV[,2:ncol(df_with_ICV)],digits=3)
df_without_ICV[,2:ncol(df_without_ICV)] <- signif(df_without_ICV[,2:ncol(df_without_ICV)],digits=3)

write_xlsx(df_with_ICV,"/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/ANOVA_pvals_with_ICV.xlsx")
write_xlsx(df_without_ICV,"/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/ANOVA_pvals_without_ICV.xlsx")



### run models with euler number as covariate and without

model_vars = c("BrainTotalVol", "CortexVol", "total_area","mean_thickness")
glob = c("with_euler", "without_euler")
ps = list()
emm = list()
min_emm = list()
max_emm = list()
model_e = list()

for (i in seq(1,length(model_vars))){
  for (g in seq(1,2)){
    if (glob[g] == "with_euler"){
      f = paste(model_vars[i],"~","+","group*sex","+","age","+","TotalEulerNumber","+","site")
    }
    else{
      f = paste(model_vars[i],"~","+","group*sex","+","age","+","site")
    }
    model_e[[glob[g]]][[i]] = lm(f,data=datab)
    xvars = attributes(anova(model_e[[glob[g]]][[i]]))$row.names
    
    if (anova(model_e[[glob[g]]][[i]])$"Pr(>F)"[xvars=="group:sex"] > 0.05){
      model_e[[glob[g]]][[i]] = update(model_e[[glob[g]]][[i]],~.-group:sex)
    } 
    
    emm[[g]] = emmeans(model_e[[glob[g]]][[i]],specs = "group",by="sex")
    ps[[g]] = pwpp(emm[[g]],adjust="none",sort = FALSE) +
      #aes(y = factor("group", level = c("BP","K","SZ"))) +
      labs(x="Uncorrected P-value") +
      ggtitle(paste(glob[g])) +
      geom_vline(xintercept = 0.05,linetype="dashed") +
      coord_flip()
    
    min_emm[[g]] = min(summary(emm[[g]])$emmean)
    max_emm[[g]] = max(summary(emm[[g]])$emmean)
  } #g
  
  for (g in seq(1,2)){
    ps[[g+2]] = plot(emm[[g]]) +
      aes(color=group) +
      facet_grid(cols =vars(sex)) +
      scale_x_continuous(limits = c(min(unlist(min_emm)), max(unlist(max_emm)))) +
      scale_y_discrete(limits=c("BP","K","SZ")) +
      coord_flip()
  }
  
  top_title = paste(model_vars[i]," models","without ICV")
  ga=grid.arrange(grobs=ps, ncol=2,top=textGrob(top_title,gp=gpar(fontsize=20)))
  ggsave(paste("eulernumber_plots/",model_vars[i],"_group_diff_pvalues_Euler",".png",sep=""),ga,width = 10,height = 10)
}




####### LSmeans estimates on the variance plot #########

y_vars = c("BrainTotalVol", "CortexVol", "total_area", "mean_thickness")
highriskG = c("BP","SZ")
glob = c("with_eICV", "without_eICV")
ss = c("sex0","sex1")
sex_name = c("Female","Male")
pd <- position_dodge(0.05) 

pls = list()


dftm = datab %>%
  pivot_longer(cols = y_vars,
               names_to = "NAME_measure",
               values_to = "VALUE_measure",
               values_drop_na = FALSE)
dftm = as.data.frame(dftm)

dftm$NORM_measure = dftm$VALUE_measure
dftm$NORM_measure_K = NA #dftm$VALUE_measure
dftm$NORM_measure_K_mt = NA #dftm$VALUE_measure


dfc_with_ICV = data.frame()
dfc_without_ICV = data.frame()

for (g in seq(1,length(glob))){
  
for (i in seq(1,length(y_vars))){
  avg_K = mean(dftm$VALUE_measure[ dftm$NAME_measure  == y_vars[i] & dftm$group  == "K"])

  if (y_vars[i] == "mean_thickness"){
    dftm$NORM_measure_K_mt[ dftm$NAME_measure  == y_vars[i] ] = ( dftm$VALUE_measure[ dftm$NAME_measure  == y_vars[i] ] - avg_K ) / avg_K
    #ls = lsmeans(model[[glob[g]]][[i]],pairwise~"group",adjust="none") ## THIS NEEDS TO BE FIXED
    }
  else{
    dftm$NORM_measure_K[ dftm$NAME_measure  == y_vars[i] ] = ( dftm$VALUE_measure[ dftm$NAME_measure  == y_vars[i] ] - avg_K ) / avg_K
    ls = lsmeans(model[[glob[g]]][[i]],pairwise~"group",by="sex",adjust="none")
    }
  
  
  lss = lsmeans(model[[glob[g]]][[i]],pairwise~"group",by="sex",adjust="none")
  cc = lss$contrasts
  

  if (GS_pvals[[glob[g]]][[y_vars[i]]] > 0.05){
    sign_GS = 0
    lss = lsmeans(model[[glob[g]]][[i]],pairwise~"group",adjust="none")
    cc = lss$contrasts
    
    #raw contrasts
    BP_diff = summary(cc)$estimate[1]
    BP_diff_tratio = summary(cc)$t.ratio[1]
    BP_diff_pval = summary(cc)$p.value[1]
    
    SZ_diff = summary(cc)$estimate[3]
    SZ_diff_tratio = summary(cc)$t.ratio[3]
    SZ_diff_pval = summary(cc)$p.value[3]
    
    rw = list(y_vars[i], sign_GS, 
              BP_diff, BP_diff_tratio, BP_diff_pval,
              SZ_diff, SZ_diff_tratio, SZ_diff_pval,
              as.numeric(NA), as.numeric(NA), as.numeric(NA),
              as.numeric(NA), as.numeric(NA), as.numeric(NA),
              as.numeric(NA), as.numeric(NA), as.numeric(NA), 
              as.numeric(NA), as.numeric(NA), as.numeric(NA))
    
  } 
  else {
    sign_GS = 1
    lss = lsmeans(model[[glob[g]]][[i]],pairwise~"group",by="sex",adjust="none")
    cc = lss$contrasts
    
    #raw contrasts
    sex0_BP_diff = summary(cc)$estimate[summary(cc)$sex == 0][1]
    sex0_BP_diff_tratio = summary(cc)$t.ratio[summary(cc)$sex == 0][1]
    sex0_BP_diff_pval = summary(cc)$p.value[summary(cc)$sex == 0][1]
    
    sex0_SZ_diff = summary(cc)$estimate[summary(cc)$sex == 0][3]
    sex0_SZ_diff_tratio = summary(cc)$t.ratio[summary(cc)$sex == 0][3]
    sex0_SZ_diff_pval = summary(cc)$p.value[summary(cc)$sex == 0][3]
    
    sex1_BP_diff = summary(cc)$estimate[summary(cc)$sex == 1][1]
    sex1_BP_diff_tratio = summary(cc)$t.ratio[summary(cc)$sex == 1][1]
    sex1_BP_diff_pval = summary(cc)$p.value[summary(cc)$sex == 1][1]
    
    sex1_SZ_diff = summary(cc)$estimate[summary(cc)$sex == 1][3]
    sex1_SZ_diff_tratio = summary(cc)$t.ratio[summary(cc)$sex == 1][3]
    sex1_SZ_diff_pval = summary(cc)$p.value[summary(cc)$sex == 1][3]
    
    rw = list(y_vars[i], sign_GS, 
              as.numeric(NA), as.numeric(NA), as.numeric(NA),
              as.numeric(NA), as.numeric(NA), as.numeric(NA),
              sex0_BP_diff, sex0_BP_diff_tratio, sex0_BP_diff_pval,
              sex0_SZ_diff, sex0_SZ_diff_tratio, sex0_SZ_diff_pval,
              sex1_BP_diff, sex1_BP_diff_tratio, sex1_BP_diff_pval,
              sex1_SZ_diff, sex1_SZ_diff_tratio, sex1_SZ_diff_pval)
    
  } 
  
  
  ls = lsmeans(model[[glob[g]]][[i]],pairwise~"group",by="sex",adjust="none")
  c = ls$contrasts
  
  K_emm = summary(ls)$lsmeans$lsmean[summary(ls$lsmeans)$sex == 0][2]
  sex1_K_emm = summary(ls)$lsmeans$lsmean[summary(ls$lsmeans)$sex == 1][2]
  
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

  
  if (glob[g] == "with_eICV"){
    dfc_with_ICV = rbindlist(list(dfc_with_ICV, rw))
  }
  else{
    dfc_without_ICV = rbindlist(list(dfc_without_ICV, rw))
  }
  
  
  
  ###
  dftm$lsmean_sex0[dftm$NAME_measure  == y_vars[i] & dftm$group == "BP" & dftm$sex == 0] = BP_diff_P
  dftm$eb_max_sex0[dftm$NAME_measure  == y_vars[i] & dftm$group == "BP" & dftm$sex == 0] = BP_diff_UCL_P
  dftm$eb_min_sex0[dftm$NAME_measure  == y_vars[i] & dftm$group == "BP" & dftm$sex == 0] = BP_diff_LCL_P
  
  dftm$lsmean_sex1[dftm$NAME_measure  == y_vars[i] & dftm$group == "BP" & dftm$sex == 1] = sex1_BP_diff_P
  dftm$eb_max_sex1[dftm$NAME_measure  == y_vars[i] & dftm$group == "BP" & dftm$sex == 1] = sex1_BP_diff_UCL_P
  dftm$eb_min_sex1[dftm$NAME_measure  == y_vars[i] & dftm$group == "BP" & dftm$sex == 1] = sex1_BP_diff_LCL_P
  
  
  ###
  dftm$lsmean_sex0[dftm$NAME_measure  == y_vars[i] & dftm$group == "SZ" & dftm$sex == 0] = SZ_diff_P
  dftm$eb_max_sex0[dftm$NAME_measure  == y_vars[i] & dftm$group == "SZ" & dftm$sex == 0] = SZ_diff_UCL_P
  dftm$eb_min_sex0[dftm$NAME_measure  == y_vars[i] & dftm$group == "SZ" & dftm$sex == 0] = SZ_diff_LCL_P
  
  dftm$lsmean_sex1[dftm$NAME_measure  == y_vars[i] & dftm$group == "SZ" & dftm$sex == 1] = sex1_SZ_diff_P
  dftm$eb_max_sex1[dftm$NAME_measure  == y_vars[i] & dftm$group == "SZ" & dftm$sex == 1] = sex1_SZ_diff_UCL_P
  dftm$eb_min_sex1[dftm$NAME_measure  == y_vars[i] & dftm$group == "SZ" & dftm$sex == 1] = sex1_SZ_diff_LCL_P
  
}
  
  if (glob[g] == "with_eICV"){
    dfc_plot_with_eICV = dftm
  }
  else{
    dfc_plot_without_eICV = dftm
  }
  
}


col_names = c("Model_yvar", "Significant_GS_interaction",
              "Contrast_BP-K","tratio_BP-K","pval_BP-K",
              "Contrast_SZ-K","tratio_SZ-K","pval_SZ-K",
              "Sex0_Contrast_BP-K","Sex0_tratio_BP-K","Sex0_pval_BP-K",
              "Sex0_Contrast_SZ-K","Sex0_tratio_SZ-K","Sex0_pval_SZ-K",
              "Sex1_Contrast_BP-K","Sex1_tratio_BP-K","Sex1_pval_BP-K",
              "Sex1_Contrast_SZ-K","Sex1_tratio_SZ-K","Sex1_pval_SZ-K")
names(dfc_with_ICV)<-col_names
names(dfc_without_ICV)<-col_names

dfc_with_ICV[,2:ncol(dfc_with_ICV)] <- signif(dfc_with_ICV[,2:ncol(dfc_with_ICV)],digits=3)
dfc_without_ICV[,2:ncol(dfc_without_ICV)] <- signif(dfc_without_ICV[,2:ncol(dfc_without_ICV)],digits=3)

#write_xlsx(dfc_with_ICV,"/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/Model_contrasts_with_ICV.xlsx")
#write_xlsx(dfc_without_ICV,"/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/Model_contrasts_without_ICV.xlsx")




#### REDO CONTRAST XLSX with male / female split

model_vars = c("BrainTotalVol", "CortexVol", "total_area", "mean_thickness", "eICV_samseg")

data_sex1 = datab[c(datab$sex == 1),]
data_sex0 = datab[c(datab$sex == 0),]

model_sex = list()

dfc_with_ICV = data.frame()
dfc_without_ICV = data.frame()

for (i in seq(1,length(model_vars))){
  if (model_vars[i] == "eICV_samseg"){
    glob = c("without_eICV")
  }
  else{
    glob = c("with_eICV", "without_eICV")
  }
  for (g in seq(1,length(glob))){
    
    if (glob[g] == "with_eICV"){
      f = paste(model_vars[i],"~","+","group*sex","+","age","+","site","+","TotalEulerNumber","+","eICV_samseg")
    }
    else{
      f = paste(model_vars[i],"~","+","group*sex","+","age","+","site","+","TotalEulerNumber")
    }
    
    model[[glob[g]]][[i]] = lm(f,data=datab)
    xvars = attributes(anova(model[[glob[g]]][[i]]))$row.names
    
    if (anova(model[[glob[g]]][[i]])$"Pr(>F)"[xvars=="group:sex"] > 0.05){
      
      model[[glob[g]]][[i]] = update(model[[glob[g]]][[i]],~.-group:sex)
      
      sign_GS = 0
      lss = lsmeans(model[[glob[g]]][[i]],pairwise~"group",adjust="none")
      cc = lss$contrasts
      
      #raw contrasts
      BP_diff = summary(cc)$estimate[1]
      BP_diff_tratio = summary(cc)$t.ratio[1]
      BP_diff_pval = summary(cc)$p.value[1]
      
      SZ_diff = summary(cc)$estimate[3]
      SZ_diff_tratio = summary(cc)$t.ratio[3]
      SZ_diff_pval = summary(cc)$p.value[3]
      
      rw = list(y_vars[i], sign_GS, 
                BP_diff, BP_diff_tratio, BP_diff_pval,
                SZ_diff, SZ_diff_tratio, SZ_diff_pval,
                as.numeric(NA), as.numeric(NA), as.numeric(NA),
                as.numeric(NA), as.numeric(NA), as.numeric(NA),
                as.numeric(NA), as.numeric(NA), as.numeric(NA), 
                as.numeric(NA), as.numeric(NA), as.numeric(NA))
      
    } 
    else {
      
      sign_GS = 1
      
      rw = list(model_vars[i], sign_GS, 
                as.numeric(NA), as.numeric(NA), as.numeric(NA),
                as.numeric(NA), as.numeric(NA), as.numeric(NA))
                
      for (s in seq(1,2)){

        if (glob[g] == "with_eICV"){
          f = paste(model_vars[i],"~","+","group","+","age","+","site","+","TotalEulerNumber","+","eICV_samseg")
        }
        else{
          f = paste(model_vars[i],"~","+","group","+","age","+","site","+","TotalEulerNumber")
        }
        
        if (s-1 == 0){
          model_sex[[glob[g]]][[model_vars[i]]][[s]] = lm(f,data=data_sex0)
        }
        else{
          model_sex[[glob[g]]][[model_vars[i]]][[s]] = lm(f,data=data_sex1)
        }
        
        lss = lsmeans(model_sex[[glob[g]]][[model_vars[i]]][[s]],pairwise~"group",adjust="none")
        cc = lss$contrasts
        
        #raw contrasts
        BP_diff = summary(cc)$estimate[1]
        BP_diff_tratio = summary(cc)$t.ratio[1]
        BP_diff_pval = summary(cc)$p.value[1]
        
        SZ_diff = summary(cc)$estimate[3]
        SZ_diff_tratio = summary(cc)$t.ratio[3]
        SZ_diff_pval = summary(cc)$p.value[3]
        
        rw = append(rw,list(BP_diff, BP_diff_tratio, BP_diff_pval,
                       SZ_diff, SZ_diff_tratio, SZ_diff_pval))
        
      } #end for s
    } #end else for GS sig 1
    
    if (glob[g] == "with_eICV"){
      dfc_with_ICV = rbindlist(list(dfc_with_ICV, rw))
    }
    else{
      dfc_without_ICV = rbindlist(list(dfc_without_ICV, rw))
    }
    
  }
  
}


col_names = c("Model_yvar", "Significant_GS_interaction",
              "Contrast_BP-K","tratio_BP-K","pval_BP-K",
              "Contrast_SZ-K","tratio_SZ-K","pval_SZ-K",
              "Sex0_Contrast_BP-K","Sex0_tratio_BP-K","Sex0_pval_BP-K",
              "Sex0_Contrast_SZ-K","Sex0_tratio_SZ-K","Sex0_pval_SZ-K",
              "Sex1_Contrast_BP-K","Sex1_tratio_BP-K","Sex1_pval_BP-K",
              "Sex1_Contrast_SZ-K","Sex1_tratio_SZ-K","Sex1_pval_SZ-K")
names(dfc_with_ICV)<-col_names
names(dfc_without_ICV)<-col_names

dfc_with_ICV[,2:ncol(dfc_with_ICV)] <- signif(dfc_with_ICV[,2:ncol(dfc_with_ICV)],digits=3)
dfc_without_ICV[,2:ncol(dfc_without_ICV)] <- signif(dfc_without_ICV[,2:ncol(dfc_without_ICV)],digits=3)

write_xlsx(dfc_with_ICV,"/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/Model_contrasts_with_ICV.xlsx")
write_xlsx(dfc_without_ICV,"/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/Model_contrasts_without_ICV.xlsx")






### PLOT THE LSMEANS 


for (g in seq(1,length(glob))){
  DF = get(paste("dfc_plot_",glob[g],sep = ""))
  
for (h in seq(1,2)){
  pls[[h]]= with(DF[!is.na(dftm$NORM_measure_K), ],
                    ggplot() + 
                     aes(x = factor(NAME_measure), y = 100*NORM_measure_K, color = group) +
                     geom_violin(position = "identity",alpha=0.3) +
                     geom_jitter(width = 0.3, size=0.1) + 
                     geom_hline(yintercept = 0,linetype="dashed") + 
                     
                     geom_point(position=pd, aes_string(x = "NAME_measure", y = paste("lsmean",ss[h],sep = "_"), color="group", group="group"), size=2) +
                     #geom_line(position=pd, aes_string(x = "NAME_measure", y = paste("lsmean",ss[h],sep = "_"), color="group", group="group")) +
                      
                     stat_summary(fun.y = mean,geom = "line",position=pd, aes_string(x = "NAME_measure", y = paste("lsmean",ss[h],sep = "_"), color="group", group="group")) + 
                   
                     geom_errorbar(position=pd, width=0.2, aes_string(group="group",color="group",ymin=paste("eb_min_",ss[h],sep=""), ymax=paste("eb_max_",ss[h],sep="")) ) +

                     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
                     labs(y = "Difference from control [%]", x = "brain measure") +
                     ggtitle(paste(sex_name[h])) +
                     
                     annotate("text", size=2.5, x = c(1,2,3), y=21, label = c(paste("G:S pval =",round(GS_pvals[[glob[g]]][[y_vars[1]]],digits=3)),
                                                                    paste("G:S pval =",round(GS_pvals[[glob[g]]][[y_vars[2]]],digits=3)),
                                                                    paste("G:S pval =",round(GS_pvals[[glob[g]]][[y_vars[3]]],digits=3))) )
  )
  
  pls[[3]]= with(DF[!is.na(dftm$NORM_measure_K_mt), ],
                 ggplot() + 
                   aes(x = factor(NAME_measure), y = 100*NORM_measure_K_mt, color = group) +
                   geom_violin(position = "identity",alpha=0.3) +
                   geom_jitter(width = 0.3, size=0.1) + 
                   geom_hline(yintercept = 0,linetype="dashed") + 
                   
                   geom_point(position=pd, aes_string(x = "NAME_measure", y = paste("lsmean",ss[h],sep = "_"), color=group, group=group), size=2) +
                   geom_errorbar(position=pd, aes_string(width=0.1,group="group",color="group",ymin=paste("eb_min_",ss[h],sep=""), ymax=paste("eb_max_",ss[h],sep="")) ) +

                   #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
                   labs(y = "Difference from control [%]", x = "brain measure") +
                   ggtitle("Both sex") + 
                 
                    annotate("text",size=2.5, x = 1, y=7.9, label = c(paste("G:S pval =",round(GS_pvals[[glob[g]]][[y_vars[4]]],digits=3))))
  )
  pls[[4]] = ggplot() + theme_void()
  
} #h = for each highrisk group in subplot

top_title = paste("LSmeans contrasts of model",glob[g])
ga=grid.arrange(grobs=pls,ncols=2, top=textGrob(top_title,gp=gpar(fontsize=20)))
ggsave(paste("LSmeans_contrasts_on_datapoints_",glob[g],".png",sep=""),ga,width = 10,height = 10)

} #g







### model diagnostics 
par(mar=c(1,1,1,1))
par(mfrow=c(length(model_vars),3),mai=c(0.3,0.3,0.3,0.3))
for (i in seq(1,length(model_vars))){
  mixed_res = rstudent(model[[i]])
  qqnorm(mixed_res,main = NULL)
  qqline(mixed_res,main = NULL)
  plot(mixed_res ~ fitted(model[[i]]),xlab="Fitted",ylab="Standardized residuals")
  plot(cooks.distance(model[[i]]), type = "p", pch = 20,ylab="Cooks distance")
  
  mtext(paste("Model diagnostics for",model_vars[[i]]), side = 3, line = -12.9*(i-1)-2, outer = TRUE)
}
par(mfrow=c(1,1))





###################################
#######     ARCHIVE     ###########
###################################


#ps[[g]] = plot(PH,comparisons=TRUE,xlab="LSmean estiamte") +
#          ggtitle(paste(glob[g])) +
#          scale_x_continuous(limits = c(min(summary(PH)$lsmeans$lsmean), max(summary(PH)$lsmeans$lsmean)))

#plot group:age interation
par(mfrow=c(1,5))
for (i in seq(1,length(model))){
  y_var = paste(model_vars[[i]])
  plot(datab$age[datab$group=="K"],datab[[y_var]][datab$group=="K"],col=group_color[1],ylab = y_var, xlab="age")
  points(datab$age[datab$group=="SZ"],datab[[y_var]][datab$group=="SZ"],col=group_color[2])
  points(datab$age[datab$group=="BP"],datab[[y_var]][datab$group=="BP"],col=group_color[3])
  
  legend("topright", legend=c(as.character(unique(datab$group))),
         col=group_color, pch=c(1,1,1), cex=1)
  
  B <- fixef(model[[i]])
  with(datab, {
    abline(a=B["groupK"], b=B["groupK:age"], lty=1, col="green")
    abline(a=B["groupSZ"], b=B["groupSZ:age"], lty=2, col="blue")
    abline(a=B["groupBP"], b=B["groupBP:age"], lty=3, col="red")
  })
}


modelt = lmer(total_area ~ group*sex + age + eICV_samseg + (1| site),data=datab)
anova(modelt)
modelt = update(modelt,~.-group:sex)
anova(modelt)
ls = lsmeans(modelt,pairwise~"group")
ls = lsmeans(modelt,pairwise~"group",by="sex")
plot(ls,comparisons=TRUE)



modelt = lmer(total_area ~ group*sex + age + (1| site),data=datab)
anova(modelt)
modelt = update(modelt,~.-group:sex)
anova(modelt)
ls = lsmeans(modelt,pairwise~"group")
ls = lsmeans(modelt,pairwise~"group",by="sex")

plot(ls,comparisons=TRUE)

emm2 = emmeans(modelt,specs = "group")
plot(emm2,comparisons = TRUE)




library(multcomp)
library(emmeans)
library(multcompView)

emm = emmeans(model[[glob[g]]][[i]],specs = "group",by="sex")

pairs(emm,
      alpha=0.05,
      Letters=letters,
      adjust="bonferroni")

pairs(emm,
      alpha=0.05,
      Letters=letters,
      adjust="tukey")





####### LSmeans estimates on the variance plot #########

y_vars = c("BrainTotalVol", "CortexVol", "total_area", "mean_thickness")
highriskG = c("BP","SZ")
glob = c("with_eICV", "without_eICV")

pls = list()

dftm = datab %>%
  pivot_longer(cols = y_vars,
               names_to = "NAME_measure",
               values_to = "VALUE_measure",
               values_drop_na = FALSE)
dftm = as.data.frame(dftm)

dftm$NORM_measure = dftm$VALUE_measure
dftm$NORM_measure_K = NA #dftm$VALUE_measure
dftm$NORM_measure_K_mt = NA #dftm$VALUE_measure


for (g in seq(1,length(glob))){
  
  for (i in seq(1,length(y_vars))){
    avg_K = mean(dftm$VALUE_measure[ dftm$NAME_measure  == y_vars[i] & dftm$group  == "K"])
    
    if (y_vars[i] == "mean_thickness"){
      dftm$NORM_measure_K_mt[ dftm$NAME_measure  == y_vars[i] ] = ( dftm$VALUE_measure[ dftm$NAME_measure  == y_vars[i] ] - avg_K ) / avg_K
    }
    else{
      dftm$NORM_measure_K[ dftm$NAME_measure  == y_vars[i] ] = ( dftm$VALUE_measure[ dftm$NAME_measure  == y_vars[i] ] - avg_K ) / avg_K
    }
    
    
    ls = lsmeans(model[[glob[g]]][[i]],pairwise~"group",by="sex",adjust="none")
    c = ls$contrasts
    K_emm = summary(ls)$lsmeans$lsmean[summary(ls$lsmeans)$sex == 0][2]
    sex1_K_emm = summary(ls)$lsmeans$lsmean[summary(ls$lsmeans)$sex == 1][2]
    
    
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
    
    #BP
    dftm$lsmean_BP[dftm$NAME_measure  == y_vars[i]] = BP_diff_P
    dftm$lsmean_BP[dftm$NAME_measure  == y_vars[i] & dftm$sex == 1] = sex1_BP_diff_P
    
    dftm$eb_max_BP[dftm$NAME_measure  == y_vars[i]] = BP_diff_UCL_P
    dftm$eb_min_BP[dftm$NAME_measure  == y_vars[i]] = BP_diff_LCL_P
    
    dftm$eb_max_BP[dftm$NAME_measure  == y_vars[i] & dftm$sex == 1] = sex1_BP_diff_UCL_P
    dftm$eb_min_BP[dftm$NAME_measure  == y_vars[i] & dftm$sex == 1] = sex1_BP_diff_LCL_P
    
    #SZ
    dftm$lsmean_SZ[dftm$NAME_measure  == y_vars[i]] = SZ_diff_P
    dftm$lsmean_SZ[dftm$NAME_measure  == y_vars[i] & dftm$sex == 1] = sex1_SZ_diff_P
    
    dftm$eb_max_SZ[dftm$NAME_measure  == y_vars[i]] = SZ_diff_UCL_P
    dftm$eb_min_SZ[dftm$NAME_measure  == y_vars[i]] = SZ_diff_LCL_P
    
    dftm$eb_max_SZ[dftm$NAME_measure  == y_vars[i] & dftm$sex == 1] = sex1_SZ_diff_UCL_P
    dftm$eb_min_SZ[dftm$NAME_measure  == y_vars[i] & dftm$sex == 1] = sex1_SZ_diff_LCL_P
    
  }
  
  
  for (h in seq(1,2)){
    pls[[2*h-1]]= with(dftm[!is.na(dftm$NORM_measure_K), ],
                       ggplot() + 
                         aes(x = factor(NAME_measure), y = 100*NORM_measure_K, color = group) +
                         geom_violin(position = "identity",alpha=0.3) +
                         geom_jitter(width = 0.3, size=0.1) + 
                         geom_hline(yintercept = 0,linetype="dashed") + 
                         
                         ggnewscale::new_scale_colour() +
                         
                         geom_point(aes_string(x = "NAME_measure", y = paste("lsmean",highriskG[h],sep = "_"), color=sex, group=sex), size=2) +
                         scale_colour_manual("Contrasts per sex", values = c("red", "blue")) + 
                         
                         geom_errorbar(aes_string(width=0.05,group="sex",color="sex",ymin=paste("eb_min_",highriskG[h],sep=""), ymax=paste("eb_max_",highriskG[h],sep="")) ) +
                         geom_line(aes_string(x = "NAME_measure", y = paste("lsmean",highriskG[h],sep = "_"), color="sex", group="sex")) +
                         
                         theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
                         labs(y = "Difference from control [%]", x = "brain measure") +
                         ggtitle(paste(highriskG[h])) +
                         
                         annotate("text", size=2.5, x = c(1,2,3), y=21, label = c(paste("G:S pval =",round(GS_pvals[[glob[g]]][[y_vars[1]]],digits=3)),
                                                                                  paste("G:S pval =",round(GS_pvals[[glob[g]]][[y_vars[2]]],digits=3)),
                                                                                  paste("G:S pval =",round(GS_pvals[[glob[g]]][[y_vars[3]]],digits=3))) )
    )
    
    pls[[2*h]]= with(dftm[!is.na(dftm$NORM_measure_K_mt), ],
                     ggplot() + 
                       aes(x = factor(NAME_measure), y = 100*NORM_measure_K_mt, color = group) +
                       geom_violin(position = "identity",alpha=0.3) +
                       geom_jitter(width = 0.3, size=0.1) + 
                       geom_hline(yintercept = 0,linetype="dashed") + 
                       
                       ggnewscale::new_scale_colour() +
                       
                       geom_point(aes_string(x = "NAME_measure", y = paste("lsmean",highriskG[h],sep = "_"), color=sex, group=sex), size=2) +
                       scale_colour_manual("Contrasts per sex", values = c("red", "blue")) + 
                       
                       geom_errorbar(aes_string(width=0.05,group="sex",color="sex",ymin=paste("eb_min_",highriskG[h],sep=""), ymax=paste("eb_max_",highriskG[h],sep="")) ) +
                       geom_line(aes_string(x = "NAME_measure", y = paste("lsmean",highriskG[h],sep = "_"), color="sex", group="sex")) +
                       
                       #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
                       labs(y = "Difference from control [%]", x = "brain measure") +
                       ggtitle(paste(highriskG[h])) + 
                       
                       annotate("text",size=2.5, x = 1, y=7.9, label = c(paste("G:S pval =",round(GS_pvals[[glob[g]]][[y_vars[4]]],digits=3))))
    )
    
  } #h = for each highrisk group in subplot
  
  top_title = paste("LSmeans contrasts of model",glob[g])
  ga=grid.arrange(grobs=pls,ncols=2, top=textGrob(top_title,gp=gpar(fontsize=20)))
  #ggsave(paste("LSmeans_contrasts_on_datapoints_",glob[g],"_alt.png",sep=""),ga,width = 10,height = 10)
  
} #g
