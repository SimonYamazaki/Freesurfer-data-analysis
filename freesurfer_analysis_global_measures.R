
setwd("/mnt/projects/VIA11/FREESURFER/Stats/Data")

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
#library(forcats)
#library(ggpcp)
#library(stringr)
#library(Cairo)
#library(grDevices)
#library(stringr)
#library(see)
#library(lmerTest)
#library(lme4)
#install.packages("devtools")
#remotes::install_github("yaweige/ggpcp", build_vignettes = TRUE)
#library(ggpcp)

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


#make new variables with shorter and contained names
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
dft$NORM_measure_K = NA
dft$NORM_measure_K_mt = NA


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

model_bvol_glob = lm(BrainTotalVol ~ group*sex + age + TotalEulerNumber + site + eICV_samseg, data=datab)
model_bvol_glob = update(model_bvol_glob,~.-group:sex)

anova(model_bvol_glob)
Anova(model_bvol_glob,type="II")
Anova(model_bvol_glob,type="III")
drop1(model_bvol_glob, test="F")

lsmeans(model_bvol_glob,pairwise~"group",adjust="none")




### run models with eICV_samseg as covariate and without
# for ANOVA tables with GS interaction

model_vars = c("BrainTotalVol", "CortexVol", "total_area","mean_thickness", "eICV_samseg")
GS_pvals = list()
ps = list()
emm = list()
min_emm = list()
max_emm = list()
model = list()

DF = data.frame()

for (i in seq(1,length(model_vars))){
  if (model_vars[i] == "eICV_samseg"){
    glob = c("without_eICV")
  }
  else{
    glob = c("without_eICV","with_eICV")
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
      }
      else{
        ICV_F = NA
        ICV_pv = NA
      }       

      if (m==1){
        rw = list(model_vars[i], 
                  group_F, group_pv, sex_F, sex_pv,
                  age_F, age_pv, site_F, site_pv, 
                  EulerNumber_F, EulerNumber_pv, 
                  ICV_F, ICV_pv, GS_F, GS_pvals[[glob[g]]][[model_vars[i]]],
                  sign_GS,m,g-1)
      }
      else{
        rw = list(model_vars[i], 
                  group_F, group_pv, sex_F, sex_pv, 
                  age_F, age_pv, site_F, site_pv, 
                  EulerNumber_F, EulerNumber_pv, 
                  ICV_F, ICV_pv, NA,   NA,
                  sign_GS,m-2,g-1)
      }
      
      DF = rbindlist(list(DF, rw))
      
    } # for m 
  } #g
} #end i

col_names = c("Model_yvar",
                  "Group_Fval","Group_pval","Sex_Fval","Sex_pval",
                  "Age_Fval","Age_pval","Site_Fval","Site_pval",
                  "EulerNumber_Fval","Eulernumber_pval",
                  "ICV_Fval","ICV_pval","Group_sex_Fval","Group_sex_pval",
                  "Significant_GS_interaction","GS_in_model","global_var_in_model")

names(DF)<-col_names

#df_with_ICV[,2:ncol(df_with_ICV)] <- signif(df_with_ICV[,2:ncol(df_with_ICV)],digits=3)
#df_without_ICV[,2:ncol(df_without_ICV)] <- signif(df_without_ICV[,2:ncol(df_without_ICV)],digits=3)

DF_xlsx_glob0 = DF[DF$global_var_in_model == 0, ]
DF_xlsx_glob1 = DF[DF$global_var_in_model == 1, ]

write_xlsx(DF_xlsx_glob1,"/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/global_variables/globvar_GS_ANOVA_pvals_with_ICV.xlsx")
write_xlsx(DF_xlsx_glob0,"/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/global_variables/globvar_GS_ANOVA_pvals_without_ICV.xlsx")




#### CONTRAST XLSX with male / female split
# and built models[sex][glob][yvar]

model_vars = c("BrainTotalVol", "CortexVol", "total_area", "mean_thickness", "eICV_samseg")
glob = c("with_eICV", "without_eICV")

sex = c("female","male","both")

data_sex1 = datab[c(datab$sex == 1),]
data_sex0 = datab[c(datab$sex == 0),]
dataf = list(data_sex0, data_sex1, datab)

models = list()

DFc = data.frame()

for (k in seq(1,length(model_vars))){
  print(model_vars[k])
  for (s in seq(0,2)){
    ss = s+1
    df_sex = dataf[[ss]]
    
    if (model_vars[k] == "mean_thickness" & sex[ss] != "both"){
      next
    }
    if (model_vars[k] != "mean_thickness" & sex[ss] == "both"){
      next
    }
    print(sex[ss])
    for (g in seq(0,1)){
      gg = g+1
      #print(glob[gg])
      
      if (glob[gg] == "without_eICV"){
        #no global measure model
        if (model_vars[k] == "mean_thickness"){
          f = paste(model_vars[k],"~","group","+","sex","+","age","+","site","+","TotalEulerNumber")
        }
        else {
          f = paste(model_vars[k],"~","group","+","age","+","site","+","TotalEulerNumber")
        }
        models[[sex[ss]]][[glob[gg]]][[model_vars[k]]] = lm(f,data=df_sex)
        print(glob[gg])
        model_ana = models[[sex[ss]]][[glob[gg]]][[model_vars[k]]]
        xvars = attributes(Anova(model_ana,type = "III"))$row.names
      }
      
      else {
        #global measure model
        if (model_vars[k] == "mean_thickness"){
          f2 = paste(model_vars[k],"~","group","+","sex","+","age","+","site","+","TotalEulerNumber","+","eICV_samseg")
        }
        else {
          f2 = paste(model_vars[k],"~","group","+","age","+","site","+","TotalEulerNumber","+","eICV_samseg")
        }        
        models[[sex[ss]]][[glob[gg]]][[model_vars[k]]] = lm(f2,data=df_sex)
        model_ana = models[[sex[ss]]][[glob[gg]]][[model_vars[k]]]
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
      
      
      #rows for contrast xlsx table
      rwc_xlsx = list(model_vars[k],K_emm,
                      BP_diff,BP_diff_tratio,BP_diff_pv,
                      SZ_diff,SZ_diff_tratio,SZ_diff_pv,
                      g,s)
      
      DFc = rbindlist(list(DFc, rwc_xlsx))
      
    } #g 
  } #s
} #k


col_names = c("Model_yvar", "K_LSmean",
              "Contrast_BP-K","tratio_BP-K","pval_BP-K",
              "Contrast_SZ-K","tratio_SZ-K","pval_SZ-K",
              "ICV_in_model","sex")
names(DFc)<-col_names


DF_xlsx_glob0 = DFc[DFc$ICV_in_model == 0, ]
DF_xlsx_glob1 = DFc[DFc$ICV_in_model == 1, ]

write_xlsx(DF_xlsx_glob1,"/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/global_variables/globvar_Model_contrasts_with_ICV.xlsx")
write_xlsx(DF_xlsx_glob0,"/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/global_variables/globvar_Model_contrasts_without_ICV.xlsx")



### plot ppwp
for (i in seq(4,length(model_vars))){
  ps = list()
  emm = list()
  min_emm = list()
  max_emm = list()
  print(model_vars[i])
  
  if (model_vars[i] == "eICV_samseg"){
    glob = c("without_eICV")
    sexplot = c("female","male")
  }
  else if (model_vars[i] == "mean_thickness"){
    glob = c("with_eICV", "without_eICV")
    sexplot = c("both")
  }
  else{
    glob = c("with_eICV", "without_eICV")
    sexplot = c("female","male")
  }
    
    for (g in seq(1,length(glob))){
      for (s in seq(1,length(sexplot))){
        
      emm[[s]] = emmeans(models[[sexplot[s]]][[glob[g]]][[model_vars[i]]],specs = "group")
      plot_idx = s+length(glob)*g-length(glob)
      ps[[plot_idx]] = pwpp(emm[[s]],adjust="none",sort = FALSE) +
        #aes(y = factor("group", level = c("BP","K","SZ"))) +
        labs(x="Uncorrected P-value") +
        ggtitle(paste(glob[g])) +
        geom_vline(xintercept = 0.05,linetype="dashed") +
        coord_flip()
      
      min_emm[[plot_idx]] = min(summary(emm[[s]])$emmean)
      max_emm[[plot_idx]] = max(summary(emm[[s]])$emmean)
      } #end g
    } #end s
    
    for (g in seq(1,length(glob))){
      for (s in seq(1,length(sexplot))){
        emm[[s]] = emmeans(models[[sexplot[s]]][[glob[g]]][[model_vars[i]]],specs = "group")
        plot_idx = s+length(glob)*g-length(glob)
        ps[[plot_idx +4]] = plot(emm[[s]]) +
          aes(color=group) +
          #facet_grid(cols =vars(sex)) +
          ggtitle(paste(sex[s])) +
          scale_x_continuous(limits = c(min(unlist(min_emm)), max(unlist(max_emm)))) +
          scale_y_discrete(limits=c("BP","K","SZ")) +
          coord_flip()
    } #end g
  } #end s
  ps<-ps[!sapply(ps,is.null)]
  top_title = paste(model_vars[i]," models")
  ga=grid.arrange(grobs=ps, ncol=min(length(glob),length(sexplot))*2,top=textGrob(top_title,gp=gpar(fontsize=20)))
  ggsave(paste(model_vars[i],"_group_diff_pvalues_ICV",".png",sep=""),ga,width = 10,height = 10)
} #end i







####### LSmeans estimates on the variance plot preparation #########

y_vars = c("BrainTotalVol", "CortexVol", "total_area","mean_thickness")
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
      sexplot = c("both")
    }
    else{
      dftm$NORM_measure_K[ dftm$NAME_measure  == y_vars[i] ] = ( dftm$VALUE_measure[ dftm$NAME_measure  == y_vars[i] ] - avg_K ) / avg_K
      sexplot = c("female","male")
      }
    
    for (s in seq(1,length(sexplot))){
      
    if (sexplot[s] == "female"){
      ls = lsmeans(models[[sexplot[s]]][[glob[g]]][[y_vars[i]]],pairwise~"group",adjust="none")
      c = ls$contrasts
      
      K_emm = summary(ls)$lsmeans$lsmean[2]
  
      #percent difference
      BP_diff_P = 100*summary(c)$estimate[1] /K_emm
      BP_diff_LCL_P = 100*confint(c)$lower.CL[1] /K_emm
      BP_diff_UCL_P = 100*confint(c)$upper.CL[1] /K_emm
      
      SZ_diff_P = -100*summary(c)$estimate[3] /K_emm
      SZ_diff_LCL_P = -100*confint(c)$lower.CL[3] /K_emm
      SZ_diff_UCL_P = -100*confint(c)$upper.CL[3] /K_emm
      
      ###
      dftm$lsmean_sex0[dftm$NAME_measure  == y_vars[i] & dftm$group == "BP" & dftm$sex == 0] = BP_diff_P
      dftm$eb_max_sex0[dftm$NAME_measure  == y_vars[i] & dftm$group == "BP" & dftm$sex == 0] = BP_diff_UCL_P
      dftm$eb_min_sex0[dftm$NAME_measure  == y_vars[i] & dftm$group == "BP" & dftm$sex == 0] = BP_diff_LCL_P
      
      ###
      dftm$lsmean_sex0[dftm$NAME_measure  == y_vars[i] & dftm$group == "SZ" & dftm$sex == 0] = SZ_diff_P
      dftm$eb_max_sex0[dftm$NAME_measure  == y_vars[i] & dftm$group == "SZ" & dftm$sex == 0] = SZ_diff_UCL_P
      dftm$eb_min_sex0[dftm$NAME_measure  == y_vars[i] & dftm$group == "SZ" & dftm$sex == 0] = SZ_diff_LCL_P
    }
    else if (sexplot[s] == "male"){
      ls = lsmeans(models[[sexplot[s]]][[glob[g]]][[model_vars[i]]],pairwise~"group",adjust="none")
      c = ls$contrasts
      
      sex1_K_emm = summary(ls)$lsmeans$lsmean[2]
      
      sex1_BP_diff_P = 100*summary(c)$estimate[1] /sex1_K_emm
      sex1_BP_diff_LCL_P = 100*confint(c)$lower.CL[1] /sex1_K_emm
      sex1_BP_diff_UCL_P = 100*confint(c)$upper.CL[1] /sex1_K_emm
      
      sex1_SZ_diff_P = -100*summary(c)$estimate[3] /sex1_K_emm
      sex1_SZ_diff_LCL_P = -100*confint(c)$lower.CL[3] /sex1_K_emm
      sex1_SZ_diff_UCL_P = -100*confint(c)$upper.CL[3] /sex1_K_emm
      
      dftm$lsmean_sex1[dftm$NAME_measure  == y_vars[i] & dftm$group == "BP" & dftm$sex == 1] = sex1_BP_diff_P
      dftm$eb_max_sex1[dftm$NAME_measure  == y_vars[i] & dftm$group == "BP" & dftm$sex == 1] = sex1_BP_diff_UCL_P
      dftm$eb_min_sex1[dftm$NAME_measure  == y_vars[i] & dftm$group == "BP" & dftm$sex == 1] = sex1_BP_diff_LCL_P
  
      dftm$lsmean_sex1[dftm$NAME_measure  == y_vars[i] & dftm$group == "SZ" & dftm$sex == 1] = sex1_SZ_diff_P
      dftm$eb_max_sex1[dftm$NAME_measure  == y_vars[i] & dftm$group == "SZ" & dftm$sex == 1] = sex1_SZ_diff_UCL_P
      dftm$eb_min_sex1[dftm$NAME_measure  == y_vars[i] & dftm$group == "SZ" & dftm$sex == 1] = sex1_SZ_diff_LCL_P
    }
      else{
        ls = lsmeans(models[[sexplot[s]]][[glob[g]]][[y_vars[i]]],pairwise~"group",adjust="none")
        c = ls$contrasts
        
        K_emm = summary(ls)$lsmeans$lsmean[2]
        
        #percent difference
        BP_diff_P = 100*summary(c)$estimate[1] /K_emm
        BP_diff_LCL_P = 100*confint(c)$lower.CL[1] /K_emm
        BP_diff_UCL_P = 100*confint(c)$upper.CL[1] /K_emm
        
        SZ_diff_P = -100*summary(c)$estimate[3] /K_emm
        SZ_diff_LCL_P = -100*confint(c)$lower.CL[3] /K_emm
        SZ_diff_UCL_P = -100*confint(c)$upper.CL[3] /K_emm
        
        ###
        dftm$lsmean[dftm$NAME_measure  == y_vars[i] & dftm$group == "BP"] = BP_diff_P
        dftm$eb_max[dftm$NAME_measure  == y_vars[i] & dftm$group == "BP"] = BP_diff_UCL_P
        dftm$eb_min[dftm$NAME_measure  == y_vars[i] & dftm$group == "BP"] = BP_diff_LCL_P
        
        ###
        dftm$lsmean[dftm$NAME_measure  == y_vars[i] & dftm$group == "SZ"] = SZ_diff_P
        dftm$eb_max[dftm$NAME_measure  == y_vars[i] & dftm$group == "SZ"] = SZ_diff_UCL_P
        dftm$eb_min[dftm$NAME_measure  == y_vars[i] & dftm$group == "SZ"] = SZ_diff_LCL_P
    } #end if sex
    } #end for sex
  } #end for i 
  
  if (glob[g] == "with_eICV"){
    dfc_plot_with_eICV = dftm
  }
  else{
    dfc_plot_without_eICV = dftm
  }
  
}



### PLOT THE LSMEANS 


for (g in seq(1,length(glob))){
  DF = get(paste("dfc_plot_",glob[g],sep = ""))
  
  for (s in seq(1,2)){
    pls[[s]]= with(DF[!is.na(dftm$NORM_measure_K), ],
                   ggplot() + 
                     aes(x = factor(NAME_measure), y = 100*NORM_measure_K, color = group) +
                     geom_violin(position = "identity",alpha=0.3) +
                     geom_jitter(width = 0.3, size=0.1) + 
                     geom_hline(yintercept = 0,linetype="dashed") + 
                     
                     geom_point(position=pd, aes_string(x = "NAME_measure", y = paste("lsmean",ss[s],sep = "_"), color="group", group="group"), size=2) +

                     stat_summary(fun = mean,geom = "line",position=pd, aes_string(x = "NAME_measure", y = paste("lsmean",ss[s],sep = "_"), color="group", group="group")) + 
                     
                     geom_errorbar(position=pd, width=0.2, aes_string(group="group",color="group",ymin=paste("eb_min_",ss[s],sep=""), ymax=paste("eb_max_",ss[s],sep="")) ) +
                     
                     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
                     labs(y = "Difference from control [%]", x = "brain measure") +
                     ggtitle(paste(sex_name[s]))
                     
                     #annotate("text", size=2.5, x = c(1,2,3), y=21, label = c(paste("G:S pval =",round(GS_pvals[[glob[g]]][[y_vars[1]]],digits=3)),
                     #                                                        paste("G:S pval =",round(GS_pvals[[glob[g]]][[y_vars[2]]],digits=3)),
                     #                                                         paste("G:S pval =",round(GS_pvals[[glob[g]]][[y_vars[3]]],digits=3))) )
    )
  } #s = for each sex

    pls[[3]]= with(DF[!is.na(dftm$NORM_measure_K_mt), ],
                   ggplot() + 
                     aes(x = factor(NAME_measure), y = 100*NORM_measure_K_mt, color = group) +
                     geom_violin(position = "identity",alpha=0.3) +
                     geom_jitter(width = 0.3, size=0.1) + 
                     geom_hline(yintercept = 0,linetype="dashed") + 
                     
                     geom_point(position=pd, aes_string(x = "NAME_measure", y = "lsmean", color=group, group=group), size=2) +
                     geom_errorbar(position=pd, aes_string(width=0.1,group="group",color="group",ymin="eb_min", ymax="eb_max") ) +
                     
                     #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
                     labs(y = "Difference from control [%]", x = "brain measure") +
                     ggtitle("Both sex") 
                     
                     #annotate("text",size=2.5, x = 1, y=7.9, label = c(paste("G:S pval =",round(GS_pvals[[glob[g]]][[y_vars[4]]],digits=3))))
    )
    pls[[4]] = ggplot() + theme_void()
    
  top_title = paste("LSmeans contrasts of model",glob[g])
  ga=grid.arrange(grobs=pls,ncols=2, top=textGrob(top_title,gp=gpar(fontsize=20)))
  ggsave(paste("LSmeans_contrasts_on_datapoints_",glob[g],".png",sep=""),ga,width = 10,height = 10)
  
} #g








