
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


#make new groups according to axis1
data_no_na = datab[!is.na(datab$diag), ]
data_SZ = data_no_na
data_SZ = data_SZ[!data_SZ$group == "BP", ]

data_SZ$new_axis1_group = "K"
data_SZ$new_axis1_group[data_SZ$diag == 0 & data_SZ$group == "SZ"] = "SZ_no_axis1"
data_SZ$new_axis1_group[data_SZ$diag == 1 & data_SZ$group == "SZ"] = "SZ_with_axis1"

data_BP = data_no_na
data_BP = data_BP[!data_BP$group == "SZ", ]

data_BP$new_axis1_group = "K"
data_BP$new_axis1_group[data_BP$diag == 0 & data_BP$group == "BP"] = "BP_no_axis1"
data_BP$new_axis1_group[data_BP$diag == 1 & data_BP$group == "BP"] = "BP_with_axis1"


#see how many subjects are in each group
sum(data_no_na$group == "K")
sum(data_no_na$group == "BP")
sum(data_no_na$group == "SZ")

sum(data_BP$new_axis1_group == "K")
sum(data_BP$new_axis1_group == "BP_no_axis1")
sum(data_BP$new_axis1_group == "BP_with_axis1")

sum(data_SZ$new_axis1_group == "K")
sum(data_SZ$new_axis1_group == "SZ_no_axis1")
sum(data_SZ$new_axis1_group == "SZ_with_axis1")


setwd("/mnt/projects/VIA11/FREESURFER/Stats/Plots/global_measures/axis1")



##############################
# whole brain measure models #
##############################

#WITHOUT global variable models 
model_bvol = lm(BrainTotalVol ~ new_axis1_group*sex + age + site, data=data_SZ)
anova(model_bvol)
summary(model_bvol)

model_cvol = lm(CortexVol ~ new_axis1_group*sex + age + site, data=data_SZ)
anova(model_cvol)
#lsmeans(model_cvol,pairwise~"group", by = "sex")

model_area = lm(total_area ~ new_axis1_group*sex + age + site, data=data_SZ)
anova(model_area)

model_mt = lm(mean_thickness ~ new_axis1_group*sex  + age + site, data=data_SZ)
anova(model_mt)
model_mt = update(model_mt,~.-new_axis1_group:sex)
anova(model_mt)
#lsmeans(model_mt,pairwise~"group")

model_icv = lm(eICV_samseg ~ new_axis1_group*sex + age + site, data=data_SZ)
anova(model_icv)



#WITH global variable covariate 
model_bvol_glob = lm(BrainTotalVol ~ new_axis1_group*sex + age + eICV_samseg + site, data=data_SZ)
anova(model_bvol_glob)
model_bvol_glob = update(model_bvol_glob,~.-new_axis1_group:sex)
anova(model_bvol_glob)

model_cvol_glob = lm(CortexVol ~ new_axis1_group*sex + age + eICV_samseg + site, data=data_SZ)
anova(model_cvol_glob)
model_cvol_glob = update(model_cvol_glob,~.-new_axis1_group:sex)
anova(model_cvol_glob)

model_area_glob = lm(total_area ~ new_axis1_group*sex + age + eICV_samseg + site, data=data_SZ)
anova(model_area_glob)
model_area_glob = update(model_area_glob,~.-new_axis1_group:sex)
anova(model_area_glob)

model_mt_glob = lm(mean_thickness ~ new_axis1_group*sex + age + eICV_samseg + site, data=data_SZ)
anova(model_mt_glob)
model_mt_glob = update(model_mt_glob,~.-new_axis1_group:sex)
anova(model_mt_glob)
#lsmeans(model_mt_glob,pairwise~"group")





### run models with eICV_samseg as covariate and without

model_vars = c("BrainTotalVol", "CortexVol", "total_area","mean_thickness")
glob = c("with_eICV", "without_eICV")
groups_df = c("BP", "SZ")
GS_pvals = list()
ps = list()
emm = list()
min_emm = list()
max_emm = list()
model = list()

for (g2 in seq(1,length(groups_df))){
for (i in seq(1,length(model_vars))){
  for (g in seq(1,2)){
    if (glob[g] == "with_eICV"){
      f = paste(model_vars[i],"~","+","new_axis1_group*sex","+","age","+","site","+","TotalEulerNumber","+","eICV_samseg")
    }
    else{
      f = paste(model_vars[i],"~","+","new_axis1_group*sex","+","age","+","site","+","TotalEulerNumber")
    }
    model[[glob[g]]][[groups_df[g2]]][[i]] = lm(f,data=get(paste("data_",groups_df[g2],sep="")))
    model_ana = model[[glob[g]]][[groups_df[g2]]][[i]]
    xvars = attributes(anova(model_ana))$row.names
    
    GS_pvals[[glob[g]]][[groups_df[g2]]][[model_vars[i]]] = anova(model_ana)$"Pr(>F)"[xvars=="new_axis1_group:sex"]
    
    if (anova(model_ana)$"Pr(>F)"[xvars=="new_axis1_group:sex"] > 0.05){
      model[[glob[g]]][[groups_df[g2]]][[i]] = update(model[[glob[g]]][[groups_df[g2]]][[i]],~.-new_axis1_group:sex)
      model_ana = update(model_ana,~.-group:sex)
    } 
    
    emm[[g]] = emmeans(model_ana,specs = "new_axis1_group",by="sex")
    ps[[g]] = pwpp(emm[[g]],adjust="none",sort = FALSE) +
      labs(x="Uncorrected P-value") +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      ggtitle(paste(glob[g])) +
      geom_vline(xintercept = 0.05,linetype="dashed") +
      coord_flip()
    
    min_emm[[g]] = min(summary(emm[[g]])$emmean)
    max_emm[[g]] = max(summary(emm[[g]])$emmean)
  } #g
  
  for (g in seq(1,2)){
    ps[[g+2]] = plot(emm[[g]],comparisons=TRUE,adjust="none") +
      aes(color=group) +
      facet_grid(cols =vars(sex)) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      
      scale_x_continuous(limits = c(min(unlist(min_emm)), max(unlist(max_emm)))) +
      
      {if (groups_df[g2] == "BP") scale_y_discrete(limits=c(paste(groups_df[g2],"_no_axis1",sep = ""),paste(groups_df[g2],"_with_axis1",sep = ""),"K"))} +
      {if (groups_df[g2] == "SZ") scale_y_discrete(limits=c("K",paste(groups_df[g2],"_no_axis1",sep = ""),paste(groups_df[g2],"_with_axis1",sep = "")))} +
      
      coord_flip()
  }
  
  top_title = paste(model_vars[i]," models for ",groups_df[g2])
  ga=grid.arrange(grobs=ps, ncol=2,top=textGrob(top_title,gp=gpar(fontsize=20)))
  ggsave(paste(model_vars[i],"_group_diff_pvalues_ICV_axis1_",groups_df[g2],".png",sep=""),ga,width = 10,height = 10)
} #i
  } #g2






####### LSmeans estimates on the variance plot #########

y_vars = c("BrainTotalVol", "CortexVol", "total_area", "mean_thickness")
glob = c("with_eICV", "without_eICV")
groups_df = c("BP", "SZ")

axis1 = c("no_axis1","with_axis1")

pls = list()

g = 2

#for (g in seq(1, 2)){

for (g2 in seq(1,length(groups_df))){
  
dftm = get(paste("data_",groups_df[g2],sep="")) %>%
  pivot_longer(cols = y_vars,
               names_to = "NAME_measure",
               values_to = "VALUE_measure",
               values_drop_na = FALSE)
dftm = as.data.frame(dftm)

dftm$NORM_measure = dftm$VALUE_measure
dftm$NORM_measure_K = NA #dftm$VALUE_measure
dftm$NORM_measure_K_mt = NA #dftm$VALUE_measure

  
  for (i in seq(1,length(y_vars))){
    avg_K = mean(dftm$VALUE_measure[ dftm$NAME_measure  == y_vars[i] & dftm$group  == "K"])
    
    if (y_vars[i] == "mean_thickness"){
      dftm$NORM_measure_K_mt[ dftm$NAME_measure  == y_vars[i] ] = ( dftm$VALUE_measure[ dftm$NAME_measure  == y_vars[i] ] - avg_K ) / avg_K
    }
    else{
      dftm$NORM_measure_K[ dftm$NAME_measure  == y_vars[i] ] = ( dftm$VALUE_measure[ dftm$NAME_measure  == y_vars[i] ] - avg_K ) / avg_K
    }
    
    ls = lsmeans(model[[glob[g]]][[groups_df[g2]]][[i]],pairwise~"new_axis1_group",by="sex",adjust="none")
    c = ls$contrasts
    K_emm = summary(ls)$lsmeans$lsmean[summary(ls$lsmeans)$sex == 0][2]
    sex1_K_emm = summary(ls)$lsmeans$lsmean[summary(ls$lsmeans)$sex == 1][2]
    
    
    #percent difference
    diff_noaxis1_P = -100*summary(c)$estimate[summary(c)$sex == 0][1] /K_emm
    diff_noaxis1_LCL_P = -100*confint(c)$lower.CL[confint(c)$sex == 0][1] /K_emm
    diff_noaxis1_UCL_P = -100*confint(c)$upper.CL[confint(c)$sex == 0][1] /K_emm
    
    diff_axis1_P = -100*summary(c)$estimate[summary(c)$sex == 0][2] /K_emm
    diff_axis1_LCL_P = -100*confint(c)$lower.CL[confint(c)$sex == 0][2] /K_emm
    diff_axis1_UCL_P = -100*confint(c)$upper.CL[confint(c)$sex == 0][2] /K_emm
    
    sex1_axis1_diff_P = -100*summary(c)$estimate[summary(c)$sex == 1][2] /sex1_K_emm
    sex1_axis1_diff_LCL_P = -100*confint(c)$lower.CL[confint(c)$sex == 1][2] /sex1_K_emm
    sex1_axis1_diff_UCL_P = -100*confint(c)$upper.CL[confint(c)$sex == 1][2] /sex1_K_emm
    
    sex1_noaxis1_diff_P = -100*summary(c)$estimate[summary(c)$sex == 1][1] /sex1_K_emm
    sex1_noaxis1_diff_LCL_P = -100*confint(c)$lower.CL[confint(c)$sex == 1][1] /sex1_K_emm
    sex1_noaxis1_diff_UCL_P = -100*confint(c)$upper.CL[confint(c)$sex == 1][1] /sex1_K_emm
    
    
    #for the particular dataframe
    dftm$lsmean[dftm$NAME_measure  == y_vars[i]] = diff_P
    dftm$lsmean[dftm$NAME_measure  == y_vars[i] & dftm$sex == 1] = sex1_diff_P
    
    dftm$eb_max[dftm$NAME_measure  == y_vars[i]] = diff_UCL_P
    dftm$eb_min[dftm$NAME_measure  == y_vars[i]] = diff_LCL_P
    
    dftm$eb_max[dftm$NAME_measure  == y_vars[i] & dftm$sex == 1] = sex1_diff_UCL_P
    dftm$eb_min[dftm$NAME_measure  == y_vars[i] & dftm$sex == 1] = sex1_diff_LCL_P
    
    #
    dftm$lsmean[dftm$NAME_measure  == y_vars[i]] = diff_P
    dftm$lsmean[dftm$NAME_measure  == y_vars[i] & dftm$sex == 1] = sex1_diff_P
    
    dftm$eb_max[dftm$NAME_measure  == y_vars[i]] = diff_UCL_P
    dftm$eb_min[dftm$NAME_measure  == y_vars[i]] = diff_LCL_P
    
    dftm$eb_max[dftm$NAME_measure  == y_vars[i] & dftm$sex == 1] = sex1_diff_UCL_P
    dftm$eb_min[dftm$NAME_measure  == y_vars[i] & dftm$sex == 1] = sex1_diff_LCL_P
    
  }

  for (a1 in seq(1,2)){
    
    
    
    pls[[2*a1-1]]= with(dftm[!is.na(dftm$NORM_measure_K), ],
                       ggplot() + 
                         aes(x = factor(NAME_measure), y = 100*NORM_measure_K, color = new_axis1_group) +
                         geom_violin(position = "identity",alpha=0.3) +
                         geom_jitter(width = 0.3, size=0.1) + 
                         geom_hline(yintercept = 0,linetype="dashed") + 
                         
                         ggnewscale::new_scale_colour() +
                         
                         geom_point(aes_string(x = "NAME_measure", y = "lsmean", color=sex, group=sex), size=2) +
                         scale_colour_manual("Contrasts per sex", values = c("red", "blue")) + 
                         
                         geom_errorbar(aes_string(width=0.05,group="sex",color="sex",ymin="eb_min", ymax="eb_max") ) +
                         geom_line(aes_string(x = "NAME_measure", y = "lsmean", color="sex", group="sex")) +
                         
                         theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
                         labs(y = "Difference from control [%]", x = "brain measure") +
                         ggtitle(paste(axis1[a1])) +
                         
                         annotate("text", size=2.5, x = c(1,2,3), y=21, label = c(paste("G:S pval =",round(GS_pvals[[glob[g]]][[groups_df[g2]]][[y_vars[1]]],digits=3)),
                                                                        paste("G:S pval =",round(GS_pvals[[glob[g]]][[groups_df[g2]]][[y_vars[2]]],digits=3)),
                                                                        paste("G:S pval =",round(GS_pvals[[glob[g]]][[groups_df[g2]]][[y_vars[3]]],digits=3))) )
    )
    
    pls[[2*a1]]= with(dftm[!is.na(dftm$NORM_measure_K_mt), ],
                     ggplot() + 
                       aes(x = factor(NAME_measure), y = 100*NORM_measure_K_mt, color = new_axis1_group) +
                       geom_violin(position = "identity",alpha=0.3) +
                       geom_jitter(width = 0.3, size=0.1) + 
                       geom_hline(yintercept = 0,linetype="dashed") + 
                       
                       ggnewscale::new_scale_colour() +
                       
                       geom_point(aes_string(x = "NAME_measure", y ="lsmean", color="sex", group="sex"), size=2) +
                       scale_colour_manual("Contrasts per sex", values = c("red", "blue")) + 
                       
                       geom_errorbar(aes_string(width=0.05,group="sex",color="sex",ymin="eb_min", ymax="eb_max") ) +
                       geom_line(aes_string(x = "NAME_measure", y = "lsmean", color="sex", group="sex")) +
                       
                       labs(y = "Difference from control [%]", x = "brain measure") +
                       ggtitle(paste(axis1[a1])) + 
                       
                       annotate("text", size=2.5, x = 1, y=7.5, label = c(paste("G:S pval =",round(GS_pvals[[glob[g]]][[groups_df[g2]]][[y_vars[4]]],digits=3))))
    )
  } #a1
top_title = paste("LSmeans contrasts of model",glob[g],groups_df[g2])
ga=grid.arrange(grobs=pls,ncols=2, top=textGrob(top_title,gp=gpar(fontsize=20)))
ggsave(paste("LSmeans_contrasts_on_datapoints_axis1",glob[g],groups_df[g2],".png",sep=""),ga,width = 10,height = 10)

} #g2


#} #g






