

### THIS IS AN ARCHIVE SCRIPT

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

data_csv_filtered <- data_csv_filtered[c(data_csv_filtered$Include_FS_studies_euler_outliers_excluded == 1),]
data_csv_filtered <- data_csv_filtered[!is.na(data_csv_filtered$Include_FS_studies_euler_outliers_excluded),]


#tell r which variables are factors
datab = data_csv_filtered
datab$sex = as.factor(datab$Sex_child)
datab$group = as.factor(datab$HighRiskStatus)
datab$site = as.factor(datab$MR_Site)
datab$diag = as.factor(datab$Axis.1_diag_v11)
datab$age = as.numeric(datab$MRI_age)


setwd("/mnt/projects/VIA11/FREESURFER/Stats/Plots/global_measures")



##############################
# whole brain measure models #
##############################


### test some things 

#WITH global variable covariate 
data_BP_K = datab[!datab$group=="SZ",]
model_bvol_glob = lm(BrainTotalVol ~ group*sex + age + TotalEulerNumber + eICV_samseg + site, data=data_BP_K)
Anova(model_bvol_glob,type="III")
#drop1(model_bvol_glob, test="F")
model_bvol_glob = update(model_bvol_glob,~.-group:sex)
Anova(model_bvol_glob,type="III")
lsmeans(model_bvol_glob,pairwise~"group",adjust="none")
#lsmeans(model_bvol_glob,pairwise~"group",by="sex",adjust="none")

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




### Run models with eICV_samseg as covariate and without

model_vars = c("BrainTotalVol", "CortexVol", "total_area","mean_thickness", "eICV_samseg")
GS_pvals = list()
ps = list()
emm = list()
min_emm = list()
max_emm = list()
model = list()

DF = data.frame()


for (i in seq(1,length(model_vars))){
  
  #Indicate which models should have ICV as covariate 
  if (model_vars[i] == "eICV_samseg"){
    glob = c("without_eICV")
  }
  else{
    glob = c("with_eICV", "without_eICV")
  }

  for (g in seq(1,length(glob))){
    
    #Define models with ICV and without 
    if (glob[g] == "with_eICV"){
      f = paste(model_vars[i],"~","+","group*sex","+","age","+","site","+","TotalEulerNumber","+","eICV_samseg")
    }
    else{
      f = paste(model_vars[i],"~","+","group*sex","+","age","+","site","+","TotalEulerNumber")
    }
    #Run the actual lm model 
    model[[glob[g]]][[model_vars[i]]] = lm(f,data=datab)
    xvars = attributes(Anova(model[[glob[g]]][[model_vars[i]]],type = "III"))$row.names
    
    
    if (Anova(model[[glob[g]]][[model_vars[i]]],type = "III")$"Pr(>F)"[xvars=="group:sex"] > 0.05){
      #Note the significance of Group/sex interaction
      sign_GS = 0
      
      #Save the model which includes the Group/sex interaction
      model_gs = model[[glob[g]]][[model_vars[i]]]
      
      #Compute G/S interaction statistics
      GS_pvals[[glob[g]]][[model_vars[i]]] = Anova(model_gs,type = "III")$"Pr(>F)"[xvars=="group:sex"]
      GS_F = Anova(model_gs,type = "III")$"F value"[xvars=="group:sex"]
      
      #Remove the interaction from the main model object, i.e. model[glob][model_var]
      model[[glob[g]]][[model_vars[i]]] = update(model[[glob[g]]][[model_vars[i]]],~.-group:sex)
      
      #save both model with G/S interaction and model without in a list
      models = list(model_gs, model[[glob[g]]][[model_vars[i]]])
      model_type = list("GS","No_interactions")
    }
    else{
      #if the G/S interaction is significant, just note the G/S interaction statistics
      GS_pvals[[glob[g]]][[model_vars[i]]] = Anova(model[[glob[g]]][[model_vars[i]]],type = "III")$"Pr(>F)"[xvars=="group:sex"]
      GS_F = Anova(model[[glob[g]]][[model_vars[i]]],type = "III")$"F value"[xvars=="group:sex"]
      sign_GS = 1
      models = list(model[[glob[g]]][[model_vars[i]]]) # only add the model with G/S interaction to the list
      model_type = list("GS")
    }
    
    #loop over either one or two models (the first model is always the model with interaction)
    for (m in seq(1,length(models))){
      #get the model
      mm = models[[m]]
      subs = length(mm$residuals)
      
      #compute anova statistics for the model
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
      
      #the first index in m is always the model with a group/sex interaction
      if (model_type[[m]]=="GS"){
        GS_in = 1
        GS_pv = GS_pvals[[glob[g]]][[model_vars[i]]]
      }
      else if (model_type[[m]]=="No_interactions"){ 
        GS_F = NA
        GS_pv = NA
        GS_in = 0
      }
      

      #append the row to the dataframe 
      DF = rbindlist(list(DF, as.list(c(model_vars[i],
                                        group_F, group_pv, sex_F, sex_pv,
                                        age_F, age_pv, site_F, site_pv, 
                                        EulerNumber_F, EulerNumber_pv, 
                                        ICV_F, ICV_pv,
                                        GS_F, GS_pv, 
                                        sign_GS, GS_in, 
                                        g-1, model_type[[m]], subs))))

    } # end for m in models
    
    #
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


#define the coloum names of the dataframe which is converted to an excel 
col_names = c("Model_yvar",
              "Group_Fval","Group_pval","Sex_Fval","Sex_pval",
              "Age_Fval","Age_pval","Site_Fval","Site_pval",
              "EulerNumber_Fval","Eulernumber_pval",
              "ICV_Fval","ICV_pval",
              "Group_sex_Fval","Group_sex_pval",
              "Significant_GS_interaction","GS_in_model",
              "ICV_in_model","model_type","n_subjects")

names(DF)<-col_names

DF_xlsx_glob0 = DF[DF$ICV_in_model == 0, ]
DF_xlsx_glob1 = DF[DF$ICV_in_model == 1, ]

df_with_ICV[,2:ncol(df_with_ICV)] <- signif(df_with_ICV[,2:ncol(df_with_ICV)],digits=3)
df_without_ICV[,2:ncol(df_without_ICV)] <- signif(df_without_ICV[,2:ncol(df_without_ICV)],digits=3)

write_xlsx(df_with_ICV,"/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/globvar_ANOVA_pvals_with_ICV.xlsx")
write_xlsx(df_without_ICV,"/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/globvar_ANOVA_pvals_without_ICV.xlsx")











###################################
#######     ARCHIVE     ###########
###################################

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




