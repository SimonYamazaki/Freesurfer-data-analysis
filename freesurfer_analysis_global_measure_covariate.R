
setwd("/mnt/projects/VIA11/FREESURFER/Stats/Data")

#load data
data_csv <- read.table("VIA11_allkey_160621_FreeSurfer_pruned_20211210.csv", header = TRUE, sep = ",", dec = ".")

#inspect the head of data and summary
head(data_csv)
summary(data_csv)

#filter the data with include variable
data_csv_filtered <- na.omit(data_csv[c(data_csv$Include_FS_studies == 1),])

#tell r which variables are factors
datab = data_csv_filtered
datab$sex = as.factor(datab$Sex_child)
datab$group = as.factor(datab$HighRiskStatus)
datab$site = as.factor(datab$MR_Site)
datab$diag = as.factor(datab$Axis.1_diag_v11)
datab$age = as.numeric(datab$MRI_age)


#4 global variable models

library(lmerTest)
library(lme4)
model_bvol = lmer(BrainTotalVol ~ group*sex + age + eICV_samseg + (1| site), data=datab)
ranova(model_bvol)
anova(model_bvol)

model_cvol = lmer(CortexVol ~ group*sex + age + eICV_samseg + (1| site), data=datab)
ranova(model_cvol)
anova(model_cvol)

model_area = lmer(total_area ~ group*sex + age + eICV_samseg + (1| site), data=datab)
ranova(model_area)
anova(model_area)

model_mt = lmer(mean_thickness ~ group + sex + age + eICV_samseg + (1| site), data=datab)
ranova(model_mt)
anova(model_mt)



### run 4 models with eICV_samseg as covariate
model_vars = c("BrainTotalVol", "CortexVol", "total_area")
model = list()
for (i in seq(1,length(model_vars))){
  f = paste(model_vars[i],"~","-1","+","group*sex","+","age","+","eICV_samseg","+","(1| site)")
  model[[i]] = lmer(f,data=datab)
}
model[[i+1]] = lmer(mean_thickness ~ -1 + group + sex + age + eICV_samseg + (1| site),data=datab)
model_vars = c(model_vars,"mean_thickness")


#####################
#run multiple models GROUP:AGE INTERACTION

model_vars = c("BrainTotalVol", "CortexVol", "total_area", "eICV_samseg")
model = list()
for (i in seq(1,length(model_vars))){
  f = paste(model_vars[i],"~","-1","+","group*sex","+","group:age","+","(1| site)")
  model[[i]] = lmer(f,data=datab)
}
model[[i+1]] = lmer(mean_thickness ~ -1 + group + sex + group:age + (1| site),data=datab)
model_vars = c(model_vars,"mean_thickness")


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


### model diagnostics 
par(mfrow=c(1,3))
for (i in seq(1,length(model_vars))){
  mixed_res = rstudent(model[[i]])
  qqnorm(mixed_res,main = NULL)
  qqline(mixed_res,main = NULL)
  plot(mixed_res ~ fitted(model[[i]]),xlab="Fitted",ylab="Standardized residuals")
  plot(cooks.distance(model[[i]]), type = "p", pch = 20,ylab="Cooks distance")
  
  mtext(paste("Model diagnostics for",model_vars[[i]]), side = 3, line = -2, outer = TRUE)
}


######################################
#    run models on each region
######################################

# extract model yvars 
library(tidyr)
library(dplyr)
library(stringr)

col_names = names(datab)
col_names_list = list()
model_names = list()
h = c("lh","rh")
m = c("area","thickness","volume")
for (i in seq(1,2)){
  for (j in seq(1,3)){
    name = paste(h[i],m[j],sep="_")
    col_names_list[[name]] = col_names
    col_names_list[[name]] = col_names_list[[name]][grepl(paste("^",h[i],sep = ""), col_names_list[[name]])]
    col_names_list[[name]] = col_names_list[[name]][grepl(paste(m[j],'$',sep = ""), col_names_list[[name]])]
    col_names_list[[name]] = col_names_list[[name]][!grepl(paste("Area_area",'$',sep = ""), col_names_list[[name]])] #remove WhiteSurfaceArea
    model_names = c(model_names, col_names_list[[name]])
  }
}

col_names_list



###### make inference with models

models = list()

mm = lmer(rh_insula_volume ~ -1 + group*sex + age + total_area +(1| site),data=datab)
anova(mm)

for (i in seq(1,length(h))){
  for (j in seq(1,length(m))){
    model_yvars = col_names_list[[paste(h[i],m[j],sep = "_")]] 
    for (k in seq(1,length(model_yvars))){

      if (m[j] == "area"){
        f = paste(model_yvars[k],"~","-1","+","group*sex","+","age","+","total_area","+","(1| site)")
      }
      else if (m[j] == "thickness"){
        f = paste(model_yvars[k],"~","-1","+","group*sex","+","age","+","mean_thickness","+","(1| site)")
      }
      else if (m[j] == "volume"){
        f = paste(model_yvars[k],"~","-1","+","group*sex","+","age","+","BrainTotalVol","+","(1| site)")
      }
      
      models[[h[i]]][[m[i]]][[model_yvars[k]]] = lmer(f,data=datab)
      pv_group_sex = anova(models[[h[i]]][[m[i]]][[model_yvars[k]]])$"Pr(>F)"[5]
      
      if (anova(models[[h[i]]][[m[i]]][[model_yvars[k]]])$"Pr(>F)"[5] > 0.05){
        models[[h[i]]][[m[i]]][[model_yvars[k]]] = update(models[[h[i]]][[m[i]]][[model_yvars[k]]],~.-group:sex)
        ls = lsmeans(models[[h[i]]][[m[i]]][[model_yvars[k]]],pairwise~"group")
        c = ls$contrasts
        
        K_emm = summary(ls)$lsmeans$lsmean[2]
        sex1_K_emm = "n/a"
        
        BP_diff = summary(c)$estimate[1] /K_emm
        BP_diff_LCL = confint(c)$lower.CL[1] /K_emm
        BP_diff_UCL = confint(c)$upper.CL[1] /K_emm
        SZ_diff = -summary(c)$estimate[3] /K_emm
        SZ_diff_LCL = -confint(c)$lower.CL[3] /K_emm
        SZ_diff_UCL = -confint(c)$upper.CL[3] /K_emm
        
        sex1_BP_diff = "n/a"
        sex1_BP_diff_LCL = "n/a"
        sex1_BP_diff_UCL = "n/a"
        
        sex1_SZ_diff = "n/a"
        sex1_SZ_diff_LCL = "n/a"
        sex1_SZ_diff_UCL = "n/a"
        
      }
      else{
        ls = lsmeans(models[[h[i]]][[m[i]]][[model_yvars[k]]],pairwise~"group", by = "sex")
        c = ls$contrasts
        
        K_emm = summary(ls)$lsmeans$lsmean[summary(ls$lsmeans)$sex == 0][2]
        sex1_K_emm = summary(ls)$lsmeans$lsmean[summary(ls$lsmeans)$sex == 1][2]
        
        BP_diff = summary(c)$estimate[summary(c)$sex == 0][1] /K_emm
        BP_diff_LCL = confint(c)$lower.CL[confint(c)$sex == 0][1] /K_emm
        BP_diff_UCL = confint(c)$upper.CL[confint(c)$sex == 0][1] /K_emm
        
        SZ_diff = -summary(c)$estimate[summary(c)$sex == 0][3] /K_emm
        SZ_diff_LCL = -confint(c)$lower.CL[confint(c)$sex == 0][3] /K_emm
        SZ_diff_UCL = -confint(c)$upper.CL[confint(c)$sex == 0][3] /K_emm
        
        sex1_BP_diff = summary(c)$estimate[summary(c)$sex == 1][1] /sex1_K_emm
        sex1_BP_diff_LCL = confint(c)$lower.CL[confint(c)$sex == 1][1] /sex1_K_emm
        sex1_BP_diff_UCL = confint(c)$upper.CL[confint(c)$sex == 1][1] /sex1_K_emm
        
        sex1_SZ_diff = -summary(c)$estimate[summary(c)$sex == 1][3] /sex1_K_emm
        sex1_SZ_diff_LCL = -confint(c)$lower.CL[confint(c)$sex == 1][3] /sex1_K_emm
        sex1_SZ_diff_UCL = -confint(c)$upper.CL[confint(c)$sex == 1][3] /sex1_K_emm
        
      }
      
      
      pv_group = anova(models[[h[i]]][[m[i]]][[model_yvars[k]]])$"Pr(>F)"[1]
      
      
      if (i==1 && j==1 && k==1){
        df<-data.frame(model_yvars[k],pv_group,pv_group_sex, K_emm, sex1_K_emm, BP_diff, BP_diff_LCL, BP_diff_UCL, SZ_diff, SZ_diff_LCL, SZ_diff_UCL, sex1_BP_diff, sex1_BP_diff_LCL, sex1_BP_diff_UCL, sex1_SZ_diff, sex1_SZ_diff_LCL, sex1_SZ_diff_UCL)
        names(df)<-c("model_yvar","Group_p_value","Group_sex_p_value","K_emm", "sex1_K_emm","BP_diff","BP_diff_LCL","BP_diff_UCL","SZ_diff","SZ_diff_LCL","SZ_diff_UCL","sex1_BP_diff","sex1_BP_diff_LCL","sex1_BP_diff_UCL","sex1_SZ_diff","sex1_SZ_diff_LCL","sex1_SZ_diff_UCL")
        
        dfs<-data.frame(model_yvars[k], K_emm, sex1_K_emm, BP_diff, BP_diff_LCL, BP_diff_UCL, SZ_diff, SZ_diff_LCL, SZ_diff_UCL,0)
        names(dfs)<-c("model_yvar","K_emm", "sex1_K_emm","BP_diff","BP_diff_LCL","BP_diff_UCL","SZ_diff","SZ_diff_LCL","SZ_diff_UCL","sex")
        dfs[nrow(dfs) + 1,] = c(model_yvars[k], K_emm, sex1_K_emm, sex1_BP_diff, sex1_BP_diff_LCL, sex1_BP_diff_UCL, sex1_SZ_diff, sex1_SZ_diff_LCL, sex1_SZ_diff_UCL,1)
      }
      else{
        df[nrow(df) + 1,] = c(model_yvars[k],pv_group,pv_group_sex, K_emm, sex1_K_emm, BP_diff, BP_diff_LCL, BP_diff_UCL, SZ_diff, SZ_diff_LCL, SZ_diff_UCL, sex1_BP_diff, sex1_BP_diff_LCL, sex1_BP_diff_UCL, sex1_SZ_diff, sex1_SZ_diff_LCL, sex1_SZ_diff_UCL)
        dfs[nrow(dfs) + 1,] = c(model_yvars[k], K_emm, sex1_K_emm, BP_diff, BP_diff_LCL, BP_diff_UCL, SZ_diff, SZ_diff_LCL, SZ_diff_UCL,0)
        dfs[nrow(dfs) + 1,] = c(model_yvars[k], K_emm, sex1_K_emm, sex1_BP_diff, sex1_BP_diff_LCL, sex1_BP_diff_UCL, sex1_SZ_diff, sex1_SZ_diff_LCL, sex1_SZ_diff_UCL,1)
        
      }
      
      
    } #k
  }
}

#plot anova group effect
p = list()

for (i in seq(1,3)){
  for (j in seq(1,2)){
    df1 = df[grepl(paste("^",h[j],sep = ""), df$model_yvar),]
    df1 = df1[grepl(paste("_",m[i],sep = ""), df1$model_yvar),]
    
    p[[i*2+j-2]] =ggplot(df1, aes(as.factor( model_yvar ), as.numeric(Group_p_value))) + 
      geom_point() + 
      geom_hline(yintercept = 0.05) + 
      labs(y = "Group anova p-value", x = "brain region") +
      ggtitle(paste(m[i],h[j]) ) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  }
}

library(gridExtra)
grid.arrange(grobs=p)


#group sex interaction
p = list()
for (i in seq(1,3)){
  for (j in seq(1,2)){
    df1 = df[grepl(paste("^",h[j],sep = ""), df$model_yvar),]
    df1 = df1[grepl(paste("_",m[i],sep = ""), df1$model_yvar),]
    
    p[[i*2+j-2]] =ggplot(df1, aes(as.factor( model_yvar ), as.numeric(Group_sex_p_value))) + 
      geom_point() + 
      geom_hline(yintercept = 0.05) + 
      labs(y = "G/S anova p-value", x = "brain region") +
      ggtitle(paste(m[i],h[j]) ) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  }
}

library(gridExtra)
grid.arrange(grobs=p)


###### mean difference plot
