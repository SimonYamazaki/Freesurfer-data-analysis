
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

#tell r which variables are factors
datab = data_csv_filtered
datab$sex = as.factor(datab$Sex_child)
datab$group = as.factor(datab$HighRiskStatus)
datab$site = as.factor(datab$MR_Site)
datab$diag = as.factor(datab$Axis.1_diag_v11)
datab$age = as.numeric(datab$MRI_age)




### whole brain measure models 

#WITH global variable covariate 
model_bvol_glob = lm(BrainTotalVol ~ group*sex + age + eICV_samseg + site, data=datab)
anova(model_bvol_glob)

model_cvol_glob = lm(CortexVol ~ group*sex + age + eICV_samseg + site, data=datab)
anova(model_cvol_glob)

model_area_glob = lm(total_area ~ group*sex + age + eICV_samseg + site, data=datab)
anova(model_area_glob)

model_mt_glob = lm(mean_thickness ~ group + sex + age + eICV_samseg + site, data=datab)
anova(model_mt_glob)


#WITHOUT global variable models 
model_bvol = lm(BrainTotalVol ~ group*sex + age + site, data=datab)
anova(model_bvol)
summary(model_bvol)

model_cvol = lm(CortexVol ~ group*sex + age + site, data=datab)
anova(model_cvol)

model_area = lm(total_area ~ group*sex + age + site, data=datab)
anova(model_area)

model_mt = lm(mean_thickness ~ group + sex  + age + site, data=datab)
anova(model_mt)

model_icv = lm(eICV_samseg ~ group*sex + age + site, data=datab)
anova(model_icv)


### run 4 models with eICV_samseg as covariate
model_vars = c("BrainTotalVol", "CortexVol", "total_area")
model = list()
for (i in seq(1,length(model_vars))){
  f = paste(model_vars[i],"~","-1","+","group*sex","+","age","+","eICV_samseg","+","site")
  model[[i]] = lm(f,data=datab)
}
model[[i+1]] = lm(mean_thickness ~ -1 + group + sex + age + eICV_samseg + site,data=datab)
model_vars = c(model_vars,"mean_thickness")


#####################
#run multiple models to be able to plot GROUP:AGE INTERACTION

model_vars = c("BrainTotalVol", "CortexVol", "total_area", "eICV_samseg")
model = list()
for (i in seq(1,length(model_vars))){
  f = paste(model_vars[i],"~","-1","+","group*sex","+","group:age","+","site")
  model[[i]] = lm(f,data=datab)
}
model[[i+1]] = lm(mean_thickness ~ -1 + group + sex + group:age + site,data=datab)
model_vars = c(model_vars,"mean_thickness")



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


######################################
#    run models on each region
######################################

# extract model yvars 

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

setwd("/mnt/projects/VIA11/FREESURFER/Stats/Plots/site_fixed")


###### make inference with models
models = list()
models_glob = list()

DF = data.frame()
DFs = data.frame()
DF_glob = data.frame()
DFs_glob = data.frame()

for (i in seq(1,2)){
  for (j in seq(1,3)){
    print(paste(h[i],m[j]))
    model_yvars = col_names_list[[paste(h[i],m[j],sep = "_")]] 
    for (k in seq(1,length(model_yvars))){
      for (mi in seq(1,2)){
        
        if (mi == 1){
          #no global measure model
          f = paste(model_yvars[k],"~","-1","+","group*sex","+","age","+","site")
          models[[h[i]]][[m[j]]][[model_yvars[k]]] = lm(f,data=datab)
          
          model_ana = models[[h[i]]][[m[j]]][[model_yvars[k]]]
          xvars = attributes(anova(model_ana))$row.names
          pv_group_sex = anova(model_ana)$"Pr(>F)"[xvars=="group:sex"]
          
          if (anova(model_ana)$"Pr(>F)"[xvars=="group:sex"] > 0.05){
            model_ana = update(model_ana,~.-group:sex)
          } 
          models[[h[i]]][[m[j]]][[model_yvars[k]]] = model_ana
          
        }
        else {
          #global measure model
          if (m[j] == "area"){
            f2 = paste(model_yvars[k],"~","-1","+","group*sex","+","age","+","total_area","+","site")
          }
          else if (m[j] == "thickness"){
            f2 = paste(model_yvars[k],"~","-1","+","group*sex","+","age","+","mean_thickness","+","site")
          }
          else if (m[j] == "volume"){
            f2 = paste(model_yvars[k],"~","-1","+","group*sex","+","age","+","BrainTotalVol","+","site")
          }
          models_glob[[h[i]]][[m[j]]][[model_yvars[k]]] = lm(f2,data=datab)
          
          model_ana = models_glob[[h[i]]][[m[j]]][[model_yvars[k]]]
          xvars = attributes(anova(model_ana))$row.names
          pv_group_sex = anova(model_ana)$"Pr(>F)"[xvars=="group:sex"]
          
          if (anova(model_ana)$"Pr(>F)"[xvars=="group:sex"] > 0.05){
            model_ana = update(model_ana,~.-group:sex)
          } 
          models_glob[[h[i]]][[m[j]]][[model_yvars[k]]] = model_ana
        }
        
        
        #use the above defined model_ana 
        
        ls = lsmeans(model_ana,pairwise~"group", by = "sex")
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
        
        rw = list(model_yvars[k],pv_group,pv_group_sex, K_emm, sex1_K_emm,
                  BP_diff, BP_diff_LCL, BP_diff_UCL, 
                  SZ_diff, SZ_diff_LCL, SZ_diff_UCL, 
                  sex1_BP_diff, sex1_BP_diff_LCL, sex1_BP_diff_UCL, 
                  sex1_SZ_diff, sex1_SZ_diff_LCL, sex1_SZ_diff_UCL,
                  BP_diff_P, BP_diff_LCL_P, BP_diff_UCL_P, 
                  SZ_diff_P, SZ_diff_LCL_P, SZ_diff_UCL_P, 
                  sex1_BP_diff_P, sex1_BP_diff_LCL_P, sex1_BP_diff_UCL_P, 
                  sex1_SZ_diff_P, sex1_SZ_diff_LCL_P, sex1_SZ_diff_UCL_P,mi)
        
        rws0 = list(model_yvars[k], K_emm, sex1_K_emm, 
                    BP_diff, BP_diff_LCL, BP_diff_UCL, 
                    SZ_diff, SZ_diff_LCL, SZ_diff_UCL,
                    BP_diff_P, BP_diff_LCL_P, BP_diff_UCL_P, 
                    SZ_diff_P, SZ_diff_LCL_P, SZ_diff_UCL_P,0)
        
        rws1 = list(model_yvars[k], K_emm, sex1_K_emm, 
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
          model_diff_pv = anova(models_glob[[h[i]]][[m[j]]][[model_yvars[k]]],models[[h[i]]][[m[j]]][[model_yvars[k]]])$"Pr(>F)"[2]
          rw = append(rw,model_diff_pv)
          
          DF_glob = rbindlist(list(DF_glob, rw))
          DFs_glob = rbindlist(list(DFs_glob, rws0))
          DFs_glob = rbindlist(list(DFs_glob, rws1))
        }
        
      } #mi 
      
    } #k
  } #j
} #i

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
names(DFs)<-c("model_yvar","K_emm", "sex1_K_emm",
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
                  "sex1_SZ_diff_P", "sex1_SZ_diff_LCL_P", "sex1_SZ_diff_UCL_P","model_type","model_diff_pv")
names(DFs_glob)<-c("model_yvar","K_emm", "sex1_K_emm",
                   "BP_diff","BP_diff_LCL","BP_diff_UCL",
                   "SZ_diff","SZ_diff_LCL","SZ_diff_UCL",
                   "BP_diff_P","BP_diff_LCL_P","BP_diff_UCL_P",
                   "SZ_diff_P","SZ_diff_LCL_P","SZ_diff_UCL_P","sex")

DFs$sex = as.factor(DFs$sex)
DFs_glob$sex = as.factor(DFs_glob$sex)



#group sex interaction
pgs = list()
df = DF_glob

for (i in seq(1,3)){
  for (j in seq(1,2)){
    df1 = df[grepl(paste("^",h[j],sep = ""), df$model_yvar),]
    df1 = df1[grepl(paste("_",m[i],sep = ""), df1$model_yvar),]
    
    pgs[[i*2+j-2]] =ggplot(df1, aes(as.factor( model_yvar ), as.numeric(Group_sex_p_value))) + 
      geom_point() + 
      geom_hline(yintercept = 0.05) + 
      labs(y = "G/S anova p-value", x = "brain region") +
      ggtitle(paste(m[i],h[j]) ) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  }
}

grid.arrange(grobs=pgs)



###### mean difference plot
groups = c("BP","SZ")
P = c("_P","")
Plab = c("Difference from control [%]","Orignial units")
h = c("lh","rh")
m = c("area","thickness","volume")
glob = c("Model WITHOUT global var","Model WITH global var")

datafs = list(DFs,DFs_glob)
sp=list()

for (k in seq(1,2)){
  dfs = datafs[[k]]
  
  g = 1
  pp = 1
  i = 1
  j = 3
  
  df1 = dfs[grepl(paste("^",h[i],sep = ""), dfs$model_yvar),]
  df1 = df1[grepl(paste("_",m[j],sep = ""), df1$model_yvar),]
  
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
  
  top_title = paste(h[i],m[j],"for group",groups[g])
}
lay <- rbind(c(1,1,1,1,1,1,1,2,2,2,2),
             c(1,1,1,1,1,1,1,2,2,2,2),
             c(1,1,1,1,1,1,1,2,2,2,2))
ps=grid.arrange(grobs=sp, layout_matrix=lay, top=textGrob(top_title,gp=gpar(fontsize=20)))
ggsave(paste(groups[g],"_",h[i],"_",m[j],"_mean_difference",".png",sep=""),ps,width = 10,height = 12)





