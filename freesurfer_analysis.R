
setwd("/mnt/projects/VIA11/FREESURFER/Stats/Data")

#load packages
library(data.table)
library(ggplot2)
library(ggdist)
library(gridExtra)
library(tidyr)
library(lsmeans)
library(grid)

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

#initial data exploration
group_color = c("green", "blue", "red")

y_vars = c("BrainTotalVol", "CortexVol", "total_area", "mean_thickness", "eICV_samseg")

setwd("/mnt/projects/VIA11/FREESURFER/Stats/Plots/global_measures")

p_inter = list()

for (i in seq(1,length(y_vars))){
p_inter[[i]]=with(datab,
  ggplot() +
  aes_string(x = "group", color = "sex", group = "sex", y = y_vars[i]) +
  geom_jitter(width = 0.1, size=0.1) + 
  stat_summary(fun = mean, geom = "point",size=3) +
  stat_summary(fun = mean, geom = "point",size=3,pch=21,colour="black") +
  stat_summary(fun = mean, geom = "line") +
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE))
)
}

ps=grid.arrange(grobs=p_inter)
ggsave(paste("group:sex_interaction_check",".png",sep=""),ps,width = 10,height = 10)


### whole brain measure models 

#WITH global variable covariate 
model_bvol_glob = lmer(BrainTotalVol ~ group*sex + age + eICV_samseg + (1| site), data=datab)
ranova(model_bvol_glob)
anova(model_bvol_glob)

model_cvol_glob = lmer(CortexVol ~ group*sex + age + eICV_samseg + (1| site), data=datab)
ranova(model_cvol_glob)
anova(model_cvol_glob)

model_area_glob = lmer(total_area ~ group*sex + age + eICV_samseg + (1| site), data=datab)
ranova(model_area_glob)
anova(model_area_glob)

model_mt_glob = lmer(mean_thickness ~ group + sex + age + eICV_samseg + (1| site), data=datab)
ranova(model_mt_glob)
anova(model_mt_glob)


#WITHOUT global variable models 
model_bvol = lmer(BrainTotalVol ~ group*sex + age + (1| site), data=datab)
ranova(model_bvol)
anova(model_bvol)
summary(model_bvol)

model_cvol = lmer(CortexVol ~ group*sex + age + (1| site), data=datab)
model_cvol = lm(CortexVol ~ group*sex + age + site, data=datab)
ranova(model_cvol)
anova(model_cvol)

model_area = lmer(total_area ~ group*sex + age + (1| site), data=datab)
model_area = lm(total_area ~ group*sex + age + site, data=datab)
ranova(model_area)
anova(model_area)

model_mt = lmer(mean_thickness ~ group + sex  + age + (1| site), data=datab)
model_mt = lm(mean_thickness ~ group + sex  + age + site, data=datab)
ranova(model_mt)
anova(model_mt)

model_icv = lmer(eICV_samseg ~ group*sex + age + (1| site), data=datab)
ranova(model_icv)
anova(model_icv)


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
#run multiple models to be able to plot GROUP:AGE INTERACTION

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

setwd("/mnt/projects/VIA11/FREESURFER/Stats/Plots/lateral")

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


#plot variance in each region for each measure

mm = "lh_area"
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
              #aes(x = factor(NAME_region_measure_h), color = group, group = group, y = 100*NORM_region_measure_h) +
              #see::geom_violinhalf(aes(group = NAME_region_measure_h)) + 
              geom_violin(position = "identity",alpha=0.3) +
              #ggdist::stat_halfeye(adjust = .5, width = .7, .width = 0, justification = -.2, point_colour = NA) +
              geom_jitter(width = 0.4, size=0.1) + 
              geom_hline(yintercept = 0) + 
              stat_summary(fun = mean, geom = "point",size=3,aes(colour = group)) +
              #stat_summary(fun = mean, geom = "point",size=3,pch=21) +
              theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
              labs(y = "Percent difference from mean [%]", x = "brain region") +
              ggtitle(paste(mm))
)
ps=grid.arrange(grobs=p_var)
ggsave(paste("original_unit_all_region_variance",mm,".png",sep=""),ps,width = 14,height = 10)

ps=grid.arrange(p_var[[2]])
ggsave(paste("all_region_variance",mm,".png",sep=""),ps,width = 14,height = 6)


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


setwd("/mnt/projects/VIA11/FREESURFER/Stats/Plots/lateral/site_random")


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
        f = paste(model_yvars[k],"~","-1","+","group*sex","+","age","+","(1| site)")
        models[[h[i]]][[m[j]]][[model_yvars[k]]] = lmer(f,data=datab,REML=FALSE)
        
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
          f2 = paste(model_yvars[k],"~","-1","+","group*sex","+","age","+","total_area","+","(1| site)")
        }
        else if (m[j] == "thickness"){
          f2 = paste(model_yvars[k],"~","-1","+","group*sex","+","age","+","mean_thickness","+","(1| site)")
        }
        else if (m[j] == "volume"){
          f2 = paste(model_yvars[k],"~","-1","+","group*sex","+","age","+","BrainTotalVol","+","(1| site)")
        }
        models_glob[[h[i]]][[m[j]]][[model_yvars[k]]] = lmer(f2,data=datab,REML=FALSE)
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
        model_LRT_pv = anova(models_glob[[h[i]]][[m[j]]][[model_yvars[k]]],models[[h[i]]][[m[j]]][[model_yvars[k]]],refit=FALSE)$"Pr(>Chisq)"[2]
        rw = append(rw,model_LRT_pv)
        
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
glob = c("Model WITHOUT global var","Model WITH global var")

DF_list = list()
DF_list[[1]] = DF
DF_list[[2]] = DF_glob


for (g in seq(1,2)){
  pgs = list()
  df = DF_list[[g]]
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
  top_title = paste(glob[g])
  ps=grid.arrange(grobs=pgs,top=textGrob(top_title,gp=gpar(fontsize=20)))
  ggsave(paste(glob[g],"group:sex_interaction_pvalue",".png",sep=""),ps,width = 10,height = 12)
}


#sanity check above results with manual models 
model_test1 = lmer(rh_parahippocampal_area ~ group*sex + age + (1| site), data=datab,REML=FALSE)
DF$Group_sex_p_value[grepl(paste("rh_parahippocampal_area"), DF$model_yvar)]
anova(model_test1)

model_test2 = lmer(rh_parahippocampal_area ~ group*sex + age + total_area + (1| site), data=datab,REML=FALSE)
DF_glob$Group_sex_p_value[grepl(paste("rh_parahippocampal_area"), DF_glob$model_yvar)]
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
  
  for (i in seq(1,2)){
    for (j in seq(1,3)){
  
      for (k in seq(1,2)){
        dfs = datafs[[k]]
        
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

      }
      top_title = paste(h[i],m[j],"for group",groups[g])
      lay <- rbind(c(1,1,1,1,1,1,1,2,2,2,2),
                   c(1,1,1,1,1,1,1,2,2,2,2),
                   c(1,1,1,1,1,1,1,2,2,2,2))
      ps=grid.arrange(grobs=sp, layout_matrix=lay, top=textGrob(top_title,gp=gpar(fontsize=20)))
      ggsave(paste(groups[g],"_",h[i],"_",m[j],"_mean_difference",".png",sep=""),ps,width = 10,height = 12)
  }
}
}



#sanity check
model_test3 = lmer(lh_parahippocampal_volume ~ group + sex + age + BrainTotalVol + (1| site), data=datab)
lsmeans(model_test3,pairwise~"group")
model_bvol2 = lmer(mean_thickness ~ -1 + group + sex + age + (1| site), data=datab)
anova(model_bvol2)
lsmeans(model_bvol2,pairwise~"group",by="sex")
plot(lsmeans(model_bvol2,pairwise~"group",by="sex"),comparison=TRUE)



#### plot evt regulær model vs glob model, som parallel coordinate??
DF$model_diff_pv <- NA
DF$glob <- 0
DF_glob$glob <- 1
new <- rbind(DF, DF_glob)
new$model_type[new$model_type == 1] = "no glob"
new$model_type[new$model_type == 2] = "with glob"

new$significant_GS[new$Group_sex_p_value < 0.05] = "G:S"
new$significant_GS[new$Group_sex_p_value > 0.05] = "NO G:S"


pcp=ggplot(new, aes(vars = vars(c(33,6,12,9,15,30)))) + #remember to change back to 32
  geom_pcp_box(boxwidth=0.1, fill="grey") +
  geom_pcp(boxwidth=0.1, aes(colour = factor(model_type))) +
  geom_pcp_text(boxwidth=0.1) +
  scale_colour_manual(values=c("darkorange", "steelblue")) + 
  guides(colour=guide_legend(override.aes = list(alpha=1)))

ps=grid.arrange(pcp)
ggsave(paste("pcp",".png",sep=""),ps,width = 10,height = 6)



#significance of model type
pm = list()
df = DF_glob

for (i in seq(1,3)){
  for (j in seq(1,2)){
    df1 = df[grepl(paste("^",h[j],sep = ""), df$model_yvar),]
    df1 = df1[grepl(paste("_",m[i],sep = ""), df1$model_yvar),]
    
    pm[[i*2+j-2]] =ggplot(df1, aes(x=as.factor( model_yvar ), y=as.numeric(model_diff_pv))) + 
      geom_point() + 
      #geom_hline(yintercept = 0.05) + 
      labs(y = "LRT p-value", x = "brain region") +
      ggtitle(paste(m[i],h[j]) ) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  }
}

ps=grid.arrange(grobs=pm)
ggsave(paste("model_comparison_pvalue",".png",sep=""),ps,width = 10,height = 12)






#############  MANOVA   ##############

# Global measures # 
y_vars = c("CortexVol", "total_area", "mean_thickness")#, "eICV_samseg")
#y_vars = col_names_list[[paste("lh","area",sep = "_")]] 
y_vars_idx = (names(datab) %in% y_vars == TRUE)
my_vars <- data.matrix(datab[,c(y_vars_idx)])

mmodel = manova(my_vars ~ group*sex + age + site, data = datab)
anova(mmodel)
mmodel = manova(my_vars ~ group + sex + age + site, data = datab)
anova(mmodel)
summary.aov(mmodel)

PH = lsmeans(mmodel,pairwise~"group", by = "sex")
plot(PH,comparison=TRUE,xlab="Multivariate lsmean")




# 34 regions combined for each measure# 
## ALL 3
y_vars = c(col_names_list[[paste("lh","area",sep = "_")]],col_names_list[[paste("lh","thickness",sep = "_")]],col_names_list[[paste("lh","volume",sep = "_")]])
y_vars_idx = (names(datab) %in% y_vars == TRUE)
my_vars <- data.matrix(datab[,c(y_vars_idx)])

mmodel1 = manova(my_vars ~ group*sex + age + site, data = datab)
anova(mmodel)
mmodel1 = manova(my_vars ~ group + sex + age + site, data = datab)
anova(mmodel)


## AREA
y_vars = col_names_list[[paste("lh","area",sep = "_")]] 
y_vars_idx = (names(datab) %in% y_vars == TRUE)
my_vars <- data.matrix(datab[,c(y_vars_idx)])

mmodel = manova(my_vars ~ group*sex + age + site, data = datab)
anova(mmodel)
mmodel = manova(my_vars ~ group + sex + age + site, data = datab)
anova(mmodel)

## VOLUME
y_vars = col_names_list[[paste("lh","volume",sep = "_")]] 
y_vars_idx = (names(datab) %in% y_vars == TRUE)
my_vars <- data.matrix(datab[,c(y_vars_idx)])

mmodel = manova(my_vars ~ group*sex + age + site, data = datab)
anova(mmodel)
mmodel = manova(my_vars ~ group + sex + age + site, data = datab)
anova(mmodel)

## THICKNESS
y_vars = col_names_list[[paste("lh","thickness",sep = "_")]] 
y_vars_idx = (names(datab) %in% y_vars == TRUE)
my_vars <- data.matrix(datab[,c(y_vars_idx)])

mmodel = manova(my_vars ~ group*sex + age + site, data = datab)
anova(mmodel)
mmodel = manova(my_vars ~ group + sex + age + site, data = datab)
anova(mmodel)



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



#col_names = names(datab)
#col_names = col_names[col_names %in% model_vars == FALSE]
#col_names_rest = col_names[!grepl("^lh", col_names)]
#col_names_rest = col_names_rest[!grepl("^rh", col_names_rest)]
#col_names_rest = col_names[col_names %in% model_vars == FALSE]

#nam = paste("sex1_",groups[g],"_diff",sep = "")
#df1[[nam]] = as.numeric(df1[[nam]])
#nam = paste("sex1_",groups[g],"_diff_LCL",sep = "")
#df1[[nam]] = as.numeric(df1[[nam]])
#nam = paste("sex1_",groups[g],"_diff_UCL",sep = "")
#df1[[nam]] = as.numeric(df1[[nam]])
