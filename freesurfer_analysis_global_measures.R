
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
    
#    if (y_vars[i] == "mean_thickness"){
#      dft$NORM_measure_K_mt[ dft$NAME_measure  == y_vars[i] ] = ( dft$VALUE_measure[ dft$NAME_measure  == y_vars[i] ] - avg_K ) / avg_K
#    }
#    else{
#      dft$NORM_measure_K[ dft$NAME_measure  == y_vars[i] ] = ( dft$VALUE_measure[ dft$NAME_measure  == y_vars[i] ] - avg_K ) / avg_K
#    }
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
model_bvol_glob = lm(BrainTotalVol ~ group*sex + age + eICV_samseg + site, data=datab)
anova(model_bvol_glob)
model_bvol_glob = update(model_bvol_glob,~.-group:sex)
anova(model_bvol_glob)

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

model_vars = c("BrainTotalVol", "CortexVol", "total_area","mean_thickness")
glob = c("with_eICV", "without_eICV")
GS_pvals = list()
ps = list()
emm = list()
min_emm = list()
max_emm = list()
model = list()

for (i in seq(1,length(model_vars))){
  for (g in seq(1,2)){
    if (glob[g] == "with_eICV"){
      f = paste(model_vars[i],"~","+","group*sex","+","age","+","site","+","TotalEulerNumber","+","eICV_samseg")
    }
    else{
      f = paste(model_vars[i],"~","+","group*sex","+","age","+","site","+","TotalEulerNumber")
    }
    model[[glob[g]]][[i]] = lm(f,data=datab)
    xvars = attributes(anova(model[[glob[g]]][[i]]))$row.names
    
    GS_pvals[[glob[g]]][[model_vars[i]]] = anova(model[[glob[g]]][[i]])$"Pr(>F)"[xvars=="group:sex"]
    
    if (anova(model[[glob[g]]][[i]])$"Pr(>F)"[xvars=="group:sex"] > 0.05){
      model[[glob[g]]][[i]] = update(model[[glob[g]]][[i]],~.-group:sex)
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








####### ARCHIVE 

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

