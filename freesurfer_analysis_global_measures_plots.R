

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







### run and plot models with euler number as covariate and without

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
