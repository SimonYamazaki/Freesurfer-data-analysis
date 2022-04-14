#author: Simon Yamazaki Jensen

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


#This script is intended for computation of statistics on global brain measures
#in multiple groups and comparing to a control group. The script generates:


#####
# - ANOVA tables with models that include a group/sex interaction 
#   an extra row of a model without the interaction is included if it turned out insignificant
#   saved in an excel sheet with a covariate and an excel sheet without

#save paths:
GS_ANOVA_with_cov = "/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/global_variables/globvar_GS_ANOVA_pvals_with_ICV.xlsx"
GS_ANOVA_without_cov = "/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/global_variables/globvar_GS_ANOVA_pvals_without_ICV.xlsx"

#####
# - An excel sheet with model relevant contrasts for each of the models in the ANOVA tables
#   saved in an excel sheet with a covariate and an excel sheet without

#save paths:
contrast_with_cov = "/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/global_variables/globvar_Model_contrasts_with_ICV.xlsx"
contrast_without_cov = "/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/global_variables/globvar_Model_contrasts_without_ICV.xlsx"

#####
# - Variance plots

#####
# - A violin plot of each measure with the individual data points as a percent difference 
#   from the control group mean as scatters. 
#   LSmeans based contrast from each group compared to the control are plotted on top
#   These contrast plots show percent difference from the control LSmean 
#   The plot is split into the two sexes for certain variables and not for others

#prefix on violin plot
LSmeans_prefix = "LSmeans_contrasts_on_datapoints_"

#####
# - Plots of LSmean estimates in each of the groups with and without a global covariate 
#   one plot for each global variable
#   this plot also include a visualisation of the p-values testing a significant contrast

#file name postfix
ppwp_postfix = "_group_diff_pvalues_ICV"

#save folder for plots:
save_folder = "/mnt/projects/VIA11/FREESURFER/Stats/Plots/global_measures"


#data path
data_path = "/mnt/projects/VIA11/FREESURFER/Stats/Data/VIA11_allkey_160621_FreeSurfer_pruned_20220126.csv"

#load data
data_csv <- read.table(data_path, header = TRUE, sep = ",", dec = ".")



#load data
data_csv <- read.table("VIA11_allkey_160621_FreeSurfer_pruned_20220126.csv", header = TRUE, sep = ",", dec = ".")

#inspect the head of data and summary
head(data_csv)
summary(data_csv)

#filter the data with include variable
# - extract rows with 1 in Include_FS_studies coloumn
data_csv_filtered <- data_csv[c(data_csv$Include_FS_studies == 1),]
data_csv_filtered <- data_csv_filtered[!is.na(data_csv_filtered$Include_FS_studies),]

# - extract rows with 1 in Include_FS_studies_euler_outliers_excluded
data_csv_filtered <- data_csv[c(data_csv$Include_FS_studies_euler_outliers_excluded == 1),]
data_csv_filtered <- data_csv_filtered[!is.na(data_csv_filtered$Include_FS_studies_euler_outliers_excluded),]


#make new variables with shorter and contained names
# - tell r which variables are factors
datab = data_csv_filtered
datab$sex = as.factor(datab$Sex_child)
datab$group = as.factor(datab$HighRiskStatus)
datab$site = as.factor(datab$MR_Site)
datab$diag = as.factor(datab$Axis.1_diag_v11)
datab$age = as.numeric(datab$MRI_age)

#add a string encoding of sex for later convenience 
datab$sexs = datab$Sex_child
datab$sexs[datab$sexs == 0] = "female"
datab$sexs[datab$sexs == 1] = "male"

#set another working directory so that files are saved in this folder
setwd(save_folder)

#specify dependent variable (global variable) column names, i.e. left hand side of statistical models
y_vars = c("BrainTotalVol", "CortexVol", "total_area", "mean_thickness", "eICV_samseg")


#### plot variance in each region for each measure ####

# reformat the dataframe for plotting
dft = datab %>%
  pivot_longer(cols = y_vars,
               names_to = "NAME_measure",
               values_to = "VALUE_measure",
               values_drop_na = FALSE)
dft = as.data.frame(dft)

#initialize a coloumn of normalized (percent difference) global variables 
dft$NORM_measure = dft$VALUE_measure
dft$NORM_measure_K = NA
dft$NORM_measure_K_mt = NA

#compute the normalized global variables based on the global variable average
for (i in seq(1,length(y_vars))){
    avg = mean(dft$VALUE_measure[ dft$NAME_measure  == y_vars[i]])
    avg_K = mean(dft$VALUE_measure[ dft$NAME_measure  == y_vars[i] & dft$group  == "K"])
    
    dft$NORM_measure[ dft$NAME_measure  == y_vars[i] ] = ( dft$VALUE_measure[ dft$NAME_measure  == y_vars[i] ] - avg) / avg
    
    dft$NORM_measure_K[ dft$NAME_measure  == y_vars[i] ] = ( dft$VALUE_measure[ dft$NAME_measure  == y_vars[i] ] - avg_K ) / avg_K

}
  
p_var = list()

#Variance plot in original units
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


#plot means of the each group in each global measure
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


## plot raw contrast means 
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

ps=grid.arrange(grobs=p_var)




##############################
# whole brain measure models #
##############################

#investigate what kind of models to run 
model_bvol_glob = lm(BrainTotalVol ~ group*sex + age + TotalEulerNumber + site + eICV_samseg, data=datab)
model_bvol_glob = update(model_bvol_glob,~.-group:sex)
anova(model_bvol_glob)
Anova(model_bvol_glob,type="II")
Anova(model_bvol_glob,type="III")
drop1(model_bvol_glob, test="F")
lsmeans(model_bvol_glob,pairwise~"group",adjust="none")



### run models with eICV_samseg as covariate and without
# generate ANOVA tables with group/sex interaction
# each row is a model 

#which global variables to model and compute statistics for
model_vars = c("BrainTotalVol", "CortexVol", "total_area","mean_thickness", "eICV_samseg")

GS_pvals = list()
ps = list()
emm = list()
min_emm = list()
max_emm = list()
model = list()

DF = data.frame()

#loop that defines models 
for (i in seq(1,length(model_vars))){
  if (model_vars[i] == "eICV_samseg"){    #only run models without_eICV for the model on eICV_samseg
    glob = c("without_eICV")
  }
  else{
    glob = c("without_eICV","with_eICV") #run both for all other variables
  }

  for (g in seq(1,length(glob))){
    if (glob[g] == "with_eICV"){
      f = paste(model_vars[i],"~","+","group*sex","+","age","+","site","+","TotalEulerNumber","+","eICV_samseg") #model with ICV
    }
    else{
      f = paste(model_vars[i],"~","+","group*sex","+","age","+","site","+","TotalEulerNumber") #model without ICV
    }
    #run the actual model
    model[[glob[g]]][[i]] = lm(f,data=datab)  
    
    #extract the variables that are modelling the global variable
    xvars = attributes(Anova(model[[glob[g]]][[i]],type = "III"))$row.names
      
    #if the group/sex interaction is insignificant 
    if (Anova(model[[glob[g]]][[i]],type = "III")$"Pr(>F)"[xvars=="group:sex"] > 0.05){
      #indicate the significance of the group/sex interation
      sign_GS = 0
      
      #save the model with group/sex interaction included
      model_gs = model[[glob[g]]][[i]]
      
      #save the anova F and p-value for the group/sex interaction
      GS_pvals = Anova(model_gs,type = "III")$"Pr(>F)"[xvars=="group:sex"]
      GS_F = Anova(model_gs,type = "III")$"F value"[xvars=="group:sex"]
      
      #save the model without the group/sex interaction
      model[[glob[g]]][[i]] = update(model[[glob[g]]][[i]],~.-group:sex)
      
      #save both models in a list
      models = list(model_gs,model[[glob[g]]][[i]])
      
    }
    else{    #if the group/sex interaction is significant 
      
      #indicate the significance of the group/sex interation
      sign_GS = 1
      
      #save the anova F and p-value for the group/sex interaction
      GS_pvals = Anova(model[[glob[g]]][[i]],type = "III")$"Pr(>F)"[xvars=="group:sex"]
      GS_F = Anova(model[[glob[g]]][[i]],type = "III")$"F value"[xvars=="group:sex"]
      
      #save the model with the group/sex interaction
      models = list(model[[glob[g]]][[i]])
    }
    
    # loop that makes rows in the ANOVA table excel 
    # makes two rows if the group/sex interaction is insignificant 
    # one row with the interaction included and one row with out
    for (m in seq(1,length(models))){
      mm = models[[m]]
      
      #extract relevant statistics for each variable in the model 
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
      
      #the first index in m is always the with a group/sex interaction
      if (m==1){
        rw = list(model_vars[i], 
                  group_F, group_pv, sex_F, sex_pv,
                  age_F, age_pv, site_F, site_pv, 
                  EulerNumber_F, EulerNumber_pv, 
                  ICV_F, ICV_pv, GS_F, GS_pvals,
                  sign_GS,m,g-1)
      }
      else{ #if m is 2 then the model is does not include the interaction
        rw = list(model_vars[i], 
                  group_F, group_pv, sex_F, sex_pv, 
                  age_F, age_pv, site_F, site_pv, 
                  EulerNumber_F, EulerNumber_pv, 
                  ICV_F, ICV_pv, NA,   NA,
                  sign_GS,m-2,g-1)
      }
      
      #append the row to the dataframe 
      DF = rbindlist(list(DF, rw))
      
    } # for m 
  } #g
} #end i

#define the coloum names of the dataframe which is converted to an excel 
col_names = c("Model_yvar",
                  "Group_Fval","Group_pval","Sex_Fval","Sex_pval",
                  "Age_Fval","Age_pval","Site_Fval","Site_pval",
                  "EulerNumber_Fval","Eulernumber_pval",
                  "ICV_Fval","ICV_pval","Group_sex_Fval","Group_sex_pval",
                  "Significant_GS_interaction","GS_in_model","ICV_in_model")
names(DF)<-col_names

#A potential way of reducing significant digits 
#df_with_ICV[,2:ncol(df_with_ICV)] <- signif(df_with_ICV[,2:ncol(df_with_ICV)],digits=3)
#df_without_ICV[,2:ncol(df_without_ICV)] <- signif(df_without_ICV[,2:ncol(df_without_ICV)],digits=3)

#split up the dataframe in models with ICV and without
DF_xlsx_ICV0 = DF[DF$ICV_in_model == 0, ]
DF_xlsx_ICV1 = DF[DF$ICV_in_model == 1, ]

#save the two dataframes in the a specified save path 
write_xlsx(DF_xlsx_ICV1,GS_ANOVA_with_cov)
write_xlsx(DF_xlsx_ICV0,GS_ANOVA_without_cov)




#### CONTRAST XLSX with male / female split
# and built models[sex][glob][yvar]

#define some relevant naming 
model_vars = c("BrainTotalVol", "CortexVol", "total_area", "mean_thickness", "eICV_samseg")
glob = c("with_eICV", "without_eICV")
sex = c("female","male","both") #3 sexes are defined, 0=female, 1=male, 2=both

#split up the data into male and female 
data_sex1 = datab[c(datab$sex == 1),]
data_sex0 = datab[c(datab$sex == 0),]
dataf = list(data_sex0, data_sex1, datab)

models = list()

DFc = data.frame()


# loop that runs the model on each global variable, each sex with and without ICV
for (k in seq(1,length(model_vars))){
  print(model_vars[k])
  for (ss in seq(1,3)){ #use the 3 defined sexes to know what data to use
    
    # define an temporary sex indicator to index
    print(sex[ss])
    
    #extract the dataframe containing the relevant sex data
    df_sex = dataf[[ss]]
    
    #only continue with both sex data when the global variable is mean_thickness 
    if (model_vars[k] == "mean_thickness" & sex[ss] != "both"){
      next
    }
    if (model_vars[k] != "mean_thickness" & sex[ss] == "both"){
      next
    }
    
    
    for (gg in seq(1,2)){
      
      if (glob[gg] == "without_eICV"){ #run models without ICV
        #define the models 
        if (model_vars[k] == "mean_thickness"){ #define a specific model for the mean_thickness variable
          f = paste(model_vars[k],"~","group","+","sex","+","age","+","site","+","TotalEulerNumber")
        }
        else {
          f = paste(model_vars[k],"~","group","+","age","+","site","+","TotalEulerNumber")
        }
        #run the model
        models[[sex[ss]]][[glob[gg]]][[model_vars[k]]] = lm(f,data=df_sex)
        print(glob[gg])
        
        #save the current model into a temporary variable to make code less cluttered
        model_ana = models[[sex[ss]]][[glob[gg]]][[model_vars[k]]]
        xvars = attributes(Anova(model_ana,type = "III"))$row.names
      }
      
      else { #run models with ICV
        #define models
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
      #compute lsmeans pairwise for each group with no correction on p-values
      ls = lsmeans(model_ana,pairwise~"group", adjust="none")
      
      #compute contrast statistics
      c = ls$contrasts
      
      #extract the lsmean of the control group
      K_emm = summary(ls)$lsmeans$lsmean[2]
      
      #estimates and statistics on raw contrasts
      #for bipolar group
      BP_diff = summary(c)$estimate[1]
      BP_diff_pv = summary(c)$"p.value"[1]
      BP_diff_tratio = summary(c)$"t.ratio"[1]
      BP_diff_LCL = confint(c)$lower.CL[1]
      BP_diff_UCL = confint(c)$upper.CL[1]
      
      #for Skizo 
      SZ_diff = -summary(c)$estimate[3]
      SZ_diff_pv = summary(c)$"p.value"[3]
      SZ_diff_tratio = summary(c)$"t.ratio"[3]
      SZ_diff_LCL = -confint(c)$lower.CL[3]
      SZ_diff_UCL = -confint(c)$upper.CL[3]
      
      
      #compute percent difference from control
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
                      gg-1,ss-1)
      
      DFc = rbindlist(list(DFc, rwc_xlsx))
      
    } #g 
  } #s
} #k

#same excel saving procedure as above 
col_names = c("Model_yvar", "K_LSmean",
              "Contrast_BP-K","tratio_BP-K","pval_BP-K",
              "Contrast_SZ-K","tratio_SZ-K","pval_SZ-K",
              "ICV_in_model","sex")
names(DFc)<-col_names

DF_xlsx_glob0 = DFc[DFc$ICV_in_model == 0, ]
DF_xlsx_glob1 = DFc[DFc$ICV_in_model == 1, ]

write_xlsx(DF_xlsx_glob1,contrast_with_cov)
write_xlsx(DF_xlsx_glob0,contrast_without_cov)



### plot ppwp (contrast p-value plot)
# and the LSmeans 

for (i in seq(1,length(model_vars))){
  ps = list()
  emm = list()
  min_emm = list()
  max_emm = list()
  print(model_vars[i])
  
  #define what models LSmeans should be plotted for 
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
    
    # loop to ppwp
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
    
    # plot LSmeans
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
  
  #save the plot here
  ggsave(paste(model_vars[i],ppwp_postfix,".png",sep=""),ga,width = 10,height = 10)
} #end i







####### LSmeans estimates on the variance plot preparation #########

#define relevant variables and names
y_vars = c("BrainTotalVol", "CortexVol", "total_area","mean_thickness")
highriskG = c("BP","SZ")
glob = c("with_eICV", "without_eICV")
sex_name = c("Female","Male")
pd <- position_dodge(0.05) 

pls = list()

#reformat data
dftm = datab %>%
  pivot_longer(cols = y_vars,
               names_to = "NAME_measure",
               values_to = "VALUE_measure",
               values_drop_na = FALSE)
dftm = as.data.frame(dftm)

#initialize new coloumns for the normalized data
dftm$NORM_measure = dftm$VALUE_measure
dftm$NORM_measure_K = NA  #percent difference from control
dftm$NORM_measure_K_mt = NA #percent difference from control only for mean_thickness

dfc_with_ICV = data.frame()
dfc_without_ICV = data.frame()


#loop that defines a new dataframe for plotting 

for (g in seq(1,length(glob))){
  
  for (i in seq(1,length(y_vars))){
    #the average value of a global variable for control group
    avg_K = mean(dftm$VALUE_measure[ dftm$NAME_measure  == y_vars[i] & dftm$group  == "K"])
    
    if (y_vars[i] == "mean_thickness"){
      #fill the coloumn of percent difference from control for mean_thickness
      #all coloumn entries are identical 
      dftm$NORM_measure_K_mt[ dftm$NAME_measure  == y_vars[i] ] = ( dftm$VALUE_measure[ dftm$NAME_measure  == y_vars[i] ] - avg_K ) / avg_K
      sexplot = c("both") #a variable describing the sexes to be plotted
    }
    else{
      #fill the coloumn of percent difference from control for the remainding global variables
      #all coloumn entries are identical 
      dftm$NORM_measure_K[ dftm$NAME_measure  == y_vars[i] ] = ( dftm$VALUE_measure[ dftm$NAME_measure  == y_vars[i] ] - avg_K ) / avg_K
      sexplot = c("female","male") #a variable describing the sexes to be plotted
      }
    
    
    for (s in seq(1,length(sexplot))){
      
      ls = lsmeans(models[[sexplot[s]]][[glob[g]]][[model_vars[i]]],pairwise~"group",adjust="none")
      c = ls$contrasts
      K_emm = summary(ls)$lsmeans$lsmean[2]
      
      #percent difference
      BP_diff_P = 100*summary(c)$estimate[1] /K_emm
      BP_diff_LCL_P = 100*confint(c)$lower.CL[1] /K_emm
      BP_diff_UCL_P = 100*confint(c)$upper.CL[1] /K_emm
      
      SZ_diff_P = -100*summary(c)$estimate[3] /K_emm
      SZ_diff_LCL_P = -100*confint(c)$lower.CL[3] /K_emm
      SZ_diff_UCL_P = -100*confint(c)$upper.CL[3] /K_emm
      
      
      if (sexplot[s] == "female" || sexplot[s] == "male"){

        ###
        dftm$lsmean_sex[dftm$NAME_measure  == y_vars[i] & dftm$group == "BP" & dftm$sexs == sexplot[s]] = BP_diff_P
        dftm$eb_max_sex[dftm$NAME_measure  == y_vars[i] & dftm$group == "BP" & dftm$sexs == sexplot[s]] = BP_diff_UCL_P
        dftm$eb_min_sex[dftm$NAME_measure  == y_vars[i] & dftm$group == "BP" & dftm$sexs == sexplot[s]] = BP_diff_LCL_P

        ###
        dftm$lsmean_sex[dftm$NAME_measure  == y_vars[i] & dftm$group == "SZ" & dftm$sexs == sexplot[s]] = SZ_diff_P
        dftm$eb_max_sex[dftm$NAME_measure  == y_vars[i] & dftm$group == "SZ" & dftm$sexs == sexplot[s]] = SZ_diff_UCL_P
        dftm$eb_min_sex[dftm$NAME_measure  == y_vars[i] & dftm$group == "SZ" & dftm$sexs == sexplot[s]] = SZ_diff_LCL_P

      }
      else{
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



### PLOT THE LSMEANS on the variance plot

#loop that plots the lsmeans 

for (g in seq(1,length(glob))){
  DF = get(paste("dfc_plot_",glob[g],sep = ""))
  
  for (s in seq(1,2)){
    pls[[s]]= with(DF[!is.na(DF$NORM_measure_K) & DF$sex == s-1, ],
                   ggplot() + 
                     aes(x = factor(NAME_measure), y = 100*NORM_measure_K, color = group) +
                     geom_violin(position = "identity",alpha=0.3) +
                     geom_jitter(width = 0.3, size=0.1) + 
                     geom_hline(yintercept = 0,linetype="dashed") + 
                     geom_point(position=pd, aes_string(x = "NAME_measure", y = "lsmean_sex", color="group", group="group"), size=2) +
                     stat_summary(fun = mean,geom = "line",position=pd, aes_string(x = "NAME_measure", y = "lsmean_sex", color="group", group="group")) + 
                     geom_errorbar(position=pd, width=0.2, aes_string(group="group",color="group",ymin="eb_min_sex", ymax="eb_max_sex") ) +
                     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
                     labs(y = "Difference from control [%]", x = "brain measure") +
                     ggtitle(paste(sex_name[s]))
                     
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
                     
    )
    pls[[4]] = ggplot() + theme_void()
    
  top_title = paste("LSmeans contrasts of model",glob[g])
  ga=grid.arrange(grobs=pls,ncols=2, top=textGrob(top_title,gp=gpar(fontsize=20)))
  
  #save the plot here
  ggsave(paste(LSmeans_prefix,glob[g],".png",sep=""),ga,width = 10,height = 10)
  
} #g



