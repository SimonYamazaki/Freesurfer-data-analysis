

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
library(rlist)

#This script is intended for computation of statistics on behavioural measures
#based on GLMs on a single group at a time (one of: SZ, BP, Control)

# models with one of the following behavioural measures as response variables is produced: 
#("CBCL_ext_cg_v11","CBCL_int_cg_v11","CBCL_totsc_cg_v11","CGASx_v11")
#with the following global brain measures as covariates: 
#("BrainTotalVol", "CortexVol", "total_area", "mean_thickness", "eICV_samseg")

#the following covariates are included in all models by default
#("age","site","TotalEulerNumber")


#returns: 
#relevant statistics from GLMs run for each combination of response variable and global brain measures and saved into an excel sheets 
#files generated are: 
# - excel sheet with variable level GLM effects without global covariate (eICV_samseg)
# - excel sheet with variable level GLM effect sizes
#For file names and directories saved into, refer to the variables run_group, save_folder, model_ANOVA_without_glob and effect_sizes_path
#data is loaded from the path specified in data_path


#which group data should the model use
run_group = "BP"

#save folder for tables and plots:
save_folder = paste("/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/global_variables/behav/behav_response/",run_group,"/",sep="")
#dir.create(file.path(dirname(save_folder_plot), basename(save_folder_plot)))

#create the folders if they dont exist
dir.create(file.path(dirname(save_folder), basename(save_folder)))

#includes one sheet with models, one row for both sex, and one row for each sex
#this sheet does not show group differences -> no contrasts either
model_ANOVA_without_glob= "behavvar_ANOVA_pvals_without_glob.xlsx"
effect_sizes_path = "ANOVA_effect_sizes.xlsx"


#data path
data_path = "/mnt/projects/VIA11/FREESURFER/Stats/Data/VIA11_allkey_160621_FreeSurfer_pruned_20220509.csv"

#load data
data_csv <- read.table(data_path, header = TRUE, sep = ",", dec = ".")


#inspect the head of data and summary
head(data_csv)
summary(data_csv)



## Preprocess 

#filter the data with include variable
# - extract rows with 1 in Include_FS_studies coloumn
data_csv_filtered <- data_csv[c(data_csv$Include_FS_studies == 1),]
data_csv_filtered <- data_csv_filtered[!is.na(data_csv_filtered$Include_FS_studies),]

# - extract rows with 1 in Include_FS_studies_euler_outliers_sibpairs_out
data_csv_filtered <- data_csv_filtered[c(data_csv_filtered$Include_FS_studies_euler_outliers_sibpairs_out == 1),]
data_csv_filtered <- data_csv_filtered[!is.na(data_csv_filtered$Include_FS_studies_euler_outliers_sibpairs_out),]

# filter a single groups data
data_csv_filtered <- data_csv_filtered[c(data_csv_filtered$HighRiskStatus_v11 == run_group),]
data_csv_filtered <- data_csv_filtered[!is.na(data_csv_filtered$HighRiskStatus_v11),]


#make new variables with shorter and contained names
# - tell r which variables are factors
datab = data_csv_filtered
datab$sex = as.factor(datab$Sex_child)
datab$site = as.factor(datab$MRI_site_v11)
datab$diag = as.factor(datab$ksads_any_diag_excl_elim_lft_v11)
datab$age = as.numeric(datab$MRI_age_v11)


#add a string encoding of sex for later convenience 
datab$sexs = datab$Sex_child
datab$sexs[datab$sexs == 0] = "female"
datab$sexs[datab$sexs == 1] = "male"




###### for loop that generates tables of global brain measure as y-variable 
#with each of the behavioural variables as covariates


#specify dependent variable (global variable) column names, i.e. left hand side of statistical models
model_vars = c("CBCL_ext_cg_v11","CBCL_int_cg_v11","CBCL_totsc_cg_v11","CGASx_v11")
brain_vars = c("BrainTotalVol", "CortexVol", "total_area", "mean_thickness", "eICV_samseg")

#covariates that are included in all models 
#encoded as a formula string
general_cov = c("age","site","TotalEulerNumber")
general_cov_formula = paste(general_cov,collapse = " + ")



#models for each sex

#### ANOVA table with male / female split
# and built models[sex][glob][yvar]

#define some relevant naming 
glob = c("without_eICV")
sex = c("female","male","both") #3 sexes are defined, 0=female, 1=male, 2=both

#split up the data into male and female 
data_sex1 = datab[c(datab$sex == 1),]
data_sex0 = datab[c(datab$sex == 0),]
dataf = list(data_sex0, data_sex1, datab)

models = list()

DF_xlsx = data.frame()

# loop that runs the model on each global variable, each sex with and without ICV
for (k in seq(1,length(model_vars))){
  print(model_vars[k])
  
  #loop that adds and analyze models WITH behav
  for (b in seq(1,length(brain_vars))){
    print(brain_vars[b])
    brain = brain_vars[b]
    
    for (ss in seq(1,3)){ #use the 3 defined sexes to know what data to use
      
      #extract the dataframe containing the relevant sex data
      df_sex = dataf[[ss]]
      
      
      for (gg in seq(1,1)){
        f = paste(model_vars[k],"~",brain,"+",general_cov_formula)
        
        if (sex[ss] == "both"){ #define a specific model for the case of both sex
          f = paste(f,"+ sex")
        }
        
        #run the model
        models[[sex[ss]]][[glob[gg]]][[brain_vars[b]]][[model_vars[k]]] = lm(f,data=df_sex)
        
        #save the current model into a temporary variable to make code less cluttered
        model_ana = models[[sex[ss]]][[glob[gg]]][[brain_vars[b]]][[model_vars[k]]]
        xvars = attributes(Anova(model_ana,type = "III"))$row.names
        
        
        
        #ANOVA sex separated
        mm=model_ana
        subs = length(model_ana$residuals)
        
        brain_F = Anova(mm,type = "III")$"F value"[xvars==brain]
        age_F = Anova(mm,type = "III")$"F value"[xvars=="age"]
        site_F = Anova(mm,type = "III")$"F value"[xvars=="site"]
        EulerNumber_F = Anova(mm,type = "III")$"F value"[xvars=="TotalEulerNumber"]
        
        brain_pv = Anova(mm,type = "III")$"Pr(>F)"[xvars==brain]
        age_pv = Anova(mm,type = "III")$"Pr(>F)"[xvars=="age"]
        site_pv = Anova(mm,type = "III")$"Pr(>F)"[xvars=="site"]
        EulerNumber_pv = Anova(mm,type = "III")$"Pr(>F)"[xvars=="TotalEulerNumber"]
        
        if (glob[gg] == "without_eICV"){
          #no global measure model
          glob_var_F = NA
          glob_var_pv = NA
        }       

        
        if (sex[ss] == "both"){
          sex_F = Anova(mm,type = "III")$"F value"[xvars=="sex"]
          sex_pv = Anova(mm,type = "III")$"Pr(>F)"[xvars=="sex"]
        }
        else {
          sex_F = NA
          sex_pv = NA
        }
        
        model_type = run_group
        
        
        #rows for ANOVA xlsx table
        rw_xlsx = list(model_vars[k], brain, 
                       sex_F, sex_pv,
                       age_F, age_pv, site_F, site_pv, 
                       EulerNumber_F, EulerNumber_pv, 
                       brain_F, brain_pv,
                       glob_var_F, glob_var_pv, "eICV_samseg",
                       gg-1,ss-1,model_type,subs)
        
        DF_xlsx = rbindlist(list(DF_xlsx, rw_xlsx))
        
      } # glob 
    } #sex 
  } #behav
} #model_vars

names(DF_xlsx)<-c("Model_yvar","Brain_var",
                  "Sex_Fval","Sex_pval",
                  "Age_Fval","Age_pval","Site_Fval","Site_pval",
                  "EulerNumber_Fval","Eulernumber_pval",
                  "Brain_Fval","Brain_pval",
                  "global_var_F","global_var_pv","global_var_name",
                  "global_var_in_model","sex","model_type","n_subjects")


DF_xlsx_glob0 = DF_xlsx[DF_xlsx$global_var_in_model == 0, ]
write_xlsx(DF_xlsx_glob0,paste(save_folder,model_ANOVA_without_glob,sep = ""))

#DF_xlsx_glob1 = DF_xlsx[DF_xlsx$global_var_in_model == 1, ]
#write_xlsx(DF_xlsx_glob1,paste(save_folder,model_ANOVA_with_glob,sep = ""))




#### Effect size ####


library(lsr)


DF = data.frame()


#all non-changing variable names
all_var_names = c("sex",general_cov)

glob = c("without_eICV")

for (k in seq(1,length(model_vars))){
  for (b in seq(1,length(brain_vars))){
    for (g in seq(1,length(glob))){
      for (s in seq(1,length(sex))){
        
        model_type = run_group
        model_eff = models[[sex[s]]][[glob[g]]][[brain_vars[b]]][[model_vars[k]]]
        
        
        #get model terms 
        xvars = attributes(etaSquared(model_eff, type = 3, anova = T))$dimnames[[1]]
        
        #get etaSquared outputs
        yvars = attributes(etaSquared(model_eff, type = 3, anova = T))$dimnames[[2]]
        
        eta_vals = as.numeric(etaSquared(model_eff, type = 3, anova = T)[,yvars="eta.sq.part"])
        
        all_col_names = c(all_var_names, brain_vars[b])
        
        var_idx = all_col_names %in% xvars
        var_idx2 = xvars %in% all_col_names
        
        row = rep(NA,length(all_col_names))
        
        xt = xvars[var_idx2]
        at = all_col_names[var_idx]
        et = eta_vals[var_idx2]
        
        row[var_idx] = et[as.numeric(sapply(at, function(x) { grep(x, xt) }))]
        
        
        
        #rows for effect size xlsx table
        DF = rbindlist(list(DF, as.list(c(model_vars[k],brain_vars[b],row,
                                          g-1,sex[s],model_type))))
        
      } #end s
    } #end g
  }
} #end i


all_col_names = c(all_var_names, "brain")

paES_col_names = paste(all_col_names, "_par_eta_sq",sep="")


names(DF) = c("model_yvar","Brain_var",paES_col_names,
              "globvar_in_model","sex","model_type")

write_xlsx(DF,paste(save_folder,effect_sizes_path,sep = ""))






# test some of the models that are to be run 


#model_eff = lm(BrainTotalVol ~ group*sex + age + site + TotalEulerNumber, data=datab)
model_eff = lm(BrainTotalVol ~ age + site + CBCL_ext_cg_v11 + TotalEulerNumber, data=dataf[[1]])
#emm_eff = emmeans(model_eff,specs="group")

#cohens d
#eff_size(emm_eff, sigma=sigma(model_eff), edf=df.residual(model_eff))

etaSquared(model_eff, type = 3, anova = T)



#with global covariate 
model_bvol_glob = lm(CBCL_ext_cg_v11 ~  BrainTotalVol + sex + age + TotalEulerNumber + site, data=datab)
Anova(model_bvol_glob,type="III")
etaSquared(model_bvol_glob, type = 3, anova = T)


model_bvol_glob = lm(BrainTotalVol ~  CBCL_ext_cg_v11 + age + TotalEulerNumber + site, data=dataf[[1]])
Anova(model_bvol_glob,type="III")
etaSquared(model_bvol_glob, type = 3, anova = T)


