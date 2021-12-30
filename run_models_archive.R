###### make inference with models
models = list()
models_glob = list()

DF = data.frame()
DFs = data.frame()
DF_glob = data.frame()
DFs_glob = data.frame()

for (i in seq(1,2)){
  for (j in seq(1,3)){
    model_yvars = col_names_list[[paste(h[i],m[j],sep = "_")]] 
    for (k in seq(1,length(model_yvars))){
      for (mi in seq(1,2)){
      
      if (mi == 1){
        #no global measure model
        f = paste(model_yvars[k],"~","-1","+","group*sex","+","age","+","(1| site)")
        models[[h[i]]][[m[j]]][[model_yvars[k]]] = lmer(f,data=datab)
        xvars = attributes(anova(models[[h[i]]][[m[j]]][[model_yvars[k]]]))$row.names
        model_ana = models[[h[i]]][[m[j]]][[model_yvars[k]]]
        
        if (anova(models[[h[i]]][[m[j]]][[model_yvars[k]]])$"Pr(>F)"[xvars=="group:sex"] > 0.05){
          models[[h[i]]][[m[j]]][[model_yvars[k]]] = update(models[[h[i]]][[m[j]]][[model_yvars[k]]],~.-group:sex)
        } 

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
        models_glob[[h[i]]][[m[j]]][[model_yvars[k]]] = lmer(f2,data=datab)
        xvars = attributes(anova(models_glob[[h[i]]][[m[j]]][[model_yvars[k]]]))$row.names
        model_ana = models_glob[[h[i]]][[m[j]]][[model_yvars[k]]]
        
        if (anova(models_glob[[h[i]]][[m[j]]][[model_yvars[k]]])$"Pr(>F)"[xvars=="group:sex"] > 0.05){
          models_glob[[h[i]]][[m[j]]][[model_yvars[k]]] = update(models_glob[[h[i]]][[m[j]]][[model_yvars[k]]],~.-group:sex)
        }  
        
      }
        
      xvars = attributes(anova(model_ana))$row.names
      pv_group_sex = anova(model_ana)$"Pr(>F)"[xvars=="group:sex"]
      
      if (anova(model_ana)$"Pr(>F)"[xvars=="group:sex"] > 0.05){
        
        model_ana = update(model_ana,~.-group:sex)
        
        ls = lsmeans(model_ana,pairwise~"group")
        c = ls$contrasts
        
        K_emm = summary(ls)$lsmeans$lsmean[2]
        sex1_K_emm = "n/a"
        
        #raw contrasts
        BP_diff = summary(c)$estimate[1]
        BP_diff_LCL = confint(c)$lower.CL[1]
        BP_diff_UCL = confint(c)$upper.CL[1]
        SZ_diff = -summary(c)$estimate[3]
        SZ_diff_LCL = -confint(c)$lower.CL[3]
        SZ_diff_UCL = -confint(c)$upper.CL[3]
        
        #percent differece
        BP_diff_P = 100*summary(c)$estimate[1] /K_emm
        BP_diff_LCL_P = 100*confint(c)$lower.CL[1] /K_emm
        BP_diff_UCL_P = 100*confint(c)$upper.CL[1] /K_emm
        SZ_diff_P = -100*summary(c)$estimate[3] /K_emm
        SZ_diff_LCL_P = -100*confint(c)$lower.CL[3] /K_emm
        SZ_diff_UCL_P = -100*confint(c)$upper.CL[3] /K_emm
        
        sex1_BP_diff = "n/a"
        sex1_BP_diff_LCL = "n/a"
        sex1_BP_diff_UCL = "n/a"
        
        sex1_SZ_diff = "n/a"
        sex1_SZ_diff_LCL = "n/a"
        sex1_SZ_diff_UCL = "n/a"
        
        sex1_BP_diff_P = "n/a"
        sex1_BP_diff_LCL_P = "n/a"
        sex1_BP_diff_UCL_P = "n/a"
        
        sex1_SZ_diff_P = "n/a"
        sex1_SZ_diff_LCL_P = "n/a"
        sex1_SZ_diff_UCL_P = "n/a"
        
      }
      else{
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
        
      }
      
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
        model_diff_pv = anova(models_glob[[h[i]]][[m[j]]][[model_yvars[k]]],models[[h[i]]][[m[j]]][[model_yvars[k]]])$"Pr(>Chisq)"[2]
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









###### make inference with models

models = list()
models_glob = list()

for (i in seq(1,2)){
  for (j in seq(1,3)){
    model_yvars = col_names_list[[paste(h[i],m[j],sep = "_")]] 
    for (k in seq(1,length(model_yvars))){
      
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
      
      models_glob[[h[i]]][[m[i]]][[model_yvars[k]]] = lmer(f2,data=datab)
      xvars = attributes(anova(models_glob[[h[i]]][[m[i]]][[model_yvars[k]]]))$row.names
      
      if (anova(models_glob[[h[i]]][[m[i]]][[model_yvars[k]]])$"Pr(>F)"[xvars=="group:sex"] > 0.05){
        models_glob[[h[i]]][[m[i]]][[model_yvars[k]]] = update(models_glob[[h[i]]][[m[i]]][[model_yvars[k]]],~.-group:sex)
      }
      
      #no global measure model
      f = paste(model_yvars[k],"~","-1","+","group*sex","+","age","+","(1| site)")
      models[[h[i]]][[m[i]]][[model_yvars[k]]] = lmer(f,data=datab)
      xvars = attributes(anova(models[[h[i]]][[m[i]]][[model_yvars[k]]]))$row.names
      
      pv_group_sex = anova(models[[h[i]]][[m[i]]][[model_yvars[k]]])$"Pr(>F)"[xvars=="group:sex"]
      
      if (anova(models[[h[i]]][[m[i]]][[model_yvars[k]]])$"Pr(>F)"[xvars=="group:sex"] > 0.05){
        models[[h[i]]][[m[i]]][[model_yvars[k]]] = update(models[[h[i]]][[m[i]]][[model_yvars[k]]],~.-group:sex)
        ls = lsmeans(models[[h[i]]][[m[i]]][[model_yvars[k]]],pairwise~"group")
        c = ls$contrasts
        
        K_emm = summary(ls)$lsmeans$lsmean[2]
        sex1_K_emm = "n/a"
        
        #raw contrasts
        BP_diff = summary(c)$estimate[1]
        BP_diff_LCL = confint(c)$lower.CL[1]
        BP_diff_UCL = confint(c)$upper.CL[1]
        SZ_diff = -summary(c)$estimate[3]
        SZ_diff_LCL = -confint(c)$lower.CL[3]
        SZ_diff_UCL = -confint(c)$upper.CL[3]
        
        #percent differece
        BP_diff_P = summary(c)$estimate[1] /K_emm
        BP_diff_LCL_P = confint(c)$lower.CL[1] /K_emm
        BP_diff_UCL_P = confint(c)$upper.CL[1] /K_emm
        SZ_diff_P = -summary(c)$estimate[3] /K_emm
        SZ_diff_LCL_P = -confint(c)$lower.CL[3] /K_emm
        SZ_diff_UCL_P = -confint(c)$upper.CL[3] /K_emm
        
        sex1_BP_diff = "n/a"
        sex1_BP_diff_LCL = "n/a"
        sex1_BP_diff_UCL = "n/a"
        
        sex1_SZ_diff = "n/a"
        sex1_SZ_diff_LCL = "n/a"
        sex1_SZ_diff_UCL = "n/a"
        
        sex1_BP_diff_P = "n/a"
        sex1_BP_diff_LCL_P = "n/a"
        sex1_BP_diff_UCL_P = "n/a"
        
        sex1_SZ_diff_P = "n/a"
        sex1_SZ_diff_LCL_P = "n/a"
        sex1_SZ_diff_UCL_P = "n/a"
        
      }
      else{
        ls = lsmeans(models[[h[i]]][[m[i]]][[model_yvars[k]]],pairwise~"group", by = "sex")
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
        BP_diff_P = summary(c)$estimate[summary(c)$sex == 0][1] /K_emm
        BP_diff_LCL_P = confint(c)$lower.CL[confint(c)$sex == 0][1] /K_emm
        BP_diff_UCL_P = confint(c)$upper.CL[confint(c)$sex == 0][1] /K_emm
        
        SZ_diff_P = -summary(c)$estimate[summary(c)$sex == 0][3] /K_emm
        SZ_diff_LCL_P = -confint(c)$lower.CL[confint(c)$sex == 0][3] /K_emm
        SZ_diff_UCL_P = -confint(c)$upper.CL[confint(c)$sex == 0][3] /K_emm
        
        sex1_BP_diff_P = summary(c)$estimate[summary(c)$sex == 1][1] /sex1_K_emm
        sex1_BP_diff_LCL_P = confint(c)$lower.CL[confint(c)$sex == 1][1] /sex1_K_emm
        sex1_BP_diff_UCL_P = confint(c)$upper.CL[confint(c)$sex == 1][1] /sex1_K_emm
        
        sex1_SZ_diff_P = -summary(c)$estimate[summary(c)$sex == 1][3] /sex1_K_emm
        sex1_SZ_diff_LCL_P = -confint(c)$lower.CL[confint(c)$sex == 1][3] /sex1_K_emm
        sex1_SZ_diff_UCL_P = -confint(c)$upper.CL[confint(c)$sex == 1][3] /sex1_K_emm
        
      }
      
      model_diff_pv = anova(models_glob[[h[i]]][[m[i]]][[model_yvars[k]]],models[[h[i]]][[m[i]]][[model_yvars[k]]])$"Pr(>Chisq)"[2]
      
      pv_group = anova(models[[h[i]]][[m[i]]][[model_yvars[k]]])$"Pr(>F)"[1]
      
      
      if (i==1 && j==1 && k==1){
        df<-data.frame(model_yvars[k],model_diff_pv,pv_group,pv_group_sex, K_emm, sex1_K_emm,
                       BP_diff, BP_diff_LCL, BP_diff_UCL, 
                       SZ_diff, SZ_diff_LCL, SZ_diff_UCL, 
                       sex1_BP_diff, sex1_BP_diff_LCL, sex1_BP_diff_UCL, 
                       sex1_SZ_diff, sex1_SZ_diff_LCL, sex1_SZ_diff_UCL,
                       BP_diff_P, BP_diff_LCL_P, BP_diff_UCL_P, 
                       SZ_diff_P, SZ_diff_LCL_P, SZ_diff_UCL_P, 
                       sex1_BP_diff_P, sex1_BP_diff_LCL_P, sex1_BP_diff_UCL_P, 
                       sex1_SZ_diff_P, sex1_SZ_diff_LCL_P, sex1_SZ_diff_UCL_P)
        names(df)<-c("model_yvar","model_diff_pv","Group_p_value","Group_sex_p_value","K_emm", "sex1_K_emm",
                     "BP_diff","BP_diff_LCL","BP_diff_UCL",
                     "SZ_diff","SZ_diff_LCL","SZ_diff_UCL",
                     "sex1_BP_diff","sex1_BP_diff_LCL","sex1_BP_diff_UCL",
                     "sex1_SZ_diff","sex1_SZ_diff_LCL","sex1_SZ_diff_UCL",
                     "BP_diff_P", "BP_diff_LCL_P", "BP_diff_UCL_P", 
                     "SZ_diff_P", "SZ_diff_LCL_P", "SZ_diff_UCL_P", 
                     "sex1_BP_diff_P", "sex1_BP_diff_LCL_P", "sex1_BP_diff_UCL_P", 
                     "sex1_SZ_diff_P", "sex1_SZ_diff_LCL_P", "sex1_SZ_diff_UCL_P")
        
        dfs<-data.frame(model_yvars[k], K_emm, sex1_K_emm, 
                        BP_diff, BP_diff_LCL, BP_diff_UCL, 
                        SZ_diff, SZ_diff_LCL, SZ_diff_UCL,
                        BP_diff_P, BP_diff_LCL_P, BP_diff_UCL_P, 
                        SZ_diff_P, SZ_diff_LCL_P, SZ_diff_UCL_P,0)
        names(dfs)<-c("model_yvar","K_emm", "sex1_K_emm",
                      "BP_diff","BP_diff_LCL","BP_diff_UCL",
                      "SZ_diff","SZ_diff_LCL","SZ_diff_UCL",
                      "BP_diff_P","BP_diff_LCL_P","BP_diff_UCL_P",
                      "SZ_diff_P","SZ_diff_LCL_P","SZ_diff_UCL_P","sex")
        dfs[nrow(dfs) + 1,] = c(model_yvars[k], K_emm, sex1_K_emm, 
                                sex1_BP_diff, sex1_BP_diff_LCL, sex1_BP_diff_UCL,
                                sex1_SZ_diff, sex1_SZ_diff_LCL, sex1_SZ_diff_UCL,
                                sex1_BP_diff_P, sex1_BP_diff_LCL_P, sex1_BP_diff_UCL_P,
                                sex1_SZ_diff_P, sex1_SZ_diff_LCL_P, sex1_SZ_diff_UCL_P,1)
      }
      else{
        df[nrow(df) + 1,] = c(model_yvars[k],model_diff_pv,pv_group,pv_group_sex, K_emm, sex1_K_emm,
                              BP_diff, BP_diff_LCL, BP_diff_UCL, 
                              SZ_diff, SZ_diff_LCL, SZ_diff_UCL, 
                              sex1_BP_diff, sex1_BP_diff_LCL, sex1_BP_diff_UCL, 
                              sex1_SZ_diff, sex1_SZ_diff_LCL, sex1_SZ_diff_UCL,
                              BP_diff_P, BP_diff_LCL_P, BP_diff_UCL_P, 
                              SZ_diff_P, SZ_diff_LCL_P, SZ_diff_UCL_P, 
                              sex1_BP_diff_P, sex1_BP_diff_LCL_P, sex1_BP_diff_UCL_P, 
                              sex1_SZ_diff_P, sex1_SZ_diff_LCL_P, sex1_SZ_diff_UCL_P)
        
        
        
        dfs[nrow(dfs) + 1,] = c(model_yvars[k], K_emm, sex1_K_emm, 
                                sex1_BP_diff, sex1_BP_diff_LCL, sex1_BP_diff_UCL,
                                sex1_SZ_diff, sex1_SZ_diff_LCL, sex1_SZ_diff_UCL,
                                sex1_BP_diff_P, sex1_BP_diff_LCL_P, sex1_BP_diff_UCL_P,
                                sex1_SZ_diff_P, sex1_SZ_diff_LCL_P, sex1_SZ_diff_UCL_P,0)
        
        dfs[nrow(dfs) + 1,] = c(model_yvars[k], K_emm, sex1_K_emm, 
                                sex1_BP_diff, sex1_BP_diff_LCL, sex1_BP_diff_UCL,
                                sex1_SZ_diff, sex1_SZ_diff_LCL, sex1_SZ_diff_UCL,
                                sex1_BP_diff_P, sex1_BP_diff_LCL_P, sex1_BP_diff_UCL_P,
                                sex1_SZ_diff_P, sex1_SZ_diff_LCL_P, sex1_SZ_diff_UCL_P,1)
        
      }
      
      
    } #k
  }
}
