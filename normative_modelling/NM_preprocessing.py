#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 19 10:11:29 2022

@author: simonyj
"""


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.model_selection import train_test_split, StratifiedShuffleSplit
from pcntoolkit.normative import estimate, evaluate
#from pcntoolkit.utils import create_bspline_basis, compute_MSLL
import pcntoolkit as pcn


data_path = "/mnt/projects/VIA11/FREESURFER/Stats/Data/VIA11_allkey_160621_FreeSurfer_pruned_20220509.csv"

save_data_dir = "/mnt/projects/VIA11/FREESURFER/Stats/Normative_modelling/data/"


#initial filtering
df = pd.read_csv(data_path) 
df = df.loc[ df['Include_FS_studies_euler_outliers_sibpairs_out'] == 1]
df = df.loc[ ~df['ksads_any_diag_excl_elim_lft_v11'].isna() ] #remove the 4 subs with nan in axis1

#extract the columns you need
new_df = df[["famlbnr","HighRiskStatus_v11","MRI_age_v11","ksads_any_diag_excl_elim_lft_v11","MRI_site_v11","Sex_child_v11","TotalEulerNumber","eICV_samseg"]]

#extract control group data for train/test splitting
K_df = new_df.loc[ new_df["HighRiskStatus_v11"].isin(["K"]) ]

#save the response variable
K_df_Y = K_df["eICV_samseg"]

#remove the response variable from covariates 
#remove info about group
K_df = K_df.drop(["HighRiskStatus_v11","eICV_samseg"],axis=1)

#define dataframes to stratify by
#the data to be stratified is only the control group data
stratify_by = ["Sex_child_v11","MRI_site_v11","ksads_any_diag_excl_elim_lft_v11"]
id_df = K_df["famlbnr"]
y_df = K_df[stratify_by]

#define train/test split with stratification
#here id_df is the indexes that should be stratified 
#y_df are the variables to stratify by
seed = 999
sss = StratifiedShuffleSplit(n_splits=1, test_size=0.2, random_state=seed)
train_idx, test_idx = next(iter(sss.split(id_df.to_numpy(), y_df.to_numpy())))


# Train DATA 


#setup the train data for modelling (from the control subjects dataframe)
K_df_X = K_df.drop(["famlbnr"],axis=1) # we dont need the via IDs for modelling
trainX = K_df_X.iloc[train_idx] #extract the training data
trainX['MRI_site_v11'] -= 1 #make site 0/1 binary instead of 1/2
trainY = K_df_Y.iloc[train_idx] #get the response variables for the training data


#save data for use in other scripts
trainX_path = save_data_dir + "trainX.txt"
trainY_path = save_data_dir + "trainY.txt"
trainX.to_csv(trainX_path,sep = ' ',header = True,index = True)
trainY.to_csv(trainY_path,sep = ' ',header = True,index = True)

#save data for modelling 
trainX_path = save_data_dir + "covariate_normsample.txt"
trainY_path = save_data_dir + "y_normsample.txt"

trainX.to_csv(trainX_path,sep = ' ',header = False,index = False)
trainY.to_csv(trainY_path,sep = ' ',header = False,index = False)



# Test DATA 

#extract test/target group data
test_df = new_df.loc[ ~new_df["HighRiskStatus_v11"].isin(["K"]) ]
test_Y = test_df["eICV_samseg"]
K_testY = K_df_Y.iloc[test_idx]

#make a whole test group both with target group and control
testY = pd.concat([test_Y,K_testY],axis=0)

test_df_HRS = test_df.drop(["eICV_samseg","famlbnr"],axis=1)
test_df = test_df.drop(["HighRiskStatus_v11","eICV_samseg","famlbnr"],axis=1)

K_testX = K_df_X.iloc[test_idx]
testX = pd.concat([test_df, K_testX],axis=0)
testX['MRI_site_v11'] -= 1
K_testX['MRI_site_v11'] -= 1 #get binary variable as [0,1] instead of [1,2]
test_df['MRI_site_v11'] -= 1



##save data for use in other scripts
K_testX_path = save_data_dir + "KtestX.txt"
K_testY_path = save_data_dir + "KtestY.txt"
K_testX.to_csv(K_testX_path,sep = ' ',header = True,index = True)
K_testY.to_csv(K_testY_path,sep = ' ',header = True,index = True)

SZBP_testX_path = save_data_dir + "SZBP_testX.txt"
SZBP_testY_path = save_data_dir + "SZBP_testY.txt"
test_df.to_csv(SZBP_testX_path,sep = ' ',header = True,index = True)
test_Y.to_csv(SZBP_testY_path,sep = ' ',header = True,index = True)

SZBP_testX_path = save_data_dir + "SZBP_testX_HRS.txt"
SZBP_testY_path = save_data_dir + "SZBP_testY_HRS.txt"
test_df_HRS.to_csv(SZBP_testX_path,sep = ' ',header = True,index = True)
test_Y.to_csv(SZBP_testY_path,sep = ' ',header = True,index = True)


testX_path = save_data_dir + "testX.txt"
testY_path = save_data_dir + "testY.txt"
testX.to_csv(testX_path,sep = ' ',header = True,index = True)
testY.to_csv(testY_path,sep = ' ',header = True, index = True)




#%% Multivariate data 

#extract the columns you need
all_cols = list(df.columns)

idx_vol = np.array([x.endswith("_volume") for x in all_cols])
cols = np.array(all_cols)[idx_vol]

m_df = df[cols]
X = m_df.to_numpy()

from sklearn.preprocessing import StandardScaler
pca_scaler = StandardScaler()
pca_scaler.fit(X)
PC1_norm = pca_scaler.transform(X)


from sklearn.decomposition import PCA
pca = PCA(n_components=1)
PC1 = pca.fit_transform(PC1_norm)

pca_df = pd.DataFrame(PC1, columns = ['PCA_volume'])



#train data

target_idx = np.invert(np.isin(np.arange(len(pca_df)), train_idx))

m_trainY = pca_df.loc[train_idx]
m_testY = pca_df.loc[target_idx]

m_trainY_path = save_data_dir + "m_trainY.txt"
m_testY_path = save_data_dir + "m_testY.txt"

m_trainY.to_csv(m_trainY_path,sep = ' ',header = True,index = True)
m_testY.to_csv(m_testY_path,sep = ' ',header = True,index = True)


m_BPY = pca_df.loc[ df['HighRiskStatus_v11'].isin(["BP"]).to_numpy() ]
m_SZY = pca_df.loc[ df['HighRiskStatus_v11'].isin(["SZ"]).to_numpy() ]
m_valY =  pca_df.loc[ test_idx ]

m_BPY_path = save_data_dir + "m_BPY.txt"
m_SZY_path = save_data_dir + "m_SZY.txt"
m_valY_path = save_data_dir + "m_valY.txt"

m_BPY.to_csv(m_BPY_path,sep = ' ', header = True, index = True)
m_SZY.to_csv(m_SZY_path,sep = ' ', header = True, index = True)
m_valY.to_csv(m_valY_path,sep = ' ', header = True, index = True)


#test data is identical to the data saved earlier 



#%%  test that stratification worked

fig, ax =  plt.subplots(2, len(stratify_by), figsize=(15, 8))

tt = [f"Train data (N control={len(trainX)})", f"Test data (N control={len(K_testX)})"]

for jj,dt in enumerate([trainX,K_testX]):
    ax[jj,0].set_ylabel(tt[jj])
    for ii,sb in enumerate(stratify_by):
        vals = list(dt[sb].unique())
        valss = sorted(vals)
        counts = dt[sb].value_counts()
        
        cc = [counts[x] for x in valss] #indices in vals and cc are now synced
        valss = [str(int(x)) for x in valss]
        
        ax[jj,ii].bar( valss, cc )
        ax[jj,ii].set_title(sb)

plt.savefig('/mnt/projects/VIA11/FREESURFER/Stats/Normative_modelling/train_test_split.png')



