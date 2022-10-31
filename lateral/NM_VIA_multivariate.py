#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  3 11:22:30 2022

@author: simonyj
"""

#%%

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
#from pcntoolkit.normative import estimate, evaluate
#from pcntoolkit.utils import create_bspline_basis, compute_MSLL
import pcntoolkit as pcn

import os 

from sklearn.preprocessing import StandardScaler


#%%

m_trainY_path = '/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/global_variables/m_trainY.txt'
m_testY_path = '/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/global_variables/m_testY.txt'
trainY = pd.read_csv(m_trainY_path,sep=' ') 
testY = pd.read_csv(m_testY_path,sep=' ') 


trainX_path = '/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/global_variables/trainX.txt'
trainX = pd.read_csv(trainX_path,sep=' ') 

testX_path = '/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/global_variables/testX.txt'
testX = pd.read_csv(testX_path,sep=' ') 



K_testX_path = '/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/global_variables/KtestX.txt'
K_testY_path = '/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/global_variables/m_valY.txt'
K_testX = pd.read_csv(K_testX_path,sep=' ') 
K_testY = pd.read_csv(K_testY_path,sep=' ') 


SZBP_testX_path = '/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/global_variables/SZBP_testX.txt'
SZBP_testY_path = '/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/global_variables/SZBP_KtestY.txt'
SZBP_testX = pd.read_csv(SZBP_testX_path,sep=' ') 
SZBP_testY = pd.read_csv(SZBP_testY_path,sep=' ') 


SZBP_testX_path = '/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/global_variables/SZBP_testX_HRS.txt'
SZBP_testX_HRS = pd.read_csv(SZBP_testX_path,sep=' ') 


BP_testY_path = '/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/global_variables/m_BPY.txt'
BP_testY = pd.read_csv(BP_testY_path,sep=' ') 
BP_testX = SZBP_testX_HRS.loc[ SZBP_testX_HRS["HighRiskStatus_v11"].isin(["BP"]) ]
BP_testX = BP_testX.drop(["HighRiskStatus_v11"],axis=1)

SZ_testY_path = '/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/global_variables/m_SZY.txt'
SZ_testY = pd.read_csv(SZ_testY_path,sep=' ') 
SZ_testX = SZBP_testX_HRS.loc[ SZBP_testX_HRS["HighRiskStatus_v11"].isin(["SZ"]) ]
SZ_testX = SZ_testX.drop(["HighRiskStatus_v11"],axis=1)



#%% load data into the right format

#training data
trainX_model = trainX.iloc[:,1:]
trainY_model = trainY[['PCA_volume']]

#testing data
testX_model = testX.iloc[:,1:]
testY_model = testY[['PCA_volume']]

#BP only
testBPX_model = BP_testX.iloc[:,1:]
testBPY_model = BP_testY[['PCA_volume']]

#SZ only
testSZX_model = SZ_testX.iloc[:,1:]
testSZY_model = SZ_testY[['PCA_volume']]

#validation data
valX_model = K_testX.iloc[:,1:]
valY_model = K_testY[['PCA_volume']]



### Training 
#normalize X
train_scaler = StandardScaler()
train_scaler.fit(trainX_model.to_numpy())
X_norm = train_scaler.transform(trainX_model.to_numpy())

#normalize y
train_scalery = StandardScaler()
train_scalery.fit(trainY_model.to_numpy())
Y_norm = train_scalery.transform(trainY_model.to_numpy())


### Testing
#normalize X
test_scaler = StandardScaler()
test_scaler.fit(testX_model.to_numpy())
testX_norm = test_scaler.transform(testX_model.to_numpy())
#testX_norm = train_scaler.transform(testX_model.to_numpy())

#normalize y
test_scalery = StandardScaler()
test_scalery.fit(testY_model.to_numpy())
testY_norm = test_scalery.transform(testY_model.to_numpy())

#Validation data 
#normalize X
valX_norm = test_scaler.transform(valX_model.to_numpy())
#normalize y
valY_norm = test_scalery.transform(valY_model.to_numpy())

#Testing data
#normalize X
testBPX_norm = test_scaler.transform(testBPX_model.to_numpy())
#normalize y
testBPY_norm = test_scalery.transform(testBPY_model.to_numpy())

#normalize X
testSZX_norm = test_scaler.transform(testSZX_model.to_numpy())
#normalize y
testSZY_norm = test_scalery.transform(testSZY_model.to_numpy())





#%% Estimate BLR instead

# set this path to wherever your ROI_models folder is located (where you copied all of the covariate & response text files to in Step 4)
data_dir = '/home/simonyj/GPR_multivariate'
os.makedirs(data_dir, exist_ok=True)

file_suffix = ["train","test","val","BP","SZ"]


X_files = [X_norm, testX_norm, valX_norm, testBPX_norm, testSZX_norm]#[trainX_model.to_numpy(), testX_model.to_numpy(), valX_model.to_numpy(), test2X_model.to_numpy()]#, X_forward] #as numpy arrays 
Y_files = [Y_norm, testY_norm, valY_norm, testBPY_norm, testSZY_norm]#[trainY_model.to_numpy(), testY_model.to_numpy(), valY_model.to_numpy(), test2Y_model.to_numpy()]#, X_forward] #as numpy arrays 



for X,Y,suf in zip(X_files,Y_files,file_suffix):

    if suf != "forward":
        resp_file = os.path.join(data_dir, f"resp_{suf}.txt")
        np.savetxt(resp_file, Y)
    
    cov_file = os.path.join(data_dir, f"cov_gpr_{suf}.txt")
    np.savetxt(cov_file, X)
    

cov_file_tr = os.path.join(data_dir, 'cov_gpr_train.txt')
cov_file_te = os.path.join(data_dir, 'cov_gpr_test.txt')
cov_file_val = os.path.join(data_dir, 'cov_gpr_val.txt')
cov_file_BP = os.path.join(data_dir, 'cov_gpr_BP.txt')
cov_file_SZ = os.path.join(data_dir, 'cov_gpr_SZ.txt')


resp_file_tr = os.path.join(data_dir, 'resp_train.txt')
resp_file_te = os.path.join(data_dir, 'resp_test.txt')
resp_file_val = os.path.join(data_dir, 'resp_val.txt')
resp_file_BP = os.path.join(data_dir, 'resp_BP.txt')
resp_file_SZ = os.path.join(data_dir, 'resp_SZ.txt')



# run a model on all test
yhat_te, s2_te, nm, Z, metrics_te = pcn.normative.estimate(covfile = cov_file_tr,
                                                           respfile = resp_file_tr,
                                                           
                                                           testcov = cov_file_te,
                                                           testresp = resp_file_te, 
                                                           
                                                           optimizer = 'cg',
                                                           cvfolds = None,
                                                           alg = 'gpr',
                                                           saveoutput=False)



yhat_val, s2_val, nm, Z_val, metrics_val = pcn.normative.estimate(covfile = cov_file_tr,
                                                           respfile = resp_file_tr,
                                                           
                                                           testcov = cov_file_val,
                                                           testresp = resp_file_val, 
                                                           
                                                           optimizer = 'cg',
                                                           cvfolds = None,
                                                           alg = 'gpr',
                                                           saveoutput=False)


yhat_val, s2_val, nm, Z_BP, metrics_val = pcn.normative.estimate(covfile = cov_file_tr,
                                                           respfile = resp_file_tr,
                                                           
                                                           testcov = cov_file_BP,
                                                           testresp = resp_file_BP, 
                                                           
                                                           optimizer = 'cg',
                                                           cvfolds = None,
                                                           alg = 'gpr',
                                                           saveoutput=False)


yhat_val, s2_val, nm, Z_SZ, metrics_val = pcn.normative.estimate(covfile = cov_file_tr,
                                                           respfile = resp_file_tr,
                                                           
                                                           testcov = cov_file_SZ,
                                                           testresp = resp_file_SZ, 
                                                           
                                                           optimizer = 'cg',
                                                           cvfolds = None,
                                                           alg = 'gpr',
                                                           saveoutput=False)



#%%

plt.hist(Z,bins=10)
plt.xlim(Z.min(),Z.max())
plt.title("Multivariate target group")
plt.axvline(x=0,color='black',linestyle='--')

plt.savefig('/mnt/projects/VIA11/FREESURFER/Stats/Normative_modelling/multivariate_GP_Z_scores_plot_target.png')


#%% Histogram plots 
fig = plt.figure(figsize=(20,4))

ax = fig.add_subplot(1, 4, 1)
ax.hist(Z_val,bins=10)
ax.set_xlim(-3.5,3.5)
ax.set_title("Validation group (only control)")
ax.axvline(x=0,color='black',linestyle='--')

ax = fig.add_subplot(1, 4, 2)
ax.hist(Z,bins=10)
ax.set_xlim(-3.5,3.5)
ax.set_title("Target group")
ax.axvline(x=0,color='black',linestyle='--')

ax = fig.add_subplot(1, 4, 3)
ax.hist(Z_BP,bins=10)
ax.set_xlim(-3.5,3.5)
ax.set_title("BP group")
ax.axvline(x=0,color='black',linestyle='--')

ax = fig.add_subplot(1, 4, 4)
ax.hist(Z_SZ,bins=10)
ax.set_xlim(-3.5,3.5)
ax.set_title("SZ group")
ax.axvline(x=0,color='black',linestyle='--')

fig.suptitle('GP based Z-scores of test groups',fontsize=15)

plt.savefig('/mnt/projects/VIA11/FREESURFER/Stats/Normative_modelling/multivariate_GP_Z_scores_plot.png')



