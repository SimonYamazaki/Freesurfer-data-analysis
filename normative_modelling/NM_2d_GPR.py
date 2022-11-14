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
import seaborn as sns
#from pcntoolkit.normative import estimate, evaluate
#from pcntoolkit.utils import create_bspline_basis, compute_MSLL
import pcntoolkit as pcn

#%%

trainX_path = '/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/global_variables/trainX.txt'
trainY_path = '/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/global_variables/trainY.txt'
trainX = pd.read_csv(trainX_path,sep=' ') 
trainY = pd.read_csv(trainY_path,sep=' ') 

testX_path = '/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/global_variables/testX.txt'
testY_path = '/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/global_variables/testY.txt'
testX = pd.read_csv(testX_path,sep=' ') 
testY = pd.read_csv(testY_path,sep=' ') 


K_testX_path = '/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/global_variables/KtestX.txt'
K_testY_path = '/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/global_variables/KtestY.txt'
K_testX = pd.read_csv(K_testX_path,sep=' ') 
K_testY = pd.read_csv(K_testY_path,sep=' ') 


SZBP_testX_path = '/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/global_variables/SZBP_testX.txt'
SZBP_testY_path = '/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/global_variables/SZBP_KtestY.txt'
SZBP_testX = pd.read_csv(SZBP_testX_path,sep=' ') 
SZBP_testY = pd.read_csv(SZBP_testY_path,sep=' ') 



#%% load data into the right format

#training data
trainX_model = trainX.iloc[:,1:3]
trainY_model = trainY[['eICV_samseg']]

#testing data
testX_model = testX.iloc[:,1:3]
testY_model = testY[['eICV_samseg']]

#only BP and SZ test data
test2X_model = SZBP_testX.iloc[:,1:3]
test2Y_model = SZBP_testY[['eICV_samseg']]

#validation data
valX_model = K_testX.iloc[:,1:3]
valY_model = K_testY[['eICV_samseg']]


### Training 
#normalize X
from sklearn.preprocessing import StandardScaler
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


#normalize X
valX_norm = test_scaler.transform(valX_model.to_numpy())

#normalize y
valY_norm = test_scalery.transform(valY_model.to_numpy())


#normalize X
test2X_norm = test_scaler.transform(test2X_model.to_numpy())

#normalize y
test2Y_norm = test_scalery.transform(test2Y_model.to_numpy())




#forward model covariates
n_vals = 100
#x_age = np.linspace(min(trainX["MRI_age_v11"].min(),testX["MRI_age_v11"].min()), max(trainX["MRI_age_v11"].max(),testX["MRI_age_v11"].max()),n_vals)
x_age = np.linspace(10,14,n_vals)

x_sex = np.vstack( [np.zeros((n_vals,1)) , np.ones((n_vals,1))] )
X_forward = np.hstack( [np.tile(x_age,2)[:,np.newaxis], x_sex] )


#normalize X
X_norm_forward = train_scaler.transform(X_forward)



#%% Estimate BLR instead
import os 
from pcntoolkit.util.utils import create_bspline_basis, compute_MSLL, create_poly_basis, create_design_matrix
import sklearn


# set this path to wherever your ROI_models folder is located (where you copied all of the covariate & response text files to in Step 4)
data_dir = '/home/simonyj/GPR'
os.makedirs(data_dir, exist_ok=True)

file_suffix = ["train","test","val","BPSZ","forward"]


X_files = [X_norm, testX_norm, valX_norm, test2X_norm, X_norm_forward]#[trainX_model.to_numpy(), testX_model.to_numpy(), valX_model.to_numpy(), test2X_model.to_numpy()]#, X_forward] #as numpy arrays 
Y_files = [Y_norm, testY_norm, valY_norm, test2Y_norm, X_norm_forward]#[trainY_model.to_numpy(), testY_model.to_numpy(), valY_model.to_numpy(), test2Y_model.to_numpy()]#, X_forward] #as numpy arrays 



for X,Y,suf in zip(X_files,Y_files,file_suffix):

    if suf != "forward":
        resp_file = os.path.join(data_dir, f"resp_{suf}.txt")
        np.savetxt(resp_file, Y)
    
    cov_file = os.path.join(data_dir, f"cov_gpr_{suf}.txt")
    np.savetxt(cov_file, X)
    

cov_file_tr = os.path.join(data_dir, 'cov_gpr_train.txt')
cov_file_te = os.path.join(data_dir, 'cov_gpr_test.txt')
cov_file_val = os.path.join(data_dir, 'cov_gpr_val.txt')
cov_file_forward = os.path.join(data_dir, 'cov_gpr_forward.txt')
cov_file_BPSZ = os.path.join(data_dir, 'cov_gpr_BPSZ.txt')

resp_file_tr = os.path.join(data_dir, 'resp_train.txt')
resp_file_te = os.path.join(data_dir, 'resp_test.txt')
resp_file_val = os.path.join(data_dir, 'resp_val.txt')
resp_file_BPSZ = os.path.join(data_dir, 'resp_BPSZ.txt')



# run a model on all test
yhat_val, s2_val, nm, Z_val, metrics_te = pcn.normative.estimate(covfile = cov_file_tr,
                                                           respfile = resp_file_tr,
                                                           
                                                           testcov = cov_file_val,
                                                           testresp = resp_file_val, 
                                                           
                                                           optimizer = 'cg',
                                                           cvfolds = None,
                                                           alg = 'gpr',
                                                           #outputsuffix = '_gpr2dtest',
                                                           saveoutput=False)


yhat_BPSZ, s2_BPSZ, nm, Z_BPSZ, metrics_te = pcn.normative.estimate(covfile = cov_file_tr,
                                                           respfile = resp_file_tr,
                                                           
                                                           testcov = cov_file_BPSZ,
                                                           testresp = resp_file_BPSZ, 
                                                           
                                                           optimizer = 'cg',
                                                           cvfolds = None,
                                                           alg = 'gpr',
                                                           #outputsuffix = '_gpr2dtest',
                                                           saveoutput=False)


yhat_te, s2_te, nm, Z, metrics_te = pcn.normative.estimate(covfile = cov_file_tr,
                                                           respfile = resp_file_tr,
                                                           
                                                           testcov = cov_file_te,
                                                           testresp = resp_file_te, 
                                                           
                                                           optimizer = 'cg',
                                                           cvfolds = None,
                                                           alg = 'gpr',
                                                           #outputsuffix = '_gpr2dtest',
                                                           saveoutput=False)

pcn.normative.estimate(covfile = cov_file_tr, 
                       respfile = resp_file_tr,
                       testcov = cov_file_forward,                 
                        cvfolds = None,
                        alg = 'gpr',
                        outputsuffix = '_gpr2dtest',
                        saveoutput=True)


#%% 

yhat = pd.read_csv('/home/simonyj/yhat_gpr2dtest.txt', sep = ' ', header=None).to_numpy()
ys2 = pd.read_csv('/home/simonyj/ys2_gpr2dtest.txt', sep = ' ', header=None).to_numpy()

sd = np.sqrt(ys2)

#extract data from a particular categorical variable 
col_names = np.array(list(trainX.columns)[1:])
var = 0
var_col_idx = 1


#x_forward_sex0 = X_norm_forward[ X_forward[:,var_col_idx]==var, 0 ]
x_forward_sex0 = X_forward[ X_forward[:,var_col_idx]==var, 0 ]
yhat_sex0 = yhat[ X_forward[:,var_col_idx]==var ]

x_forward_plot = x_forward_sex0

yhat_plot = yhat_sex0
sd_plot = sd[ X_forward[:,var_col_idx]==var ]
#yhat_plot = train_scalery.inverse_transform( yhat_sex0 )
#sd_plot = np.sqrt(train_scalery.inverse_transform( ys2[ X_forward[:,var_col_idx]==var ] ))


dataX = trainX_model.to_numpy()[ trainX_model.to_numpy()[:,var_col_idx]==var, 0 ]
#norm_dataX = X_norm[trainX_model.to_numpy()[:,var_col_idx]==var, 0]
norm_dataY = Y_norm[trainX_model.to_numpy()[:,var_col_idx]==var]

val_dataX = valX_model.to_numpy()[ valX_model.to_numpy()[:,var_col_idx]==var, 0 ]
#test_norm_dataX = valX_norm[valX_model.to_numpy()[:,var_col_idx]==var, 0]
test_norm_dataY = valY_norm[valX_model.to_numpy()[:,var_col_idx]==var]

test_data2X = test2X_model.to_numpy()[ test2X_model.to_numpy()[:,var_col_idx]==var, 0 ]
#test_norm_data2X = test2X_norm[test2X_model.to_numpy()[:,var_col_idx]==var, 0]
test_norm_data2Y = test2Y_norm[test2X_model.to_numpy()[:,var_col_idx]==var]



fig = plt.figure(figsize=(10,5))
ax = fig.add_subplot(1, 1, 1)

ax.plot(x_forward_plot, yhat_plot)
ax.fill_between(x_forward_plot, (yhat_plot-2*sd_plot).squeeze(), (yhat_plot+2*sd_plot).squeeze(),alpha=0.4,label='95% (2*sd)')

ax.scatter(dataX,norm_dataY,color='red',label='K train')
ax.scatter(val_dataX,test_norm_dataY,color='green',label='K validation',marker='+')
ax.scatter(test_data2X,test_norm_data2Y,color='blue',label='SZ+BP',s=10)


ax.set_title(f"{col_names[var_col_idx]}={var}")
ax.set_xlabel("Age")
ax.set_ylabel("Brain Volume")

#ax.scatter(norm_dataX,norm_dataY,color='red',label='K train',s=10) #trainX['MRI_age_v11'],trainY['eICV_samseg']
#ax.scatter(test_norm_dataX,test_norm_dataY,color='green',label='K validation',marker='+')
#ax.scatter(test_norm_data2X,test_norm_data2Y,color='blue',label='SZ/BP test set',marker='+')

#ax.scatter(trainX['MRI_age_v11'],trainY['eICV_samseg'],color='red',label='K train')
#ax.scatter(K_testX['MRI_age_v11'],K_testY['eICV_samseg'],color='green',label='K validation',marker='+')
#ax.scatter(SZBP_testX['MRI_age_v11'],SZBP_testY['eICV_samseg'],color='blue',label='SZ+BP',s=10)

ax.legend()


#plt.savefig('/mnt/projects/VIA11/FREESURFER/Stats/Normative_modelling/GP_model_plot.png')



#%%
fig = plt.figure(figsize=(10,4))

ax = fig.add_subplot(1, 2, 1)
ax.hist(Z_BPSZ,bins=10)
ax.set_xlim(Z.min(),Z.max())
ax.set_title("Target group")

ax = fig.add_subplot(1, 2, 2)
ax.hist(Z_val,bins=10)
ax.set_xlim(Z.min(),Z.max())
ax.set_title("Validation group (only control)")

fig.suptitle('GP based Z-scores of test groups',fontsize=15)

#manual computation of Z 
#Z[ts, nz[i]] = (Ytest - Yhat[ts, nz[i]]) / np.sqrt(S2[ts, nz[i]]) 
#dobbelttjek evt med den BLR model du har

plt.savefig('/mnt/projects/VIA11/FREESURFER/Stats/Normative_modelling/GP_Z_scores_plot.png')



#%% 

from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import WhiteKernel, RBF, ConstantKernel, DotProduct

X = trainX_model.to_numpy()#[:,0][:,np.newaxis]
y = trainY_model.to_numpy() #/ trainY_model.to_numpy().std()#+ resp_mean

kernel = DotProduct() + ConstantKernel(constant_value=1e1, constant_value_bounds=(1e-3, 1e2)) * RBF(length_scale=1, length_scale_bounds=(1e-1, 1e5)) + WhiteKernel(noise_level=1, noise_level_bounds=(1e-3, 1e5))

gpr = GaussianProcessRegressor(kernel=kernel, normalize_y=True).fit(X, y) #, alpha=3e5

#m, sd = gpr.predict(X_forward[:100][:,np.newaxis], return_std=True)
m, sd = gpr.predict(X_forward[:100,:], return_std=True)


plt.plot(X_forward[:100,0],m)
plt.fill_between(X_forward[:100,0], (m.squeeze()-2*sd), (m.squeeze()+2*sd),alpha=0.5,label='95% (2*sd)')
#plt.plot(X_forward[:100,0], (m.squeeze()-200*sd))
#plt.plot(X_forward[:100,0], (m.squeeze()+200*sd))
#plt.ylim([1e6, 2e6])
plt.scatter(trainX['MRI_age_v11'],trainY['eICV_samseg'],color='red',label='K train',s=10)
#plt.scatter(trainX['MRI_age_v11'], (trainY['eICV_samseg']-resp_mean)/trainY_model.to_numpy().std(),color='red',label='K train',s=10)

#plt.scatter(testX['MRI_age_v11'],testY['eICV_samseg'],color='green',label='test',s=10)


#Xtest = testX[['MRI_age_v11','ksads_any_diag_excl_elim_lft_v11']].to_numpy()
Xtest = testX[['MRI_age_v11']].to_numpy()
Ytest = testY['eICV_samseg'].to_numpy()

#m, sd = gpr.predict(Xtest, return_std=True)
#Z_scores = (Ytest - m) / sd 

#%%
print(gpr.kernel_.get_params())


