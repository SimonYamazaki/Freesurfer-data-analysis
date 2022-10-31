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

trainX_model_path = '/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/global_variables/train_covariate_normsample.txt'
trainY_model_path = '/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/global_variables/train_y_normsample.txt'
trainX_model = trainX.iloc[:,1:3]
trainY_model = trainY[['eICV_samseg']]
#trainX_model.to_csv(trainX_model_path,sep = ' ',header = False,index = False)
#trainY_model.to_csv(trainY_model_path,sep = ' ',header = False,index = False)

testX_model_path = '/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/global_variables/test_covariate_normsample.txt'
testY_model_path = '/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/global_variables/test_y_normsample.txt'
testX_model = testX.iloc[:,1:3]
testY_model = testY[['eICV_samseg']]
#testX_model.to_csv(testX_model_path,sep = ' ',header = False,index = False)
#testY_model.to_csv(testY_model_path,sep = ' ',header = False,index = False)


#only BP and SZ test data
test2X_model = SZBP_testX.iloc[:,1:3]
test2Y_model = SZBP_testY[['eICV_samseg']]


#validation data
valX_model = K_testX.iloc[:,1:3]
valY_model = K_testY[['eICV_samseg']]



#forward model covariates
n_vals = 10
x_age = np.linspace(min(trainX["MRI_age_v11"].min(),testX["MRI_age_v11"].min()), max(trainX["MRI_age_v11"].max(),testX["MRI_age_v11"].max()),n_vals)
x_euler = np.linspace(trainX["TotalEulerNumber"].min(), trainX["TotalEulerNumber"].max(),n_vals)
x_sex = np.array([0,1])
x_site = np.array([0,1])
x_axis1 = np.array([0,1])


x_sex = np.vstack( [np.zeros((n_vals,1)) , np.ones((n_vals,1))] )
X_forward = np.hstack( [np.tile(x_age,2)[:,np.newaxis], x_sex] )




#%% Estimate BLR instead
import os 
from pcntoolkit.util.utils import create_bspline_basis, compute_MSLL, create_poly_basis, create_design_matrix
import sklearn


# set this path to wherever your ROI_models folder is located (where you copied all of the covariate & response text files to in Step 4)
data_dir = '/home/simonyj/BLR'
os.makedirs(data_dir, exist_ok=True)

file_suffix = ["train","test","val","BPSZ","forward"]


X_files = [trainX_model.to_numpy(), testX_model.to_numpy(), valX_model.to_numpy(), test2X_model.to_numpy(), X_forward] #as numpy arrays 
Y_files = [trainY_model.to_numpy(), testY_model.to_numpy(), valY_model.to_numpy(), test2Y_model.to_numpy(), X_forward] #as numpy arrays 


# Create a cubic B-spline basis (used for regression)
xmin = 10.5 # xmin & xmax are the boundaries for ages of participants in the dataset
xmax = 13

#B = create_bspline_basis(xmin, xmax)
PF = sklearn.preprocessing.PolynomialFeatures(degree=2)


for X,Y,suf in zip(X_files,Y_files,file_suffix):

    # load train & test response files
    if suf != "forward":
        resp_file = os.path.join(data_dir, f"resp_{suf}.txt")
        np.savetxt(resp_file, Y)
    
    #Phi = create_poly_basis(X[:,0],dimpoly=2)
    X_poly = PF.fit_transform(X)
    
    cov_file = os.path.join(data_dir, f"cov_poly_{suf}.txt")
    np.savetxt(cov_file, X_poly)
    
    # add intercept column
    #X_spline = np.concatenate((X, np.ones((X.shape[0],1))), axis=1)
    #np.savetxt(os.path.join(data_dir, 'cov_int_{suf}.txt'), X_spline)

    # create Bspline basis set
    #Phi = np.array([B(i) for i in X_spline[:,0]])
    #X_spline = np.concatenate((X_spline, Phi), axis=1)

    #cov_file = os.path.join(data_dir, f"cov_bspline_{suf}.txt")
    #np.savetxt(cov_file, X_spline)


cov_file_tr = os.path.join(data_dir, 'cov_poly_train.txt')
cov_file_te = os.path.join(data_dir, 'cov_poly_test.txt')
cov_file_val = os.path.join(data_dir, 'cov_poly_val.txt')
cov_file_forward = os.path.join(data_dir, 'cov_poly_forward.txt')
cov_file_BPSZ = os.path.join(data_dir, 'cov_poly_BPSZ.txt')

resp_file_tr = os.path.join(data_dir, 'resp_train.txt')
resp_file_te = os.path.join(data_dir, 'resp_test.txt')
resp_file_val = os.path.join(data_dir, 'resp_val.txt')
resp_file_BPSZ = os.path.join(data_dir, 'resp_BPSZ.txt')


#run forward model
pcn.normative.estimate(covfile = cov_file_tr,
                       respfile = resp_file_tr,
                       testcov = cov_file_forward,
                       alg = 'blr',
                       optimizer = 'powell',
                       outputsuffix = '_2dforwardblr',
                       saveoutput=True)



# run a model on all test
yhat_te, s2_te, nm, Z, metrics_te = pcn.normative.estimate(cov_file_tr,
                                             resp_file_tr,
                                             testresp=resp_file_te,
                                             testcov=cov_file_te,
                                             alg = 'blr',optimizer = 'powell',
                                             savemodel = False,
                                             saveoutput = False, standardize = False)


#run model on validation data
yhat_te, s2_te, nm, Z_val, metrics_te = pcn.normative.estimate(cov_file_tr,
                                             resp_file_tr,
                                             testresp=resp_file_val,
                                             testcov=cov_file_val,
                                             alg = 'blr',optimizer = 'powell',
                                             savemodel = False,
                                             saveoutput = False, standardize = False)


#run model on BP SZ data
yhat_te, s2_te, nm, Z_BPSZ, metrics_te = pcn.normative.estimate(cov_file_tr,
                                             resp_file_tr,
                                             testresp=resp_file_BPSZ,
                                             testcov=cov_file_BPSZ,
                                             alg = 'blr',optimizer = 'powell',
                                             savemodel = False,
                                             saveoutput = False, standardize = False)


#%% 

yhat = pd.read_csv('/home/simonyj/yhat_2dforwardblr.txt', sep = ' ', header=None).to_numpy()
ys2 = pd.read_csv('/home/simonyj/ys2_2dforwardblr.txt', sep = ' ', header=None).to_numpy()

sd = np.sqrt(ys2)

#extract data from a particular categorical variable 
var = 0
var_col_idx = 1

x_forward_sex0 = X_forward[ X_forward[:,var_col_idx]==var, 0 ] #0 is the age column
yhat_sex0 = yhat[ X_forward[:,var_col_idx]==var ]

x_forward_plot = x_forward_sex0
yhat_plot = yhat_sex0

sd_plot = sd[ X_forward[:,var_col_idx]==var ]


fig = plt.figure(figsize=(10,5))
ax = fig.add_subplot(1, 1, 1)
ax.scatter(trainX['MRI_age_v11'],trainY['eICV_samseg'],color='red',label='K train',s=10)


ax.plot(x_forward_plot, yhat_plot)
ax.fill_between(x_forward_plot, (yhat_plot-2*sd_plot).squeeze(), (yhat_plot+2*sd_plot).squeeze(),alpha=0.4,label='95% (2*sd)')
ax.scatter(trainX['MRI_age_v11'],trainY['eICV_samseg'],color='red',label='K train',s=10)
#ax.scatter(K_testX['MRI_age_v11'],K_testY['eICV_samseg'],color='green',label='K validation',marker='+')
#ax.scatter(SZBP_testX['MRI_age_v11'],SZBP_testY['eICV_samseg'],color='blue',label='SZ+BP',s=10)
ax.legend()


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

fig.suptitle('Z-scores of test groups',fontsize=15)


