#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 10 16:00:00 2022

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
import os 


#%%

save_data_dir = "/mnt/projects/VIA11/FREESURFER/Stats/Normative_modelling/data/"
blr_out_dir = "/mnt/projects/VIA11/FREESURFER/Stats/Normative_modelling/BLR/"
working_dir = blr_out_dir
os.chdir(working_dir)

trainX_path = save_data_dir + "trainX.txt"
trainY_path = save_data_dir + "trainY.txt"
trainX = pd.read_csv(trainX_path,sep=' ') 
trainY = pd.read_csv(trainY_path,sep=' ') 

testX_path = save_data_dir + "testX.txt"
testY_path = save_data_dir + "testY.txt"
testX = pd.read_csv(testX_path,sep=' ') 
testY = pd.read_csv(testY_path,sep=' ') 


K_testX_path = save_data_dir + "KtestX.txt"
K_testY_path = save_data_dir + "KtestY.txt"
K_testX = pd.read_csv(K_testX_path,sep=' ') 
K_testY = pd.read_csv(K_testY_path,sep=' ') 


SZBP_testX_path = save_data_dir + "SZBP_testX.txt"
SZBP_testY_path = save_data_dir + "SZBP_testY.txt"
SZBP_testX = pd.read_csv(SZBP_testX_path,sep=' ') 
SZBP_testY = pd.read_csv(SZBP_testY_path,sep=' ') 



#%% load data into the right format

trainX_model_path = save_data_dir + "train_covariate_normsample.txt"
trainY_model_path = save_data_dir + "train_y_normsample.txt"
trainX_model = trainX.iloc[:,1:]
trainY_model = trainY[['eICV_samseg']]

testX_model_path = save_data_dir + "test_covariate_normsample.txt"
testY_model_path = save_data_dir + "test_y_normsample.txt"
testX_model = testX.iloc[:,1:]
testY_model = testY[['eICV_samseg']]


#only BP and SZ test data
test2X_model = SZBP_testX.iloc[:,1:]
test2Y_model = SZBP_testY[['eICV_samseg']]


#validation data
valX_model = K_testX.iloc[:,1:]
valY_model = K_testY[['eICV_samseg']]



#forward model covariates
n_vals = 100

euler_mean = trainX["TotalEulerNumber"].mean()

x_age = np.linspace(10,14,n_vals)[:,np.newaxis]
x_diag = np.zeros((n_vals,1))
x_site = np.zeros((n_vals,1))
x_sex = np.zeros((n_vals,1))
x_euler = np.array([[euler_mean]]).repeat(n_vals,1).T

X_forward = np.hstack( [x_age, x_diag, x_site, x_sex, x_euler] )

#normalize X
#X_norm_forward = train_scaler.transform(X_forward)



#%% Estimate BLR instead
import os 
from pcntoolkit.util.utils import create_bspline_basis, compute_MSLL, create_poly_basis, create_design_matrix
import sklearn


# set this path to wherever your ROI_models folder is located (where you copied all of the covariate & response text files to in Step 4)
os.makedirs(blr_out_dir, exist_ok=True)

file_suffix = ["train","test","val","BPSZ","forward"]


X_files = [trainX_model.to_numpy(), testX_model.to_numpy(), valX_model.to_numpy(), test2X_model.to_numpy(), X_forward] #as numpy arrays 
Y_files = [trainY_model.to_numpy(), testY_model.to_numpy(), valY_model.to_numpy(), test2Y_model.to_numpy(), X_forward] #as numpy arrays 


poly_order = 1

PF = sklearn.preprocessing.PolynomialFeatures(degree=poly_order)


for X,Y,suf in zip(X_files,Y_files,file_suffix):

    # load train & test response files
    if suf != "forward":
        resp_file = os.path.join(blr_out_dir, f"resp_{suf}.txt")
        np.savetxt(resp_file, Y)
    
    #Phi = create_poly_basis(X[:,0],dimpoly=2)
    X_poly = PF.fit_transform(X)
    
    cov_file = os.path.join(blr_out_dir, f"cov_poly_{suf}.txt")
    np.savetxt(cov_file, X_poly)
    
    # add intercept column
    #X_spline = np.concatenate((X, np.ones((X.shape[0],1))), axis=1)
    #np.savetxt(os.path.join(data_dir, 'cov_int_{suf}.txt'), X_spline)

    # create Bspline basis set
    #Phi = np.array([B(i) for i in X_spline[:,0]])
    #X_spline = np.concatenate((X_spline, Phi), axis=1)

    #cov_file = os.path.join(data_dir, f"cov_bspline_{suf}.txt")
    #np.savetxt(cov_file, X_spline)


cov_file_tr = os.path.join(blr_out_dir, 'cov_poly_train.txt')
cov_file_te = os.path.join(blr_out_dir, 'cov_poly_test.txt')
cov_file_val = os.path.join(blr_out_dir, 'cov_poly_val.txt')
cov_file_BPSZ = os.path.join(blr_out_dir, 'cov_poly_BPSZ.txt')
cov_file_forward = os.path.join(blr_out_dir, 'cov_poly_forward.txt')


resp_file_tr = os.path.join(blr_out_dir, 'resp_train.txt')
resp_file_te = os.path.join(blr_out_dir, 'resp_test.txt')
resp_file_val = os.path.join(blr_out_dir, 'resp_val.txt')
resp_file_BPSZ = os.path.join(blr_out_dir, 'resp_BPSZ.txt')


#run forward model
pcn.normative.estimate(covfile = cov_file_tr,
                       respfile = resp_file_tr,
                       testcov = cov_file_forward,
                       alg = 'blr',optimizer = 'powell',
                       outputsuffix = '_forwardblr',
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
yhat_te, s2_te, nm, Z_val, metrics_val = pcn.normative.estimate(cov_file_tr,
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

yhat = pd.read_csv(working_dir + 'yhat_forwardblr.txt', sep = ' ', header=None).to_numpy()
ys2 = pd.read_csv(working_dir + 'ys2_forwardblr.txt', sep = ' ', header=None).to_numpy()

sd = np.sqrt(ys2)

#extract data from a particular categorical variable 
col_names = np.array(list(trainX.columns)[1:])
var = [0,0,0]
var_col_idx = [1,2,3]


x_forward_plot = X_forward[:,0]
yhat_plot = yhat
sd_plot = sd

trainXn = trainX_model.to_numpy()
valXn = valX_model.to_numpy()
test2Xn = test2X_model.to_numpy()

train_idx = ((trainXn[:,1]==0.)*1 * (trainXn[:,2]==0.)*1 * (trainXn[:,3]==0.)*1)==1
val_idx = ((valXn[:,1]==0.)*1 * (valXn[:,2]==0.)*1 * (valXn[:,3]==0.)*1)==1
test_idx = ((test2Xn[:,1]==0.)*1 * (test2Xn[:,2]==0.)*1 * (test2Xn[:,3]==0.)*1)==1


dataX = trainXn[train_idx,0]
dataY = trainY_model.to_numpy()[train_idx]

val_dataX = valXn[val_idx,0]
val_dataY = valY_model.to_numpy()[val_idx]

test_data2X = test2Xn[test_idx,0]
test_data2Y = test2Y_model.to_numpy()[test_idx]



fig = plt.figure(figsize=(10,5))
ax = fig.add_subplot(1, 1, 1)

ax.plot(x_forward_plot, yhat_plot)
ax.fill_between(x_forward_plot, (yhat_plot-2*sd_plot).squeeze(), (yhat_plot+2*sd_plot).squeeze(),alpha=0.4,label='95% (2*sd)')

ax.scatter(dataX,dataY,color='red',label='K train')
ax.scatter(val_dataX,val_dataY,color='green',label='K validation',marker='+')
ax.scatter(test_data2X,test_data2Y,color='blue',label='SZ+BP',s=10)


ax.set_title(f"{col_names[var_col_idx]}={var}")
ax.set_xlabel("Age")
ax.set_ylabel("Brain Volume")

ax.legend()


plt.savefig(f"/mnt/projects/VIA11/FREESURFER/Stats/Normative_modelling/figures/BLR_order{poly_order}_model_plot.png")




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

fig.suptitle(f"Polynomial BLR (order {poly_order}) Z-scores of test groups",fontsize=15)


plt.savefig(f"/mnt/projects/VIA11/FREESURFER/Stats/Normative_modelling/figures/BLR_order{poly_order}_Z_scores_plot.png")

