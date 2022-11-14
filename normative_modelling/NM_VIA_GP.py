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

import os 

#%%
save_data_dir = "/mnt/projects/VIA11/FREESURFER/Stats/Normative_modelling/data/"
gpr_out_dir = "/mnt/projects/VIA11/FREESURFER/Stats/Normative_modelling/GPR/"
working_dir = gpr_out_dir
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


SZBP_testX_path = save_data_dir + "SZBP_testX_HRS.txt"
SZBP_testY_path = save_data_dir + "SZBP_testY_HRS.txt"
SZBP_testX_HRS = pd.read_csv(SZBP_testX_path,sep=' ') 
SZBP_testY_HRS = pd.read_csv(SZBP_testY_path,sep=' ') 

SZ_testX = SZBP_testX_HRS.loc[ SZBP_testX_HRS["HighRiskStatus_v11"].isin(["SZ"]) ]
SZ_testX = SZ_testX.drop(["HighRiskStatus_v11"],axis=1)
SZ_testY = SZBP_testY_HRS.loc[ SZBP_testX_HRS["HighRiskStatus_v11"].isin(["SZ"]) ]

BP_testX = SZBP_testX_HRS.loc[ SZBP_testX_HRS["HighRiskStatus_v11"].isin(["BP"]) ]
BP_testX = BP_testX.drop(["HighRiskStatus_v11"],axis=1)
BP_testY = SZBP_testY_HRS.loc[ SZBP_testX_HRS["HighRiskStatus_v11"].isin(["BP"]) ]


#%% load data into the right format

#training data
trainX_model = trainX.iloc[:,1:]
trainY_model = trainY[['eICV_samseg']]

#testing data
testX_model = testX.iloc[:,1:]
testY_model = testY[['eICV_samseg']]

#only BP and SZ test data
test2X_model = SZBP_testX.iloc[:,1:]
test2Y_model = SZBP_testY[['eICV_samseg']]

#BP only
testBPX_model = BP_testX.iloc[:,1:]
testBPY_model = BP_testY[['eICV_samseg']]

#SZ only
testSZX_model = SZ_testX.iloc[:,1:]
testSZY_model = SZ_testY[['eICV_samseg']]

#validation data
valX_model = K_testX.iloc[:,1:]
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

#Validation data 
#normalize X
valX_norm = test_scaler.transform(valX_model.to_numpy())
#normalize y
valY_norm = test_scalery.transform(valY_model.to_numpy())

#Testing data
#normalize X
test2X_norm = test_scaler.transform(test2X_model.to_numpy())
#normalize y
test2Y_norm = test_scalery.transform(test2Y_model.to_numpy())

#normalize X
testBPX_norm = test_scaler.transform(testBPX_model.to_numpy())
#normalize y
testBPY_norm = test_scalery.transform(testBPY_model.to_numpy())

#normalize X
testSZX_norm = test_scaler.transform(testSZX_model.to_numpy())
#normalize y
testSZY_norm = test_scalery.transform(testSZY_model.to_numpy())




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
X_norm_forward = train_scaler.transform(X_forward)



#%% Estimate BLR instead
import os 
from pcntoolkit.util.utils import create_bspline_basis, compute_MSLL, create_poly_basis, create_design_matrix
import sklearn


# set this path to wherever your ROI_models folder is located (where you copied all of the covariate & response text files to in Step 4)
os.makedirs(gpr_out_dir, exist_ok=True)

file_suffix = ["train","test","val","BPSZ","BP","SZ","forward"]


X_files = [X_norm, testX_norm, valX_norm, test2X_norm, testBPX_norm, testSZX_norm, X_norm_forward]#[trainX_model.to_numpy(), testX_model.to_numpy(), valX_model.to_numpy(), test2X_model.to_numpy()]#, X_forward] #as numpy arrays 
Y_files = [Y_norm, testY_norm, valY_norm, test2Y_norm, testBPY_norm, testSZY_norm, X_norm_forward]#[trainY_model.to_numpy(), testY_model.to_numpy(), valY_model.to_numpy(), test2Y_model.to_numpy()]#, X_forward] #as numpy arrays 



for X,Y,suf in zip(X_files,Y_files,file_suffix):

    if suf != "forward":
        resp_file = os.path.join(gpr_out_dir, f"resp_{suf}.txt")
        np.savetxt(resp_file, Y)
    
    cov_file = os.path.join(gpr_out_dir, f"cov_gpr_{suf}.txt")
    np.savetxt(cov_file, X)
    

# save the file paths as string variables
cov_file_tr = os.path.join(gpr_out_dir, 'cov_gpr_train.txt')
cov_file_te = os.path.join(gpr_out_dir, 'cov_gpr_test.txt')
cov_file_val = os.path.join(gpr_out_dir, 'cov_gpr_val.txt')
cov_file_forward = os.path.join(gpr_out_dir, 'cov_gpr_forward.txt')
cov_file_BPSZ = os.path.join(gpr_out_dir, 'cov_gpr_BPSZ.txt')
cov_file_BP = os.path.join(gpr_out_dir, 'cov_gpr_BP.txt')
cov_file_SZ = os.path.join(gpr_out_dir, 'cov_gpr_SZ.txt')


resp_file_tr = os.path.join(gpr_out_dir, 'resp_train.txt')
resp_file_te = os.path.join(gpr_out_dir, 'resp_test.txt')
resp_file_val = os.path.join(gpr_out_dir, 'resp_val.txt')
resp_file_BPSZ = os.path.join(gpr_out_dir, 'resp_BPSZ.txt')
resp_file_BP = os.path.join(gpr_out_dir, 'resp_BP.txt')
resp_file_SZ = os.path.join(gpr_out_dir, 'resp_SZ.txt')



# All these functions saves in the current working dir

# run a model on all test
yhat_val, s2_val, nm, Z_val, metrics_val = pcn.normative.estimate(covfile = cov_file_tr,
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



#Only BP
yhat_BP, s2_BP, nm, Z_BP, metrics_te = pcn.normative.estimate(covfile = cov_file_tr,
                                                           respfile = resp_file_tr,
                                                           
                                                           testcov = cov_file_BP,
                                                           testresp = resp_file_BP, 
                                                           
                                                           optimizer = 'cg',
                                                           cvfolds = None,
                                                           alg = 'gpr',
                                                           #outputsuffix = '_gpr2dtest',
                                                           saveoutput=False)

#only SZ
yhat_SZ, s2_SZ, nm, Z_SZ, metrics_te = pcn.normative.estimate(covfile = cov_file_tr,
                                                           respfile = resp_file_tr,
                                                           
                                                           testcov = cov_file_SZ,
                                                           testresp = resp_file_SZ, 
                                                           
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

yhat = pd.read_csv(working_dir + 'yhat_gpr2dtest.txt', sep = ' ', header=None).to_numpy()
ys2 = pd.read_csv(working_dir + 'ys2_gpr2dtest.txt', sep = ' ', header=None).to_numpy()

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
norm_dataY = Y_norm[train_idx]

val_dataX = valXn[val_idx,0]
test_norm_dataY = valY_norm[val_idx]

test_data2X = test2Xn[test_idx,0]
test_norm_data2Y = test2Y_norm[test_idx]



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


plt.savefig('/mnt/projects/VIA11/FREESURFER/Stats/Normative_modelling/figures/GP_model_plot.png')

#%%

plt.hist(Z,bins=10)
plt.xlim(Z.min(),Z.max())
plt.title("Univariate target group")
plt.axvline(x=0,color='black',linestyle='--')


#%% Histogram plots 
fig = plt.figure(figsize=(20,4))

ax = fig.add_subplot(1, 4, 1)
ax.hist(Z_val,bins=10)
ax.set_xlim(-3.5,3.5)
ax.set_title("Validation group (only control)")
ax.axvline(x=0,color='black',linestyle='--')

ax = fig.add_subplot(1, 4, 2)
ax.hist(Z_BPSZ,bins=10)
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

#manual computation of Z 
#Z[ts, nz[i]] = (Ytest - Yhat[ts, nz[i]]) / np.sqrt(S2[ts, nz[i]]) 
#dobbelttjek evt med den BLR model du har

plt.savefig('/mnt/projects/VIA11/FREESURFER/Stats/Normative_modelling/figures/GP_Z_scores_plot.png')




#%% Custom GP

from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import WhiteKernel, RBF, ConstantKernel, DotProduct

X = trainX_model.to_numpy()#[:,0][:,np.newaxis]
y = trainY_model.to_numpy() #/ trainY_model.to_numpy().std()#+ resp_mean

kernel = DotProduct() + ConstantKernel(constant_value=1e1, constant_value_bounds=(1e-3, 1e2)) * RBF(length_scale=1, length_scale_bounds=(1e-1, 1e5)) + WhiteKernel(noise_level=1, noise_level_bounds=(1e-3, 1e5))

gpr = GaussianProcessRegressor(kernel=kernel, normalize_y=True).fit(X, y) #, alpha=3e5

print(gpr.kernel_.get_params())

#m, sd = gpr.predict(X_forward[:100][:,np.newaxis], return_std=True)
#m, sd = gpr.predict(X_forward, return_std=True)

#plt.plot(X_forward[:100,0],m)
#plt.fill_between(X_forward[:100,0], (m.squeeze()-2*sd), (m.squeeze()+2*sd),alpha=0.5,label='95% (2*sd)')
#plt.scatter(trainX['MRI_age_v11'],trainY['eICV_samseg'],color='red',label='K train',s=10)

#plt.scatter(trainX['MRI_age_v11'], (trainY['eICV_samseg']-resp_mean)/trainY_model.to_numpy().std(),color='red',label='K train',s=10)
#plt.scatter(testX['MRI_age_v11'],testY['eICV_samseg'],color='green',label='test',s=10)
#plt.ylim([1e6, 2e6])

valXn = valX_model.to_numpy()
valYn = valY_model.to_numpy().squeeze()
m, sd = gpr.predict(valXn, return_std=True)
Z_scores_val = (valYn - m.squeeze()) / sd 

test2Xn = test2X_model.to_numpy()
test2Yn = test2Y_model.to_numpy().squeeze()
m, sd = gpr.predict(test2Xn, return_std=True)
Z_scores_BPSZ = (test2Yn - m.squeeze()) / sd 

Xtest = testX_model.to_numpy()
Ytest = testY_model.to_numpy().squeeze()
m, sd = gpr.predict(Xtest, return_std=True)
Z_scores = (Ytest - m.squeeze()) / sd 


#%%

fig = plt.figure(figsize=(10,4))

ax = fig.add_subplot(1, 2, 1)
ax.hist(Z_scores_BPSZ,bins=10)
ax.set_xlim(Z_scores.min(),Z_scores.max())
ax.set_title("Target group")

ax = fig.add_subplot(1, 2, 2)
ax.hist(Z_scores_val,bins=10)
ax.set_xlim(Z_scores.min(),Z_scores.max())
ax.set_title("Validation group (only control)")

fig.suptitle('GP based Z-scores of test groups',fontsize=15)

#manual computation of Z 
#Z[ts, nz[i]] = (Ytest - Yhat[ts, nz[i]]) / np.sqrt(S2[ts, nz[i]]) 
#dobbelttjek evt med den BLR model du har

plt.savefig('/mnt/projects/VIA11/FREESURFER/Stats/Normative_modelling/figures/custum_GP_Z_scores_plot.png')







