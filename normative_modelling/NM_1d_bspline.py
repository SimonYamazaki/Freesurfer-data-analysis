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



#%% test on less covariates

trainX_2d_path = '/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/global_variables/train_covariate_normsample.txt'
trainY_2d_path = '/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/global_variables/train_y_normsample.txt'

#trainX_2d = trainX[['MRI_age_v11','Sex_child_v11']]
#trainX_2d.Sex_child_v11 = trainX_2d.Sex_child_v11.astype(int)
#trainX_2d = pd.get_dummies(trainX_2d, columns=['Sex_child_v11'])
#trainY_2d = trainY.copy()

trainX_2d_loaded = trainX[['MRI_age_v11']]
trainY_2d = trainY[['eICV_samseg']].copy()
trainX_2d = trainX_2d_loaded.copy()
#trainX_2d = (trainX_2d - trainX_2d.mean())/trainX_2d.std()
trainX_2d.to_csv(trainX_2d_path,sep = ' ',header = False,index = False)
trainY_2d.to_csv(trainY_2d_path,sep = ' ',header = False,index = False)


testX_2d_loaded = testX[['MRI_age_v11']]
testY_2d = testY[['eICV_samseg']].copy()
testX_2d = testX_2d_loaded.copy()
#testX_2d.to_csv(testX_2d_path,sep = ' ',header = False,index = False)
#testY_2d.to_csv(testY_2d_path,sep = ' ',header = False,index = False)


#compute forward model 
covariate_forwardmodel = {'age': [11, 11.3, 11.5, 11.8, 12, 12.5, 12.8,
                                  11, 11.3, 11.5, 11.8, 12, 12.5, 12.8],
                           'sex0': [0, 0, 0, 0, 0, 0, 0,
                                    1, 1, 1, 1, 1, 1, 1],
                           'sex1': [1, 1, 1, 1, 1, 1, 1,
                                    0, 0, 0, 0, 0, 0, 0]}

covariate_forwardmodel = {'age': [11, 11.3, 11.5, 11.8, 12, 12.5, 12.8]}
covariate_forwardmodel = pd.DataFrame(data=covariate_forwardmodel)


#xf = np.array([11.3, 11.5, 11.8, 12, 12.3, 12.5, 12.8])
#xf = (xf-trainX_2d_loaded.to_numpy().mean()) / trainX_2d_loaded.to_numpy().std()
#covariate_forwardmodel = pd.DataFrame(xf)


plot_test_path = '/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/global_variables/covariate_forwardmodel.txt'
covariate_forwardmodel.to_csv(plot_test_path,
                           sep = ' ',
                           header = False,
                           index = False)

#%% modelling

pcn.normative.estimate(covfile = trainX_2d_path,
                       respfile = trainY_2d_path,
                       testcov = plot_test_path,
                       cvfolds = None,
                       alg = 'gpr',
                       outputsuffix = '_2dforward')

pcn.normative.estimate(covfile = trainX_2d_path,
                        respfile = trainY_2d_path,             
                        cvfolds = 2,
                        alg = 'gpr',
                        outputsuffix = '_2fold')


#%% Estimate BLR instead

import os 
from pcntoolkit.util.utils import create_bspline_basis, compute_MSLL

# set this path to wherever your ROI_models folder is located (where you copied all of the covariate & response text files to in Step 4)
data_dir = '/home/simonyj/BLR'
os.makedirs(data_dir, exist_ok=True)

file_suffix = ["train","test","val","forward"]

valX_2d = K_testX["MRI_age_v11"].to_numpy()[:,np.newaxis]
valY_2d = K_testY["eICV_samseg"].to_numpy()

x_forward = [11, 11.3, 11.5, 11.8, 12, 12.3, 12.5, 12.8, 13]
X_forward = np.array(x_forward)[:,np.newaxis]

X_files = [trainX_2d.to_numpy(), testX_2d.to_numpy(), valX_2d, X_forward] #as numpy arrays 
Y_files = [trainY_2d.to_numpy(), testY_2d.to_numpy(), valY_2d, X_forward] #as numpy arrays 


# Create a cubic B-spline basis (used for regression)
xmin = 10.5 # xmin & xmax are the boundaries for ages of participants in the dataset
xmax = 13

B = create_bspline_basis(xmin, xmax)


for X,Y,suf in zip(X_files,Y_files,file_suffix):

    # load train & test response files
    if suf != "forward":
        resp_file = os.path.join(data_dir, f"resp_{suf}.txt")
        np.savetxt(resp_file, Y)
    
    # add intercept column
    X_spline = np.concatenate((X, np.ones((X.shape[0],1))), axis=1)
    #np.savetxt(os.path.join(data_dir, 'cov_int_{suf}.txt'), X_spline)

    # create Bspline basis set
    Phi = np.array([B(i) for i in X_spline[:,0]])
    X_spline = np.concatenate((X_spline, Phi), axis=1)

    cov_file = os.path.join(data_dir, f"cov_bspline_{suf}.txt")
    np.savetxt(cov_file, X_spline)


cov_file_tr = os.path.join(data_dir, 'cov_bspline_train.txt')
cov_file_te = os.path.join(data_dir, 'cov_bspline_test.txt')
cov_file_val = os.path.join(data_dir, 'cov_bspline_val.txt')
cov_file_forward = os.path.join(data_dir, 'cov_bspline_forward.txt')

resp_file_tr = os.path.join(data_dir, 'resp_train.txt')
resp_file_te = os.path.join(data_dir, 'resp_test.txt')
resp_file_val = os.path.join(data_dir, 'resp_val.txt')


#run forward model
pcn.normative.estimate(covfile = cov_file_tr,
                       respfile = resp_file_tr,
                       testcov = cov_file_forward,
                       alg = 'blr',
                       optimizer = 'powell',
                       outputsuffix = '_2dforwardblr')


# run a model on all test
yhat_te, s2_te, nm, Z, metrics_te = pcn.normative.estimate(cov_file_tr,
                                             resp_file_tr,
                                             testresp=resp_file_te,
                                             testcov=cov_file_te,
                                             alg = 'blr',
                                             optimizer = 'powell',
                                             savemodel = False,
                                             saveoutput = False,
                                             standardize = False)


#run model on validation data
yhat_te, s2_te, nm, Z_val, metrics_te = pcn.normative.estimate(cov_file_tr,
                                             resp_file_tr,
                                             testresp=resp_file_val,
                                             testcov=cov_file_val,
                                             alg = 'blr',
                                             optimizer = 'powell',
                                             savemodel = False,
                                             saveoutput = False,
                                             standardize = False)

#%% 

yhat = pd.read_csv('/home/simonyj/yhat_2dforwardblr.txt', sep = ' ', header=None).to_numpy()
ys2 = pd.read_csv('/home/simonyj/ys2_2dforwardblr.txt', sep = ' ', header=None).to_numpy()

sd = np.sqrt(ys2)
#n=3
#z=1.96
#sd = z*np.power(ys2/n,.5)

fig = plt.figure(figsize=(10,5))
ax = fig.add_subplot(1, 1, 1)


ax.plot(x_forward, yhat)
ax.fill_between(x_forward, (yhat-2*sd).squeeze(), (yhat+2*sd).squeeze(),alpha=0.4,label='95% (2*sd)')
ax.scatter(trainX['MRI_age_v11'],trainY['eICV_samseg'],color='red',label='K train',s=10)
ax.scatter(K_testX['MRI_age_v11'],K_testY['eICV_samseg'],color='green',label='K validation',marker='+')
ax.scatter(SZBP_testX['MRI_age_v11'],SZBP_testY['eICV_samseg'],color='blue',label='SZ+BP',s=10)
ax.legend()


#%%
fig = plt.figure(figsize=(10,4))

ax = fig.add_subplot(1, 2, 1)
ax.hist(Z,bins=10)
ax.set_xlim(Z.min(),Z.max())
ax.set_title("Entire test group")

ax = fig.add_subplot(1, 2, 2)
ax.hist(Z_val,bins=10)
ax.set_xlim(Z.min(),Z.max())
ax.set_title("Validation group (only control)")

fig.suptitle('Z-scores of test groups',fontsize=15)




#%% plot from tutorial 


def confidence_interval(s2,x,z,n_responses_var=1):
  CI=np.zeros((len(x_forward),n_responses_var))
  for i,xdot in enumerate(x_forward):
    ci_inx = np.where(np.logical_and(x>=xdot-0.3, x<=xdot+0.3))
    S2=s2[ci_inx]
    S_hat=np.mean(S2,axis=0)
    n=S2.shape[0]
    CI[i,:]=z*np.power(S_hat/n,.5)
  return CI



feature_names=['eICV_samseg']
sex_covariates=['Female','Male']

# Creating plots for Female and male
for i,sex in enumerate(sex_covariates):
    #forward model data
    forward_yhat = pd.read_csv('/home/simonyj/yhat_2dforward.txt', sep = ' ', header=None)
    yhat_forward = forward_yhat.values
    yhat_forward = yhat_forward[7*i:7*(i+1)] 
    x_forward = [11.3, 11.5, 11.8, 12, 12.3, 12.5, 12.8]

    # Find the index of the data exclusively for one sex. Female:0, Male: 1
    inx = np.where(trainX_2d[f"Sex_child_v11_{i}"]==i)[0]
    #x=trainX.values[inx,1]
    x = trainX_2d.values[inx,0]

    # actual data
    y = pd.read_csv(trainY_2d_path, sep = ' ', header=None)
    y = y.values[inx]
    
    # confidence Interval yhat+ z *(std/n^.5) 
    #s2 = pd.read_csv('/home/simonyj/ys2_2fold.txt', sep = ' ', header=None) #brug s2 fra estimeret træningsmodel data, til at finde s2 værdier omkring x_forward værdier
    #s2 = s2.values[inx]

    s2 = pd.read_csv('/home/simonyj/ys2_2dforward.txt', sep = ' ', header=None)
    s2 = s2.values[7*i:7*(i+1)]
    
    CI_95 = confidence_interval(s2,np.array(x_forward),1.96)
    CI_99 = confidence_interval(s2,np.array(x_forward),2.58)


    # Creat a trejactroy for each point
    for j,name in enumerate(feature_names):
         fig=plt.figure()
         ax=fig.add_subplot(111)
         ax.plot(x_forward,yhat_forward[:,j], linewidth=4, label='Normative trejactory')


         ax.plot(x_forward,CI_95[:,j]+yhat_forward[:,j], linewidth=2,linestyle='--',c='g', label='95% confidence interval')
         ax.plot(x_forward,-CI_95[:,j]+yhat_forward[:,j], linewidth=2,linestyle='--',c='g')

         ax.plot(x_forward,CI_99[:,j]+yhat_forward[:,j], linewidth=1,linestyle='--',c='k', label='99% confidence interval')
         ax.plot(x_forward,-CI_99[:,j]+yhat_forward[:,j], linewidth=1,linestyle='--',c='k')

         ax.scatter(x,y[:,j],c='r', label=name)
         plt.legend(loc='upper left')
         plt.title('Normative trejectory of ' +name+' in '+sex+' cohort')
         plt.show()
         plt.close()
