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
trainX_path = '/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/global_variables/trainX.txt'
trainY_path = '/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/global_variables/trainY.txt'

trainX.to_csv(trainX_path,
                sep = ' ',
                header = True,
                index = True)

trainY.to_csv(trainY_path,
                sep = ' ',
                header = True,
                index = True)

#save data for modelling 
trainX_path = '/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/global_variables/covariate_normsample.txt'
trainY_path = '/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/global_variables/y_normsample.txt'

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
K_testX_path = '/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/global_variables/KtestX.txt'
K_testY_path = '/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/global_variables/KtestY.txt'
K_testX.to_csv(K_testX_path,sep = ' ',header = True,index = True)
K_testY.to_csv(K_testY_path,sep = ' ',header = True,index = True)

SZBP_testX_path = '/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/global_variables/SZBP_testX.txt'
SZBP_testY_path = '/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/global_variables/SZBP_KtestY.txt'
test_df.to_csv(SZBP_testX_path,sep = ' ',header = True,index = True)
test_Y.to_csv(SZBP_testY_path,sep = ' ',header = True,index = True)

SZBP_testX_path = '/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/global_variables/SZBP_testX_HRS.txt'
SZBP_testY_path = '/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/global_variables/SZBP_testY_HRS.txt'
test_df_HRS.to_csv(SZBP_testX_path,sep = ' ',header = True,index = True)
test_Y.to_csv(SZBP_testY_path,sep = ' ',header = True,index = True)


testX_path = '/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/global_variables/testX.txt'
testY_path = '/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/global_variables/testY.txt'
testX.to_csv(testX_path,sep = ' ',header = True,index = True)
testY.to_csv(testY_path,sep = ' ',header = True, index = True)


#save for this script
testX_path = '/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/global_variables/test_covariate_normsample.txt'
testY_path = '/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/global_variables/test_y_normsample.txt'
testX.to_csv(testX_path,sep = ' ',header = False,index = False)
testY.to_csv(testY_path,sep = ' ',header = False,index = False)


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

m_trainY_path = '/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/global_variables/m_trainY.txt'
m_testY_path = '/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/global_variables/m_testY.txt'

m_trainY.to_csv(m_trainY_path,
                sep = ' ',
                header = True,
                index = True)

m_testY.to_csv(m_testY_path,
                sep = ' ',
                header = True,
                index = True)



m_BPY = pca_df.loc[ df['HighRiskStatus_v11'].isin(["BP"]).to_numpy() ]
m_SZY = pca_df.loc[ df['HighRiskStatus_v11'].isin(["SZ"]).to_numpy() ]
m_valY =  pca_df.loc[ test_idx ]

m_BPY_path = '/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/global_variables/m_BPY.txt'
m_SZY_path = '/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/global_variables/m_SZY.txt'
m_valY_path = '/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/global_variables/m_valY.txt'

m_BPY.to_csv(m_BPY_path,sep = ' ', header = True, index = True)
m_SZY.to_csv(m_SZY_path,sep = ' ', header = True, index = True)
m_valY.to_csv(m_valY_path,sep = ' ', header = True, index = True)



#%%
#test data 

K_testX_path = '/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/global_variables/KtestX.txt'
K_testY_path = '/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/global_variables/KtestY.txt'
K_testX.to_csv(K_testX_path,sep = ' ',header = True,index = True)
K_testY.to_csv(K_testY_path,sep = ' ',header = True,index = True)

SZBP_testX_path = '/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/global_variables/SZBP_testX.txt'
SZBP_testY_path = '/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/global_variables/SZBP_KtestY.txt'
test_df.to_csv(SZBP_testX_path,sep = ' ',header = True,index = True)
test_Y.to_csv(SZBP_testY_path,sep = ' ',header = True,index = True)

SZBP_testX_path = '/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/global_variables/SZBP_testX_HRS.txt'
SZBP_testY_path = '/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/global_variables/SZBP_testY_HRS.txt'
test_df_HRS.to_csv(SZBP_testX_path,sep = ' ',header = True,index = True)
test_Y.to_csv(SZBP_testY_path,sep = ' ',header = True,index = True)


testX_path = '/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/global_variables/testX.txt'
testY_path = '/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/global_variables/testY.txt'
testX.to_csv(testX_path,sep = ' ',header = True,index = True)
testY.to_csv(testY_path,sep = ' ',header = True, index = True)





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


#%% First attempt to run a model on training and test data

yhat_te, s2_te, nm, Z, metrics_te = pcn.normative.estimate(covfile = trainX_path,
                                                           respfile = trainY_path,
                                                           
                                                           testresp = testY_path, 
                                                           testcov = testX_path,
                                                           
                                                           cvfolds = None,
                                                           alg = 'gpr',
                                                           outputsuffix = '_initialTest',
                                                           saveoutput=True)



#%% Make a forward model (predictions?)


#testX.MRI_age_v11.min(), 11.540041, 11.750856, 12.030116, 12.210815, 12.369610, testX.MRI_age_v11.max()

e = -70
covariate_forwardmodel = {'age': [11, 11.3, 11.5, 11.8, 12, 12.5, 12.8,
                                  11, 11.3, 11.5, 11.8, 12, 12.5, 12.8
                                  #testX.MRI_age_v11.min(), 11.540041, 11.750856, 12.030116, 12.210815, 12.369610, testX.MRI_age_v11.max(),
                                  #testX.MRI_age_v11.min(), 11.540041, 11.750856, 12.030116, 12.210815, 12.369610, testX.MRI_age_v11.max()
                                  ],
                           'axis1': [0, 0, 0, 0, 0, 0, 0,
                                     0, 0, 0, 0, 0, 0, 0],
                           'site': [0, 0, 0, 0, 0, 0, 0,
                                    0, 0, 0, 0, 0, 0, 0],
                          'sex': [0, 0, 0, 0, 0, 0, 0,
                                  0, 0, 0, 0, 0, 0, 0],
                          'euler': [e, e, e, e, e, e, e,
                                    e, e, e, e, e, e, e],}

covariate_forwardmodel = pd.DataFrame(data=covariate_forwardmodel)

plot_test_path = '/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/global_variables/covariate_forwardmodel.txt'
covariate_forwardmodel.to_csv(plot_test_path,
                           sep = ' ',
                           header = False,
                           index = False)

# estimate forward model
pcn.normative.estimate(covfile = trainX_path,
                       respfile = trainY_path,
                       testcov = plot_test_path,
                       cvfolds = None,
                       alg = 'gpr',
                       outputsuffix = '_forward')


#yhat, s2 = predict(cov_file_dummy)


#%% Make plots of the forward model 
# try to replicate the tutorial


# confidence interval calculation at x_forward
def confidence_interval(s2,x,z,n_responses_var=1):
  CI=np.zeros((len(x_forward),n_responses_var))
  for i,xdot in enumerate(x_forward):
    #ci_inx=np.isin(x,xdot)
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
    forward_yhat = pd.read_csv('/home/simonyj/yhat_forward.txt', sep = ' ', header=None)
    yhat_forward=forward_yhat.values
    yhat_forward=yhat_forward[7*i:7*(i+1)]
    x_forward=[11.3, 11.5, 11.8, 12, 12.3, 12.5, 12.8]

    # Find the index of the data exclusively for one sex. Female:0, Male: 1
    inx=np.where(trainX.Sex_child_v11==i)[0]
    #x=trainX.values[inx,1]
    x=trainX.values[inx,0]

    # actual data
    y = pd.read_csv(trainY_path, sep = ' ', header=None)
    y=y.values[inx]
    
    # confidence Interval yhat+ z *(std/n^.5)-->.95 % CI:z=1.96, 99% CI:z=2.58
    s2= pd.read_csv('/home/simonyj/ys2_initialTest.txt', sep = ' ', header=None)
    s2=s2.values[inx]

    CI_95=confidence_interval(s2,x,1.96)
    CI_99=confidence_interval(s2,x,2.58)


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


#%% Try to make other variance intervals 


# confidence interval calculation at x_forward
def variance_interval(s2,x,n_responses_var=1):
  CI=np.zeros((len(x_forward),n_responses_var))
  for i,xdot in enumerate(x_forward):
    ci_inx = np.where(np.logical_and(x>=xdot-0.3, x<=xdot+0.3))
    S2=s2[ci_inx]
    S_hat=np.mean(S2,axis=0)
    CI[i,:]=S_hat#np.sqrt(S_hat)
  return CI


feature_names=['eICV_samseg']
sex_covariates=['Female','Male']

# Creating plots for Female and male
for i,sex in enumerate(sex_covariates):
    #forward model data
    forward_yhat = pd.read_csv('/home/simonyj/yhat_forward.txt', sep = ' ', header=None)
    yhat_forward=forward_yhat.values
    yhat_forward=yhat_forward[7*i:7*(i+1)]
    x_forward=[11.3, 11.5, 11.8, 12, 12.3, 12.5, 12.8]

    # Find the index of the data exclusively for one sex. Female:0, Male: 1
    inx=np.where(trainX.Sex_child_v11==i)[0]
    #x=trainX.values[inx,1]
    x=trainX.values[inx,0]

    # actual data
    y = pd.read_csv(trainY_path, sep = ' ', header=None)
    y=y.values[inx]
    
    # confidence Interval yhat+ z *(std/n^.5)-->.95 % CI:z=1.96, 99% CI:z=2.58
    s2= pd.read_csv('/home/simonyj/ys2_initialTest.txt', sep = ' ', header=None)
    s2=s2.values[inx]

    VI=variance_interval(s2,x)


    # Creat a trejactroy for each point
    for j,name in enumerate(feature_names):
         fig=plt.figure()
         ax=fig.add_subplot(111)
         ax.plot(x_forward,yhat_forward[:,j], linewidth=4, label='Normative trejactory')

         ax.plot(x_forward,VI[:,j]+yhat_forward[:,j], linewidth=2,linestyle='--',c='black', label='Variance')
         ax.plot(x_forward,-VI[:,j]+yhat_forward[:,j], linewidth=2,linestyle='--',c='black')

         ax.scatter(x,y[:,j],c='r', label=name)
         plt.legend(loc='upper left')
         plt.title('Normative trejectory of ' +name+' in '+sex+' cohort')
         plt.show()
         plt.close()
         
         
         
         
#%% test on less covariates

trainX_2d_path = '/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/global_variables/train_covariate_normsample.txt'
trainY_2d_path = '/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/global_variables/train_y_normsample.txt'

trainX_2d = trainX[['MRI_age_v11','Sex_child_v11']]
trainY_2d = trainY.copy()

trainX_2d.Sex_child_v11 = trainX_2d.Sex_child_v11.astype(int)

trainX_2d = pd.get_dummies(trainX_2d, columns=['Sex_child_v11'])



trainX_2d.to_csv(trainX_2d_path,
                sep = ' ',
                header = False,
                index = False)

trainY_2d.to_csv(trainY_2d_path,
                sep = ' ',
                header = False,
                index = False)


#compute forward model 
covariate_forwardmodel = {'age': [11, 11.3, 11.5, 11.8, 12, 12.5, 12.8,
                                  11, 11.3, 11.5, 11.8, 12, 12.5, 12.8],
                           'sex0': [0, 0, 0, 0, 0, 0, 0,
                                    1, 1, 1, 1, 1, 1, 1],
                           'sex1': [1, 1, 1, 1, 1, 1, 1,
                                    0, 0, 0, 0, 0, 0, 0]}
covariate_forwardmodel = pd.DataFrame(data=covariate_forwardmodel)

plot_test_path = '/mnt/projects/VIA11/FREESURFER/Stats/Model_tables/global_variables/covariate_forwardmodel.txt'
covariate_forwardmodel.to_csv(plot_test_path,
                           sep = ' ',
                           header = False,
                           index = False)


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

# estimate forward model
pcn.normative.estimate(covfile = trainX_2d_path,
                       respfile = trainY_2d_path,
                       testcov = plot_test_path,
                       cvfolds = None,
                       alg = 'gpr',
                       outputsuffix = '_2dforward')


#yhat_te, s2_te, nm, Z, metrics_te = pcn.normative.estimate(covfile = trainX_2d_path,
#                                                           respfile = trainY_2d_path,
#                                                           
#                                                           testresp = testY_2d_path, 
#                                                           testcov = testX_2d_path,
#                                                           
#                                                           cvfolds = None,
#                                                           alg = 'gpr',
#                                                           outputsuffix = '_2dtest',
#                                                           saveoutput=False)

#%%


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
