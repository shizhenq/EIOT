    
"""

    
Needs eiot and eiot_packages both available at:
    https://github.com/salvadorgarciamunoz/eiot
    
"""

import scipy.io as spio
import numpy as np
import eiot
import eiot_extras as ee
import matplotlib.pyplot as plt


#LOAD THE DATA FROM A MATLAB FILE
NIRData      = spio.loadmat('RESULTS_B100_EIOT_vs_CLS_20200312.mat')
nir_spectra  = np.array(NIRData['SPEC_CAL_OLD'])
Ck           = np.array(NIRData['Y_CAL_OLD_1'])
wavenumbers  = np.array(NIRData['WL_AXIS'])
S_I         = np.array(NIRData['RES_mean'])
spec_test   = np.array(NIRData['SPEC_B100_COMBINE'])
HPLC        =np.array(NIRData['HPLC_B100_COMBINE'])
rhat_matlab       =np.array(NIRData['r_hat_Pete'])
#dose_source  = np.array(NIRData['dose_source'])

# PRE-PROCESS SPECTRA 
nir_spectra_use = ee.snv(nir_spectra)
#nir_spectra_2_use,M = ee.savgol(5,1,2,nir_spectra)


# Divide the set into Calibration and Validation taking one in every two samples
#nir_spectra_2_use_cal = nir_spectra_2_use[::2,:]
#nir_spectra_2_use_val = nir_spectra_2_use[1:nir_spectra_2_use.shape[0]:2,:]
#Ck_cal                = Ck[::2,:]
#Ck_val                = Ck[1:Ck.shape[0]:2,:]
#dose_source_cal       = dose_source[::2,:]
#dose_source_val       = dose_source[1:dose_source.shape[0]:2,:]


# Build  Unsupervised EIOT Model and plot lambdas
eiot_obj = eiot.build(nir_spectra_use,Ck)

# Add batch-by-batch non-chemical interference
eiot_obj['S_I']=S_I
S_E=np.vstack((eiot_obj['S_hat'],eiot_obj['S_I']))
eiot_obj['S_E']=S_E
eiot_obj['num_e_sI']=11
eiot_obj['num_sI']=11
eiot_obj['pyo_M']   = np.arange(1,eiot_obj['S_I'].shape[0]+1)
eiot_obj['pyo_M'] = eiot_obj['pyo_M'].tolist()
eiot_obj['pyo_S_I'] = ee.np2D2pyomo(eiot_obj['S_I'])
eiot_obj['pyo_Me']= np.arange(eiot_obj['num_sI']-eiot_obj['num_e_sI']+1,eiot_obj['num_sI']+1)
eiot_obj['pyo_Me']=eiot_obj['pyo_Me'].tolist()
eiot_obj['abs_max_exc_ri']=2
rhat=[]
ri_hat=[]
for i in range(spec_test.shape[0]):
    pred_unsup  = eiot.calc(ee.snv(spec_test[i,:]),eiot_obj)
    rhat.append(pred_unsup['r_hat'])
    ri_hat.append(pred_unsup['r_I_hat'])
ri_hat=np.array(ri_hat)
rhat=np.array(rhat)

A=list(range(178))
B=list(range(11))
fig,ax = plt.subplots()
im=ax.pcolormesh(B,A,ri_hat, cmap='viridis') 
fig.colorbar(im, orientation='vertical')  
ax.set_title('rI_hat')
ax.set_xlabel('batch #')
ax.set_ylabel('sample #')
fig.show()

fig1,ax1=plt.subplots()
ax1.plot(range(178),rhat[:,0],label='pyEIOT w exclusive augmentation')
ax1.plot(range(178),rhat_matlab[:,0],label='mEIOT w exclusive augmentation')
ax1.plot(range(178),HPLC/100*0.4,label='HPLC')
ax1.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
           ncol=2, mode="expand", borderaxespad=0.)
ax1.set_xlabel('sample #')
ax1.set_ylabel('API conc (w/w)')      

#Plot the Lambdas using MATPLOTLIB
#fig,ax=plt.subplots()
#ax.plot(list(range(1,11)),eiot_obj['lambdas'][0:10],'ob')
#ax.set_title('Lambda plot for Unsupervised EIOT')
#ax.set_ylabel('Eigenvalues of $\epsilon_{ch}$')
#plt.show()
#print('Lambdas :' + str(eiot_obj['lambdas'][0:7]))


#Buid EIOT objects with increreasing number of NCI and calculate RMSE
#rmse_vs_nci = []
#for nci in [0,1,2,3,4,5,6,7]:
#    eiot_obj    = eiot.build(nir_spectra_2_use_cal,Ck_cal,num_si_u=nci)
#    pred_unsup  = eiot.calc(nir_spectra_2_use_val,eiot_obj)
#    rmse_unsup  = np.sqrt(np.mean((Ck_val[:,0] - pred_unsup['r_hat'][:,0])**2))
#    rmse_vs_nci.append(rmse_unsup)

#PLOT RMSE vs # of NCI
#fig,ax=plt.subplots()
#ax.plot([0,1,2,3,4,5,6,7],rmse_vs_nci,'ob')
#ax.set_title('RMSE vs # of NCI Unsupervised EIOT')
#ax.set_ylabel('RMSE')
#ax.set_xlabel('# NCI')
#plt.show()

# Build  Supervised EIOT Model and plot lambdas
#eiot_obj_S = eiot.build(nir_spectra_2_use_cal,Ck_cal,R_ik=dose_source_cal)
#pred_sup_ps = eiot.calc(nir_spectra_2_use_val,eiot_obj_S)
#pred_sup_as = eiot.calc(nir_spectra_2_use_val,eiot_obj_S,r_ik=dose_source_val)
#    
#fig,ax=plt.subplots()
#ax.set_title('Lambda plot for Supervised EIOT')
#ax.plot(list(range(1,11)),eiot_obj_S['lambdas'][0:10],'ob')
#ax.set_ylabel('Eigenvalues of $\epsilon_{ch}$')
#plt.show()
#print('Lambdas :' + str(eiot_obj_S['lambdas'][0:7]))

#PLOT RMSE vs # of NCI
#rmse_vs_nci_ps = []
#rmse_vs_nci_as = []
#for nci in [0,1,2,3,4,5,6,7]:
#    eiot_obj_S  = eiot.build(nir_spectra_2_use_cal,Ck_cal,R_ik=dose_source_cal,num_si_u=nci)
#    pred_sup_ps = eiot.calc(nir_spectra_2_use_val,eiot_obj_S)
#    pred_sup_as = eiot.calc(nir_spectra_2_use_val,eiot_obj_S,r_ik=dose_source_val)
#    rmse_ps     = np.sqrt(np.mean((Ck_val[:,0] - pred_sup_ps['r_hat'][:,0])**2))
#    rmse_as     = np.sqrt(np.mean((Ck_val[:,0] - pred_sup_as['r_hat'][:,0])**2))
#    rmse_vs_nci_ps.append(rmse_ps)
#    rmse_vs_nci_as.append(rmse_as)
#
#fig,ax=plt.subplots()
#ax.plot(list(range(0,8)),rmse_vs_nci_ps,'ob',label='Passive Supervision')
#ax.plot(list(range(0,8)),rmse_vs_nci_as,'or',label='Active Supervision')
#ax.set_title('RMSE vs # of NCI Supervised EIOT')
#ax.set_ylabel('RMSE')
#ax.set_xlabel('# NCI')
#ax.legend()
#plt.show()
#

    

# Build  Unsupervised EIOT Model with 1 NCI
#eiot_obj = eiot.build(nir_spectra_2_use_cal,Ck_cal,num_si_u=1)
#print("Lambda threshold for EIOT Unsup = "+ str(eiot_obj['lambdas']))

# Build  Supervised EIOT Model with 1 NCI 
#eiot_obj_S = eiot.build(nir_spectra_2_use_cal,Ck_cal,R_ik=dose_source_cal,num_si_u=1)
#print("Lambda threshold for EIOT Sup = "+ str(eiot_obj_S['lambdas']))

#Predict validation data w/ Unsup EIOT
#print('Making predictions of Validation set using Unsupervised object')
#pred_unsup=eiot.calc(nir_spectra_2_use_val,eiot_obj)

#Plot obs vs Pred.
#fig,ax=plt.subplots()
#ax.set_title('Pred. vs Obs - EIOT Unsupervised')
#ax.plot(Ck_val[:,0],pred_unsup['r_hat'][:,0],'ob')
#ax.set_xlabel('Observed r[API]')
#ax.set_ylabel('Predicted r[API]')
#plt.show()

#Predict validation data w/ Supervised EIOT using Passive Supervision
#print('Making predictions of Validation set using Supervised object w PS')
#pred_sup_ps=eiot.calc(nir_spectra_2_use_val,eiot_obj_S)

#Plot obs vs Pred.
#fig,ax=plt.subplots()
#ax.set_title('Pred. vs Obs - EIOT Supervised w PS')
#ax.plot(Ck_val[:,0],pred_sup_ps['r_hat'][:,0],'ob')
#ax.set_xlabel('Observed r[API]')
#ax.set_ylabel('Predicted r[API]')
#plt.show()
#
##Predict validation data w/ Supervised EIOT and ACTIVE Supervision
#print('Making predictions of Validation set using Supervised object w AS')
#pred_sup_as=eiot.calc(nir_spectra_2_use_val,eiot_obj_S,r_ik=dose_source_val)
#
##Plot obs vs Pred.
#fig,ax=plt.subplots()
#ax.set_title('Pred. vs Obs - EIOT Supervised w AS')
#ax.plot(Ck_val[:,0],pred_sup_as['r_hat'][:,0],'ob')
#ax.set_xlabel('Observed r[API]')
#ax.set_ylabel('Predicted r[API]')
#plt.show()
#
##Calculate RMSE for API prediction
#rmse_unsup  = np.sqrt(np.mean((Ck_val[:,0] - pred_unsup['r_hat'][:,0])**2))
#rmse_sup_ps = np.sqrt(np.mean((Ck_val[:,0] - pred_sup_ps['r_hat'][:,0])**2))
#rmse_sup_as = np.sqrt(np.mean((Ck_val[:,0] - pred_sup_as['r_hat'][:,0])**2))
#
#print('RMSE EIOT Unsupervised   :' + str (rmse_unsup) )
#print('RMSE EIOT Supervised w PS:' + str (rmse_sup_ps) )
#print('RMSE EIOT Supervised w AS:' + str (rmse_sup_as) )
#
#
#
#TSS       = np.sum(Ck_val[:,0]**2)
#RSS_unsup = np.sum((Ck_val[:,0]  - pred_unsup['r_hat'][:,0] )**2)
#RSS_supPS = np.sum((Ck_val[:,0]  - pred_sup_ps['r_hat'][:,0])**2)
#RSS_supAS = np.sum((Ck_val[:,0]  - pred_sup_as['r_hat'][:,0])**2)
#
#
#R2_unsup= 1 - RSS_unsup/TSS
#R2_supPS= 1 - RSS_supPS/TSS
#R2_supAS= 1 - RSS_supAS/TSS
#
#print('R2Y EIOT Unsupervised   :' + str (R2_unsup) )
#print('R2Y EIOT Supervised w PS:' + str (R2_supPS) )
#print('R2Y EIOT Supervised w AS:' + str (R2_supAS) )
#
#
