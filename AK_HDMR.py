# -*- coding: utf-8 -*-
"""
Created on Thu Apr 20 11:07:32 2017

@author: nlelievre
"""
print('Importation of library')
import time

import numpy as np
import scipy as sp
import openturns as ot
import Functions_AK as Fct_AK
import matplotlib.pyplot as plt
import cma as cma


debut = time.time()
""" Data """

print('Definition of parameters')
Pseuil = 0.999

nMC = 100000
Seuil = 0
cmin = 0.05
nini = 4
Umin = 0.01

#Loi = [1,1]
#Moy = [0,0]
#Stdev = [1,1]
dim = 12
Loi = 2*np.ones(dim)
Moy = 1*np.ones(dim)
Stdev = 0.2*np.ones(dim)


#Loi = [2,2,2,2]
#Moy = [300,130,100,50]
#Stdev = [20,8,5,2.5]
nva = np.size(Loi)
print('Definition of the first DoE')
DOEs = np.random.normal(0,1,(nini,nva))
#DOE_y = Fct_AK.Perf(DOE_u, Seuil, Loi, Moy, Stdev)
#DOE_y = DOE_y.reshape((nini,1))

S = (1 / (DOEs.shape[0] - 1)) * np.dot(DOEs.transpose(),DOEs)
eigenvalues_desor, eigenvectors_desor = np.linalg.eigh(S)
eigenvalues_argOrder = np.argsort(eigenvalues_desor,0)
eigenvalues_argOrder = eigenvalues_argOrder[::-1]
eigenvalues= np.sort(eigenvalues_desor,0)[::-1]
eigenvectors = eigenvectors_desor[:,eigenvalues_argOrder]

def pass_standard_pca(X_stand):
    X_pca = np.dot(eigenvectors.transpose(), X_stand.transpose()).transpose()
    return X_pca

def pass_pca_standard(X_pca):
    X_stand = np.dot( sp.linalg.inv(eigenvectors.transpose()) ,X_pca.transpose()).transpose()
    return X_stand

def _Perf(X_pca, Seuil, Loi, Moy, Stdev):
    X_stand = pass_pca_standard(X_pca)
    y = Fct_AK.Perf(X_stand, Seuil, Loi, Moy, Stdev)
    return y

DOEp = pass_standard_pca(DOEs)

DOE_u = np.zeros((nini,nva,nva))
DOE_y = np.zeros((nini,1,nva))
g0 = _Perf(np.zeros((1,nva)), Seuil, Loi, Moy, Stdev)

for idoe in np.arange(nva):
    DOE_u[:,idoe,idoe] = DOEp[:,idoe]
    DOE_y[:,0,idoe] = _Perf(DOE_u[:,:,idoe], Seuil, Loi, Moy, Stdev)

stop
print('Evaluation of the true failure probability')
Base = np.random.normal(0,1,(nMC,nva))
REP = Fct_AK.Perf(Base, Seuil, Loi, Moy, Stdev)
Pf_vrai = REP[REP<=0].size/nMC

ite =0
crit_OK =0
new_u = []
new_y = []
 	
ienrich = nini*np.ones((nva,))

print('Starting analysis')
while crit_OK == 0:
    ite += 1
    
    print('iteration')
    print(ite)
    
    if ite > 1:
        posU = np.where(np.all(Base==new_u,axis=1))[0][0]
        posv = np.argmax(var_Gnva[posU,:,:])
        ienrich[posv] += 1
        neue_u = new_u[0,posv]
        if ienrich[posv] < DOE_u.shape[0]:
            DOE_u[int(ienrich[posv])-1,posv,posv] = neue_u
            DOE_y[int(ienrich[posv])-1,0,posv] = _Perf( DOE_u[int(ienrich[posv])-1,:,posv].reshape((1,DOE_u[int(ienrich[posv])-1,:,posv].size)) , Seuil, Loi, Moy, Stdev)            
        else:
            doe_uadd = np.zeros((1,DOE_u.shape[1],DOE_u.shape[2]))
            doe_yadd = np.zeros((1,DOE_y.shape[1],DOE_y.shape[2]))
            DOE_u = np.concatenate((DOE_u,doe_uadd))
            DOE_y = np.concatenate((DOE_y,doe_yadd))                 
            DOE_u[int(ienrich[posv])-1,posv,posv] = neue_u
            DOE_y[int(ienrich[posv])-1,0,posv] = _Perf( DOE_u[int(ienrich[posv])-1,:,posv].reshape((1,DOE_u[int(ienrich[posv])-1,:,posv].size)) , Seuil, Loi, Moy, Stdev)
            
            
    
    basis = ot.ConstantBasisFactory(1).build()
    covarianceModel = ot.SquaredExponential(1)
    mu_Gnva = np.zeros((nMC,1,nva))
    var_Gnva = np.zeros((nMC,1,nva))
    for inva in np.arange(nva):
        inputSample = ot.Sample(DOE_u[0:int(ienrich[inva]),inva,inva].reshape((DOE_u[0:int(ienrich[inva]),inva,inva].size,1)))
        outputSample = ot.Sample(DOE_y[0:int(ienrich[inva]),0,inva].reshape((DOE_y[0:int(ienrich[inva]),0,inva].size,1)))
        algo  = ot.KrigingAlgorithm(inputSample, outputSample, covarianceModel, basis)
        algo.run()
        result = algo.getResult()
        metamodel = result.getMetaModel()
        mu_Gnva[:,:,inva] = np.array(metamodel(Base[:,inva].reshape((Base[:,inva].size,1)))) - g0
        for id_point in np.arange(nMC):
            var_Gnva[id_point,0,inva] = np.array(result.getConditionalCovariance(np.array([Base[id_point,inva]])))
        
    """                                 R                                   """
    
    mu_G = g0+np.sum(mu_Gnva,2)
    var_G = np.sum(var_Gnva,2)
    var_G[var_G<np.finfo(float).eps] = np.finfo(float).eps
    U = np.abs(mu_G) / np.sqrt(var_G)
    pos_minU= np.argmin(U)
    minU = U[pos_minU]
    Pf = mu_G[mu_G<0].size/nMC
    if Pf != 0 :
        Cov_Pf = np.sqrt((1-Pf)/(nMC*Pf))
    else:
        Cov_Pf = np.Inf
    
    Pinf=mu_G[np.logical_and(mu_G<=0,U>2)].size/nMC;
    Psup=1-mu_G[np.logical_and(mu_G>=0,U>2)].size/nMC;
    
    if Pinf == 0:
        crit = 1
    else:
        crit = abs((Pinf-Psup)/Pinf)
        
    if crit <= Umin:
        crit_OK = 1
    else:
        new_u = Base[pos_minU,:]
        new_u = new_u.reshape((1,dim))
        new_y = Fct_AK.Perf(new_u, Seuil, Loi, Moy, Stdev)

    print('crit')
    print(crit)
    print('Pf')
    print(Pf)
    print('Pf_vrai')
    print(Pf_vrai)
    print('Cov_Pf')
    print(Cov_Pf)

fin = time.time()-debut

file = open('HDMR_test_python.txt','w')
file.write('Time\n')
file.write(str(fin))
file.write('\n')
file.write('Pf\n')
file.write(str(Pf))
file.write('\n')
file.write('Cov(pf)\n')
file.write(str(Cov_Pf))
file.write('\n')
file.write('Itération\n')
file.write(str(ite))
file.write('\n')
file.write('N calls\n')
file.write(str(nini + ite -1))
file.write('\n')
file.close()




