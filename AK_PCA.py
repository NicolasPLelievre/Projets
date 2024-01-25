# Script permettant le calcul de probabilités de défaillance couplant les méthodes AK et PCA pour l'analyse de fiabilité en grande dimension.

# -*- coding: utf-8 -*-
"""
Created on Thu Apr 20 11:07:32 2017

@author: nlelievre
"""
print('Importation of library')
import time

import numpy as np
import openturns as ot
#import scipy as sp
#from scipy.stats import norm
import Functions_AK as Fct_AK
#import matplotlib.pyplot as plt
import cma as cma
#from mpl_toolkits.mplot3d import Axes3D



""" Data """

print('Definition of parameters')
Pseuil = 0.999

nMC = 10000
Seuil = 0
cmin = 0.05
nini = 20

#Loi = [1,1]
#Moy = [0,0]
#Stdev = [1,1]
dim = 20
Loi = 2*np.ones(dim)
Moy = 1*np.ones(dim)
Stdev = 0.2*np.ones(dim)


#Loi = [2,2,2,2]
#Moy = [300,130,100,50]
#Stdev = [20,8,5,2.5]
nva = np.size(Loi)
print('Definition of the first DoE')
DOE_u = np.random.normal(0,1,(nini,nva))
DOE_y = Fct_AK.Perf(DOE_u, Seuil, Loi, Moy, Stdev)
DOE_y = DOE_y.reshape((nini,1))

print('Evaluation of the true failure probability')
Base = np.random.normal(0,1,(nMC,nva))
REP = Fct_AK.Perf(Base, Seuil, Loi, Moy, Stdev)
Pf_vrai = REP[REP<=0].size/nMC
ite =0
crit_OK =0
new_u = []
new_y = []
 	
time.sleep(1)
print('Starting analysis')
while crit_OK == 0:
    ite += 1
    
    print('iteration')
    print(ite)
    time.sleep(1)
    if ite > 1:
        DOE_u = np.concatenate((DOE_u,new_u))
        DOE_y = np.concatenate((DOE_y,new_y))
    
    
    inputSample = ot.Sample(DOE_u)
    outputSample = ot.Sample(DOE_y)
    basis = ot.ConstantBasisFactory(nva).build()
    covarianceModel = ot.SquaredExponential(nva)
    algo  = ot.KrigingAlgorithm(inputSample, outputSample, covarianceModel, basis)
    
#    algo.run()
#    result1 = algo.getResult()
    LogLikelihood = algo.getReducedLogLikelihoodFunction()
    size = 100
    
    
    distribution = ot.ComposedDistribution([ot.Uniform(1e-10, 100)] * dim)
    inputDesign = ot.SobolIndicesAlgorithmImplementation.Generate(distribution, size, True)
    outputDesign = LogLikelihood(inputDesign)
    stop
    
    out = np.array(outputDesign)
    minout = np.zeros((100,1))
    inp = np.zeros((100,dim))
    kk = 0
    for k in np.arange(100):
        minout[kk] = np.argmin(out)
        out[np.int(minout[kk])] = 0
        inp[kk,:] = np.array(inputDesign[np.int(minout[kk])])
        kk += 1
        
    S = (1 / (inp.shape[0] - 1)) * np.dot(inp.transpose(),inp)
    eigenvalues, eigenvectors = np.linalg.eig(S)
    
    jstop = 0
    jit = 0
    while jstop == 0:
        pctg = np.sum(eigenvalues.real[0:1+jit])/np.sum(eigenvalues.real)
        if pctg < 0.9:
            jit += 1
        else:
            jstop = 1
    if jit == 0:
        jit = 1
    eigenvectors_red = eigenvectors[:,0:jit+1]
    
    
    optim_bounds = [1e-10, 100]
    
    def F(theta_red):
        theta = np.dot( np.linalg.inv(eigenvectors.transpose())[:,0:jit+1] ,theta_red.transpose()).transpose()
        if theta.shape[0] == theta.size:
            th = theta
            if th[th<optim_bounds[0]].size != 0 or th[th>optim_bounds[1]].size != 0 :
                FctLogLikelihood = np.Inf
            else:
                FctLogLikelihood = -float(np.array(LogLikelihood(th)))
        else:       
            FctLogLikelihood = np.zeros((theta.shape[0]))
            for i in np.arange(theta.shape[0]):
                th = theta[i,:]
                if th[th<optim_bounds[0]].size != 0 or th[th>optim_bounds[1]].size != 0 :
                    FctLogLikelihood[i] = np.Inf
                else:
                    FctLogLikelihood[i] = -float(np.array(LogLikelihood(th)))
            
        return FctLogLikelihood
    
    
    theta0 = np.array(inputDesign[np.int(minout[0])]).reshape((1,dim))
    theta0_red = np.dot(eigenvectors_red.transpose(), theta0.transpose()).transpose()
    theta0_red = np.ndarray.tolist(theta0_red)[0]
    sigma0 = 33
    es = cma.CMAEvolutionStrategy(theta0_red,sigma0,{'verb_disp': 0})
    res = es.optimize(F).result()
#    cma.pprint(es.result())
    theta_red_o = es.result()[0]
    print('dimension of the optimization')
    print(theta_red_o.size)
    time.sleep(1)
    theta_o = np.dot( np.linalg.inv(eigenvectors.transpose())[:,0:jit+1] ,theta_red_o.transpose()).transpose()
    
    
    
    covarianceModel = ot.SquaredExponential(nva)
    covarianceModel.setScale(theta_o)
    algo  = ot.KrigingAlgorithm(inputSample, outputSample, covarianceModel, basis)
    algo.setOptimizeParameters(False)
    algo.run()
    result = algo.getResult()
    metamodel = result.getMetaModel()
    stop
    mu_G = np.array(metamodel(Base))
    var_G = np.zeros_like(mu_G)
    for id_point in np.arange(nMC):
        var_G[id_point] = np.array(result.getConditionalCovariance(Base[id_point]))[0]
    var_G[var_G<np.finfo(float).eps] = np.finfo(float).eps
    U = np.abs(mu_G) / np.sqrt(var_G)
    pos_minU= np.argmin(U)
    minU = U[pos_minU]
    Pf = mu_G[mu_G<0].size/nMC
    if Pf != 0 :
        Cov_Pf = np.sqrt((1-Pf)/(nMC*Pf))
    else:
        Cov_Pf = np.Inf
    
    if minU > 2:
        crit_OK = 1
    else:
        new_u = Base[pos_minU,:]
        new_u = new_u.reshape((1,dim))
        new_y = Fct_AK.Perf(new_u, Seuil, Loi, Moy, Stdev)
    print('Pf')
    print(Pf)
    print('Cov_Pf')
    print(Cov_Pf)
    time.sleep(1)        
    




