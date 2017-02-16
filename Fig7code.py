# This software is provided as supplementary material for the following publication:
#
# Eyherabide HG, Disambiguating the role of noise correlations when decoding neural
# populations together, arXiv (2017), https://arxiv.org/abs/1608.05501v2.
#
# Should you use this code, I kindly request you to cite the aforementioened publication.
#
# DESCRIPTION:
#
# Computes the descriptive and communication information losses in Figure 7a, 7b and 7c 
# in the limit for large numbers of independent information streams using all neurons in
# all populations.
#
# DEPENDENCIES:
#
# The software requires the packages
#
# - scipy (https://www.scipy.org/)
# - vegas (https://pypi.python.org/pypi/vegas)
# - json
# - math 
# - numpy
#
# ARGUMENTS AND VARIABLES:
#
# In the codes below, I have consistently employed the following variable names for
# represent the same magnitudes
#
#   - q     denotes the probability that the stimulus is composed of a box.
#   - a     denotes the probability that the response [2,2] is elicited by stimuli
#           composed of boxes.
#   - rho1 and rho2      denote the correlation coefficients of the responses elicited by
#                        stimuli composed of boxes and circles, respectively.
#   - dipmv and dipsd    denote the mean values and the standard deviations, respectively, 
#                        of the descriptive information loss caused by joint NI decoders. 
#   - dilmv and dilsd    denote the mean values and the standard deviations, respectively, 
#                        of the communication information loss caused by joint NI decoders. 
#   - infomv and infosd  denote the mean values and the standard deviations, respectively, 
#                        of the total transmitted information. 
#   - samplesize         denotes the number of samples used to compute the integrals.
#
# EXAMPLE:
#
# import Fig7code as F7c
#
# res7a = F7c.resultsFig7a	
# res7b = F7c.resultsFig7b	
# res7c = F7c.resultsFig7c	
#
# VERSION CONTROL
# 
# V1.000 Hugo Gabriel Eyherabide (10 Feb 2017)
# 
# Should you find bugs, please contact Hugo Gabriel Eyherabide (neuralinfo@eyherabidehg.com)
#
# LICENSE
# 
# Copyright (c) 2017, Hugo Gabriel Eyherabide 
# All rights reserved.
# 
# Redistribution and use in source and binary forms, with or without modification,
# are permitted provided that the following conditions are met:
# 
# 1.  Redistributions of source code must retain the above copyright notice, 
#     this list of conditions and the following disclaimer.
# 
# 2.  Redistributions in binary form must reproduce the above copyright notice, 
#     this list of conditions and the following disclaimer in the documentation 
#     and/or other materials provided with the distribution.
# 
# 3.  Neither the name of the copyright holder nor the names of its contributors 
#     may be used to endorse or promote products derived from this software 
#     without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED 
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. 
# IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, 
# INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT 
# NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR 
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
# WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY 
# OF SUCH DAMAGE.
    

from scipy.optimize import minimize_scalar as minimize
from scipy.integrate import dblquad, quad
import vegas
import math as m
import json
import numpy

# Integrand for computing communication information loss in Figure 7a
def dinidlintFig7a(q,a,theta):
    return (q*a*m.log(1+(1-q)/q*((1-a)/a)**theta))

# Integrand for computing the transmitted information in Figure 7a
def infointFig7a(q):
    return (-q*m.log(q)-(1-q)*m.log(1-q))

# Descriptive information loss in Figure 7a
def dinidFig7a(amax,opt={'xtol':1E-4}):
    return dblquad(lambda y,x:dinidlintFig7a(x,y,1),0.05,0.95,lambda x:0.05,lambda x: amax,epsabs = 1E-6, epsrel=1E-3)
    
# Communication information loss in Figure 7a
def dinidlFig7a(amax,opt={'xtol':1E-4}):
    theta = minimize(lambda theta: dblquad(lambda y,x:dinidlintFig7a(x,y,theta),0.05,0.95,lambda x:0.05,lambda x: amax,epsabs = 1E-6, epsrel=1E-3)[0], method='brent',options=opt)
    return dblquad(lambda y,x:dinidlintFig7a(x,y,theta.x),0.05,0.95,lambda x:0.05,lambda x: amax,epsabs = 1E-6, epsrel=1E-3)
    
    
# Total transmitted information in Figure 7a
def infoFig7a():
    return quad(infointFig7a,.05,.95,epsabs = 1E-6, epsrel=1E-3)



# Compute the descriptive and communication losses for large number of independent
# information streams in Figure 7a
def resultsFig7a():

    # amax denotes the maximum value of the interval from which the probability
    # of the response [2,2] given that the stimulus was a box is chosen for each
    # independent information stream
    amax = [ind*.025+.05 for ind in range(1,100) if ind<0.91/0.025]
    alen = len(amax)
    
    # data will contain the results. 
    datazero = lambda: [0 for ind in range(0,alen)];
    data = {'amax': datazero(),'dipmv':datazero(),'dipsd':datazero(),
            'dilmv':datazero(),'dilsd':datazero(),'infomv':datazero(),'infosd':datazero()}   
                
    # Computes the transmitted information, which is independent of amax
    aux = infoFig7a()
    data['infomv'][0] = aux[0]/0.9
    data['infosd'][0] = aux[1]/0.9

    for ind in range(0,alen):
        amaxnow = amax[ind]
        data['amax'][ind] = amaxnow  
        data['infomv'][ind] = data['infomv'][0]
        data['infosd'][ind] = data['infosd'][0]
        
        probArea = 0.9*(amaxnow-0.05)
        
        # Computes the descriptive information loss
        aux = dinidFig7a(amaxnow)
        data['dipmv'][ind] = aux[0]/probArea
        data['dipsd'][ind] = aux[1]/probArea

        # Computes the communication information loss
        aux = dinidlFig7a(amaxnow)
        data['dilmv'][ind] = aux[0]/probArea
        data['dilsd'][ind] = aux[1]/probArea
        
        print([amax,data['dipmv'][ind],data['dilmv'][ind],data['infomv'][ind]])
    
    return data



# Integrand for computing the descriptive and the communication information loss in Figure 7b
# Recall that the former is equal to the latter with theta=1, and that in Figure 7b, the
# correlation coefficients of the responses associated with boxes and circles are the same
def dinidlintFig7b(data,theta):
    data = numpy.append(data,data[3])
    return dinidlintFig7c(data,theta)

# Integrand for computing the transmitted information in Figure 7b
def infointFig7b(data):
    data = numpy.append(data,data[3])
    return infointFig7c(data)

# Descriptive information loss in Figure 7b
def dinidFig7b(samplesize,amax,rhomax,opt={'xtol':1E-4}):
    # The integration is performed 10 times in order to train the integrator
    # and then 10 times more in order to compute the actual values.
    # Check the documentaion of vegas for more information.
    integ = vegas.Integrator([[-5,5],[-5,5],[0.05,amax],[-.95,rhomax]])
    integ(lambda data: dinidlintFig7b(data,1), nitn=10,neval=samplesize)
    return integ(lambda data: dinidlintFig7b(data,1), nitn=10,neval=samplesize)

    
# Communication information loss in Figure 7b
def dinidlFig7b(samplesize,amax,rhomax,opt={'xtol':1E-4}):
    # The integration is performed 10 times in order to train the integrator
    # and then 10 times more in order to compute the actual values.
    # Check the documentaion of vegas for more information.
    integ = vegas.Integrator([[-5,5],[-5,5],[0.05,amax],[-.95,rhomax]])
    integ(lambda data: dinidlintFig7b(data,1), nitn=10,neval=samplesize)
    theta = minimize(lambda theta: integ(lambda data: dinidlintFig7b(data,theta), nitn=10,neval=samplesize).mean,method='brent',options=opt)
    return integ(lambda data: dinidlintFig7b(data,theta.x), nitn=10,neval=samplesize)
    

# Total transmitted information in Figure 7b
def infoFig7b(samplesize,amax,rhomax,opt={'xtol':1E-4}):
    # The integration is performed 10 times in order to train the integrator
    # and then 10 times more in order to compute the actual values.
    # Check the documentaion of vegas for more information.
    integ = vegas.Integrator([[-5,5],[-5,5],[0.05,amax],[-.95,rhomax]])
    integ(infointFig7b, nitn=10,neval=samplesize)
    return integ(infointFig7b, nitn=10,neval=samplesize)


# Compute the descriptive and communication losses for large number of independent
# information streams in Figure 7b
def resultsFig7b():

    # rhomax denotes the maximum value of the interval from which the correlation coefficients
    # are chosen for each independent information stream
    rhomax = [ind*.025-.95 for ind in range(1,100) if ind<1.91/0.025]
    rholen = len(rhomax)
    
    # data will contain the results. 
    datazero = lambda: [0 for ind in range(0,rholen)];
    data = {'rhomax': datazero(),'dipmv':datazero(),'dipsd':datazero(),
            'dilmv':datazero(),'dilsd':datazero(),'infomv':datazero(),'infosd':datazero()}   


    for ind in range(0,rholen):
        rhomaxnow = rhomax[ind]
        data['rhomax'][ind] = rhomaxnow  
                
        # Computes the descriptive information loss
        aux = dinidFig7b(100000,0.95,rhomaxnow)
        data['dipmv'][ind] = aux.mean
        data['dipsd'][ind] = aux.sdev

        # Computes the communication information loss
        aux = dinidlFig7b(100000,0.95,rhomaxnow)
        data['dilmv'][ind] = aux.mean
        data['dilsd'][ind] = aux.sdev
                   
        # Computes the transmitted information
        aux = infoFig7b(100000,0.95,rhomaxnow)
        data['infomv'][ind] = aux.mean
        data['infosd'][ind] = aux.sdev
        
        print([rhomaxnow,data['dipmv'][ind],data['dilmv'][ind],data['infomv'][ind]])
  

    return data



# Integrand for computing the descriptive and the communication information loss in Figure 7c
# Recall that the former is equal to the latter with theta=1.
def dinidlintFig7c(data,theta):
    x,y,q,rho1,rho2 = data    

    x1 = x+1; y1 = y+1; xpy1 = x1**2+y1**2
    x2 = x-1; y2 = y-1; xpy2 = x2**2+y2**2
    det1 = 1-rho1**2
    det11 = -0.5/det1;
    det12 = rho1/det1;
    det13 = det1**0.5;
    det2 = 1-rho2**2
    det21 = -0.5/det2;
    det22 = rho2/det2;
    det23 = det2**0.5;
    
    k = 0.5/m.pi; k1 = q*k; k2 = (1-q)*k
    
    prs1 = k1/det13*m.exp(det11*xpy1+det12*x1*y1)
    prs2 = k2/det23*m.exp(det21*xpy2+det22*x2*y2)
    pr = prs1 + prs2
    if pr>0: ps1dr = prs1/pr; ps2dr = prs2/pr

    pnis1dr = k1*m.exp(-0.5*theta*xpy1)
    pnis2dr = k2*m.exp(-0.5*theta*xpy2)
    pnisum = pnis1dr + pnis2dr
    if pnisum>0: pnis1dr = pnis1dr/pnisum; pnis2dr = pnis2dr/pnisum
    
    intval = 0
    if prs1>0: intval = intval + prs1*m.log(ps1dr/pnis1dr)
    if prs2>0: intval = intval + prs2*m.log(ps2dr/pnis2dr)
    return intval

# Integrand for computing the transmitted information in Figure 7c
def infointFig7c(data):
    x,y,q,rho1,rho2 = data

    x1 = x+1; y1 = y+1; xpy1 = x1**2+y1**2
    x2 = x-1; y2 = y-1; xpy2 = x2**2+y2**2
    det1 = 1-rho1**2
    det2 = 1-rho2**2
    k = 0.5/m.pi; k1 = q*k; k2 = (1-q)*k
    
    prs1 = k1/det1**0.5*m.exp((-0.5*xpy1+rho1*x1*y1)/det1)
    prs2 = k2/det2**0.5*m.exp((-0.5*xpy2+rho2*x2*y2)/det2)
    pr = prs1 + prs2
    if pr>0: ps1dr = prs1/pr; ps2dr = prs2/pr

    intval = 0
    if prs1>0: intval = intval + prs1*m.log(ps1dr/q)
    if prs2>0: intval = intval + prs2*m.log(ps2dr/(1-q))
    return intval


# Descriptive information loss in Figure 7c
def dinidFig7c(samplesize,amax,rhomax,opt={'xtol':1E-4}):
    # The integration is performed 10 times in order to train the integrator
    # and then 10 times more in order to compute the actual values.
    # Check the documentaion of vegas for more information.
    integ = vegas.Integrator([[-5,5],[-5,5],[0.05,amax],[-.95,rhomax],[-.95,rhomax]])
    integ(lambda data: dinidlintFig7c(data,1), nitn=10,neval=samplesize)
    return integ(lambda data: dinidlintFig7c(data,1), nitn=10,neval=samplesize)

    
# Communication information loss in Figure 7c
def dinidlFig7c(samplesize,amax,rhomax,opt={'xtol':1E-4}):
    # The integration is performed 10 times in order to train the integrator
    # and then 10 times more in order to compute the actual values.
    # Check the documentaion of vegas for more information.
    integ = vegas.Integrator([[-5,5],[-5,5],[0.05,amax],[-.95,rhomax],[-.95,rhomax]])
    integ(lambda data: dinidlintFig7c(data,1), nitn=10,neval=samplesize)
    theta = minimize(lambda theta: integ(lambda data: dinidlintFig7c(data,theta), nitn=10,neval=samplesize).mean,method='brent',options=opt)
    return integ(lambda data: dinidlintFig7c(data,theta.x), nitn=10,neval=samplesize)
    

# Total transmitted information in Figure 7c
def infoFig7c(samplesize,amax,rhomax,opt={'xtol':1E-4}):
    # The integration is performed 10 times in order to train the integrator
    # and then 10 times more in order to compute the actual values.
    # Check the documentaion of vegas for more information.
    integ = vegas.Integrator([[-5,5],[-5,5],[0.05,amax],[-.95,rhomax],[-.95,rhomax]])
    integ(infointFig7c, nitn=10,neval=samplesize)
    return integ(infointFig7c, nitn=10,neval=samplesize)


# Compute the descriptive and communication losses for large number of independent
# information streams in Figure 7c
def resultsFig7c():

    # rhomax denotes the maximum value of the interval from which the correlation coefficients
    # are chosen for each independent information stream
    rhomax = [ind*.025-.95 for ind in range(1,100) if ind<1.91/0.025]
    rholen = len(rhomax)
    
    # data will contain the results. 
    datazero = lambda: [0 for ind in range(0,rholen)];
    data = {'rhomax': datazero(),'dipmv':datazero(),'dipsd':datazero(),
            'dilmv':datazero(),'dilsd':datazero(),'infomv':datazero(),'infosd':datazero()}   


    for ind in range(0,rholen):
        rhomaxnow = rhomax[ind]
        data['rhomax'][ind] = rhomaxnow  
                
        # Computes the descriptive information loss
        aux = dinidFig7c(100000,0.95,rhomaxnow)
        data['dipmv'][ind] = aux.mean
        data['dipsd'][ind] = aux.sdev

        # Computes the communication information loss
        aux = dinidlFig7c(100000,0.95,rhomaxnow)
        data['dilmv'][ind] = aux.mean
        data['dilsd'][ind] = aux.sdev
                   
        # Computes the transmitted information
        aux = infoFig7c(100000,0.95,rhomaxnow)
        data['infomv'][ind] = aux.mean
        data['infosd'][ind] = aux.sdev
        
        print([rhomaxnow,data['dipmv'][ind],data['dilmv'][ind],data['infomv'][ind]])
  
    return data

# Reads results in json format
def readJson(filename):
    with open(filename,'r') as infile:
        data = json.load(infile)
    return data
    
# Saves results in json format
def saveJson(filename,data):
    with open(filename,'w') as outfile:
        json.dump(data,outfile)

    
    
