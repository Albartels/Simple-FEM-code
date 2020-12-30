#!/usr/bin/env python3
# -*- coding: utf-8 -*-
###############################################################################
# Linear elastic Hook material depending on lame constants mu and lambda using
# the input Youngs-modulus and Poisson ratio. Material response is defined in 
# R3 (tdm=3). For the definition of the 4th
# order stiffness tensor for-loops (Einstein summation) should be avoided, to 
# save compuation time. 
#
# coded by abartels
###############################################################################

import numpy as np


def hook_material(tdm,mat_parameter,eps):
    
    # identity matrix
    II2=np.identity(tdm)
    
    # initialization of sigma and material stiffness
    sig=np.zeros((tdm,tdm))
    EE=np.zeros((tdm,tdm,tdm,tdm))
    
    # material parameter
    mu1=mat_parameter[0]/2.0/(1.0+mat_parameter[1]) #d(1)/2.0d0/(1.0d0+d(2))
    lambda1=mat_parameter[0]*mat_parameter[1] \
        /(1.0+mat_parameter[1])/(1.0-2.0*mat_parameter[1]) #d(1)*d(2)/(1.0d0+d(2))/(1.0d0-2.0d0*d(2)) 
    
    # based on the energy function: 
    # \psi = \mu * \varepsilon:\varepsilon+1/2*\lambda*tr(\varepsilon)^2
    # calculate sigma and EE
    sig=2.0*mu1*eps+lambda1*np.trace(eps)*II2
    EE=mu1*(np.einsum('ik,jl->ijkl',II2,II2)+np.einsum('jk,il->ijkl',II2,II2))\
        +lambda1*np.einsum('ij,kl->ijkl',II2,II2)
    #for ii in range(tdm):
    #    for jj in range(tdm):
    #        for kk in range(tdm):
    #            for ll in range(tdm):
    #                EE[ii,jj,kk,ll]=mu1*(II2[ii,kk]*II2[jj,ll]  \
    #                      + II2[jj,kk]*II2[ii,ll])              \
    #                      + lambda1*II2[ii,jj]*II2[kk,ll] 
    
    
    return (sig,EE)