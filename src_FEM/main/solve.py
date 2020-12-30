#!/usr/bin/env python3
# -*- coding: utf-8 -*-
###############################################################################
# Solving the linearized system of equations deleting the dirichlet related 
# dofs. Solving the Dofs of free nodes.
# Using here the linalg.solve function of numpy
#
# coded by abartels
###############################################################################

import numpy as np


def solve(ndf,nnp,iter1,numfreeDofs,freeDofs,tol,rsd,K,u):
    
    
    # 2-norm of residuum
    rsn=np.linalg.norm(rsd)
    
    if rsn>=tol:
        
        # define solution vector of free DOFs
        #duF=np.zeros(numfreeDofs)
        
        # define corresponding rhs
        rsdF=rsd[freeDofs]
        
        # reducing stiffness matrix to free Dofs
        KFF=K[np.ix_(freeDofs,freeDofs)]
                
        # solve system of equations
        duF=np.linalg.solve(KFF,rsdF)
               
        # update displacemnet vector
        u[freeDofs,0]+=duF  
    
    return (rsn,u)