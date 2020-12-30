#!/usr/bin/env python3
# -*- coding: utf-8 -*-
###############################################################################
# Linear elastic element formulation for 3d and 2d plane strain 
# In the folder/modules element the functions ansatz_function and hook_material
# are called. For loops should be avoided and vector based programming should 
# be used instead to reduce the assembly time.
#
# coded by abartels
###############################################################################



import numpy as np
import time
import sys

# self defined modules
from src_FEM.element.ansatz_function import *
from src_FEM.element.hook_material import *

Ke4= lambda Ga,EE,Gb,dV: np.einsum('k,ikjl,l->ij',Ga,EE,Gb)*dV

def element_linear_elastic(e,ndm,ndf,nen,nqp,nst,mat,b_vec,
    Xe,ue,statVar,statVar0):
    
    # initializing vectors and stiffness matrix
    finte=np.zeros(ndf*nen)
    fvole=np.zeros(ndf*nen)
    fdyne=np.zeros(ndf*nen)
    Ke=np.zeros((ndf*nen,ndf*nen))
    
    # tensor dimension
    tdm=3
    
    # some history variables
    totdV=0.0
    hstress=np.zeros(tdm*tdm)
    
    # current displacements for small strain setting
    # xe=Xe+ue for finite strain
    xe=ue[0:ndm,0:nen]
    
    #start_gauss=time.time()
    # gauss points for complete element
    w1,q1=gaussintegrationdata(ndm,nqp)
    
    # shape functions for whole element
    shp,dshp=shape_function(ndm,nen,nqp,q1) # dshp(nen,ndm)

    for q in range(nqp):
        
        # calculate the jacobian
        Jac=np.tensordot(Xe,dshp[range(nqp*q,nen+nqp*q)],1)
        
        # calculate the determinate of jacobian
        detJac=np.linalg.det(Jac)
        
        if detJac<=0.0:
            sys.exit('Error: detJac <= 0 -->> element distortion')
            
        # calculate Volume in every q and sum up
        dV=detJac*w1[q]
        totdV+=dV
        
        # invert Jac
        invJac=np.linalg.inv(Jac)
        
        # gradient of shape functions G=dN/ds * ds/dX
        G=np.tensordot(dshp[range(nqp*q,nqp*q+nen)],invJac,1)

        # calculate gradient of deformation
        h1=xe.dot(G)
        # calculate strain: here linearized theory
        eps=np.zeros((tdm,tdm))
        eps[0:ndm,0:ndm]=0.5*(np.transpose(h1)+h1)
        
        # call linear elastic material behavior: hooks law
        [sig,EE]=hook_material(tdm,mat.Parameter,eps)

        # calculate the stress in element and store it
        if tdm*tdm<nst+1:
            hstress+=sig.reshape([tdm*tdm])*dV
        
        # store history in statVar: stress
        statVar[:tdm*tdm]=hstress[:tdm*tdm]/totdV
        # store the element volume
        if nst>=tdm*tdm+1:
            statVar[tdm*tdm]=totdV
            
        # assembly of element rhs and stiffness Kab
        # loop over entry aa
        for aa in range(nen):
            # fdyne is in this case zero
            fvole[range(aa*ndm,(aa+1)*ndm)]+=shp[aa+nqp*q]*b_vec*dV
            finte[range(aa*ndm,(aa+1)*ndm)]+=\
                -np.tensordot(sig[:ndm,:ndm],G[aa],1)*dV
                    
            # loop over entry bb
            for bb in range(nen):
                # stiffness matrix KAB
                Ke[np.ix_(range(aa*ndm,(aa+1)*ndm),range(bb*ndm,(bb+1)*ndm))]\
                    +=Ke4(G[aa],EE[:ndm,:ndm,:ndm,:ndm],G[bb],dV)
                                    
    return (fvole,fdyne,finte,Ke)
    