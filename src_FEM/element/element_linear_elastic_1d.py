# -*- coding: utf-8 -*-
###############################################################################
# Linear elastic element formulation for 1D 
# In the folder/modules element the functions ansatz_function 
# is called. For loops should be avoided and vector based programming should 
# be used instead to reduce the assembly time.
#
# coded by abartels
###############################################################################


import numpy as np
import time
import sys

# self defined modules
from src_FEM.element.ansatz_function import *

def element_linear_elastic_1d(e,ndm,ndf,nen,nqp,nst,mat,b_vec,
    Xe,ue,statVar,statVar0):
    
    if ndm>1:
        sys.exit('Error: element formulation is restricted to 1d')
    
    # initializing vectors and stiffness matrix
    finte=np.zeros(ndf*nen)
    fvole=np.zeros(ndf*nen)
    fdyne=np.zeros(ndf*nen)
    Ke=np.zeros((ndf*nen,ndf*nen))
    
    
    # some history variables
    totdV=0.0
    hstress=0.0
    
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
        Jac=np.array(Xe).dot(dshp[range(nqp*q,nen+nqp*q)])
        # calculate the determinate of jacobian
        detJac=np.linalg.det(Jac)
        
        if detJac<=0.0:
            sys.exit('Error: detJac <= 0 -->> element distortion')
        
        # calculate Volume in every q and sum up
        dx=detJac*w1[q]
        dV=mat.Parameter[1]*dx
        #print('Area',mat.Parameter[1])
        totdV+=dV
        
        # invert Jac
        invJac=1.0/detJac
        
        # gradient of shape functions G=dN/ds * ds/dX
        G=dshp[range(nqp*q,nqp*q+nen)]*invJac

        # calculate gradient of deformation
        h1=xe.dot(G)
        # calculate strain: here linearized theory
        eps=h1

        
        # call linear elastic material behavior: hooks law
        sig=mat.Parameter[0]*eps
        EE=mat.Parameter[0]
        
        # calculate the stress in element and store it
        if 1<nst+1:
            hstress+=sig*dV
        
        # store history in statVar: stress
        statVar[0]=hstress/totdV
        # store the element volume
        if nst>=2:
            statVar[1]=totdV
 
        # assembly of element rhs and stiffness Kab           
        # loop over entry aa
        for aa in range(nen):
            # fdyne is in this case zero
            fvole[range(aa*ndm,(aa+1)*ndm)]+=shp[aa+nqp*q]*b_vec*dV
            finte[range(aa*ndm,(aa+1)*ndm)]+=-sig[0,0]*G[aa]*dV
                    
            # loop over entry bb
            for bb in range(nen):
                # stiffness matrix KAB
                Ke[np.ix_(range(aa*ndm,(aa+1)*ndm),range(bb*ndm,(bb+1)*ndm))]\
                    +=G[aa]*EE*G[bb]*dV

    return (fvole,fdyne,finte,Ke)
    