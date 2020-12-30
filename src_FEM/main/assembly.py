#!/usr/bin/env python3
# -*- coding: utf-8 -*-
###############################################################################
# Assembly of all element rhs (finte,fvole,fdyne) and stiffnesses (Ke).
# New element routines should be added here using a new Element_type
#
# coded by abartels
###############################################################################

import numpy as np

# self defined elements
from src_FEM.element.element_linear_elastic import  *
from src_FEM.element.element_linear_elastic_1d import  *

def assembly(ndm,ndf,nen,nst,nnp,nel,nqp,elements,b_vec,u,
                 fvol,fdyn,fint,K):
    
    # initialization of force vectors and matrices on element level
    finte=np.zeros(ndf*nen)
    fvole=np.zeros(ndf*nen)
    fdyne=np.zeros(ndf*nen)
    Ke=np.zeros((ndf*nen,ndf*nen))
    
    # loop over all elements
    for e in range(nel):
        
        # gdof list: list of global degree's of freedom for element e
        # already stored in elements[e].edof -->> now in vector
        gdof=np.reshape(elements[e].edof,(ndf*nen),order='F')
        
        # element solution vector in ue [ndf]x[nen]
        ue=u[np.array(elements[e].edof).reshape(ndf,nen),0].reshape([ndf,nen])
        
        # calculate element stiffness depending on the material law,
        # element type
        if elements[e].mat.Element_type==1:
            
            [fvole,fdyne,finte,Ke]=element_linear_elastic(
                    e,ndm,ndf,nen,nqp,nst,
                    elements[e].mat,b_vec,elements[e].Xe,
                    ue,elements[e].statVar,elements[e].statVar0)
            
        elif elements[e].mat.Element_type==2:
            #print('element type:', elements[e].mat.Element_type)
            
            [fvole,fdyne,finte,Ke]=element_linear_elastic_1d(
                    e,ndm,ndf,nen,nqp,nst,
                    elements[e].mat,b_vec,elements[e].Xe,
                    ue,elements[e].statVar,elements[e].statVar0)

        else:
            print('element type:',' not defined - check input!')
            
            
        # assembling of element force vectors
        fint[gdof]+=finte
        fvol[gdof]+=fvole
        fdyn[gdof]+=fdyne
        
        # assembly of element stiffness Ke into global stiffness K
        K[np.ix_(gdof,gdof)]+=Ke
    
    return 0