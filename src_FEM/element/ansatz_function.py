#!/usr/bin/env python3
# -*- coding: utf-8 -*-
###############################################################################
# The function gaussintegrationdata() defines the Gauss points in the local 
# coordinates q and the weighting factors w for the element integration.
# The function shape_function() defines the shape functiosn and its derivatives
# w.r.t the local coordinates q of the different element types 
#
# coded by abartels
###############################################################################


import numpy as np
import sys
from numpy import sqrt


def gaussintegrationdata(ndm,nqp):
    
    # initialization of weighting factors and Gauss points
    w=np.zeros(nqp)
    q=np.zeros((nqp,ndm))
    
    if ndm==1:
        if nqp==1:
            w[0]=2.0
            q[0][0]=0.0
        else:
            sys.exit('so far no further implementation for 1D elements')
    elif ndm==2:
        if nqp==1:
            w[0]=4.0;
            q[0][0]=0.0; q[0][1]=0.0; 
        elif nqp==4:
            w[0]=1.0; w[1]=1.0; w[2]=1.0; w[3]=1.0;
            a1=1.0/sqrt(3.0);
            q[0][0]=-a1; q[0][1]=-a1;
            q[1][0]= a1; q[1][1]=-a1;
            q[2][0]= a1; q[2][1]= a1;
            q[3][0]=-a1; q[3][1]= a1;
        else:
            sys.exit('2D error: so far no implementation for nqp:',str(nqp))
    elif ndm==3:
        if nqp==8:
            w[0]=1.0; w[1]=1.0; w[2]=1.0; w[3]=1.0;
            w[4]=1.0; w[5]=1.0; w[6]=1.0; w[7]=1.0;
            a1=1.0/sqrt(3.0);
            q[0][0]=-a1; q[0][1]=-a1; q[0][2]=-a1;
            q[1][0]= a1; q[1][1]=-a1; q[1][2]=-a1;
            q[2][0]= a1; q[2][1]= a1; q[2][2]=-a1;
            q[3][0]=-a1; q[3][1]= a1; q[3][2]=-a1;
            q[4][0]=-a1; q[4][1]=-a1; q[4][2]= a1;
            q[5][0]= a1; q[5][1]=-a1; q[5][2]= a1;
            q[6][0]= a1; q[6][1]= a1; q[6][2]= a1;
            q[7][0]=-a1; q[7][1]= a1; q[7][2]= a1;
        else:
            sys.exit('so far no further implementation for 3D elements')
    else:
        sys.exit('fatal error ndm:',str(ndm))
    
    
    return (w,q)


def shape_function(ndm,nen,nqp,q):
    
    # initialization of shape function and its derivative
    shp=np.zeros(nen*nqp)
    dshp=np.zeros((nen*nqp,ndm))
    
    if ndm==1:
        if nen==2:
            for ii in range(nqp):
                shp[0+nqp*ii]=0.5*(1.0-q[ii][0]);
                shp[1+nqp*ii]=0.5*(1.0+q[ii][0]);
    			# derivative of shp w.r.t general. coordinate
                dshp[0+nqp*ii][0]=-0.5;
                dshp[1+nqp*ii][0]=0.5;
        else:
            sys.exit('so far no further implementation for 1D elements')
    elif ndm==2:
        if nen==4:
            for ii in range(nqp):
                shp[0+nqp*ii]=0.25*(1.0-q[ii][0])*(1.0-q[ii][1]);
                shp[1+nqp*ii]=0.25*(1.0+q[ii][0])*(1.0-q[ii][1]);
                shp[2+nqp*ii]=0.25*(1.0+q[ii][0])*(1.0+q[ii][1]);
                shp[3+nqp*ii]=0.25*(1.0-q[ii][0])*(1.0+q[ii][1]);
                # derivative of shp w.r.t general. coordinate
                dshp[0+nqp*ii][0]=-0.25*(1.0-q[ii][1]);
                dshp[0+nqp*ii][1]=-0.25*(1.0-q[ii][0]);
                dshp[1+nqp*ii][0]= 0.25*(1.0-q[ii][1]);
                dshp[1+nqp*ii][1]=-0.25*(1.0+q[ii][0]);
                dshp[2+nqp*ii][0]= 0.25*(1.0+q[ii][1]);
                dshp[2+nqp*ii][1]= 0.25*(1.0+q[ii][0]);
                dshp[3+nqp*ii][0]=-0.25*(1.0+q[ii][1]);
                dshp[3+nqp*ii][1]= 0.25*(1.0-q[ii][0]);
        else:
            sys.exit('2D error: so far no implementation for nen:',str(nen))
    elif ndm==3:
        if nen==8:
            for ii in range(nqp):
                shp[0+nqp*ii]=0.125*(1.0-q[ii][0])*(1.0-q[ii][1])*(1.0-q[ii][2]);
                shp[1+nqp*ii]=0.125*(1.0+q[ii][0])*(1.0-q[ii][1])*(1.0-q[ii][2]);
                shp[2+nqp*ii]=0.125*(1.0+q[ii][0])*(1.0+q[ii][1])*(1.0-q[ii][2]);
                shp[3+nqp*ii]=0.125*(1.0-q[ii][0])*(1.0+q[ii][1])*(1.0-q[ii][2]);
                shp[4+nqp*ii]=0.125*(1.0-q[ii][0])*(1.0-q[ii][1])*(1.0+q[ii][2]);
                shp[5+nqp*ii]=0.125*(1.0+q[ii][0])*(1.0-q[ii][1])*(1.0+q[ii][2]);
                shp[6+nqp*ii]=0.125*(1.0+q[ii][0])*(1.0+q[ii][1])*(1.0+q[ii][2]);
                shp[7+nqp*ii]=0.125*(1.0-q[ii][0])*(1.0+q[ii][1])*(1.0+q[ii][2]);

			    # derivatives of shape functions w.r.t. local coordinates
                dshp[0+nqp*ii][0]=-0.125*(1.0-q[ii][1])*(1.0-q[ii][2]);
                dshp[0+nqp*ii][1]=-0.125*(1.0-q[ii][0])*(1.0-q[ii][2]);
                dshp[0+nqp*ii][2]=-0.125*(1.0-q[ii][0])*(1.0-q[ii][1]);
                
                dshp[1+nqp*ii][0]= 0.125*(1.0-q[ii][1])*(1.0-q[ii][2]);
                dshp[1+nqp*ii][1]=-0.125*(1.0+q[ii][0])*(1.0-q[ii][2]);
                dshp[1+nqp*ii][2]=-0.125*(1.0+q[ii][0])*(1.0-q[ii][1]);
                
                dshp[2+nqp*ii][0]= 0.125*(1.0+q[ii][1])*(1.0-q[ii][2]);
                dshp[2+nqp*ii][1]= 0.125*(1.0+q[ii][0])*(1.0-q[ii][2]);
                dshp[2+nqp*ii][2]=-0.125*(1.0+q[ii][0])*(1.0+q[ii][1]);
                
                dshp[3+nqp*ii][0]=-0.125*(1.0+q[ii][1])*(1.0-q[ii][2]);
                dshp[3+nqp*ii][1]= 0.125*(1.0-q[ii][0])*(1.0-q[ii][2]);
                dshp[3+nqp*ii][2]=-0.125*(1.0-q[ii][0])*(1.0+q[ii][1]);
                
                dshp[4+nqp*ii][0]=-0.125*(1.0-q[ii][1])*(1.0+q[ii][2]);
                dshp[4+nqp*ii][1]=-0.125*(1.0-q[ii][0])*(1.0+q[ii][2]);
                dshp[4+nqp*ii][2]= 0.125*(1.0-q[ii][0])*(1.0-q[ii][1]);
                
                dshp[5+nqp*ii][0]= 0.125*(1.0-q[ii][1])*(1.0+q[ii][2]);
                dshp[5+nqp*ii][1]=-0.125*(1.0+q[ii][0])*(1.0+q[ii][2]);
                dshp[5+nqp*ii][2]= 0.125*(1.0+q[ii][0])*(1.0-q[ii][1]);
                
                dshp[6+nqp*ii][0]= 0.125*(1.0+q[ii][1])*(1.0+q[ii][2]);
                dshp[6+nqp*ii][1]= 0.125*(1.0+q[ii][0])*(1.0+q[ii][2]);
                dshp[6+nqp*ii][2]= 0.125*(1.0+q[ii][0])*(1.0+q[ii][1]);
                
                dshp[7+nqp*ii][0]=-0.125*(1.0+q[ii][1])*(1.0+q[ii][2]);
                dshp[7+nqp*ii][1]= 0.125*(1.0-q[ii][0])*(1.0+q[ii][2]);
                dshp[7+nqp*ii][2]= 0.125*(1.0-q[ii][0])*(1.0+q[ii][1]);
        else:
            sys.exit('so far no further implementation for 3D elements')
    else:
        sys.exit('fatal error ndm:',str(ndm))    
    
    
    return (shp,dshp)




