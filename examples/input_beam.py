#!/usr/bin/env python3
# -*- coding: utf-8 -*-
###############################################################################
# Input example: 3D beam example: linear elastic material behaviour
# 8 node elements in 3d (linear interpolation)
#
# numpy arrays should be used
# 
# output:
# name: name of input. is used for generation of output
# ndf: number degrees of freedom
# nmat: number of materials
# nst: number of state variables
# nqp: number of quadrature points
# ndm: number of dimension
# nnp: number of node points
# nel: number of elements
# nen: number of elment nodes
# X1: Coordinates of Nodes in array form: 1st is node 1, 2nd is node 2,...
# connectivity: array of element related nodes
# material_number: array of Element_type (see assembly.py) which are assigned 
# drlt: array including dirichlet boundary condition information
# drltLoad: array describing the dirichlet boundary condition value of drlt
# neum: similar to drlt, only Neumann boundary condition
# neumLoad: similar to drltLoad, only Neumann boundary condition
# b_vec: body force vector: optional (default is zero)
# loadcurve: array describing the load-time relation
# solveparameter: [dt,stopTime,maxiter,maxstep,tol]=[time increment,stop time, 
#   maximal number of iterations before stop of execution, 
#   maximum number of steps before stop of execution, 
#   solution precision tolerance]
#
# coded by abartels
###############################################################################

import numpy as np
import pandas as pd
import os
import sys

get_path=os.getcwd()    # path of FEM.py, if file is loaded elsewhere the path 
                        # location needs to be defined --> .../examples/.. .csv

# 1 sided clamped beam
def input_3d_beam():
    name='beam_'
    #self defined variables
    nmat=2; ndf=3; nparts=2; nst=10; nqp=8;
    
    #coordinates: import using pandas
    X1=pd.read_csv(get_path+'/examples/beam_3d_coordinates.csv',header=0,
        names=['NodeID','x','y','z'])
    X1=X1.values[:,1:]
    [nnp,ndm]=np.shape(X1)

    #connectivity: import using pandas
    connectivity=pd.read_csv(get_path+'/examples/beam_3d_elements.csv',
        header=0,names=['Id','mat','1','2','3','4','5','6','7','8'])
    connectivity=connectivity.values[:,2:].astype(int)
    [nel,nen]=np.shape(connectivity)
    # assigning material number
    material_number=[1]*nel
    #material_number[2]=2
    
    
    # defining material parameter
    material=[[0 for i in range(4)] for j in range(nmat)]
    material[0][0]=1 # material 1
    material[0][1]=1
    material[0][2:4]=[1000,0.3]
    material[1][0]=1 # material 2
    material[1][1]=2
    material[1][2:4]=[1500,0.2]
    
    
    # dirichlet boundary conditons: pandas import
    drlt_l=pd.read_csv(get_path+'/examples/beam_3d_drlt.csv',header=0,
        names=['node','LC','dof','value'])
    drlt=drlt_l.values[:,:3].astype(int)
    drltLoad=drlt_l.values[:,3]
    
    # neumenn boundary conditions: neum and neumLoad need to be defined
    # although zero
    neum=np.zeros([1,3]).astype(int)
    neumLoad=np.zeros(1)
    #neum[0]=[7,0,1]
    #neumLoad[0]=100.0
    
    # body force vector
    b_vec=np.zeros(3)
    
    # load curve
    loadcurve=np.array([[0.0,1.0],[0.0,1.0]]) # first row time, second row slope

    # solution properties
    dt=0.1
    stopTime=0.3
    maxiter=10
    maxstep=2
    tol=10e-8
    
    
    solveparameter=[dt,stopTime,maxiter,maxstep,tol]
    
    return (name,nparts,nmat,ndf,ndm,nst,nnp,nel,nen,nqp,X1,connectivity,
            material_number,material,
            drlt,drltLoad,neum,neumLoad,b_vec,loadcurve,
            solveparameter)