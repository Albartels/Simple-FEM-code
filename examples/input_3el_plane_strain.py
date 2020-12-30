#!/usr/bin/env python3
# -*- coding: utf-8 -*-
###############################################################################
# Input example: 3 element test with Neumann boundary condition
# numpy arrays should be used
# Problem definition in Plane strain: 4 noded elements in 2d
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

def input_3el_plane_strain():
    name='L_shape_test_'
    #self defined variables
    nmat=1; ndf=2; nparts=2; nst=10; nqp=4;
    
    #coordinates
    # x-direction, y-direction
    X1=np.array([[0.0,0.0], #first node
                [1.0,0.0],  # second node
                [2.0,0.0],
                [0.0,1.0],
                [1.2,1.4],
                [2.0,1.0],
                [0.0,2.0],
                [1.0,2.0]])
    
    [nnp,ndm]=np.shape(X1)
    
    #connectivity
    connectivity=np.array([[1,2,5,4],
                  [2,3,6,5],
                  [4,5,8,7]])
    
    [nel,nen]=np.shape(connectivity)
    material_number=[1]*nel
    
    
    # defining material parameter
    material=[[0 for i in range(4)] for j in range(nmat)]
    material[0][0]=1 # material 1
    material[0][1]=1
    material[0][2:4]=[1000,0.3]  
    
    # dirichlet boundary conditons
    # drlt: node,LC,DOF
    drlt=np.zeros([4,3]).astype(int)
    drlt[0]=[1,0,1]
    drlt[1]=[1,0,2]
    drlt[2]=[3,0,1]
    drlt[3]=[3,0,2]
    drltLoad=np.zeros(4) # float type
    #drltLoad[1]=4.0
    

    # neumenn boundary conditions
    # neum: node,LC,DOF
    neum=np.zeros([1,3]).astype(int) # integer type
    neumLoad=np.zeros(1)
    neum[0]=[7,0,1]
    neumLoad[0]=100.0 # float type
    
    # body force vector
    b_vec=np.zeros(2)
    
    # load curve
    loadcurve=np.array([[0.0,1.0],[0.0,1.0]]) # first row time, second row slope

    
    dt=0.1
    stopTime=1
    maxiter=5
    maxstep=2
    tol=10e-8
    
    solveparameter=[dt,stopTime,maxiter,maxstep,tol]
    
    return (name,nparts,nmat,ndf,ndm,nst,nnp,nel,nen,nqp,X1,connectivity,
            material_number,material,
            drlt,drltLoad,neum,neumLoad,b_vec,loadcurve,
            solveparameter)