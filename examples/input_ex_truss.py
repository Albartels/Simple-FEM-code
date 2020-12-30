#!/usr/bin/env python3
# -*- coding: utf-8 -*-
###############################################################################
# Input example: 1D truss example
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


def input_ex_1d_truss():
    name='truss_ex_'
    #self defined variables
    nmat=2; ndf=1; nparts=1; nst=2; nqp=1;
    
    #coordinates
    el=10
    X1=np.linspace(0,10,el)
    X1=X1.reshape([10,1])
    (nnp,ndm)=np.shape(X1)

    #connectivity
    connectivity=np.array([[ii,ii+1] for ii in range(1,el)])
    
    [nel,nen]=np.shape(connectivity)
    # assigning material number
    material_number=[2]*nel
    
    
    # defining material parameter
    material=[[0 for i in range(4)] for j in range(nmat)]
    material[0][0]=1 # material 1
    material[0][1]=2
    material[0][2:4]=[1000,0.5] # e-modul and cross section
    material[1][0]=1 # material 2
    material[1][1]=2
    material[1][2:4]=[1500,0.7]
    
    
    # dirichlet boundary conditons
    drlt=np.zeros([1,3]).astype(int)
    drlt[0]=[1,0,1]
    drltLoad=np.zeros(1)
    #drltLoad[1]=4.0
    

    # neumenn boundary conditions
    neum=np.zeros([1,3]).astype(int)
    neumLoad=np.zeros(1)
    neum[0]=[10,0,1]
    neumLoad[0]=100.0

    # body force vector
    b_vec=np.zeros(1)
    
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
