#!/usr/bin/env python3
# -*- coding: utf-8 -*-
###############################################################################
# In this function the input data is sorted. Imported for the structs module
# the structs solveParam, mat and elements are filled. 
#
# coded by abartels
###############################################################################

from src_FEM.main.structs import *

def input_data(nmat,nel,nen,ndm,ndf,nst,nnp,
        material,material_number,connectivity,X1,
        drlt,neum,
        solveparameter):
    
    # include solveparameter
    solveParam=solve_param()
    solveParam.dt=solveparameter[0]
    solveParam.stopTime=solveparameter[1]
    solveParam.maxstep=solveparameter[3]
    solveParam.maxiter=solveparameter[2]
    solveParam.tol=solveparameter[4]

    # filling material parameter to mat
    mat=[mat_param() for i in range(nmat)]

    for m in range(nmat):
        mat[m].Mat_number=material[m][0]
        mat[m].Element_type=material[m][1]
        mat[m].Parameter=material[m][2:]
    


    # multiplying element class
    elements=[element() for i in range(nel)]

    # filling elements as struct with connectivity, elements dof's, 
    # stateVar, statVar0
    for e in range(nel):
        elements[e].mat=mat[material_number[e]-1]
        #
        elements[e].cn=[0]*nen
        for ii in range(nen):
            elements[e].cn[ii]=connectivity[e][ii]
            #
        elements[e].edof=[[0 for i in range(nen)] for j in range(ndf)]
        for ii in range(nen):
            for jj in range(ndf):
                elements[e].edof[jj][ii]=ndf*connectivity[e][ii]+jj-ndf

        # coordinates
        elements[e].Xe= [[0 for i in range(nen)] for j in range(ndm)]    
        for ii in range(nen):
            for jj in range(ndm):
                elements[e].Xe[jj][ii]=X1[connectivity[e][ii]-1][jj]
                    #
        # initializing state variables
        elements[e].statVar=[0 for i in range(nst)]
        elements[e].statVar0=[0 for i in range(nst)] 
    
    # define vector with all dofs
    numallDofs=nnp*ndf
    allDofs=list(range(numallDofs))

    numdrltDofs=len(drlt)
    drltDofs=[0]*numdrltDofs
    for ii in range(numdrltDofs):
        drltDofs[ii]=(drlt[ii][0]-1)*ndf+drlt[ii][2]-1
      
    numneumDofs=len(neum)
    neumDofs=[0]*numneumDofs
    for ii in range(numneumDofs):
        neumDofs[ii]=(neum[ii][0]-1)*ndf+neum[ii][2]-1
            
    
    numfreeDofs=numallDofs-numdrltDofs
    freeDofs=list(set(allDofs)-set(drltDofs))
    
    
    return (elements,solveParam,numallDofs,allDofs,numdrltDofs,drltDofs,
            numneumDofs,neumDofs,numfreeDofs,freeDofs)

