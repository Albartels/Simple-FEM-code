#!/usr/bin/env python3
# -*- coding: utf-8 -*-
###############################################################################
# Simple FE-code for developing material routines (without contact)
# The code is written using free accessible python libaries 
# (numpy,time,sys,vtk), where vtk is not default installed usually.
# For the output .vtu files are generated which can be viewed using paraview.
#
# coded by abartels
###############################################################################

# for installation of vtk
#$ conda install -c anaconda vtk

import numpy as np
import time


# self defined routines and modules
from src_FEM.main.structs import *
from src_FEM.main.get_input import *
from src_FEM.main.help_functions import *
from src_FEM.main.assembly import *
from src_FEM.main.solve import *
from src_FEM.main.output_vtk import *

# boundary value problems
from examples.input_ex_truss import *
from examples.input_3el_plane_strain import *
from examples.input_beam import *


###############################################################################
# start of defining input
###############################################################################
start_time=time.time()

# example 1: of a 1d truss: linear elastic material response
# 2 noded line elements
#[io_file_name,nparts,nmat,ndf,ndm,nst,nnp,nel,nen,nqp,X1,connectivity,
#    material_number,material,
#    drlt,drltLoad,neum,neumLoad,b_vec,loadcurve,
#    solveparameter]=input_ex_1d_truss()

# example 2: of a 3 element test in 2d - plane strain: linear elastic material response
# 4-noded linear elements
[io_file_name,nparts,nmat,ndf,ndm,nst,nnp,nel,nen,nqp,X1,connectivity,
    material_number,material,
    drlt,drltLoad,neum,neumLoad,b_vec,loadcurve,
    solveparameter]=input_3el_plane_strain()   

# example 3: of 3d beam: linear elastic material response
# 8 noded linear 3d elements
#[io_file_name,nparts,nmat,ndf,ndm,nst,nnp,nel,nen,nqp,X1,connectivity,
#    material_number,material,
#    drlt,drltLoad,neum,neumLoad,b_vec,loadcurve,
#    solveparameter]=input_3d_beam()

###############################################################################
# end of defining input
###############################################################################



###############################################################################
# getting input data
###############################################################################

[elements,solveParam,numallDofs,allDofs,numdrltDofs,drltDofs,
    numneumDofs,neumDofs,numfreeDofs,freeDofs] = input_data(nmat,nel,
    nen,ndm,ndf,nst,nnp,
    material,material_number,connectivity,X1,
    drlt,neum,
    solveparameter)
            
###############################################################################
# end getting input
###############################################################################
pre_time=time.time()
solve_time0=pre_time
print('Preprocessing Time: %12.8fs \n' %(pre_time-start_time))

###############################################################################
# start analysis
###############################################################################

# creating displacement vector
u=np.zeros((ndf*nnp,2))

tot_time=0.0
step=0

while (tot_time<solveParam.stopTime and step<=(solveParam.maxstep-1)):
    ###
    # updating time and step
    tot_time+=solveParam.dt
    step+=1
    print('Load Time: %6.4f   Step: %2i' %(tot_time,step))
    print('Load Factor:', scale_factor(loadcurve, tot_time))
    
    # define rhs: forces due to neumann load fsur
    fsur=np.zeros(ndf*nnp)
    fsur[neumDofs]=neumLoad*scale_factor(loadcurve, tot_time)

    # insert initial values of solution vector: dirichlet boundary
    u[drltDofs,0]=drltLoad*scale_factor(loadcurve, tot_time)  #for current time  
    
    # initializing force vectors: internal, volume, dynamic
    fint=np.zeros(nnp*ndf)
    fvol=np.zeros(nnp*ndf)
    fdyn=np.zeros(nnp*ndf)
    
    # initialize stiffness matrix
    K=np.zeros((nnp*ndf,nnp*ndf))
    
    # define residual: also for dummy case + iteration counter
    rsn=1.0
    iter1=0
    
    # start of residual iteration -->> till solution
    while (rsn>solveParam.tol and (iter1<=solveParam.maxiter)):
        
        #######################################################################
        # start of assembling
        #######################################################################
        s_assembly_time=time.time()
        
        # initializing force vectors: internal, volume, dynamic
        fint=np.zeros(nnp*ndf)
        fvol=np.zeros(nnp*ndf)
        fdyn=np.zeros(nnp*ndf)
    
        # initialize stiffness matrix
        K=np.zeros((nnp*ndf,nnp*ndf))
        assembly(ndm,ndf,nen,nst,nnp,nel,nqp,elements,b_vec,u,
                 fvol,fdyn,fint,K)
        
        e_assembly_time=time.time()                
        #######################################################################
        # end of assembling
        #######################################################################
    

        #######################################################################
        # start of solution procedure
        #######################################################################
        s_solving_time=time.time()
        
        # define residuums vector
        rsd=np.zeros(ndf*nnp)
        rsd[freeDofs]=fint[freeDofs]+fvol[freeDofs] \
            +fsur[freeDofs]+fdyn[freeDofs]
        
        [rsn,u]=solve(ndf,nnp,iter1,numfreeDofs,freeDofs, \
              solveParam.tol,rsd,K,u)
        
        e_solving_time=time.time()
        #######################################################################
        # end of solution procedure
        #######################################################################


        # iteration counter
        print('Iter.: %2i # Res. 2-norm: %+12.7E \
# Ass. Time: %8.8fs # Sol. Time: %8.8fs'\
              %(iter1,rsn,e_assembly_time-s_assembly_time,\
                e_solving_time-s_solving_time))
        iter1+=1
        ### end of while loop
        
    # for new time step: 
    # update of difference \Delta u^{iter} = u^{iter}_{t+1} - u_{t}
    u[0:nnp*ndf,1]=u[0:nnp*ndf,0]

    # for new time step: update of state variables
    for ee in range(nel):
        elements[ee].statVar0=elements[ee].statVar
            
    # computing reaction forces
    duF=u[np.ix_(freeDofs),0]
    KDF=K[np.ix_(drltDofs,freeDofs)]
    #
    freac=np.zeros(nnp*ndf)
    freac[drltDofs]=KDF.dot(duF[0])

    solve_time=time.time()
    print('Solution time for %2i. step: %12.8fs'%(step,solve_time-solve_time0))
    solve_time0=solve_time
    ###########################################################################
    # start of writing output -->> .vtk
    ###########################################################################

    output_vtk(io_file_name,ndm,ndf,nen,nel,nnp,step,tot_time,
        X1,elements,fsur,fvol,fdyn,freac,u)

    ###########################################################################
    # end of writing output -->> .vtk
    ###########################################################################
    
    
    ###########################################################################
    # start of writing output terminal and file
    ###########################################################################
    print('\n'+'*'*72)
    print('*'*3+' D I S P L A C E M E N T S')
    print('*'*72)
    output_text='  Node'
    #output_text_coord=['%1i Coord     ' %(i+1) for i in range(ndm)]
    #output_text_displ=['%1i Displ     ' %(i+1) for i in range(ndf)]
    for ii in range(ndm):
        output_text+='      %1i Coord' %(ii+1)
    for ii in range(ndf):
        output_text+='      %1i Displ' %(ii+1) 
    print(output_text)
    for inp in range(nnp):
        output_text='%6i' %(inp+1)
        for ii in range(ndm):
            output_text+='  %+8.4E' %X1[inp,ii]
        for ii in range(ndf):
            output_text+='  %+8.4E' %u[inp*ndf+ii,0]
        print(output_text)    
    print('\n'+'*'*72)
    print('*'*3+' F O R C E S')
    print('*'*72)
    output_text='  Node'
    #output_text_coord=['%1i Coord     ' %(i+1) for i in range(ndm)]
    #output_text_displ=['%1i Displ     ' %(i+1) for i in range(ndf)]
    for ii in range(ndm):
        output_text+='      %1i Coord' %(ii+1)
    for ii in range(ndf):
        output_text+='      %1i Force' %(ii+1) 
    print(output_text)
    for inp in range(nnp):
        output_text='%6i' %(inp+1)
        for ii in range(ndm):
            output_text+='  %+8.4E' %X1[inp,ii]
        for ii in range(ndf):
            output_text+='  %+8.4E' %(fsur[inp*ndf+ii]+fvol[inp*ndf+ii]+\
                                      fdyn[inp*ndf+ii]+freac[inp*ndf+ii])
        print(output_text)    
    print('\n')
    ###########################################################################
    # end of writing output terminal and file
    ###########################################################################

    post_time=time.time()
    print('Postprocessing time for %2i. step: %12.8fs \n'%(step,post_time-solve_time0))

print('Total computation Time: %12.8fs' %(time.time()-start_time))



