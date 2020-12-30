#!/usr/bin/env python3
# -*- coding: utf-8 -*-
###############################################################################
# For this function the module vtk is required. 
# In this function the output is generated in a .vtu/.vtk format. The generated
# data is stored accordingly, here in an unstructured grid format. 
#
# coded by abartels
###############################################################################

#from vtk import *
import vtk

# self defined module
from src_FEM.main.create_unstructuredgrid_vtk import *

def output_vtk(io_file_name,ndm,ndf,nen,nel,nnp,step,time,
        X1,elements,fsur,fvol,fdyn,freac,u):
    
    # in this routine only writing data into files
    # point, elements and connectivity is managed by 
    # create_unstructuredgrid_vtk()
    
    print('creation of output file:',io_file_name+str(step)+'.vtu')
    
    
    # storing the time
    ttime=vtk.vtkDoubleArray()
    ttime.SetName('TIME')
    ttime.SetNumberOfTuples(1)
    ttime.SetTuple1(0,time)
    
    # start writing point data into vtk arrays
    # vector data: displacement: always 3D
    u_data=vtk.vtkDoubleArray()
    u_data.SetName('Displacements')
    u_data.SetNumberOfComponents(3)
    # vector data: volume, dynamic and surface forces
    fvol_data=vtk.vtkDoubleArray()
    fvol_data.SetName('F_vol')
    fvol_data.SetNumberOfComponents(3)
    fdyn_data=vtk.vtkDoubleArray()
    fdyn_data.SetName('F_dyn')
    fdyn_data.SetNumberOfComponents(3)
    fsur_data=vtk.vtkDoubleArray()    
    fsur_data.SetName('F_sur')
    fsur_data.SetNumberOfComponents(3)
    freac_data=vtk.vtkDoubleArray()    
    freac_data.SetName('F_freac')
    freac_data.SetNumberOfComponents(3)
    for inp in range(nnp):
        if (ndm==1):
            u_data.InsertNextTuple3(u[ndf*inp,0],0.0,0.0)
            fvol_data.InsertNextTuple3(fvol[ndf*inp],0.0,0.0)
            fdyn_data.InsertNextTuple3(fdyn[ndf*inp],0.0,0.0)
            fsur_data.InsertNextTuple3(fsur[ndf*inp],0.0,0.0)
            freac_data.InsertNextTuple3(freac[ndf*inp],0.0,0.0)
        if (ndm==2):
            u_data.InsertNextTuple3(u[ndf*inp,0],u[ndf*inp+1,0],0.0)
            fvol_data.InsertNextTuple3(fvol[ndf*inp],fvol[ndf*inp+1],0.0)
            fdyn_data.InsertNextTuple3(fdyn[ndf*inp],fdyn[ndf*inp+1],0.0)
            fsur_data.InsertNextTuple3(fsur[ndf*inp],fsur[ndf*inp+1],0.0)
            freac_data.InsertNextTuple3(freac[ndf*inp],freac[ndf*inp+1],0.0)
        if (ndm==3):
            u_data.InsertNextTuple3(u[ndf*inp,0],u[ndf*inp+1,0],u[ndf*inp+2,0])
            fvol_data.InsertNextTuple3(fvol[ndf*inp],fvol[ndf*inp+1],fvol[ndf*inp+2])
            fdyn_data.InsertNextTuple3(fdyn[ndf*inp],fdyn[ndf*inp+1],fdyn[ndf*inp+2])
            fsur_data.InsertNextTuple3(fsur[ndf*inp],fsur[ndf*inp+1],fsur[ndf*inp+2])
            freac_data.InsertNextTuple3(freac[ndf*inp],freac[ndf*inp+1],freac[ndf*inp+2])
    ### end point data
    
    
    # start writing cell data
    # material number
    mat_data=vtk.vtkDoubleArray()
    mat_data.SetName('Material_Element-Type')
    mat_data.SetNumberOfComponents(1)
    volume_data=vtk.vtkDoubleArray()
    volume_data.SetName('Volume_Element')
    volume_data.SetNumberOfComponents(1)
    stress_data=vtk.vtkDoubleArray()
    stress_data.SetName('Sigma')
    if ndm==1:
        stress_data.SetNumberOfComponents(1)
        for iel in range(nel):
            mat_data.InsertNextValue(elements[iel].mat.Mat_number)
            volume_data.InsertNextValue(elements[iel].statVar[1])
            stress_data.InsertNextValue(elements[iel].statVar[0])
    else:
        stress_data.SetNumberOfComponents(9)    
        for iel in range(nel):
            mat_data.InsertNextValue(elements[iel].mat.Mat_number)
            volume_data.InsertNextValue(elements[iel].statVar[9])
            stress_data.InsertNextTuple9(elements[iel].statVar[0],
                elements[iel].statVar[1],elements[iel].statVar[2],
                elements[iel].statVar[3],elements[iel].statVar[4],
                elements[iel].statVar[5],elements[iel].statVar[6],
                elements[iel].statVar[7],elements[iel].statVar[8])   
    ### end cell data
    
    
    # writing into unstructured grid
    unstructuredGrid=vtk.vtkUnstructuredGrid()
    create_unstructuredgrid_vtk(ndm,nen,nel,nnp,X1,elements,
                                unstructuredGrid)
    unstructuredGrid.GetFieldData().AddArray(ttime)
    unstructuredGrid.GetPointData().AddArray(u_data)
    unstructuredGrid.GetPointData().AddArray(fsur_data)
    unstructuredGrid.GetPointData().AddArray(fvol_data)
    unstructuredGrid.GetPointData().AddArray(fdyn_data)
    unstructuredGrid.GetPointData().AddArray(freac_data)
    unstructuredGrid.GetCellData().AddArray(mat_data)
    unstructuredGrid.GetCellData().AddArray(volume_data)
    unstructuredGrid.GetCellData().AddArray(stress_data)


    # writing to .vtk file in XML format
    writer=vtk.vtkXMLUnstructuredGridWriter()
    vtk_file_name=io_file_name+str(step)+'.vtu'
    writer.SetFileName(vtk_file_name)
    writer.SetInputData(unstructuredGrid)
    writer.SetDataModeToAscii()
    writer.Write()

    return 0