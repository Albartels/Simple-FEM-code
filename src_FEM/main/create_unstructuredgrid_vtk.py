#!/usr/bin/env python3
# -*- coding: utf-8 -*-
###############################################################################
# the function is reqired for the module output_vtk. An unstructured grid is
# defined using 3d vtkPoints and InsertNextCell and SetCell. See also vtk
# references 
#
# coded by abartels
###############################################################################

#from vtk import *
import vtk
import sys

def create_unstructuredgrid_vtk(ndm,nen,nel,nnp,X1,elements,
                                unstructuredGrid):


    ### nsert all points/coordinates: here in 3D independent on ndm
    points=vtk.vtkPoints()
    if ndm==1:
        for ii in range(nnp):
            points.InsertNextPoint(X1[ii,0],0.0,0.0)
    elif ndm==2:
        for ii in range(nnp):
            points.InsertNextPoint(X1[ii,0],X1[ii,1],0.0)
    else:
        for ii in range(nnp):
            points.InsertNextPoint(X1[ii,0],X1[ii,1],X1[ii,2])
            
    ### initialize all avaiable element types
    connectivity_line=vtk.vtkLine()
    #connectivity_triangle=vtk.vtkTriangle()
    connectivity_quad=vtk.vtkQuad()
    connectivity_hex=vtk.vtkHexahedron()
    # ToDo: other element type of vtk
    
    
    ### define the cells and write to connectivity list in vtk
    connectivity=vtk.vtkCellArray()
    for ee in range(nel):
        for jj in range(nen):
            if ndm==1 and nen==2:
                connectivity_line.GetPointIds(). \
                    SetId(jj,elements[ee].cn[jj]-1)
            #elif ndm==2 and nen==3:    
            #    connectivity_triangle.GetPointIds(). \
            #        SetId(jj,elements[ee].cn[jj]-1)
            elif ndm==2 and nen==4:
                connectivity_quad.GetPointIds(). \
                    SetId(jj,elements[ee].cn[jj]-1)
            elif ndm==3 and nen==8:
                connectivity_hex.GetPointIds(). \
                    SetId(jj,elements[ee].cn[jj]-1)
            else:
                sys.exit('so far no further implementation for this element')
        ### end element connectivity            
        ### add defined element connectivity to total connectivity
        if ndm==1 and nen==2:
            connectivity.InsertNextCell(connectivity_line)
        #elif ndm==2 and nen==3:
        #    connectivity.InsertNextCell(connectivity_triangle)            
        elif ndm==2 and nen==4:
            connectivity.InsertNextCell(connectivity_quad)
        elif ndm==3 and nen==8:
            connectivity.InsertNextCell(connectivity_hex)
        else:
            sys.exit('so far no further implementation for this element')
    ### end element loop

    ### write to unstructuredGrid --> vtk
    unstructuredGrid.SetPoints(points)
    if ndm==1 and nen==2:
        unstructuredGrid.SetCells(vtk.VTK_LINE,connectivity)
    #elif ndm==2 and nen==3:
    #    unstructuredGrid.SetCells(vtk.VTK_TRIANGLE,connectivity)
    elif ndm==2 and nen==4:
        unstructuredGrid.SetCells(vtk.VTK_QUAD,connectivity)
    elif ndm==3 and nen==8:
        unstructuredGrid.SetCells(vtk.VTK_HEXAHEDRON,connectivity)
    else:
        sys.exit('so far no further implementation for this element')
            
    return 0