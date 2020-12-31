# Simple-FEM-code

### A simple Finite Element Code for testing material routines and material models
The purpose of this Finite Element program is a have a simple Python-based enviroment for testing and demonstrating material models. The code is written in a generalized form such that even coupled problems like thermo-mechanics or magneto-mechanics can be included. However, also non-linearities in material evaluation laws can be included. This code is therefore not restricted to pure static problems. It can be extended also to dynamics or other field problems. In this version the FEM program accounts on simple Hookian material behavior in small strain theory.

### Restrictions of the Finite Element Code
Although its generalized program setup, the material behavior is restricted to small strain theory and Hookian material behavior. Contact mechanics is not included in this code. 

### Running the code
The following Finite Element code is written in Python. For execution of the program the standard module numpy and the visualization toolkit vtk are required. Under Pyton 3.6 the code was developed and can be started easily in Anaconda-Spyder. Simply copy the files and start to run the file FEM.py. The output is a vtk format, which can be read with the visualization program Paraview.

### Examples: 
In the main file FEM.py three boundary value problems can be tested. A 1D truss, a 2D plane strain example of a L-shape (3 elements) and a 3D beam, which is clamped on one side. Of course, it is possible to generate own boundary value problems. Only the input convention (see examples) should be guaranteed, such that the code runs successfully. 
