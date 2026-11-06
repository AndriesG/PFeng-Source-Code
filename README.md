# Simple guide to the files that were changed/optimised.

Most files of relevance are in MiniFE/ref/src

## *.yaml files

These files contain the output of the several runs of the program

## matvecnew.hpp & driver.hpp

An additional file that implements the optimised Functor passed to cg_solve.<br>
driver.hpp then creates and passes the new Functor to cg_solve in lines 284-285. The old code is commented out in lines 281-283.

## SparseMatrix_functions.hpp
this file contained the modified access code for the use of the ELL Matrix implementation
