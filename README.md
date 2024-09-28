# Forward adjoint model
![forward adjoint_picture](https://github.com/BIOPTIMOD/Forward_Adjoint/blob/main/DOC/PICTURES/jmse-09-00176-g001.png)
## What is the Forwar_adjoint model?
Forward_adjoint is a forwrad and adjoint models to solve the propagation of light along the water column. The model is used with the code [ogstm](https://github.com/inogs/ogstm) and also with the [FABM tool](https://github.com/inogs/bfmforfabm). 
The model is based on fortran and also coded in python.
The forward model has been tested in several publications ([Lazzari et al., 2021](https://www.mdpi.com/2077-1312/9/2/176), [Alvarez et al., 2022](https://www.sciencedirect.com/science/article/pii/S0079661122000507?via%3Dihub), [Alvarez et al., 2023](https://bg.copernicus.org/articles/20/4591/2023/)). In the forward mode the model is used to derive light propagation influenced by optical constituents in the water column.

The adjoint model has been tested at the BOUSSOLE site and ca be used to retrieve optical constituents in sea waters such as chlorophyll, detritus or colored dissolved organic carbon, [Lazzari et al., 2024](https://rdcu.be/dVi82).

# How to download

To download the code cloning the repository and compile the library do the following:

git clone git@github.com:BIOPTIMOD/Forward_Adjoint.git

cd Forward_Adjoint/src

git checkout tags/release-1.0

make clean

make lib

The default compiler is gfortran, but editing the FC variable in the makefile you can select ifort.



