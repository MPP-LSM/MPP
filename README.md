# MPP [![Build Status](https://travis-ci.org/MPP-LSM/MPP.svg?branch=master)](https://travis-ci.org/MPP-LSM/MPP)

Multi-Physics Problem (MPP) library is a standalone library that
solves biophysics problems relevant to global land surface models (LSMs).


## Installation

### Requirements

* Mac
* Git 
* Fortran compiler
* CMake
* PETSc
* Message Passing Interface


### Installation instructions

#### Install PETSc

1. Clone PETSc and check out the **supported** version.

```
git clone https://bitbucket.org/petsc/petsc petsc
cd petsc
git checkout v3.8.3
```

2. (Detailed PETSc installation instructions are available [here](http://www.mcs.anl.gov/petsc/documentation/installation.html))

```

export PETSC_DIR=$PWD
export PETSC_ARCH=gcc5-debug

./configure \
--with-cc=/opt/local/bin/mpicc-openmpi-gcc5   \
--with-cxx=/opt/local/bin/mpicxx-openmpi-gcc5 \
--with-fc=/opt/local/bin/mpif90-openmpi-gcc5  \
--with-mpiexec=/opt/local/bin/mpiexec-openmpi-gcc5 \
--with-blas-lapack-lib=/System/Library/Frameworks/Accelerate.framework/Versions/Current/Accelerate
```
