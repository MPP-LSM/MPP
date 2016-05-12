#!/bin/sh

git clone https://bitbucket.org/petsc/petsc.git

PETSC_GIT_HASH=c41c7662de68b036bda6be236f939e8b55959cb0

cd petsc

git checkout ${PETSC_GIT_HASH}

export PETSC_DIR=$PWD
export PETSC_ARCH=linux-gnu

./configure PETSC_ARCH=linux-gnu --with-mpi=1 --with-debug=$DEBUG --with-shared-libraries=1 --download-hdf5 --download-metis --download-parmetis --with-c2html=0

make 

