#!/bin/bash
# Build hipSYCL on Intel devcloud
# use
# qsub -l nodes=1:ppn=2 -d . gethipSYCL.sh 
cd $HOME
git clone --recurse-submodules https://github.com/illuhad/hipSYCL
cd hipSYCL/
git checkout fe8465c
mkdir build
cd build/
cmake -DCMAKE_INSTALL_PREFIX=$HOME/hipSYCL-install -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ ..
make
make install
cd ..
cd ..
