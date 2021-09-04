#!/bin/bash
# Build hipSYCL on Intel devcloud
# use
# qsub -l nodes=1:gen9:ppn=2 -d . gethipSYCL.sh 
cd $HOME
wget https://downloads.sourceforge.net/project/boost/boost/1.77.0/boost_1_77_0.tar.gz
tar -xf boost_1_77_0.tar.gz
cd boost_1_77_0/
./bootstrap.sh --prefix=$HOME/boost_1_77_0-install
./b2 --build-dir=$HOME/boost_1_77_0-build toolset=intel stage
./b2 install --prefix=$HOME/boost_1_77_0-install --toolset=intel --build-type=complete --layout=versioned
cd ..
git clone --recurse-submodules https://github.com/illuhad/hipSYCL
cd hipSYCL/
git checkout fe8465c
mkdir build
cd build/
cmake -DCMAKE_INSTALL_PREFIX=$HOME/hipSYCL-install -DBOOST_ROOT=$HOME/boost_1_77_0-install -DCMAKE_C_COMPILER=icc -DCMAKE_CXX_COMPILER=icpc ..
#cmake -DCMAKE_INSTALL_PREFIX=$HOME/hipSYCL-install -DBOOST_ROOT=$HOME/boost_1_77_0-install -DCMAKE_C_COMPILER=icx -DCMAKE_CXX_COMPILER=icpx ..
make
make install
cd ..
cd ..
