#!/bin/bash

# Load necessary Hyak modules
module load ssmc/miniconda/3.9
module load cuda/11.8.0
module load cmake/3.20.0
module load gcc/10.2.0

# Get HOOMD source code from Github
git clone --recursive https://github.com/the-z-lab/Dielectric.git
cd Dielectric/
git checkout d_Eq
cd ../

# Create and activate a Python virtual environment
python -m venv hoomd-venv --system-site-packages
source hoomd-venv/bin/activate

# Add the AMGX libraries and include paths
export LD_LIBRARY_PATH=/gscratch/zeelab/haoqing/HOOMD_pluginfix/AMGX/build:$LD_LIBRARY_PATH
export AMGX_INCLUDE_DIRECTORY=/gscratch/zeelab/haoqing/HOOMD_pluginfix/AMGX/base/include
export AMGX_LIBRARY=/gscratch/zeelab/haoqing/HOOMD_pluginfix/AMGX/build

# Build HOOMD-dielectric-plugin
cd Dielectric/
rm -rf build
mkdir build
cd build

cmake \
  -DAMGX_INCLUDE_DIRECTORY=$AMGX_INCLUDE_DIRECTORY \
  -DCMAKE_PREFIX_PATH="/gscratch/zeelab/haoqing/HOOMD_pluginfix/AMGX/build" \
  ..

make -j6
make install
