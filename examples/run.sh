#!/bin/bash

# Load necessary Hyak modules
module load cuda/11.8.0

# Load hoomd by activating the python virtual environment in which HOOMD was built.
# Different file path needed for different users.
source /gscratch/zeelab/haoqing/HOOMD_pluginfix/hoomd-venv/bin/activate

# Add the amgx libraries to the file path. Should consider doing this automatically during install?
export LD_LIBRARY_PATH=/gscratch/zeelab/haoqing/HOOMD_pluginfix/AMGX/build:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/gscratch/zeelab/haoqing/HOOMD_pluginfix/Dielectric/build:$LD_LIBRARY_PATH

# Run simulations
python Dielectric_FD_grouptag_nonactivegroup.py