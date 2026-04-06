#!/bin/bash

#SBATCH --mail-type=END,FAIL,REQUEUE
#SBATCH --mail-user=hz322@uw.edu

#SBATCH --job-name=FD_dielectric
#SBATCH --account=zeelab
#SBATCH --partition=gpu-l40
#SBATCH --gpus=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=10:00:00
#SBATCH --mem=30gb
#SBATCH -o %j.out

# Inputs
PHI=(0.01)
GRAD=(0.010 0.020 0.050 0.100)

# Load necessary Hyak modules
module load cuda/11.8.0
module load ssmc/miniconda/3.9
module load cmake/3.20.0
module load gcc/10.2.0

# Load hoomd by activating the python virtual environment in which HOOMD was built.
# Different file path needed for different users.
source /gscratch/zeelab/haoqing/HOOMD_pluginfix/hoomd-venv/bin/activate

# Add the amgx libraries to the file path. Should consider doing this automatically during install?
export LD_LIBRARY_PATH=/gscratch/zeelab/haoqing/HOOMD_pluginfix/AMGX/build:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/gscratch/zeelab/haoqing/HOOMD_pluginfix/Dielectric/build:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/mmfs1/sw/contrib/ssmc-src/miniconda/3.9/lib:$LD_LIBRARY_PATH
#export LD_PRELOAD=/mmfs1/sw/contrib/ssmc-src/miniconda/3.9/lib/libmkl_core.so:/mmfs1/sw/contrib/ssmc-src/miniconda/3.9/lib/libmkl_sequential.so
#export LD_PRELOAD=/mmfs1/sw/contrib/ssmc-src/miniconda/3.9/lib/libmkl_core.so:/mmfs1/sw/contrib/ssmc-src/miniconda/3.9/lib/libmkl_sequential.so:/mmfs1/sw/contrib/ssmc-src/miniconda/3.9/lib/libmkl_def.so


# so_plugin=/mmfs1/gscratch/zeelab/haoqing/HOOMD_pluginfix/hoomd-venv/lib/python3.9/site-packages/hoomd/Dielectric/_Dielectric.cpython-39-x86_64-linux-gnu.so

# echo "== NEEDED deps in the plugin =="
# objdump -p "$so_plugin" | grep NEEDED

# echo "== ldd on the plugin =="
# ldd "$so_plugin"



# Run simulations at different params
for phi in "${PHI[@]}"; do
    for G in "${GRAD[@]}"; do
        python Dielectric_FD_grouptag_nonactivegroup.py $phi $G
    done
done

