from hoomd import *
import hoomd
import hoomd.md
import hoomd.Dielectric

import numpy as np
from datetime import datetime

# Quantities to be specified.
dt = 1e-3  # time step
N = 6000  # number of particles
phi = 0.10 # volume fraction
strength = 1.  # strength of (valence 1) Coulombic interaction at contact
field = 1.0  # field strength
gradient = 0.01  # field gradient strength
T = 1  # temperature
z_p = 1 # valence of positive charges
z_m = 1 # valence of negative charges
t_rand = 10 # randomization time
t_eq = 10  # equilibration time
t_run = 100  # run time
N_image = 100  # number of output snapshots in gsd file
N_txt = 10 # number of output txt files 
error = 1e-3  # desired error tolerance
xi = 0.5  # Ewald splitting parameter

# Construct the output file name
fileprefix = 'N{}_phii{:.2f}_G{:.2f}_E{:.2f}_FD'.format(N, phi, gradient, field)

# Adjust parameters according to the temperature.  This ensures that all energies are scaled by kT,
# even if kT is not zero in the simulation.
strength = strength*T
field = field*np.sqrt(T)
gradient = gradient*np.sqrt(T)
t_rand = t_rand/T
t_eq = t_eq/T
t_run = t_run/T

# Enforce electroneutrality
Nnonactive = 2000
Nsalt = int(np.round((N-Nnonactive)/(z_p+z_m)))

# Construct arrays for the particle charges and types
q = np.sqrt(8.*np.pi*strength)  # valence 1 charge
charge = ([z_p*q]*z_m + [-z_m*q]*z_p)*Nsalt + [0]*Nnonactive
type = (['pos']*z_m + ['neg']*z_p)*Nsalt + ['nonactive']*Nnonactive
lambda_p = ([0]*z_m + [4]*z_p)*Nsalt + [-1]*Nnonactive # particle conductivity; for non-active particles, set lambda_p = -1

# Randomize particle order
ind = np.random.permutation(N)
type = np.array(type)[ind].tolist()
charge = np.array(charge)[ind].tolist()
lambda_p = np.array(lambda_p)[ind].tolist()

# Typical nondimensionalization
mass = dt  # particle mass; this ensures overdamped dynamics if a two step integrator is used
diameter = 2.  # particle diameter
radius = 1.  # particle radius; sets the length scale to be the particle radius
gamma = 1.  # particle drag coefficient; sets the time scale to be the diffusion time

# Other quantities
L = (4.*np.pi*N/(3.*phi))**(1./3.)  # box dimension

# Compute the number of time steps
N_rand = int(np.round(t_rand/dt))
N_eq = int(np.round(t_eq/dt))
N_run = int(np.round(t_run/dt))
N_imageperiod = int(np.round(N_run/N_image)) 
N_txtperiod = int(np.round(N_run/N_txt)) 

# Parameters for initializing system on a simple cubic lattice
m = int(np.ceil(N**(1./3.)))  # smallest latticle dimension that can hold all of the particles
a = L/m  # distance between lattice sites
lo = -L/2.  # used to center the lattice about the origin

# Create a snapshot of an empty system
context.initialize()
snapshot = data.make_snapshot(N=N, box=data.boxdim(L=L), particle_types=['pos','neg','nonactive'])

# Initialize the system from the snapshot
system = init.read_snapshot(snapshot)

# Initialize a simple cubic array of particles
for p in system.particles:
    (i,j,k) = (p.tag % m, p.tag/m % m, p.tag/m**2  % m)
    p.position = (lo + i*a + a/2, lo + j*a + a/2, lo + k*a + a/2)
    p.mass = mass
    p.diameter = diameter
    p.type = type[p.tag]
    p.charge = charge[p.tag]

# In the Heyes-Melrose algorithm, hard sphere interactions takes on the form of a Hookean spring with a time step dependent spring constant.
def hs_potential_BD(r,rmin,rmax,k,req):
    V = k/2.*(r - req)**2.
    F = -k*(r - req) 
    return (V, F)

# Initialize the neighbor list
nl = hoomd.md.nlist.cell()

# Create an interpolation table for the BD hard sphere interactions
table_hs_BD = hoomd.md.pair.table(width=1000, nlist=nl)
table_hs_BD.pair_coeff.set('pos','pos',func=hs_potential_BD,rmin=0.,rmax=diameter,coeff=dict(k=1./(2.*dt),req=diameter))
table_hs_BD.pair_coeff.set('pos','neg',func=hs_potential_BD,rmin=0.,rmax=diameter,coeff=dict(k=1./(2.*dt),req=diameter))
table_hs_BD.pair_coeff.set('neg','neg',func=hs_potential_BD,rmin=0.,rmax=diameter,coeff=dict(k=1./(2.*dt),req=diameter))
table_hs_BD.pair_coeff.set('pos','nonactive',func=hs_potential_BD,rmin=0.,rmax=diameter,coeff=dict(k=0,req=diameter))
table_hs_BD.pair_coeff.set('neg','nonactive',func=hs_potential_BD,rmin=0.,rmax=diameter,coeff=dict(k=0,req=diameter))
table_hs_BD.pair_coeff.set('nonactive','nonactive',func=hs_potential_BD,rmin=0.,rmax=diameter,coeff=dict(k=0,req=diameter))

# Establish Brownian dynamics integrator.
all = group.all()
hoomd.md.integrate.mode_standard(dt=dt)
bd = hoomd.md.integrate.brownian(group=all, kT=T, seed=datetime.now().microsecond)
bd.set_gamma('pos', gamma=gamma)
bd.set_gamma('neg', gamma=gamma)
bd.set_gamma('nonactive', gamma=gamma)

# Thermalize with only steric repulsions
run(N_rand)

# Turn on the electrostatic interactions.
# dipoleflag = 0: default settings; solves for the induced dipoles
# dipoleflag = 1: constant dipole; sets each dipole to the isolated particle dipole
# dipoleflag = 2: charge only; ignores dipolar interactions

# Plugin Input: (group, group_tag, conductivity, field = [0., 0., 0.], gradient = [0., 0., 0.], xi = 0.5, errortol = 1e-3, fileprefix = "", period = 0, dipoleflag = 0)
# group_tag is an array size of Ntotal, should match with tag (Ntotal) and conductivity (group_size) arrays. 
# For a non-active particleL: group_tag = -1; not included in the conductivity array

# The logic behind:
    # unsigned int idx = d_group_members[group_idx] -> looping through group_idx (make sure all particles are in the active group), getting the hoomd idx from group_idx 
    # unsigned int tag = d_tag[idx] -> getting particle tag from idx
    # int group_tag = d_group_tag[tag] -> getting group_tag from tag

    # For a specific particle in the active group, we are getting its charge and conductivity from: 
    # charge[idx], dipole[group_tag] -> getting particle charge from charge and dipole array (user input)

# Example: 
    # For a system with in total 8 particles (Ntotal=8)
    # Active group: groupA (group_size=6); for A1, lambda=0; for A2, lambda=4
    # Non-active group: groupB

    # tag array (size Ntotal=8): [0, 1, 2, 3, 4, 5, 6, 7]
    # type array (size Ntotal=8): [A1, A1, A2, B, A2, B, A1, A2]
    # group_tag array (size Ntotal=8): [0, 1, 2, -1, 3, -1, 4, 5] (Set group_tag=-1 for non-active particles; for active particles, group_tag=tag)
    # conductivity arrray (active group size group_size=6): [0, 0, 4, 4, 0, 4]

    # The particle order of tag, group_tag, conductivity should align with each other.

groupPOS = group.type(name='pos', type='pos')
groupNEG = group.type(name='neg', type='neg')
groupACTIVE = group.union(name="active", a=groupPOS, b=groupNEG)
N_active = len(groupACTIVE)

# conductivity for active particles only: remove all -1, keep order
lambda_active = [x for x in lambda_p if x != -1]

# group_tag for all particles: active particles get 0,1,2,... in order; non-active get -1
group_tag_array = []
active_id = 0
for x in lambda_p:
    if x == -1:
        group_tag_array.append(-1)
    else:
        group_tag_array.append(active_id)
        active_id += 1

dielectric = hoomd.Dielectric.compute.Dielectric(
    group=groupACTIVE,
    group_tag=group_tag_array,
    conductivity=lambda_active,
    field=[0.0, 0.0, field],
    gradient=[gradient, 0.0, 0.0],
    xi=xi,
    errortol=error,
    fileprefix=fileprefix,
    period=N_imageperiod,
    dipoleflag=0
)

# Equilibrate
run(N_eq)

# Set up sampling
hoomd.dump.gsd(filename=fileprefix+'.gsd', period=N_imageperiod, group=all, overwrite=True)

# Run the simulation
run(N_run+1)
