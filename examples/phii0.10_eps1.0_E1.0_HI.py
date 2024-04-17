from hoomd import *
import hoomd
import hoomd.md
import hoomd.Dielectric
import hoomd.PSEv1

import numpy as np
from datetime import datetime

# Quantities to be specified.
dt = 1e-3  # time step
N = 50000  # number of particles
phi = 0.10  # volume fraction
strength = 1.  # strength of (valence 1) Coulombic interaction at contact
field = 1.0  # field strength
gradient = 0.  # field gradient strength
T = 1  # temperature
lambda_p = 0. # particle conductivity
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
fileprefix = 'N{}_phii{:.2f}_eps{:.1f}_E{:.2f}'.format(N, phi, strength, E_0)

# Adjust parameters according to the temperature.  This ensures that all energies are scaled by kT,
# even if kT is not zero in the simulation.
strength = strength*T
field = field*np.sqrt(T)
gradient = gradient*np.sqrt(T)
t_rand = t_rand/T
t_eq = t_eq/T
t_run = t_run/T

# Enforce electroneutrality
Nsalt = int(np.round(N/(z_p+z_m)))
N = int(np.round(Nsalt*(z_p+z_m)))

# Construct arrays for the particle charges and types
q = np.sqrt(8.*np.pi*strength)  # valence 1 charge
charge = ([z_p*q]*z_m + [-z_m*q]*z_p)*Nsalt
type = (['pos']*z_m + ['neg']*z_p)*Nsalt

# Randomize the particle order
ind = np.random.permutation(N)
charge = np.array(charge)[ind]
type = np.array(type)[ind]

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
snapshot = data.make_snapshot(N=N, box=data.boxdim(L=L), particle_types=['pos','neg'])

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

# In the Heyes-Melrose algorithm, hard sphere interactions takes on the form of this potential when HI are included
def hs_potential_HI(r, rmin, rmax, dt, a):
    V = 8.*a/(3.*dt)*(2.*a*np.log(2.*a/r) + (r - 2.*a))
    F = 8.*a/(3.*dt*r)*(2.*a - r)
    return (V, F)

# Initialize the neighbor list
nl = hoomd.md.nlist.cell()

# Create an interpolation table for the HI hard sphere interactions
table_hs_HI = hoomd.md.pair.table(width=1000, nlist=nl)
table_hs_HI.pair_coeff.set('pos','pos',func=hs_potential_HI,rmin=0.001,rmax=diameter,coeff=dict(dt=dt,a=radius))
table_hs_HI.pair_coeff.set('pos','neg',func=hs_potential_HI,rmin=0.001,rmax=diameter,coeff=dict(dt=dt,a=radius))
table_hs_HI.pair_coeff.set('neg','neg',func=hs_potential_HI,rmin=0.001,rmax=diameter,coeff=dict(dt=dt,a=radius))

# Establish HI integrator.
all = group.all()
hoomd.md.integrate.mode_standard(dt=dt)
S = hoomd.PSEv1.integrate.PSEv1(group=all, seed=datetime.now().microsecond, T=T, xi=xi, error=error)

# Thermalize with only steric repulsions
run(N_rand)

# Turn on the electrostatic interactions.
# dipoleflag = 0: default settings; solves for the induced dipoles
# dipoleflag = 1: constant dipole; sets each dipole to the isolated particle dipole
# dipoleflag = 2: charge only; ignores dipolar interactions
dielectric = hoomd.Dielectric.compute.Dielectric(group=all, conductivity=[lambda_p]*N, field=[0.0, 0.0, field], gradient=[gradient, 0.0, 0.0],
                                                   xi=xi, errortol=error, fileprefix=fileprefix, period=N_txtperiod, dipoleflag=2)

# Equilibrate
run(N_eq)

# Set up sampling
hoomd.dump.gsd(filename=fileprefix+'.gsd', period=N_imageperiod, group=all, overwrite=True)

# Run the simulation
run(N_run+1)
