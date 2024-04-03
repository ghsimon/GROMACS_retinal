'''
This script calculates the diffusion coefficient profile from an umbrella sampling simulation.

The calculation is based on the method described in:
Gerhard Hummer 2005 New J. Phys.7 34
https://iopscience.iop.org/article/10.1088/1367-2630/7/1/034/pdf

The script reads in simulation data from a series of files named 'COLVAR' followed by a number indicating the position of the US harmonic biasing potentials. 
These files are expected to contain data from different simulations that have been performed under different biasing potentials.

The script uses a class called 'diffusion' to calculate the diffusion coefficient for each simulation. 
The 'diffusion' class takes as input the number of timesteps per simulation, the timestep, and the number of timesteps used by PLUMED to print out the data.

Usage:
    python calculate_diffusion.py
'''

# Import packages
import numpy as np
import scipy.constants as constants

###########################
###     CONSTANTS       ###
###########################

Nsims   = 125 # number of umbrellas
Nsteps  = 2000001
dt      = 0.002 #ps
stride  = 1
T       = 300 #K
kT      = T*constants.k*constants.N_A/1000 #kJ/mol
beta    = 1/kT
eqsteps = 0
xj      = np.arange(-3.1,3.11,0.05)

assert xj.size == Nsims # check if the defined positions of the umbrellas correspond to the number of umbrellas

###########################
###   DIFFUSION CLASS   ###
###########################

class diffusion():
    '''
    This class calculates diffusion coefficient profile from an umbrella sampling simulation.

    The calculation is based on the method described in:
    Gerhard Hummer 2005 New J. Phys.7 34
    https://iopscience.iop.org/article/10.1088/1367-2630/7/1/034/pdf

    Attributes:
    Nsteps (int): The number of timesteps per simulation.
    dt (float): The timestep in picoseconds.
    stride (int): The number of timesteps used by PLUMED to print out the data.

    Methods:
    compute(x): Calculates the diffusion coefficient based on the input array x.
                x is expected to be an array of positions from an umbrella sampling simulation.
    '''
    def __init__(self,
                       Nsteps,
                       dt,                                        # ps
                       stride,                                    # timesteps
                 ):
        '''
        Constructs the diffusion object with the number of timesteps per simulation (Nsteps), 
        the timestep (dt), and the number of timesteps used by PLUMED to print out the data (stride).
        '''
        self.Nsteps = Nsteps
        self.dt     = dt
        self.stride = stride
    def compute(self, x):
        '''
        Calculates the diffusion coefficient based on the input array x.
        x is expected to be an array of positions from an umbrella sampling simulation.
        The method returns the diffusion coefficient.
        '''
        varx = np.var(x, ddof = 1)
        avgx = np.mean(x)
        #implementation with np.correlate:
        C = np.correlate(x - avgx, x - avgx, mode='full')
        C = C[self.Nsteps - 1 :: -1]  / C[self.Nsteps - 1]
        print(f'C[0] = {C[0]}')
        print(f'varx = {varx}')
        h = np.nonzero(C < 0)[0][0]
        integral = np.trapz(C[:h] * self.stride * self.dt)
        print(f'integral = {integral}')
        return varx / integral

########################################
###   IMPORT DIFFUSION FROM CYTHON   ###
########################################

# To build the Cython diffusion module:
# python3 buildCythonModule.py build_ext --inplace
from diffusionC import diffusionC

#################################
###   DIFFUSION CALCULATION   ###
#################################

# Initialize arrays to store the diffusion coefficients and average positions
diff_arr = np.zeros(Nsims)
avg_arr = np.zeros(Nsims)

# Loop over the number of simulation trajectories
for rep in range(Nsims):
    # Construct the directory name based on the simulation number
    directory = "../COL/COLVAR" + "{:.2f}".format(xj[rep])
    print(directory)

    # Load the data from the COLVAR file, skipping the first row and loading only the second column
    # The second column contains the relevant coordinate data
    with open(directory) as f:
        lines = (line for line in f if not line.startswith('#'))
        r     = np.loadtxt(lines, usecols = (1), unpack = True)[eqsteps:]
    
    # Unwrap the coordinate data to remove discontinuities due to periodic boundary conditions
    x = np.unwrap(r)
    
    # Calculate the average position, adjusting for periodic boundary conditions
    avgx = np.mean(x)
    if avgx > np.pi:
        avgx -= 2*np.pi
    elif avgx < -np.pi:
        avgx += 2*np.pi
    
    # Calculate the diffusion coefficient
    #diff_arr[rep]   = diffusion(Nsteps, dt, stride).compute(x)  # using diffusion class employing np.correlate
    diff_arr[rep]   = diffusionC(Nsteps, dt, stride).compute(x)  # using diffusionC from Cython class
    
    # Store the average position
    avg_arr[rep]    = avgx

    print(f'diffusion: {diff_arr[rep]}')
    
np.save('diff_arr.npy',diff_arr)
np.save('avg_arr.npy',avg_arr)
