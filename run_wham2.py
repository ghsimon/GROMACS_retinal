'''
This script performs the Weighted Histogram Analysis Method (WHAM) on a set of simulation data.

The script reads in biasing potential data from a series of files named 'ALLCOLVAR' followed by a number. 
These files are expected to be in the current working directory and contain data from the concatenated trajectory 
reweighted under the different biasing potentials, with the number indicating the position of the umbrella/bias.

The biasing potential data is then used to perform the WHAM calculation using the `wham` function from the `wham` module. 
The result of the WHAM calculation is added to a pandas DataFrame and written to a file named "bias_multi2.dat".

This script is based on the tutorial at https://www.plumed.org/doc-v2.7/user-doc/html/masterclass-21-3.html.

Usage:
    python run_wham2.py
'''

# Package imports
import numpy as np
import plumed
import wham
import os

# Simulation constants
kB      = 0.008314462618    # Boltzmann constant in kJ/(mol*K)
T       = 300               # Temperature in K
kBT     = kB*T

# Check number of ALLCOLVAR files in the current directory
n=0
for fname in os.listdir(os.getcwd()):
    if ('ALLCOLVAR' in fname and 'bck' not in fname):
        n+=1
print('Number of files: %i'%n)

# Load the biasing potential data from the ALLCOLVAR files
i=0
print('Load col and bias')
col0        = plumed.read_as_pandas('ALLCOLVAR-3.10') # Read the first file
bias        = np.zeros((len(col0["restraint-phi.bias"]),n)) #each column of bias should contain bias for full concatenate trajectory according to one umbrella
bias[:,i]   = col0["restraint-phi.bias"][-len(bias):] # Add the biasing potential data from the first file to the array
print('ALLCOLVAR-3.10 (number %i)'%(i))
i+=1
# Repeat for the remaining files
for j in np.arange(-3.05,-0.01,0.05):
    bias[:,i]=plumed.read_as_pandas('ALLCOLVAR%.2f'%j)["restraint-phi.bias"][-len(bias):]
    print('ALLCOLVAR%.2f (number %i)'%(j,i))
    i+=1
for j in np.arange(0.0,3.14,0.05):
    bias[:,i]=plumed.read_as_pandas('ALLCOLVAR%.2f'%j)["restraint-phi.bias"][-len(bias):]
    print('ALLCOLVAR%.2f (number %i)'%(j,i))
    i+=1
print('Finished')

# Check that the number of files matches the number of biasing potential data arrays
assert n==i

# Perform the WHAM calculation
print('Do WHAM')
w = wham.wham(bias,T=kBT,verbose=True)
print('Finished')

# Add the WHAM results to the DataFrame and write it to a file
print('make colvar')
colvar = col0
colvar["logweights"]=w["logW"]
print('Write')
plumed.write_pandas(colvar,"bias_multi2.dat")
