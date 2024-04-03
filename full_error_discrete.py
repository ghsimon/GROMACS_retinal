'''
This script performs error estimation based on bootstrapping and the Weighted Histogram Analysis Method (WHAM) on a set of simulation data.

The script reads in biasing potential data from a series of files named 'ALLCOLVAR' followed by a number. 
These files are expected to be in the current working directory and contain data from the concatenated trajectory 
reweighted under the different biasing potentials, with the number indicating the position of the umbrella/bias.

The biasing potential data is then used to perform the WHAM calculation using the `wham` function from the `wham` module. 
Bootstrapping is performed by randomly selecting blocks of data from the biasing potential data and performing WHAM on these bootstrapped data sets.
Errors are calculated based on the bootstrapped data sets.

This script is based on the tutorial at https://www.plumed.org/doc-v2.7/user-doc/html/masterclass-21-3.html
or alternatively https://github.com/plumed/masterclass-21-3/blob/main/INSTRUCTIONS.md

Usage:
    python full_error_discrete.py
'''

# Package imports
import numpy as np
import plumed
import wham
import os
import subprocess
import copy
import scipy.constants as constants
import pandas as pd

# Simulation constants
kB              = constants.k
T               = 300
kBT             = kB*T*constants.N_A/1000
traj_len        = 60000 #length of each umbrella, generally 6e6/stride with 6e6 number of integration steps in mdp file

# Bootstrap and WHAM constants
NB              = 20
cycles          = 200
thresholdstring = "1e-8"
threshold       = float(thresholdstring)
filename        = "bootstrap_c%i_N%i_t%s_discrete.txt"%(cycles,NB,thresholdstring)

# Check number of ALLCOLVAR files in the current directory
n=0
for fname in os.listdir('ALLCOLVAR'):
    if ('ALLCOLVAR' in fname and 'bck' not in fname):
        n+=1
print('Number of files: %i'%n)

# Load the biasing potential data from the ALLCOLVAR files
i=0
print('Load col and bias')
col0        = plumed.read_as_pandas('ALLCOLVAR/ALLCOLVAR-3.10') # Read the first file
bias        = np.zeros((len(col0["restraint-phi.bias"]),n)) #each column of bias should contain bias for full concatenate trajectory according to one umbrella
bias[:,i]   = col0["restraint-phi.bias"][-len(bias):]  # Add the biasing potential data from the first file to the array
print('ALLCOLVAR-3.10 (number %i)'%(i))

# Repeat for the remaining files
i += 1
for j in np.arange(-3.05,-0.01,0.05):
    bias[:,i]=plumed.read_as_pandas('ALLCOLVAR/ALLCOLVAR%.2f'%j)["restraint-phi.bias"][-len(bias):]
    print('ALLCOLVAR%.2f (number %i)'%(j,i))
    i+=1
for j in np.arange(0.0,3.11,0.05):
    bias[:,i]=plumed.read_as_pandas('ALLCOLVAR/ALLCOLVAR%.2f'%j)["restraint-phi.bias"][-len(bias):]
    print('ALLCOLVAR%.2f (number %i)'%(j,i))
    i+=1
print('Finished')

# Check that the number of files matches the number of biasing potential data arrays
assert n==i

# Initiate text file 
with open('bootstrap/' + filename,"w") as f:
    print(f"""
    {cycles} BOOTSTRAP CYCLES WITH {NB} BLOCKS

    ignored points: {bias.shape[0]-n*traj_len}
    """,file=f)

# WHAM without bootstrapping
print(f'ignored points: {bias.shape[0]-n*traj_len}')
bb          = bias[-n*traj_len:,:].reshape((n,-1,n)).reshape((n,NB,traj_len//NB,n)) #split the bias in blocks per umbrella
cc          = np.array(col0.phi)[-n*traj_len:].reshape((n,-1)).reshape((n,NB,traj_len//NB)) #do the same for the phi values

tr          = cc.flatten() # phi-values of trajectory
is_in_B     = np.int_(np.logical_and(tr>-1.6,tr<1.6)) # one if trajectory point is cis and zero if trans
is_in_A     = np.int_(np.logical_or(tr<=-1.6,tr>=1.6)) # one if trajectory point is trans and zero if cis

w0          = wham.wham(bb.reshape((-1,n)),T=kBT,verbose=True,threshold=threshold) # do WHAM conversion

PB          = np.average(is_in_B,weights=np.exp(w0["logW"])) # calculate cis population according to weights calculated using WHAM
PA          = np.average(is_in_A,weights=np.exp(w0["logW"])) # calculate trans population according to weights calculated using WHAM
dF          = -kBT*np.log(np.average(is_in_B,weights=np.exp(w0["logW"]))/np.average(is_in_A,weights=np.exp(w0["logW"]))) # calculate free energy difference cis/trans
print("population trans: ",PA)
print("population cis: ",PB)
print("Delta F(cis-trans): ",dF)

# WHAM with bootstrapping
popA        = [] # list of estimated population trans for each cycle
popB        = [] # list of estimated population cis for each cycle
deltaF      = [] # list of estimated free energy difference trans vs cis for each cycle
ff          = [] # list of calculated free energy surface for each cycle

for i in range(cycles):
    # we then analyze the bootstrapped trajectories
    
    print('Bootstrap iteration %i'%i)
    c       = np.random.choice(NB,NB)  # c is an array of length NB with numbers between zero and NB-1
    # Perform WHAM on the bootstrapped bias data
    w       = wham.wham(bb[:,c,:,:].reshape((-1,n)),T=kBT,verbose=True,threshold=threshold) #wham with blocks randomized for each umbrella
    
    # Get the phi-values of the trajectory corresponding to the randomized blocks
    tr      = cc[:,c,:].flatten() # phi-values of trajectory corresponding to randomized blocks
    # Determine which points are in the cis and trans states
    is_in_B = np.int_(np.logical_and(tr>-1.6,tr<1.6)) # cis
    is_in_A = np.int_(np.logical_or(tr<=-1.6,tr>=1.6)) # trans
    
     # Calculate the populations in the cis and trans states and the free energy difference
    popA.append(np.average(is_in_A,weights=np.exp(w["logW"]))) # calculate population according to weights calculated using WHAM with bootstrapping
    popB.append(np.average(is_in_B,weights=np.exp(w["logW"]))) # calculate population according to weights calculated using WHAM with bootstrapping
    deltaF.append(-kBT*np.log(
        np.average(is_in_B,weights=np.exp(w["logW"]))/np.average(is_in_A,weights=np.exp(w["logW"]))
    ))
    print('Creating temp_bias_multi_discrete.dat')
    
    #create new plumed pandas file with bootstrapped trajectory
    colvar = pd.DataFrame(np.array([np.array(col0['time'])[:len(tr)],tr,w["logW"]]).T,columns=['time','phi','logweights'])
    plumed.write_pandas(colvar,'temp_bias_multi_discrete.dat') # temporary bias file
    
    print('Run plumed driver')
    # Run the plumed driver with the bootstrapped data
    subprocess.run('plumed driver --noatoms --plumed error_plumed_multi_discrete.dat --kt 2.4943387854',shell=True)
    # Remove files from the previous iteration
    subprocess.run('rm -f bck*',shell=True) 
    
    print('Create temp_fes_phi_catr_discrete.dat') # temporary free energy surface file
    fes_phir=plumed.read_as_pandas("temp_fes_phi_catr_discrete.dat").replace([np.inf, -np.inf], np.nan).dropna() # read in temporary free energy surface file
    ff.append(np.array(fes_phir.ffphir)) # append temporary FES of this iteration to ff list
    
    with open('bootstrap/' + filename,"a") as f:
        print(f"""
        Bootstrap iteration {i}
        popA = {np.average(is_in_A,weights=np.exp(w["logW"]))}
        popB = {np.average(is_in_B,weights=np.exp(w["logW"]))}
        deltaF = {-kBT*np.log(np.average(is_in_B,weights=np.exp(w["logW"]))/np.average(is_in_A,weights=np.exp(w["logW"])))}

        """,file=f)

# save the data
np.save('bootstrap/popA_bootstrap_%i_%i_discrete.npy'%(cycles,NB),popA)
np.save('bootstrap/popB_bootstrap_%i_%i_discrete.npy'%(cycles,NB),popB)
np.save('bootstrap/deltaF_bootstrap_%i_%i_discrete.npy'%(cycles,NB),deltaF)

# calculate errors
ff = np.array(ff)
error_ff = np.std(ff,axis=0) # calculate errors on full FES
np.save('bootstrap/error_ff_bootstrap_%i_%i_discrete.npy'%(cycles,NB),error_ff) # save

# calculate errors on population and free energy difference
error_popA      = np.std(popA)
error_popB      = np.std(popB)
error_deltaF    = np.std(deltaF)

# print and save errors
print("error population A:",error_popA)
print("error population B:",error_popB)
print("error Delta F:",error_deltaF)
with open('bootstrap/' + filename,"a") as f:
    print(f"""
    #############################################################
    From trajectory:
    P_trans = {PA}
    P_cis = {PB}
    Delta F(cis-trans) = {dF}

    From bootstrapping with {NB} bins of size {traj_len//NB}:
    error pop A = {error_popA}
    error pop B = {error_popB}
    error Delta F = {error_deltaF}
    #############################################################
    """,file=f)
