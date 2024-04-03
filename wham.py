'''
Script downloaded at https://github.com/plumed/masterclass-21-3.git
'''

import numpy as np
import sys

def wham(bias,
        *,
        frame_weight=None,
        traj_weight=None,
        T: float = 1.0,
        maxiter: int = 5000,
        threshold: float = 1e-8,
        verbose: bool = False):
    """
    Implements the Weighted Histogram Analysis Method (WHAM) for combining information from multiple simulations.

    Parameters:
    bias (numpy.ndarray): A 2D array representing the biasing potential from each simulation.
    frame_weight (numpy.ndarray, optional): Weights for each frame. Defaults to None, in which case all frames are equally weighted.
    traj_weight (numpy.ndarray, optional): Weights for each trajectory. Defaults to None, in which case all trajectories are equally weighted.
    T (float, optional): Temperature. Defaults to 1.0.
    maxiter (int, optional): Maximum number of iterations for the WHAM algorithm. Defaults to 5000.
    threshold (float, optional): Convergence threshold for the WHAM algorithm. Defaults to 1e-8.
    verbose (bool, optional): If True, print progress messages to stderr. Defaults to False.

    Returns:
    dict: A dictionary containing the log weights ('logW'), log partition functions ('logZ'), number of iterations ('nit'), and final change in partition functions ('eps').
    """
    
    # The number of frames and trajectories are determined from the shape of the 'bias' array.
    nframes = bias.shape[0]
    ntraj = bias.shape[1]

    # If no weights are provided, they are set to 1 for all frames and trajectories.
    if frame_weight is None:
        frame_weight = np.ones(nframes)
    if traj_weight is None:
        traj_weight = np.ones(ntraj)

    assert len(traj_weight) == ntraj
    assert len(frame_weight) == nframes

    # divide by T once for all
    shifted_bias = bias/T
    # track shifts
    shifts0 = np.min(shifted_bias, axis=0)
    shifted_bias -= shifts0[np.newaxis,:]
    shifts1 = np.min(shifted_bias, axis=1)
    shifted_bias -= shifts1[:,np.newaxis]

    # do exponentials only once
    expv = np.exp(-shifted_bias)

    # Initialize the partition functions Z for each trajectory to 1
    Z = np.ones(ntraj)

    # Store the old partition functions for convergence checking
    Zold = Z

    # If verbose is True, print a start message to stderr
    if verbose:
        sys.stderr.write("WHAM: start\n")
    
    # Main loop for the WHAM iterations
    for nit in range(maxiter):
        # Calculate the unnormalized weights for each frame
        weight = 1.0/np.matmul(expv, traj_weight/Z)*frame_weight

        # Update the partition functions Z
        Z = np.matmul(weight, expv)

        # Normalize the partition functions
        Z /= np.sum(Z*traj_weight)

        # Calculate the change in the partition functions for convergence checking
        eps = np.sum(np.log(Z/Zold)**2)

        # Store the current partition functions for the next iteration
        Zold = Z

        # If verbose is True and the iteration number is a multiple of 10, print a progress message to stderr
        if (verbose and nit%10==0):
            sys.stderr.write("WHAM: iteration "+str(nit)+" eps "+str(eps)+"\n")
        
        # If the change in the partition functions is below the threshold, break the loop
        if eps < threshold:
            break
    
    # Store the number of function evaluations
    nfev=nit

    # Calculate the log weights for each frame
    logW = np.log(weight) + shifts1

    # If verbose is True, print an end message to stderr
    if verbose:
        sys.stderr.write("WHAM: end")

    # Return a dictionary with the log weights, log partition functions, number of iterations, and final change in partition functions
    return {"logW":logW, "logZ":np.log(Z)-shifts0, "nit":nit, "eps":eps}
