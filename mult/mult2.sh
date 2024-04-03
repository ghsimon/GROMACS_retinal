# This script performs umbrella sampling simulations for multiple values of a torsional angle.
# It creates directories for storing the output files and performs equilibration and production runs for each angle.
# The script uses GROMACS and PLUMED for the simulations.

# The main steps in the script are as follows for each value of the dihedral angle:
# 1. Perform equilibration with weaker umbrella restraints
# 2. Perform equilibration with actual umbrella restraints
# 3. Run the production trajectory with actual umbrella restraints

# The script uses a for loop to iterate over a range of values for the dihedral angle.
# For each iteration, it performs the equilibration and production runs with the corresponding angle.

# Note: The script assumes that the necessary input files (mdp, gro, top, cpt) are already present in the working directory.

# Example usage: bash mult2.sh

# Loop over a sequence of values for the dihedral angle from 0.05 radians to 1.5 radians with a step of 0.05
for AT in $(seq 0.05 0.05 1.5)
do

# Create a directory named EQCOLVAR if it doesn't exist
mkdir -p EQCOLVAR

### EQUILIBRATE WITH WEAKER UMBRELLA RESTRAINTS
gmx grompp -f eq1_nvt.mdp -c em.gro -p topol.top -r 6eid_capped_matched.gro -o eq1_nvt$AT.tpr

# Create a PLUMED input file with a restraint on a torsion angle phi with a force constant of 100 kJ/mol/rad and a target value of AT
cat >eq1_plumed2.dat << EOF
phi: TORSION ATOMS=34,29,27,25
restraint-phi: RESTRAINT ARG=phi KAPPA=100.0 AT=$AT
PRINT STRIDE=100 ARG=phi,restraint-phi.bias FILE=EQCOLVAR/EQ1COLVAR$AT
EOF

# Run the simulation with GROMACS mdrun
gmx_mpi mdrun -plumed eq1_plumed2.dat -deffnm eq1_nvt$AT -x eq1_traj$AT.xtc

### EQUILIBRATE WITH ACTUAL UMBRELLA RESTRAINTS
gmx grompp -f eq2_nvt.mdp -c eq1_nvt$AT.gro -p topol.top -r 6eid_capped_matched.gro -o eq2_nvt$AT.tpr

# Create a PLUMED input file with a restraint on a torsion angle phi with a force constant of 750 kJ/mol/rad and a target value of AT
cat >eq2_plumed2.dat << EOF
phi: TORSION ATOMS=34,29,27,25
restraint-phi: RESTRAINT ARG=phi KAPPA=750.0 AT=$AT
PRINT STRIDE=100 ARG=phi,restraint-phi.bias FILE=EQCOLVAR/EQ2COLVAR$AT
EOF

# Run the simulation with GROMACS mdrun
gmx_mpi mdrun -plumed eq2_plumed2.dat -deffnm eq2_nvt$AT -x eq2_traj$AT.xtc

### RUN TRAJECTORY WITH ACTUAL UMBRELLA RESTRAINTS
gmx grompp -f md.mdp -c eq2_nvt$AT.gro -t eq2_nvt$AT.cpt -p topol.top -r 6eid_capped_matched.gro -o us$AT.tpr

# Create a PLUMED input file with a restraint on a torsion angle phi with a force constant of 750 kJ/mol/rad and a target value of AT
cat >plumed2.dat << EOF
phi: TORSION ATOMS=34,29,27,25
restraint-phi: RESTRAINT ARG=phi KAPPA=750.0 AT=$AT
PRINT STRIDE=1 ARG=phi,restraint-phi.bias FILE=COLVAR$AT
EOF

# Run the simulation with GROMACS mdrun
mpirun -np 1 gmx_mpi mdrun -plumed plumed2.dat -deffnm us$AT -x traj$AT.xtc

done
