# The script iterates over the range of values 'AT' for the dihedral angle 'phi' at which US trajectries have previously been generated using the mult*.sh scripts in the mult directory.
# For each value of 'AT', it generates a PLUMED input file 'plumed.dat' with the appropriate parameters.
# Plumed is then used to calculate weights for the full concatenated trajectory at the given 'AT' value.
# The weights can be used in the weighted histogram analysis method (WHAM) to calculate the free energy profile along the dihedral angle 'phi'.

# Loop over a sequence of values for the dihedral angle from 0.0 radians to 3.1 radians with a step of 0.05
for AT in $(seq 0.0 0.05 3.1)
do

# Generate the PLUMED input file 'plumed2.dat' with the appropriate parameters
cat >plumed2.dat << EOF
phi: TORSION ATOMS=34,29,27,25
restraint-phi: RESTRAINT ARG=phi KAPPA=750.0 AT=$AT
PRINT STRIDE=100 ARG=phi,restraint-phi.bias FILE=ALLCOLVAR$AT
EOF

# Run the simulation using the generated input file
plumed driver --mf_xtc alltraj.xtc --trajectory-stride=100 --plumed plumed2.dat

done
