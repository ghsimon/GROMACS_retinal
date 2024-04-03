# This script is used to clean up the directory structure for umbrella sampling simulations.
# It creates necessary directories, moves files to their respective directories, and removes temporary directories.

# Create a temporary directory
mkdir tmp

# Move all files starting with "eq" and ending with ".mdp" to the temporary directory
mv eq*.mdp tmp

# Create a directory for umbrella sampling
mkdir us

# Move all files starting with "us" to the umbrella sampling directory
mv us* us

# Create a directory for eq1
mkdir eq1

# Move all files starting with "eq1" to the eq1 directory
mv eq1* eq1

# Create a directory for eq2
mkdir eq2

# Move all files starting with "eq2" to the eq2 directory
mv eq2* eq2

# Create a directory for COL
mkdir COL

# Move all files starting with "COLV" to the COL directory
mv COLV* COL

# Move all files from the temporary directory to the current directory
mv tmp/* .

# Remove the temporary directory
rm -r tmp
