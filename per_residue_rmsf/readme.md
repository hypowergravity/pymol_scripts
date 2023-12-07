
# Script to calculate per residue RMSF in refined MD trajectory

### RMSF

The formula of the RMSF is given as below:


```math
\rho \text{RMSF}_i = \sqrt{\frac{1}{N} \sum_{j=1}^{N} (r_{i,j} - \langle r_i \rangle)^2}
```

Where 
$\rho \text{RMSF}_i$ is per residue RMSD.

and $r_{i,j}$ is atom of interest, $\langle r_i \rangle$ is mean of co-ordinates of the specified atom.

RMSF.py plots per residue RMSF for all the non hydrogen atoms using pymol. 
