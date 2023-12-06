
# Script to calculate per residue RMSF in refined MD trajectory

### RMSF

The formula of the RMSF is given as below:

$$ \rho \text{RMSF}_i = \sqrt{\frac{1}{N} \sum_{j=1}^{N} (r\_{i,j} - \langle r_i \rangle)^2} $$

RMSF.py plots per residue RMSF for all the non hydrogen atoms using pymol. 
