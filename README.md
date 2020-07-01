# periodic-Coulomb-Tree-Method 
This is the accompanying software to the paper, Periodic Coulomb Tree Method: An Alternative To Parallel Particle Mesh Ewald (J. Chem. Theory Comput. 2020, 16, 1, 7-17). The specific code for the periodic Coulomb Tree (PCT) method is in the file tree_module.f which is incorporated into DL_POLY Classic
version 1.9. This release contains the DL_POLY Classic files as well. The supporting information for the paper provides details of our implementation in 
DL_POLY Classic. If you are new to DL_POLY Classic, the user manual provides a useful startup guide. 

To use PCT instead of other methods for electrostatics in periodic boundary conditions, include the following directives in the CONTROL file:

ltree           - causes DL_POLY to compute electrostatics with PCT

taylor/tricubic - For the particle-cluster interaction, PCT will use a Taylor polynomial approximation if the key 'taylor' is included in the CONTROL 
                  file and tricubic approximation if the key 'tricubic' is included. Be sure not to include both but you must include at least one  
                  
mactheta    x   - Sets the criterion for approximating a particle-cluster interaction via a Taylor polynomial/tricubic polynomial approximation. 
                  x should be between 0(zero) and 1(one). The method is more accurate for smaller y values but also more computationally expensive.
                  
maxparnode  y   - Sets the maximum number of particles in a leaf of a tree. A small y results in a deeper tree which usually results in an increase in 
                  compute time
                  
torder      z   - Sets the order of the Taylor polynomial approximation. This is required if you choose the 'key' taylor. z is 0 or a positive integer

lattlim     w   - Sets the range of periodic images. w is 0 or a positive interger. w = 0 means no periodic images. w = 1 means include nearest 
                  neihgbors of the fundmental box as periodic images
                 
                  

                  
