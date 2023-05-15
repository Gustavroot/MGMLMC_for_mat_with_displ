# MGMLMC for a matrix with a displacement - particular case of Schwinger (1+1)-D

## Installation

Only MATLAB is needed to run these tests.

## Running

Run the script: main.m

## Matrices

The matrices used in these tests were created by first generating the Schwinger (1+1)-D gauge configurations via the code in https://github.com/Gustavroot/2p1D-Schwinger, and then the 
corresponding sparse matrices were constructed with the script in matrices/schwinggen.m.

## Next steps

 -- compute full spectrum of the finest-level D. Double-check with Jacob that this spectrum makes physical sense. NOTE : for 128^2 and m0=-0.1, I did beta=1.0, but I couldn't go up to e.g. 
beta=3.0 .. why is this? Compute and plot the 5000 smalest eigenvalues for the 256^2 lattice as well, and also double-check this spectrum with Jacob

 -- implement calling block power iteration on the MGMLMC difference-level operators

 -- implement calling Hutchinson on the MGMLMC difference-level operators

 -- compute variances of the operators involved in both Huchinson and MGMLMC, using 1000 sample size

 -- implement stopping criterium in Hutchinson based on the tolerance, and correspondingly adjust for MGMLMC
