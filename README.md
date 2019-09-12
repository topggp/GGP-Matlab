# GGP-Matlab
GGP Matlab code
This is supposed to be supplementary material for Generalized Geometry Projection Framework [1].
The proposed framework make the use of Method of Moving Asymptotes optimization solvers [2].
Moreover for Gauss quadrature we use the code in [3].
In this folder you will find:
1. **GGP.mlx**: Matlab live script that contains the main code and comments to run Genaralized Geometry Ptojection. At the first look the script may seems long. In reality many parts are copied by top88 matlab code [4]. Moreover graphics also take many lines but are not essential for the understanding of the code.
2. **GGP.m**: Matlab script that contains the main code and comments to run Genaralized Geometry Ptojection. Could be usefull for Matlab users with versions prior to R2016a.
3. **Wgp.m**: Matlab function that take as inputs the projection parameters provided in the main, the coordinates of sampling window Gauss points and reads in output smooth characteristic function and derivatives computed for each component and sampling window gauss point.
4. **lgwt.m**: Matlab function that provides Gauss point coordinates and weights for a 1D domain [3].
5. **model_updateM.m** and **model_updateV.m**: Matlab functions that takes as inputs local volume fractions and provide Finite Element Young Modulus, densities and their sensitivities. 
6. **Aggregation_Pi.m** function adopted to make the union of components ussing different methods. This also provide sensitivities with respect to entries.
7. **smooth_sat.m**: matlab function that applies smooth saturation to avoid aggregation related issues [1].
8. **mmasub.m**, **subsolv.m**, **kktcheck.m**: matlab functions for MMA optimization solver and KKT residual norm [3].
9. **norato_bar.m**: matlab function that provide the distance of a component boundary frome the component center in function of the polar coordinate angle. This is needed for component plots.

# References
[1] Coniglio, Simone et al. "Generalized Geometry Projection a unified approach for geometric feature based topology optimization."

[2] Svanberg, Krister. "MMA and GCMMA, versions September 2007." Optimization and Systems Theory 104 (2007).

[3] Von Winckel, Greg. "[Legendre-Gauss Quadrature Weights and Nodes](https://fr.mathworks.com/matlabcentral/fileexchange/4540-legendre-gauss-quadrature-weights-and-nodes)"  (2004).

[4] Andreassen, Erik, et al. "Efficient topology optimization in MATLAB using 88 lines of code." Structural and Multidisciplinary Optimization 43.1 (2011): 1-16.
