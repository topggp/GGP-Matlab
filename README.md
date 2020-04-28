# GGP-Matlab code
This is supposed to be supplementary material for Generalized Geometry Projection Framework [1].
The proposed framework uses Method of Moving Asymptotes (MMA) optimization solvers [2]. The associated MMA codes can also be downloaded from [http://www.smoptit.se](http://www.smoptit.se). 
Moreover for Gauss quadrature we use the code in [[3](https://fr.mathworks.com/matlabcentral/fileexchange/4540-legendre-gauss-quadrature-weights-and-nodes)].
Three test cases are currently available in the provided code: the short cantilever beam, the MBB beam and the L-shape beam.
Other test cases can be easily implemented as in [4]. 
In this folder you will find:
1. **GGP.mlx** and  "[**GGP.html**](https://htmlpreview.github.io/?https://github.com/topggp/GGP-Matlab/blob/master/GGP.html)" : Matlab live script and html that contain the main code and comments to run Genaralized Geometry Ptojection. At the first look the script may seems long. In reality many parts are copied by top88 matlab code [4]. Moreover graphics also take many lines but are not essential for the understanding of the code.
2. **GGP_main.m**: Matlab script that contains the main code and comments to run Genaralized Geometry Ptojection. Could be usefull for Matlab users with versions prior to R2016a.
3. **Wgp.m**: Matlab function that takes as inputs the projection parameters provided in the main file, the coordinates of sampling window Gauss points and reads in output smooth characteristic function and derivatives computed for each component and sampling window Gauss point.
4. **lgwt.m**: Matlab function that provides Gauss point coordinates and weights for a 1D domain [3].
5. **model_updateM.m** and **model_updateV.m**: Matlab functions that takes as inputs local volume fractions and provide Finite Element Young's Modulus, densities and their sensitivities. 
6. **Aggregation_Pi.m** function adopted to make the union of components using different methods. This also provides sensitivities with respect to entries.
7. **smooth_sat.m**: matlab function that applies smooth saturation to avoid aggregation related issues [1].
8. **mmasub.m**, **subsolv.m**, **kktcheck.m**: Matlab functions for MMA optimization solver and KKT residual norm [3].
9. **norato_bar.m**: Matlab function that provides the distance of a component boundary from the component center as a function of the polar coordinate angle. This is needed for component plots.

[![View GGP-Matlab on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://fr.mathworks.com/matlabcentral/fileexchange/75249-ggp-matlab)

# References
[1] Coniglio, Simone et al. "Generalized Geometry Projection: A unified approach for geometric feature based topology optimization." Archives of Computational Methods in Engineering, DOI: 10.1007/s11831-019-09362-8 (2019) 

[2] Svanberg, Krister. "MMA and GCMMA, versions September 2007." Optimization and Systems Theory 104 (2007).

[3] Von Winckel, Greg. "[Legendre-Gauss Quadrature Weights and Nodes](https://fr.mathworks.com/matlabcentral/fileexchange/4540-legendre-gauss-quadrature-weights-and-nodes)"  (2004).

[4] Andreassen, Erik, et al. "Efficient topology optimization in MATLAB using 88 lines of code." Structural and Multidisciplinary Optimization 43.1 (2011): 1-16.
