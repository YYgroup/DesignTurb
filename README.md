# DesignTurb: Vortex Construction Tool

## Description
DesignTurb offers a versatile method for constructing the vorticity of vortex tubes with customizable centerline topology, differential twist, and variable thickness. This tool is ideal for designing classical turbulence fields, conceptualized with quantum vortex tangles as the core framework and complemented with customizable vortices as modular elements. This approach allows for free and precise adjustment on the distribution and shape of elemental vortices. The code is written in FORTRAN and supports MPI parallel computing. If you are interested in using the code for your own research, please contact yyg@pku.edu.cn for more details.

## How It Works
- **Profile Shape**: Modify the profile shape of vortices via the `tpFunc` subroutine.
- **Thickness Distribution**: Adjust the thickness along the centerlines in the `sigmaFunc` subroutine.
- **Twist Distribution**: Change the twist distribution along the centerlines using the `etaFunc` subroutine.
- **Centerline Data**: Input as discrete points, exemplified in `TrefoilCenterline.dat` for a trefoil knot and `QTCenterline.dat` for entangled vortices. Entangled centerline can be obtained from superfluid simulation based on the vortex filament method (see [qvort](https://github.com/abag/qvort) on GitHub).

The code transforms discrete control points on centerlines into cubic spline curves, defined by polynomial parametric equations. It then constructs a 3D vorticity field for vortex tubes within a periodic box, translating curved cylindrical coordinates to Cartesian coordinates. The process ensures the vector field adheres to periodic boundary conditions.

### Centerline Parameters
The centerline data file should be organized in order by the following parameters:
- `nline`: Total number of centerlines
- `npointlist(i)`: Number of discrete points per centerline
- `cxall(i)`: x-coordinates of discrete points
- `cyall(i)`: y-coordinates of discrete points
- `czall(i)`: z-coordinates of discrete points

## Reference
W. Shen, J. Yao, and Y. Yang, "Weaving Classical Turbulence with Quantum Skeleton" (2024). [arXiv:2401.11149](https://arxiv.org/abs/2401.11149)
W. Shen, J. Yao, F. Hussain, and Y. Yang, “Role of internal structures within a vortex in helicity dynamics”, Journal of Fluid Mechanics, 970, A26 (2023). [JFM](https://doi.org/10.1017/jfm.2023.620)
W. Shen, J. Yao, F. Hussain, and Y. Yang, “Topological transition and helicity conversion of	vortex knots and links”, Journal of Fluid Mechanics, 943, A41 (2022). [JFM](https://doi.org/10.1017/jfm.2022.464)
S. Xiong, Y. Yang, “Construction of knotted vortex tubes with the writhe-dependent helicity”, Physics of Fluids, 31, 047101 (2019). [POF](https://doi.org/10.1063/1.5088015)
