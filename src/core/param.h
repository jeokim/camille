#ifndef CORE_PARAM_H
#define CORE_PARAM_H

//

const int TRUE = 1;
const int FALSE = 0;
const int OK = 1;
const int NOT_OK = 0;
const int MINIMUM = 0;
const int MAXIMUM = 1;

const int NONE = -1;

const int DIM_MAX = 3;

const int FIRST = 0;
const int SECOND = 1;
const int THIRD = 2;
const int FOURTH = 3;
const int FIFTH = 4;
const int SIXTH = 5;
const int SEVENTH = 6;
const int EIGHTH = 7;
const int NINTH = 8;
const int TENTH = 9;
//
const int LEFT = 0;
const int RIGHT = 1;
//
const int NEXT = 0;
const int NEXTNEXT = 1;
//
const int IN = 0;
const int OUT = 1;

const int DUMMY_INT = -1;
const double DUMMY_DOUBLE = -1E+09;
const double DUMMY_LARGE = 1E+99;
const double DUMMY_SMALL = 1E-99;

// coordinate specification
const int XI = 0;
const int ETA = 1;
const int ZETA = 2;
//
// Cartesian coordinates defined as (x, y, z)
const int XDIR = 0;
const int YDIR = 1;
const int ZDIR = 2;
//
// Cylindrical coordinates defined as (x, r, theta)
const int RDIR = 1;
//
const int CARTESIAN = 0;
const int AXISYMMETRIC = 1;

// domain periodicity
const int NONPERIODIC = 0;
const int PERIODIC_PLANE = 1;

// physical model
//
// do NOT modify the following numbers associated with those starting with
// IVAR_ and IAUX_ since that may affect other physical models; if needed, 
// define new ones
// by defaults, time-advanced solution variables are \rho, u_x, u_y, u_z, & p
const int IVAR_RHO = 0;
const int IVAR_UX = 1;
const int IVAR_UY = 2;
const int IVAR_UZ = 3;
const int IVAR_P = 4;
// in the cylindrical coordinates
const int IVAR_UR = 2;
const int IVAR_UTHETA = 3;
// if entropy replaces density
const int IVAR_S = 0;
// mixture fraction
const int IVAR_Z = 5;

// auxiliary or dependent variables which are not time advanced
const int IAUX_RHO = 0;
const int IAUX_RHO_MEAN = 1;
const int IAUX_T = 2;
const int IAUX_T_MEAN = 3;
const int IAUX_CP = 4;
// if viscous and heat transfer matters
const int IAUX_MU = 5;
const int IAUX_MU_MEAN = 6;
const int IAUX_LAMBDA = 7;
const int IAUX_LAMBDA_MEAN = 8;
const int IAUX_D = 9;
const int IAUX_D_MEAN = 10;
// dc_p/dZ, \Psi for composition noise in the linearized Euler formulation
const int IAUX_DCPDZ_LEE = 5;
const int IAUX_PSI_LEE = 6;
// dc_p/dZ, \Psi for composition noise in the linearized Navier--Stokes formulation
const int IAUX_DCPDZ_LNS = 11;
const int IAUX_PSI_LNS = 12;

// temporal discretization

// spatial discretization
const int MAX_ORDER_ACCURACY = 10;
const int STANDARD_CENTRAL = 1;
const int SYMMETRIC = 0;
const int ANTISYMMETRIC = 1;
//
const int FDO11P = 11; // Bogey & Bailly (JCP 2004)
const int SFO11P = 11; // Bogey & Bailly (JCP 2004)

// types of boundary
const int BOUNDARY_PERIODIC = 0; // if periodic (neither boundary condition nor interpolation is applied)
const int BOUNDARY_BC = 1; // if boundary condition (e.g. Dirichlet, Neumann, characteristic) is applied
const int BOUNDARY_BUFFERZONE = 2;
const int BOUNDARY_OVERSET = 3; // if overset-grid interpolation is applied
//
const int BC_DIRICHLET = 0;
const int BC_DIRICHLET_ALLZERO = 1;
const int BC_DIRICHLET_HARMONICWAVE = 2;
//
const int BC_NEUMANN = 100;
//
const int BC_WALL_SLIP_KINEMATIC = 200;
const int BC_WALL_SLIP_KINEMATIC_X = 201;
const int BC_WALL_SLIP_KINEMATIC_Y = 202;
const int BC_WALL_SLIP_KINEMATIC_Z = 203;
//
const int BC_CENTERLINE_CART_NORM2X = 900;
const int BC_CENTERLINE_CART_NORM2Y = 901;
const int BC_CENTERLINE_CART_NORM2Z = 902;
const int BC_CENTERLINE_AXISYM = 903;
//
const int SPONGE_FREUND_AMBIENT = 1000;
const int SPONGE_FREUND_DIRICHLET = 1001;
const int SPONGE_FREUND_HARMONICWAVE = 1002;

// type of data files
const int PLOT3D = 0;

const int PLOT3D_GRID = 0;
const int PLOT3D_SOLUTION = 1;
const int PLOT3D_FUNCTION = 2;

// overset

// IBLANK
const int NOTBLANKED = 0;
const int BLANKED = 1;
const int INTERPOLATED = -1;

// scheme for computing metrics

// global parameters
const int dir_other[DIM_MAX][2] = {{ETA, ZETA}, {XI, ZETA}, {XI, ETA}};
const int ij2by2[2][2] = {{0, 1}, {2, 3}};
const int ij3by3[3][3] = {{0, 1, 2}, {3, 4, 5}, {6, 7, 8}};

// miscellaneous parameters
const char file_base[] = "camille";

//

#endif
