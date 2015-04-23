#include </usr/local/include/fftw3.h>  // This is for My Mac Pro
//#include </opt/sharcnet/fftw/3.3.2/intel/include/fftw3.h> // This is for Sharcnet
//#include </usr/include/fftw3.h> // This is for use on Landua
//#include </usr/local/include/fftw3.h> // This is for elmspace2


#include <stdio.h>     //Include the standard input/output libraries
#include <iostream>  //Cout and Cin etc.
#include <fstream>
#include <complex>
#include <time.h>
#include <stdlib.h>    //Include standard fucntion libraries
#include <math.h>      //Use the math function libraries
#include "./include/smemory.hh"  //Use my custom memory handling class


using namespace std;
 

#define Nx 20
#define Ny 20
#define Nz 20

#define ChainType 8
#define Pi 3.14159

// Trigger parameters
int Iomega;
int box_min;
int box_min_xy_relax;
int box_min_xyz_relax;
int LAM, HEX;
int HOR, VER;

int Round; // This is used to check the loop condition (because I want to use previous results for Round>1)

double Numb_of_Periods;
double Lam_Period;
double Hex_Period;

double xBAir;
int Numb_of_Arms;

double global_fE, flobal_HomfE, global_xAB, global_xAAir, global_xASub, global_xBAir, global_xBSub;

// FFTW parameters
fftw_plan forward_plan, inverse_plan;
double *input_q, *transformed_q, *final_q;

// Other parameters
double h_AAir, h_BAir, h_ASub, h_BSub;
double Lx, Ly, Lz;
double pA1ave, pA2ave, pA3ave, pA4ave, pB1ave, pB2ave,pB3ave, pB4ave;





