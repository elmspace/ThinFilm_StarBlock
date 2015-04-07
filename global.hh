#include <stdio.h>     //Include the standard input/output libraries
#include <iostream>  //Cout and Cin etc.
#include <fstream>
#include <stdlib.h>    //Include standard fucntion libraries
#include <math.h>      //Use the math function libraries
#include <time.h>      //Call system time libraries to define the integer seed for random numbers
#include "./include/smemory.hh"  //Use my custom memory handling class
#include </opt/sharcnet/fftw/3.3.2/intel/include/fftw3.h>
//#include "mpi.h"     //Use this for MPI parallel implimentation later 

#define Nx 32
#define Ny 32
#define Nz 32

#define ChainType 8
#define Pi 3.14159

// Trigger parameters
int Iomega;
int box_min;

// FFTW parameters
fftw_plan forward_plan, inverse_plan;
double *input_q, *transformed_q, *final_q;

// Other parameters
double pA1ave, pA2ave, pA3ave, pA4ave, pB1ave, pB2ave,pB3ave, pB4ave;


