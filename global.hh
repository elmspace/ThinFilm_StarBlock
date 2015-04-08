//#include </usr/local/include/fftw3.h>  // This is for My Mac Pro
//#include </opt/sharcnet/fftw/3.3.2/intel/include/fftw3.h> // This is for Sharcnet
//#include </usr/include/fftw3.h> // This is for use on Landua
#include </usr/local/include/fftw3.h> // This is for elmspace2


#include <stdio.h>     //Include the standard input/output libraries
#include <iostream>  //Cout and Cin etc.
#include <fstream>
#include <stdlib.h>    //Include standard fucntion libraries
#include <math.h>      //Use the math function libraries
#include "./include/smemory.hh"  //Use my custom memory handling class

 

#define Nx 22
#define Ny 22
#define Nz 22

#define ChainType 8
#define Pi 3.14159

// Trigger parameters
int Iomega;
int box_min;

// FFTW parameters
fftw_plan forward_plan, inverse_plan;
double *input_q, *transformed_q, *final_q;

// Other parameters
double h_AAir, h_BAir, h_ASub, h_BSub;
double Lx, Ly, Lz;
double pA1ave, pA2ave, pA3ave, pA4ave, pB1ave, pB2ave,pB3ave, pB4ave;


// Used in solvediffeq.hh
double         ***wds, ***kds;


// Used in ConcAB.hh
double      ****qA1,****qA2,****qA3,****qA4;
double      ****qB1,****qB2,****qB3,****qB4;
double      ****qdagA1,****qdagA2,****qdagA3,****qdagA4;
double      ****qdagB1,****qdagB2,****qdagB3,****qdagB4;
double      ***qint;
double      ***qintA1,***qintA2,***qintA3,***qintA4;
double      ***qintB1,***qintB2,***qintB3,***qintB4;
