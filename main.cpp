
#include "./global.hh"
#include "./ABCDE/Arms.hh"
#include "./ABCDE/parametersABCD.hh"
#include "./ABCDE/WaveVectors.hh"
#include "./ABCDE/omega.hh"
#include "./ABCDE/solvediffeq.hh"
#include "./ABCDE/ConcAB.hh"
#include "./ABCDE/fEhomo.hh"
#include "./ABCDE/Incomp.hh"
#include "./ABCDE/FreeEnergy_Box_Edition.hh"
#include "./ABCDE/size_adjust.hh"
#include "./ABCDE/size_adjust_2D_xy.hh"
#include "./ABCDE/SaveData.hh"
#include "./ABCDE/FreeEnergy.hh"

using namespace std;


int main(){

  double ****w;
  double ***eta;
  double ****phi;
  double ****h;
  double *chi;
  double *f;
  double ds;
  double *Ns;
  double ***k_vector;
  double *dxyz;
  double **chiMatrix;
  //int FFTW_MEASURE;
  
  w=create_4d_double_array(ChainType,Nx,Ny,Nz,"w");
  eta=create_3d_double_array(Nx,Ny,Nz,"eta");
  phi=create_4d_double_array(ChainType,Nx,Ny,Nz,"phi");
  h=create_4d_double_array(ChainType,Nx,Ny,Nz,"h");
  chi=create_1d_double_array(1,"chi");
  f=create_1d_double_array(ChainType,"f");
  Ns=create_1d_double_array(ChainType,"Ns");
  k_vector=create_3d_double_array(Nx,Ny,Nz,"k_vector");
  dxyz=create_1d_double_array(3,"dxyz");
  chiMatrix=create_2d_double_array(ChainType,ChainType,"chiMatrix");

  long iseed;
  time_t t;
  iseed=time(&t);
  srand48(iseed);

  input_q=(double*)fftw_malloc(sizeof(double)*(Nx*Ny*Nz));
  transformed_q=(double*)fftw_malloc(sizeof(double)*(Nx*Ny*Nz));
  final_q=(double*)fftw_malloc(sizeof(double)*(Nx*Ny*Nz));


  forward_plan=fftw_plan_r2r_3d(Nx,Ny,Nz,input_q,transformed_q,FFTW_REDFT10,FFTW_REDFT10,FFTW_REDFT10,FFTW_PRESERVE_INPUT);
  inverse_plan=fftw_plan_r2r_3d(Nx,Ny,Nz,transformed_q,final_q,FFTW_REDFT01,FFTW_REDFT01,FFTW_REDFT01,FFTW_PRESERVE_INPUT);
  
  parametersAB(chi,f,ds,Ns,dxyz,chiMatrix,h);
 
  omega(w);

  FreeEnergy(w,phi,eta,Ns,ds,k_vector,chi,dxyz,chiMatrix,h);
  


  //Destroy memory allocations------------
  fftw_free(input_q);
  fftw_free(transformed_q);
  fftw_free(final_q);
  
  fftw_destroy_plan(forward_plan);
  fftw_destroy_plan(inverse_plan);

  destroy_4d_double_array(w);
  destroy_3d_double_array(eta);
  destroy_4d_double_array(phi);
  destroy_4d_double_array(h);
  destroy_1d_double_array(chi);
  destroy_1d_double_array(Ns);
  destroy_1d_double_array(f);
  destroy_3d_double_array(k_vector);
  destroy_1d_double_array(dxyz);
  destroy_2d_double_array(chiMatrix);
  //-------------------------------------

  return 0;
}
