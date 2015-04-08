
#include "./global.hh"
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
  chi=create_1d_double_array(2,"chi");
  f=create_1d_double_array(ChainType,"f");
  Ns=create_1d_double_array(ChainType,"Ns");
  k_vector=create_3d_double_array(Nx,Ny,Nz,"k_vector");
  dxyz=create_1d_double_array(3,"dxyz");
  chiMatrix=create_2d_double_array(ChainType,ChainType,"chiMatrix");

  // Used in solvediffeq.hh
  wds=create_3d_double_array(Nx,Ny,Nz,"wds");
  kds=create_3d_double_array(Nx,Ny,Nz,"kds");
  //+++++++++++++++++++++++++++++++++++++++++

  // Used in ConcAB.hh
  qA1=create_4d_double_array(Nx,Ny,Nz,((int)Ns[0]+1),"qA1");
  qA2=create_4d_double_array(Nx,Ny,Nz,((int)Ns[1]+1),"qA2");
  qA3=create_4d_double_array(Nx,Ny,Nz,((int)Ns[2]+1),"qA3");
  qA4=create_4d_double_array(Nx,Ny,Nz,((int)Ns[3]+1),"qA4");
  // --
  qB1=create_4d_double_array(Nx,Ny,Nz,((int)Ns[4]+1),"qB1");
  qB2=create_4d_double_array(Nx,Ny,Nz,((int)Ns[5]+1),"qB2");
  qB3=create_4d_double_array(Nx,Ny,Nz,((int)Ns[6]+1),"qB3");
  qB4=create_4d_double_array(Nx,Ny,Nz,((int)Ns[7]+1),"qB4");
  // --
  qdagA1=create_4d_double_array(Nx,Ny,Nz,((int)Ns[0]+1),"qdagA1");
  qdagA2=create_4d_double_array(Nx,Ny,Nz,((int)Ns[1]+1),"qdagA2");
  qdagA3=create_4d_double_array(Nx,Ny,Nz,((int)Ns[2]+1),"qdagA3");
  qdagA4=create_4d_double_array(Nx,Ny,Nz,((int)Ns[3]+1),"qdagA4");
  // --
  qdagB1=create_4d_double_array(Nx,Ny,Nz,((int)Ns[4]+1),"qdagB1");
  qdagB2=create_4d_double_array(Nx,Ny,Nz,((int)Ns[5]+1),"qdagB2");
  qdagB3=create_4d_double_array(Nx,Ny,Nz,((int)Ns[6]+1),"qdagB3");
  qdagB4=create_4d_double_array(Nx,Ny,Nz,((int)Ns[7]+1),"qdagB4");
  // --
  qint=create_3d_double_array(Nx,Ny,Nz,"qint");
  // --
  qintA1=create_3d_double_array(Nx,Ny,Nz,"qintA1");
  qintA2=create_3d_double_array(Nx,Ny,Nz,"qintA2");
  qintA3=create_3d_double_array(Nx,Ny,Nz,"qintA3");
  qintA4=create_3d_double_array(Nx,Ny,Nz,"qintA4");
  // --
  qintB1=create_3d_double_array(Nx,Ny,Nz,"qintB1");
  qintB2=create_3d_double_array(Nx,Ny,Nz,"qintB2");
  qintB3=create_3d_double_array(Nx,Ny,Nz,"qintB3");
  qintB4=create_3d_double_array(Nx,Ny,Nz,"qintB4");
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


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

  // Used in solvediffeq.hh
  destroy_3d_double_array(wds);
  destroy_3d_double_array(kds);
  //+++++++++++++++++++++++++++

  // Used in ConcAB.hh
  destroy_4d_double_array(qA1);
  destroy_4d_double_array(qA2);
  destroy_4d_double_array(qA3);
  destroy_4d_double_array(qA4);
  destroy_4d_double_array(qB1);
  destroy_4d_double_array(qB2);
  destroy_4d_double_array(qB3);
  destroy_4d_double_array(qB4);
  destroy_4d_double_array(qdagA1);
  destroy_4d_double_array(qdagA2);
  destroy_4d_double_array(qdagA3);
  destroy_4d_double_array(qdagA4);
  destroy_4d_double_array(qdagB1);  
  destroy_4d_double_array(qdagB2);  
  destroy_4d_double_array(qdagB3);  
  destroy_4d_double_array(qdagB4);
  destroy_3d_double_array(qint);
  destroy_3d_double_array(qintA1);
  destroy_3d_double_array(qintA2);
  destroy_3d_double_array(qintA3);
  destroy_3d_double_array(qintA4);
  destroy_3d_double_array(qintB1);
  destroy_3d_double_array(qintB2);
  destroy_3d_double_array(qintB3);
  destroy_3d_double_array(qintB4);
  //++++++++++++++++++++++++++++++


  return 0;
}
