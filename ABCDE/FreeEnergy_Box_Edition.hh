double FreeEnergy_Box_Edition(double ****w_temp, double ****phi, double ***eta, double *Ns, double ds, double ***k_vector, double *chi, double *dxyz_temp, double **chiMatrix){

  
  double  currentfE; 
 
  int     i,j,k,chain,ii,jj; 
  double  precision=1.0e-3; 
  double  QAB; 
  double  fES; 

  WaveVectors(k_vector,dxyz_temp);
  
  currentfE=0.0;
  fES=0.0;
  
  QAB=ConcAB(phi,w_temp,Ns,ds,k_vector,dxyz_temp);
  
  fES=log(QAB);
  currentfE=-fES;
  
  return currentfE;

  
};
