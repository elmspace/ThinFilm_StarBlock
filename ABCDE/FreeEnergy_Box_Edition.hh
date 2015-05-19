double FreeEnergy_Box_Edition(double ****w_temp, double ****phi, double ***eta, double *Ns, double ds, double ***k_vector, double *chi, double *dxyz_temp, double **chiMatrix){

  
  double  currentfE; 
 
  int     i,j,k,chain,ii,jj; 
  double  precision=1.0e-3; 
  double  QAB,QHA,QHS; 
  double  fES; 

  WaveVectors(k_vector,dxyz_temp);
  
  currentfE=0.0;
  fES=0.0;
  
  if(Numb_of_Arms==1){
    QAB=ConcAB_1Arm(phi,w_temp,Ns,ds,k_vector,dxyz_temp);
  }else if(Numb_of_Arms==2){
    QAB=ConcAB_2Arm(phi,w_temp,Ns,ds,k_vector,dxyz_temp);
  }else if(Numb_of_Arms==3){
    QAB=ConcAB_3Arm(phi,w_temp,Ns,ds,k_vector,dxyz_temp);
  }else if(Numb_of_Arms==4){
    QAB=ConcAB_4Arm(phi,w_temp,Ns,ds,k_vector,dxyz_temp);
  }
  QHA=ConcHA(phi,w_temp,Ns,ds,k_vector,dxyz_temp);
  QHS=ConcHS(phi,w_temp,Ns,ds,k_vector,dxyz_temp);
  
  //fES=(pMultiAve)*log(QAB)+(pAirAve/kappa_HA)*log(QHA)+(pSubAve/kappa_HS)*log(QHS); // This includes the Air and Substrate Homopolymers
  // I think this is the fES we want, for the film alone.
  fES=(pMultiAve)*log(QAB); 
  currentfE=-fES;
  
  return currentfE;

  
};
