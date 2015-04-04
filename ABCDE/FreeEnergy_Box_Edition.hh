double FreeEnergy_Box_Edition(double ****w_temp, double ****phi, double ***eta, double *Ns, double ds, double ***k_vector, double *chi, double *dxyz_temp, double **chiMatrix){

  
  double  currentfE, oldfE, deltafE; 
  int     maxIter=1; 
  int     i,j,k,iter,chain,ii,jj; 
  double  precision=1.0e-3; 
  double  QAB,QD,QE; 
  double  fEW, fEchi, fES; 
  double  epsilon, gamma;
  double  ***delphi;
  double  ****delW;
  double  ****newW;
  double  deltaW;
  double  fE_homo;

  delW=create_4d_double_array(ChainType,Nx,Ny,Nz,"delW");
  delphi=create_3d_double_array(Nx,Ny,Nz,"delphi");
  newW=create_4d_double_array(ChainType,Nx,Ny,Nz,"newW");

  WaveVectors(k_vector,dxyz_temp);
  
  currentfE=0.0;
  deltafE=0.0;
  
  epsilon=0.01;
  gamma=0.01;
  
  iter=0;  
  do{
    iter++;
    
    fEW=0.0;
    fEchi=0.0;
    fES=0.0;
    
    QAB=ConcAB(phi,w_temp,Ns,ds,k_vector,dxyz_temp);
    QD=ConcD(phi);
    QE=ConcE(phi);

    Incomp(eta,phi,delphi);
    
    fEW=0.0;
    fEchi=0.0;
      
    deltaW=0.0;
    
    
    for(i=0;i<Nx;i++){
      for(j=0;j<Ny;j++){
	for(k=0;k<Nz;k++){
	  for(ii=0;ii<ChainType;ii++){
	    newW[ii][i][j][k]=0.0;  
	  }
	}
      }
    }
    for(i=0;i<Nx;i++){
      for(j=0;j<Ny;j++){
	for(k=0;k<Nz;k++){
	  
	  for(ii=0;ii<ChainType;ii++){
	    for(jj=0;jj<ChainType;jj++){
	      
	      newW[ii][i][j][k]+=((chiMatrix[ii][jj]*phi[jj][i][j][k]));
	      fEchi+=phi[ii][i][j][k]*chiMatrix[ii][jj]*phi[jj][i][j][k]*dxyz_temp[0]*dxyz_temp[1]*dxyz_temp[2];
	    
	    }

	    if(ii==0){newW[ii][i][j][k]+=eta[i][j][k];}    //A1
	    if(ii==1){newW[ii][i][j][k]+=eta[i][j][k];}    //A2
	    if(ii==2){newW[ii][i][j][k]+=eta[i][j][k];}    //A3
	    if(ii==3){newW[ii][i][j][k]+=eta[i][j][k];}    //A4
	    if(ii==4){newW[ii][i][j][k]+=eta[i][j][k];}    //B1
	    if(ii==5){newW[ii][i][j][k]+=eta[i][j][k];}    //B2
	    if(ii==6){newW[ii][i][j][k]+=eta[i][j][k];}    //B3
	    if(ii==7){newW[ii][i][j][k]+=eta[i][j][k];}    //B4
	    if(ii==8){newW[ii][i][j][k]=0.0;}              //E
	    if(ii==9){newW[ii][i][j][k]=0.0;}              //D

	    fEW+=(newW[ii][i][j][k]*phi[ii][i][j][k]*dxyz_temp[0]*dxyz_temp[1]*dxyz_temp[2]);
	    delW[ii][i][j][k]=newW[ii][i][j][k]-w_temp[ii][i][j][k];
	    deltaW+=fabs(delW[ii][i][j][k]);
	  }
	}
      }
    }
    deltaW/=(Nx*Ny*Nz);
    fEchi/=(2.0*((Nx*dxyz_temp[0])*(Ny*dxyz_temp[1])*(Nz*dxyz_temp[2])));
    fEW/=(((Nx*dxyz_temp[0])*(Ny*dxyz_temp[1])*(Nz*dxyz_temp[2])));
    fES=p1ave*log(QAB);
        
    oldfE=currentfE;
    
    fE_homo=homogenousfE(chiMatrix);
    
    //currentfE=-fES-fEW+fEchi-fE_homo;
    currentfE=-fES;

    deltafE=fabs(currentfE-oldfE);
      
    for(i=0;i<Nx;i++){
      for(j=0;j<Ny;j++){
	for(k=0;k<Nz;k++){
	  for(chain=0;chain<ChainType;chain++){
	    w_temp[chain][i][j][k]+=(gamma*delW[chain][i][j][k]-epsilon*delphi[i][j][k]);
	  }
	}
      }
    }
    
   
  }while(iter<maxIter);//(deltafE>precision);

  return currentfE;

  
  destroy_3d_double_array(delphi);
  destroy_4d_double_array(delW);
  destroy_4d_double_array(newW);
  
};
