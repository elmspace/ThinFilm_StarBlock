void FreeEnergy(double ****w, double ****phi, double ***eta, double *Ns, double ds, double ***k_vector, double *chi, double *dxyz, double **chiMatrix, double ****h){

  
 
  int     maxIter=50; 
  int     msg=1;
  int     i,j,k,iter,chain,ii,jj; 
  double  currentfE, oldfE, deltafE,oldfE_iter; 
  double  precision=1.0e-2; 
  double  QAB; 
  double  fEW, fEchi, fES, fEsurf; 
  double  epsilon, gamma;
  double  ***delphi;
  double  ****delW;
  double  ****newW;
  double  deltaW;
  double  fE_homo;
  double  *avphi;

  avphi=create_1d_double_array(ChainType,"avphi");
  delW=create_4d_double_array(ChainType,Nx,Ny,Nz,"delW");
  delphi=create_3d_double_array(Nx,Ny,Nz,"delphi");
  newW=create_4d_double_array(ChainType,Nx,Ny,Nz,"newW");

  avphi[0]=pA1ave; // A1 average
  avphi[1]=pA2ave; // A2 average
  avphi[2]=pA3ave; // A3 average
  avphi[3]=pA4ave; // A4 average
  avphi[4]=pB1ave; // B1 average
  avphi[5]=pB2ave; // B2 average
  avphi[6]=pB3ave; // B3 average
  avphi[7]=pB4ave; // B4 average

  oldfE=1.0e2;
  std::ofstream outputFile("./RESULTS/fE.dat");
  do{
   
    WaveVectors(k_vector,dxyz);
    currentfE=0.0;
    deltafE=0.0;
  
    epsilon=0.025;
    gamma=0.025;
  
    iter=0;  
    
    do{
      
      iter++;
    
      // setting free energies to 0
      fEW=0.0;
      fEchi=0.0;
      fES=0.0;
      fEsurf=0.0;
      deltaW=0.0;

      QAB=ConcAB(phi,w,Ns,ds,k_vector,dxyz);

      Incomp(eta,phi,delphi);
    
      // Initializing newW
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
	  
		newW[ii][i][j][k]+=chiMatrix[ii][jj]*(phi[jj][i][j][k]-avphi[jj]);
		fEchi+=phi[ii][i][j][k]*chiMatrix[ii][jj]*phi[jj][i][j][k]*dxyz[0]*dxyz[1]*dxyz[2];

	      }

	      newW[ii][i][j][k]+=h[ii][i][j][k]+eta[i][j][k];

	      fEW+=(newW[ii][i][j][k]*phi[ii][i][j][k]*dxyz[0]*dxyz[1]*dxyz[2]);
	      fEsurf+=phi[ii][i][j][k]*h[ii][i][j][k]*dxyz[0]*dxyz[1]*dxyz[2];
	      delW[ii][i][j][k]=newW[ii][i][j][k]-w[ii][i][j][k];
	      deltaW+=fabs(delW[ii][i][j][k]);
	    }
	 
	  }
	}
      }
  
      deltaW/=(Nx*Ny*Nz);
      fEchi/=(2.0*((Nx*dxyz[0])*(Ny*dxyz[1])*(Nz*dxyz[2])));
      fEW/=(((Nx*dxyz[0])*(Ny*dxyz[1])*(Nz*dxyz[2])));
      fEsurf/=(((Nx*dxyz[0])*(Ny*dxyz[1])*(Nz*dxyz[2])));

      fES=log(QAB);   
      fE_homo=homogenousfE(chiMatrix);

      currentfE=-fES-fEW+fEchi+fEsurf;
      
      deltafE=fabs(currentfE-oldfE_iter);

      std::cout<<"Iter="<<iter<<"  fE="<<currentfE<<"  delW=" <<deltaW<<"  delfE="<<currentfE-fE_homo<<std::endl;
      oldfE_iter=currentfE;

      // Updating the new W-field
      for(i=0;i<Nx;i++){
	for(j=0;j<Ny;j++){
	  for(k=0;k<Nz;k++){
	    for(chain=0;chain<ChainType;chain++){
	      w[chain][i][j][k]+=(gamma*delW[chain][i][j][k]-epsilon*delphi[i][j][k]);
	    }
	  }
	}
      }
  
      SaveData(phi,w,dxyz);

    }while((deltaW>precision)||(iter<maxIter));
  
    outputFile <<currentfE<<" "<<fE_homo<<" "<<dxyz[0]*Nx<<" "<<dxyz[1]*Ny<<" "<<dxyz[2]*Nz<<std::endl;
  
    if(box_min_xy_relax==1){size_adjust_2D_xy(w,phi,eta,Ns,ds,k_vector,chi,dxyz,chiMatrix);}
    if(box_min_xyz_relax==1){size_adjust(w,phi,eta,Ns,ds,k_vector,chi,dxyz,chiMatrix);}

    if((oldfE<currentfE)||((abs(oldfE)-abs(currentfE))<1.0e-4)){
      msg=0;
    }
    if(msg==1){
      oldfE=currentfE;
    }
    if(box_min==0){
      msg=0;
    }
    
  }while(msg==1);

  outputFile <<"Done"<<std::endl;
  outputFile.close();

  
  destroy_3d_double_array(delphi);
  destroy_4d_double_array(delW);
  destroy_4d_double_array(newW);
  destroy_1d_double_array(avphi);
  
};
