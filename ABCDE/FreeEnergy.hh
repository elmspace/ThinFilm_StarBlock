void FreeEnergy(double ****w, double ****phi, double ***eta, double *Ns, double ds, double ***k_vector, double *chi, double *dxyz, double **chiMatrix){

  
  double  currentfE, oldfE, deltafE,oldfE_iter; 
  int     maxIter=4000; 
  int     i,j,k,iter,chain,ii,jj; 
  double  precision=1.0e-2; 
  double  QAB,QD,QE; 
  double  fEW, fEchi, fES; 
  double  epsilon, gamma;
  double  ***delphi;
  double  ****delW;
  double  ****newW;
  double  deltaW;
  double  fE_homo;
  double  box,msg;
  double  *avphi;
  double  order_parameter;


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
  avphi[8]=pEave; //  E  average
  avphi[9]=pDave; //  D  average

  msg=1.0;
  oldfE=1.0e2;
  std::ofstream outputFile("./fE.dat");
  do{
   
    WaveVectors(k_vector,dxyz);
    currentfE=0.0;
    deltafE=0.0;
  
    epsilon=0.025;
    gamma=0.025;
  
    iter=0;  
    
    do{
      
      iter++;
    
      fEW=0.0;
      fEchi=0.0;
      fES=0.0;

      
      QAB=ConcAB(phi,w,Ns,ds,k_vector,dxyz);
      
      QD=ConcD(phi);
      QE=ConcE(phi);
     
      Incomp(eta,phi,delphi);
      
      fEW=0.0;
      fEchi=0.0;
      
      deltaW=0.0;

      order_parameter=0.0;     
      // Some output for checking the progress
      std::ofstream outputFile7("./phi.dat");
      for(k=0;k<Nz;k++){
	outputFile7 <<k<<" "<<phi[0][Nx/2][Ny/2][k]+phi[1][Nx/2][Ny/2][k]+phi[2][Nx/2][Ny/2][k]+phi[3][Nx/2][Ny/2][k]<<" "<<phi[4][Nx/2][Ny/2][k]+phi[5][Nx/2][Ny/2][k]+phi[6][Nx/2][Ny/2][k]+phi[7][Nx/2][Ny/2][k]<<" "<<phi[8][Nx/2][Ny/2][k]<<" "<<phi[9][Nx/2][Ny/2][k]<<std::endl;
	order_parameter+=fabs(0.5-(phi[0][Nx/2][Ny/2][k]+phi[1][Nx/2][Ny/2][k]+phi[2][Nx/2][Ny/2][k]));
      }
      outputFile7.close();     
      order_parameter/=Nz;
      

      // Some output for checking the progress
      std::ofstream outputFile17("./phi2D.dat");
      for(i=0;i<Nx;i++){
	for(k=0;k<Nz;k++){
	  outputFile17 <<i*dxyz[0]<<" "<<k*dxyz[2]<<" "<<phi[0][i][Ny/2][k]+phi[1][i][Ny/2][k]+phi[2][i][Ny/2][k]+phi[3][i][Ny/2][k]<<" "<<phi[4][i][Ny/2][k]+phi[5][i][Ny/2][k]+phi[6][i][Ny/2][k]+phi[7][i][Ny/2][k]<<std::endl;
	}}
      outputFile17.close();     


            
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

	      if(ii==0){newW[ii][i][j][k]+=eta[i][j][k];}    //A1
	      if(ii==1){newW[ii][i][j][k]+=eta[i][j][k];}    //A2
	      if(ii==2){newW[ii][i][j][k]+=eta[i][j][k];}    //A3
	      if(ii==3){newW[ii][i][j][k]+=eta[i][j][k];}    //A4

	      if(ii==4){newW[ii][i][j][k]+=eta[i][j][k];}    //B1
	      if(ii==5){newW[ii][i][j][k]+=eta[i][j][k];}    //B2
	      if(ii==6){newW[ii][i][j][k]+=eta[i][j][k];}    //B3
	      if(ii==7){newW[ii][i][j][k]+=eta[i][j][k];}    //B4

	      if(ii==8){newW[ii][i][j][k]=0.0;}  //E
	      if(ii==9){newW[ii][i][j][k]=0.0;}  //D

	      fEW+=(newW[ii][i][j][k]*phi[ii][i][j][k]*dxyz[0]*dxyz[1]*dxyz[2]);
	      delW[ii][i][j][k]=newW[ii][i][j][k]-w[ii][i][j][k];
	      deltaW+=fabs(delW[ii][i][j][k]);
	    }
	 
	  }
	}
      }

      deltaW/=(Nx*Ny*Nz);
      fEchi/=(2.0*((Nx*dxyz[0])*(Ny*dxyz[1])*(Nz*dxyz[2])));
      fEW/=(((Nx*dxyz[0])*(Ny*dxyz[1])*(Nz*dxyz[2])));
    
      fES=p1ave*log(QAB);
   
      fE_homo=homogenousfE(chiMatrix);

      //currentfE=-fES-fEW+fEchi-fE_homo;
      currentfE=-fES-fEW+fEchi;               // taking out the homogenous fE

      deltafE=fabs(currentfE-oldfE_iter);

      std::cout<<iter<<" "<<currentfE<< " " << deltaW<<" "<<deltafE<<std::endl;
      oldfE_iter=currentfE;

      for(i=0;i<Nx;i++){
	for(j=0;j<Ny;j++){
	  for(k=0;k<Nz;k++){

	    for(chain=0;chain<ChainType;chain++){
	      w[chain][i][j][k]+=(gamma*delW[chain][i][j][k]-epsilon*delphi[i][j][k]);
	    }

	  }
	}
      }
  

    }while(deltaW>precision);
   
    //+++++++++++++++++++++++++++++ This output is setup for the matlab plotting +++++++++++++++++++
    std::ofstream outputFile7("./xyz.dat");
    for (i=0;i<Nx;i++){
      outputFile7<<i*dxyz[0]<<" "<<i*dxyz[1]<<" "<<i*dxyz[2]<<std::endl;
    }
    outputFile7.close();  
    std::ofstream outputFile8("./ABCD.dat");
    for (i=0;i<Nx;i++){
      for(j=0;j<Ny;j++){
	for(k=0;k<Nz;k++){
	  outputFile8<<phi[0][i][j][k]+phi[1][i][j][k]+phi[2][i][j][k]+phi[3][i][j][k]<<" "<<phi[4][i][j][k]+phi[5][i][j][k]+phi[6][i][j][k]+phi[7][i][j][k]<<" "<<phi[8][i][j][k]<<" "<<phi[9][i][j][k]<<std::endl;
	}}}
    outputFile8.close();
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    outputFile <<currentfE<<" "<<fE_homo<<" "<<dxyz[0]*Nx<<" "<<dxyz[1]*Ny<<" "<<dxyz[2]*Nz<<std::endl;
  
    
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // Writting to data files
    std::ofstream outputFile6("./omega.dat");
    for(i=0;i<Nx;i++){
      for(j=0;j<Ny;j++){
	for(k=0;k<Nz;k++){
	  outputFile6 <<i<<" "<<j<<" "<<k<< " "<<w[0][i][j][k]<<" "<<w[1][i][j][k]<<" "<<w[2][i][j][k]<<" "<<w[3][i][j][k]<<" "<<w[4][i][j][k]<< " "<<w[5][i][j][k]<<" "<<w[6][i][j][k]<<" "<<w[7][i][j][k]<<" "<<w[8][i][j][k]<<" "<<w[9][i][j][k]<<std::endl;
	}
      }
    }
    outputFile6.close();
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      
  

    box=size_adjust_2D_xy(w,phi,eta,Ns,ds,k_vector,chi,dxyz,chiMatrix);
 
   
    if(oldfE<currentfE){
      msg=0.0;
    }
    if(msg>0.5){
      oldfE=currentfE;
    }
    if(box_min<0.0){
      msg=0.0;
    }
    
  }while(msg>0.5);

  outputFile <<"Done"<<std::endl;
  outputFile.close();

  
  destroy_3d_double_array(delphi);
  destroy_4d_double_array(delW);
  destroy_4d_double_array(newW);
  destroy_1d_double_array(avphi);
  
};
