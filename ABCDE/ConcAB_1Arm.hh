double ConcAB_1Arm(double ****phi,double ****w,double *Ns,double ds,double ***k_vector,double *dxyz){

  int         i,j,l,s;
  double      Q; 
  double      ****qA1;
  double      ****qB1;
  double      ****qdagA1;
  double      ****qdagB1;
  double      ***qint;
  double      ***qintA1;
  double      ***qintB1;

  // Used in ConcAB.hh
  qA1=create_4d_double_array(Nx,Ny,Nz,((int)Ns[0]+1),"qA1");
  // --
  qB1=create_4d_double_array(Nx,Ny,Nz,((int)Ns[4]+1),"qB1");
  // --
  qdagA1=create_4d_double_array(Nx,Ny,Nz,((int)Ns[0]+1),"qdagA1");
  // --
  qdagB1=create_4d_double_array(Nx,Ny,Nz,((int)Ns[4]+1),"qdagB1");
  // --
  qint=create_3d_double_array(Nx,Ny,Nz,"qint");
  // --
  qintA1=create_3d_double_array(Nx,Ny,Nz,"qintA1");
  // --
  qintB1=create_3d_double_array(Nx,Ny,Nz,"qintB1");
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


  // Here is the for loop for doing the qint, setting it to 1.0
  for(i=0;i<Nx;i++){
    for(j=0;j<Ny;j++){
      for(l=0;l<Nz;l++){
	qint[i][j][l]=1.0;
      }
    }
  }

  // Here we will solve the diffusion question
  solveModDiffEqn_FFT(qB1,w[4],qint,ds,(int)Ns[4],1,k_vector,dxyz);

  // The result from the above calculation becomes qdags initial cond
  for(i=0;i<Nx;i++){
    for(j=0;j<Ny;j++){
      for(l=0;l<Nz;l++){
	qintA1[i][j][l]=qB1[i][j][l][(int)Ns[4]];
      }
    }
  }
 
 
  // Here we will solve the diffusion question
  solveModDiffEqn_FFT(qA1,w[0],qintA1,ds,(int)Ns[0],1,k_vector,dxyz);

  
  // The result from the above calculation becomes qdags initial cond
  for(i=0;i<Nx;i++){
    for(j=0;j<Ny;j++){
      for(l=0;l<Nz;l++){
	qintA1[i][j][l]=1.0;
      }
    }
  }
  

  // Here we will solve the diffusion question
  solveModDiffEqn_FFT(qdagA1,w[0],qintA1,ds,(int)Ns[0],-1,k_vector,dxyz);
  

  // The result from the above calculation becomes qdags initial cond
  for(i=0;i<Nx;i++){
    for(j=0;j<Ny;j++){
      for(l=0;l<Nz;l++){
	qintB1[i][j][l]=qdagA1[i][j][l][(int)Ns[0]];
      }
    }
  }
  
  // Here we will solve the diffusion question
  solveModDiffEqn_FFT(qdagB1,w[4],qintB1,ds,(int)Ns[4],-1,k_vector,dxyz);
 
  //*******************************************************************************************************************

  //std::cout<<"+++++"<<std::endl;
  // Here we are doing the sum to get the single chain partition function
  Q=0.0;
  for(i=0;i<Nx;i++){
    for(j=0;j<Ny;j++){
      for(l=0;l<Nz;l++){
	Q+=(qA1[i][j][l][(int)Ns[0]])*dxyz[0]*dxyz[1]*dxyz[2];
      }
    }
  }
  // Normalizing with respect to the volume of the box
  Q/=((dxyz[0]*Nx)*(dxyz[1]*Ny)*(dxyz[2]*Nz));


  
  // Here we do the concentration calculation
  for(i=0;i<Nx;i++){
    for(j=0;j<Ny;j++){
      for(l=0;l<Nz;l++){

	phi[0][i][j][l]=0.0;  //A1

	phi[4][i][j][l]=0.0;  //B1

	//A1
	for(s=0;s<(Ns[0]+1);s++){
	  if(s==0 || s==(int)Ns[0]){
	    phi[0][i][j][l]+=0.5*qA1[i][j][l][s]*qdagA1[i][j][l][(int)Ns[0]-s]*ds;
	  }else{
	    phi[0][i][j][l]+=qA1[i][j][l][s]*qdagA1[i][j][l][(int)Ns[0]-s]*ds;
	  }
	}

	//B1
	for(s=0;s<(Ns[4]+1);s++){
	  if(s==0 || s==(int)Ns[4]){
	    phi[4][i][j][l]+=0.5*qB1[i][j][l][s]*qdagB1[i][j][l][(int)Ns[4]-s]*ds;
	  }else{
	    phi[4][i][j][l]+=qB1[i][j][l][s]*qdagB1[i][j][l][(int)Ns[4]-s]*ds;
	  }
	}

	phi[0][i][j][l]*=(pMultiAve/Q);
	phi[1][i][j][l]=0.0;
	phi[2][i][j][l]=0.0;
	phi[3][i][j][l]=0.0;
	phi[4][i][j][l]*=(pMultiAve/Q);
	phi[5][i][j][l]=0.0;
	phi[6][i][j][l]=0.0;
	phi[7][i][j][l]=0.0;

      }
    }
  }
  
  // Used in ConcAB.hh
  destroy_4d_double_array(qA1);
  destroy_4d_double_array(qB1);
  destroy_4d_double_array(qdagA1);
  destroy_4d_double_array(qdagB1);  
  destroy_3d_double_array(qint);
  destroy_3d_double_array(qintA1);
  destroy_3d_double_array(qintB1);
  //++++++++++++++++++++++++++++++

  return Q;


};
