double ConcAB(double ****phi,double ****w,double *Ns,double ds,double ***k_vector,double *dxyz){

  int         i,j,l,s;
  double      Q;
 
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
  solveModDiffEqn_FFT(qB2,w[5],qint,ds,(int)Ns[5],1,k_vector,dxyz);
  solveModDiffEqn_FFT(qB3,w[6],qint,ds,(int)Ns[6],1,k_vector,dxyz);
  solveModDiffEqn_FFT(qB4,w[7],qint,ds,(int)Ns[7],1,k_vector,dxyz);


  // The result from the above calculation becomes qdags initial cond
  for(i=0;i<Nx;i++){
    for(j=0;j<Ny;j++){
      for(l=0;l<Nz;l++){
	qintA1[i][j][l]=qB1[i][j][l][(int)Ns[4]];
	qintA2[i][j][l]=qB2[i][j][l][(int)Ns[5]];
	qintA3[i][j][l]=qB3[i][j][l][(int)Ns[6]];
	qintA4[i][j][l]=qB4[i][j][l][(int)Ns[7]];
      }
    }
  }
 
 
  // Here we will solve the diffusion question
  solveModDiffEqn_FFT(qA1,w[0],qintA1,ds,(int)Ns[0],1,k_vector,dxyz);
  solveModDiffEqn_FFT(qA2,w[1],qintA2,ds,(int)Ns[1],1,k_vector,dxyz);
  solveModDiffEqn_FFT(qA3,w[2],qintA3,ds,(int)Ns[2],1,k_vector,dxyz);
  solveModDiffEqn_FFT(qA4,w[3],qintA4,ds,(int)Ns[3],1,k_vector,dxyz);

  
  // The result from the above calculation becomes qdags initial cond
  for(i=0;i<Nx;i++){
    for(j=0;j<Ny;j++){
      for(l=0;l<Nz;l++){
	qintA1[i][j][l]=qA2[i][j][l][(int)Ns[1]]*qA3[i][j][l][(int)Ns[2]]*qA4[i][j][l][(int)Ns[3]];
	qintA2[i][j][l]=qA1[i][j][l][(int)Ns[0]]*qA3[i][j][l][(int)Ns[2]]*qA4[i][j][l][(int)Ns[3]];
	qintA3[i][j][l]=qA1[i][j][l][(int)Ns[0]]*qA2[i][j][l][(int)Ns[1]]*qA4[i][j][l][(int)Ns[3]];
	qintA4[i][j][l]=qA1[i][j][l][(int)Ns[0]]*qA2[i][j][l][(int)Ns[1]]*qA3[i][j][l][(int)Ns[2]];
      }
    }
  }
  

  // Here we will solve the diffusion question
  solveModDiffEqn_FFT(qdagA1,w[0],qintA1,ds,(int)Ns[0],-1,k_vector,dxyz);
  solveModDiffEqn_FFT(qdagA2,w[1],qintA2,ds,(int)Ns[1],-1,k_vector,dxyz);
  solveModDiffEqn_FFT(qdagA3,w[2],qintA3,ds,(int)Ns[2],-1,k_vector,dxyz);
  solveModDiffEqn_FFT(qdagA4,w[3],qintA4,ds,(int)Ns[3],-1,k_vector,dxyz);
  

  // The result from the above calculation becomes qdags initial cond
  for(i=0;i<Nx;i++){
    for(j=0;j<Ny;j++){
      for(l=0;l<Nz;l++){
	qintB1[i][j][l]=qdagA1[i][j][l][(int)Ns[0]];
	qintB2[i][j][l]=qdagA2[i][j][l][(int)Ns[1]];
	qintB3[i][j][l]=qdagA3[i][j][l][(int)Ns[2]];
	qintB4[i][j][l]=qdagA4[i][j][l][(int)Ns[3]];
      }
    }
  }
  

  // Here we will solve the diffusion question
  solveModDiffEqn_FFT(qdagB1,w[4],qintB1,ds,(int)Ns[4],-1,k_vector,dxyz);
  solveModDiffEqn_FFT(qdagB2,w[5],qintB2,ds,(int)Ns[5],-1,k_vector,dxyz);
  solveModDiffEqn_FFT(qdagB3,w[6],qintB3,ds,(int)Ns[6],-1,k_vector,dxyz);
  solveModDiffEqn_FFT(qdagB4,w[7],qintB4,ds,(int)Ns[7],-1,k_vector,dxyz);
  
  //*******************************************************************************************************************

  //std::cout<<"+++++"<<std::endl;
  // Here we are doing the sum to get the single chain partition function
  Q=0.0;
  for(i=0;i<Nx;i++){
    for(j=0;j<Ny;j++){
      for(l=0;l<Nz;l++){
	Q+=(qA1[i][j][l][(int)Ns[0]]*qA2[i][j][l][(int)Ns[1]]*qA3[i][j][l][(int)Ns[2]]*qA4[i][j][l][(int)Ns[3]])*dxyz[0]*dxyz[1]*dxyz[2];
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
	phi[1][i][j][l]=0.0;  //A2
	phi[2][i][j][l]=0.0;  //A3
	phi[3][i][j][l]=0.0;  //A4

	phi[4][i][j][l]=0.0;  //B1
	phi[5][i][j][l]=0.0;  //B2
	phi[6][i][j][l]=0.0;  //B3
	phi[7][i][j][l]=0.0;  //B4

	//A1
	for(s=0;s<(Ns[0]+1);s++){
	  if(s==0 || s==(int)Ns[0]){
	    phi[0][i][j][l]+=0.5*qA1[i][j][l][s]*qdagA1[i][j][l][(int)Ns[0]-s]*ds;
	  }else{
	    phi[0][i][j][l]+=qA1[i][j][l][s]*qdagA1[i][j][l][(int)Ns[0]-s]*ds;
	  }
	}

	//A2
	for(s=0;s<(Ns[1]+1);s++){
	  if(s==0 || s==(int)Ns[1]){
	    phi[1][i][j][l]+=0.5*qA2[i][j][l][s]*qdagA2[i][j][l][(int)Ns[1]-s]*ds;
	  }else{
	    phi[1][i][j][l]+=qA2[i][j][l][s]*qdagA2[i][j][l][(int)Ns[1]-s]*ds;
	  }
	}

	//A3
	for(s=0;s<(Ns[2]+1);s++){
	  if(s==0 || s==(int)Ns[2]){
	    phi[2][i][j][l]+=0.5*qA3[i][j][l][s]*qdagA3[i][j][l][(int)Ns[2]-s]*ds;
	  }else{
	    phi[2][i][j][l]+=qA3[i][j][l][s]*qdagA3[i][j][l][(int)Ns[2]-s]*ds;
	  }
	}

	//A4
	for(s=0;s<(Ns[3]+1);s++){
	  if(s==0 || s==(int)Ns[3]){
	    phi[3][i][j][l]+=0.5*qA4[i][j][l][s]*qdagA4[i][j][l][(int)Ns[3]-s]*ds;
	  }else{
	    phi[3][i][j][l]+=qA4[i][j][l][s]*qdagA4[i][j][l][(int)Ns[3]-s]*ds;
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

	//B2
	for(s=0;s<(Ns[5]+1);s++){
	  if(s==0 || s==(int)Ns[5]){
	    phi[5][i][j][l]+=0.5*qB2[i][j][l][s]*qdagB2[i][j][l][(int)Ns[5]-s]*ds;
	  }else{
	    phi[5][i][j][l]+=qB2[i][j][l][s]*qdagB2[i][j][l][(int)Ns[5]-s]*ds;
	  }
	}

	//B3
	for(s=0;s<(Ns[6]+1);s++){
	  if(s==0 || s==(int)Ns[6]){
	    phi[6][i][j][l]+=0.5*qB3[i][j][l][s]*qdagB3[i][j][l][(int)Ns[6]-s]*ds;
	  }else{
	    phi[6][i][j][l]+=qB3[i][j][l][s]*qdagB3[i][j][l][(int)Ns[6]-s]*ds;
	  }
	}

	//B4
	for(s=0;s<(Ns[7]+1);s++){
	  if(s==0 || s==(int)Ns[7]){
	    phi[7][i][j][l]+=0.5*qB4[i][j][l][s]*qdagB4[i][j][l][(int)Ns[7]-s]*ds;
	  }else{
	    phi[7][i][j][l]+=qB4[i][j][l][s]*qdagB4[i][j][l][(int)Ns[7]-s]*ds;
	  }
	}

	phi[0][i][j][l]*=(1.0/Q);
	phi[1][i][j][l]*=(1.0/Q);
	phi[2][i][j][l]*=(1.0/Q);
	phi[3][i][j][l]*=(1.0/Q);
	phi[4][i][j][l]*=(1.0/Q);
	phi[5][i][j][l]*=(1.0/Q);
	phi[6][i][j][l]*=(1.0/Q);
	phi[7][i][j][l]*=(1.0/Q);

      }
    }
  }
  
  return Q;


};
