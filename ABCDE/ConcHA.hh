double ConcHA(double ****phi,double ****w,double *Ns,double ds,double ***k_vector,double *dxyz){

  int         i,j,l,s;
  double      Q;

  double      ****qHA;
  double      ***qintHA;


  qHA=create_4d_double_array(Nx,Ny,Nz,((int)Ns[8]+1),"qHA");

  qintHA=create_3d_double_array(Nx,Ny,Nz,"qintHA");

  //+++++++++++++++++++++++++++++++++++++++Forward++++++++++++++++++++++++++++++++++++++++++++++++++
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  for(i=0;i<Nx;i++){
    for(j=0;j<Ny;j++){
      for(l=0;l<Nz;l++){
	qintHA[i][j][l]=1.0;
      }
    }
  }
  solveModDiffEqn_FFT(qHA,w[8],qintHA,ds,(int)Ns[8],1,k_vector,dxyz);
 

  //++++++++++++++++++++++++++++++++++++++Single Chain Partition Function+++++++++++++++++++++++++++
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  Q=0.0;
  for(i=0;i<Nx;i++){
    for(j=0;j<Ny;j++){
      for(l=0;l<Nz;l++){
	Q+=(qHA[i][j][l][(int)Ns[8]])*dxyz[0]*dxyz[1]*dxyz[2];
      }
    }
  }
  // Normalizing with respect to the volume of the box
  Q/=((dxyz[0]*Nx)*(dxyz[1]*Ny)*(dxyz[2]*Nz));

  
  // Here we do the concentration calculation
  for(i=0;i<Nx;i++){
    for(j=0;j<Ny;j++){
      for(l=0;l<Nz;l++){

	phi[8][i][j][l]=0.0;  //HA

	//HA
	for(s=0;s<(Ns[8]+1);s++){
	  if(s==0 || s==(int)Ns[8]){
	    phi[8][i][j][l]+=0.5*qHA[i][j][l][s]*qHA[i][j][l][(int)Ns[8]-s]*ds;
	  }else{
	    phi[8][i][j][l]+=qHA[i][j][l][s]*qHA[i][j][l][(int)Ns[8]-s]*ds;
	  }
	}

	phi[8][i][j][l]*=(pAirAve/(Q*kappa_HA));

      }
    }
  }
  
  //clearing the memory
  destroy_4d_double_array(qHA);
  destroy_3d_double_array(qintHA);



  return Q;


};
