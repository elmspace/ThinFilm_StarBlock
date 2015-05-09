double ConcHS(double ****phi,double ****w,double *Ns,double ds,double ***k_vector,double *dxyz){

  int         i,j,l,s;
  double      Q;

  double      ****qHS;
  double      ***qintHS;


  qHS=create_4d_double_array(Nx,Ny,Nz,((int)Ns[9]+1),"qHS");

  qintHS=create_3d_double_array(Nx,Ny,Nz,"qintHS");

  //+++++++++++++++++++++++++++++++++++++++Forward++++++++++++++++++++++++++++++++++++++++++++++++++
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  for(i=0;i<Nx;i++){
    for(j=0;j<Ny;j++){
      for(l=0;l<Nz;l++){
	qintHS[i][j][l]=1.0;
      }
    }
  }
  solveModDiffEqn_FFT(qHS,w[9],qintHS,ds,(int)Ns[9],1,k_vector,dxyz);
 

  //++++++++++++++++++++++++++++++++++++++Single Chain Partition Function+++++++++++++++++++++++++++
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  Q=0.0;
  for(i=0;i<Nx;i++){
    for(j=0;j<Ny;j++){
      for(l=0;l<Nz;l++){
	Q+=(qHS[i][j][l][(int)Ns[9]])*dxyz[0]*dxyz[1]*dxyz[2];
      }
    }
  }
  // Normalizing with respect to the volume of the box
  Q/=((dxyz[0]*Nx)*(dxyz[1]*Ny)*(dxyz[2]*Nz));

  
  // Here we do the concentration calculation
  for(i=0;i<Nx;i++){
    for(j=0;j<Ny;j++){
      for(l=0;l<Nz;l++){

	phi[9][i][j][l]=0.0;  //HS

	//HS
	for(s=0;s<(Ns[9]+1);s++){
	  if(s==0 || s==(int)Ns[9]){
	    phi[9][i][j][l]+=0.5*qHS[i][j][l][s]*qHS[i][j][l][(int)Ns[9]-s]*ds;
	  }else{
	    phi[9][i][j][l]+=qHS[i][j][l][s]*qHS[i][j][l][(int)Ns[9]-s]*ds;
	  }
	}

	phi[9][i][j][l]*=(pSubAve/(Q*kappa_HA));

      }
    }
  }
  
  //clearing the memory
  destroy_4d_double_array(qHS);
  destroy_3d_double_array(qintHS);



  return Q;


};
