void parametersAB(double *chi,double *f,double &ds,double *Ns,double *dxyz,double **chiMatrix, double ****h){
  
  int i,j,k;
  int Ds;
  double xAB;

  // Minimize with respect to box size (yes=1, No=0)
  box_min=1;

  // Degree of polymerization (Each arm of the star is 100)
  Ns[0]=35;  // A1
  Ns[1]=Ns[0];  // A2
  Ns[2]=Ns[0];  // A3
  Ns[3]=Ns[0];  // A4
  Ns[4]=100-Ns[0];  // B1
  Ns[5]=100-Ns[0];  // B2
  Ns[6]=100-Ns[0];  // B3
  Ns[7]=100-Ns[0];  // B4
 
  // Total length
  Ds=Ns[0]+Ns[1]+Ns[2]+Ns[3]+Ns[4]+Ns[5]+Ns[6]+Ns[7];
 
  // Setting the generic chi parameters
  xAB=(0.14)*Ds;
  h_AAir=(0.08)*Ds;
  h_BAir=(0.08)*Ds;
  h_ASub=(0.08)*Ds;
  h_BSub=(0.06)*Ds;

 
  // 0 read
  // 1 make
  // 2 random
  Iomega=1;
  
  // For 4Arm L=2.2
  Lx=5.0;
  Ly=5.0;
  Lz=5.0;
  dxyz[0]=Lx/Nx;
  dxyz[1]=Ly/Ny;
  dxyz[2]=Lz/Nz;
 
  f[0]=Ns[0]/Ds;  // fA1
  f[1]=Ns[1]/Ds;  // fA2
  f[2]=Ns[2]/Ds;  // fA3
  f[3]=Ns[3]/Ds;  // fA4
  f[4]=Ns[4]/Ds;  // fB1
  f[5]=Ns[5]/Ds;  // fB2
  f[6]=Ns[6]/Ds;  // fB3
  f[7]=Ns[7]/Ds;  // fB4

  // Setting up the individual chi values
  chi[0]=0.0;  
  chi[1]=xAB;  
  //++++++++++++++++++++++++++++++++++++++++++++++++

  // Average Concentrations
  pA1ave=f[0]; //A1
  pA2ave=f[1]; //A2
  pA3ave=f[2]; //A3
  pA4ave=f[3]; //A4

  pB1ave=f[4]; //B1
  pB2ave=f[5]; //B2
  pB3ave=f[6]; //B3
  pB4ave=f[7]; //B4

  ds=1.0/Ds;
  

  // Setting up the chi matrix in this case 2 by 2
  
  chiMatrix[0][0]=0.0;  
  chiMatrix[0][1]=0.0;
  chiMatrix[0][2]=0.0;
  chiMatrix[0][3]=0.0;
  chiMatrix[0][4]=chi[0];
  chiMatrix[0][5]=chi[0];
  chiMatrix[0][6]=chi[0];
  chiMatrix[0][7]=chi[0];
  
  chiMatrix[1][0]=0.0;  
  chiMatrix[1][1]=0.0;
  chiMatrix[1][2]=0.0;
  chiMatrix[1][3]=0.0;
  chiMatrix[1][4]=chi[0];
  chiMatrix[1][5]=chi[0];
  chiMatrix[1][6]=chi[0];
  chiMatrix[1][7]=chi[0];

  chiMatrix[2][0]=0.0;  
  chiMatrix[2][1]=0.0;
  chiMatrix[2][2]=0.0;
  chiMatrix[2][3]=0.0;
  chiMatrix[2][4]=chi[0];
  chiMatrix[2][5]=chi[0];
  chiMatrix[2][6]=chi[0];
  chiMatrix[2][7]=chi[0];

  chiMatrix[3][0]=0.0;  
  chiMatrix[3][1]=0.0;
  chiMatrix[3][2]=0.0;
  chiMatrix[3][3]=0.0;
  chiMatrix[3][4]=chi[0];
  chiMatrix[3][5]=chi[0];
  chiMatrix[3][6]=chi[0];
  chiMatrix[3][7]=chi[0];

  chiMatrix[4][0]=chi[0];  
  chiMatrix[4][1]=chi[0];
  chiMatrix[4][2]=chi[0];
  chiMatrix[4][3]=chi[0];
  chiMatrix[4][4]=0.0;
  chiMatrix[4][5]=0.0;
  chiMatrix[4][6]=0.0;
  chiMatrix[4][7]=0.0;

  chiMatrix[5][0]=chi[0];  
  chiMatrix[5][1]=chi[0];
  chiMatrix[5][2]=chi[0];
  chiMatrix[5][3]=chi[0];
  chiMatrix[5][4]=0.0;
  chiMatrix[5][5]=0.0;
  chiMatrix[5][6]=0.0;
  chiMatrix[5][7]=0.0;

  chiMatrix[6][0]=chi[0];  
  chiMatrix[6][1]=chi[0];
  chiMatrix[6][2]=chi[0];
  chiMatrix[6][3]=chi[0];
  chiMatrix[6][4]=0.0;
  chiMatrix[6][5]=0.0;
  chiMatrix[6][6]=0.0;
  chiMatrix[6][7]=0.0;

  chiMatrix[7][0]=chi[0];  
  chiMatrix[7][1]=chi[0];
  chiMatrix[7][2]=chi[0];
  chiMatrix[7][3]=chi[0];
  chiMatrix[7][4]=0.0;
  chiMatrix[7][5]=0.0;
  chiMatrix[7][6]=0.0;
  chiMatrix[7][7]=0.0;

 
  // Setting up the surface field
  for(i=0;i<Nx;i++){
    for(j=0;j<Ny;j++){
      for(k=0;k<Nz;k++){
	
	if(k==0){ // k=0 is the substrate surface
	  h[0][i][j][k]=h_ASub;
	  h[1][i][j][k]=h_ASub;
	  h[2][i][j][k]=h_ASub;
	  h[3][i][j][k]=h_ASub;
	  h[4][i][j][k]=h_BSub;
	  h[5][i][j][k]=h_BSub;
	  h[6][i][j][k]=h_BSub;
	  h[7][i][j][k]=h_BSub;
	}else if(k==(Nz-1)){ // k=Nz-1 is the air interface
	  h[0][i][j][k]=h_AAir;
	  h[1][i][j][k]=h_AAir;
	  h[2][i][j][k]=h_AAir;
	  h[3][i][j][k]=h_AAir;
	  h[4][i][j][k]=h_BAir;
	  h[5][i][j][k]=h_BAir;
	  h[6][i][j][k]=h_BAir;
	  h[7][i][j][k]=h_BAir;
	}else{ // No surface interaction in bulk
	  h[0][i][j][k]=0.0;
	  h[1][i][j][k]=0.0;
	  h[2][i][j][k]=0.0;
	  h[3][i][j][k]=0.0;
	  h[4][i][j][k]=0.0;
	  h[5][i][j][k]=0.0;
	  h[6][i][j][k]=0.0;
	  h[7][i][j][k]=0.0;
	}

      }
    }
  }
 


};