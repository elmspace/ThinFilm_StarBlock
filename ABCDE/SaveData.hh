void SaveData(double ****phi,double ***PHI_0, double ****w, double *dxyz){

  int i, j ,k;

  //+++++++++++++++++++++++++++++ This output is setup for the matlab plotting +++++++++++++++++++
  std::ofstream outputFile7("./MATLAB/xyz.dat");
  for (i=0;i<Nx;i++){
    outputFile7<<i*dxyz[0]<<" "<<i*dxyz[1]<<" "<<i*dxyz[2]<<std::endl;
  }
  outputFile7.close();  
  std::ofstream outputFile8("./MATLAB/ABCD.dat");
  for (i=0;i<Nx;i++){
    for(j=0;j<Ny;j++){
      for(k=0;k<Nz;k++){
	outputFile8<<phi[0][i][j][k]+phi[1][i][j][k]+phi[2][i][j][k]+phi[3][i][j][k]<<" "<<phi[4][i][j][k]+phi[5][i][j][k]+phi[6][i][j][k]+phi[7][i][j][k]<<std::endl;
      }
    }
  }
  outputFile8.close();
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
    

  // ________________________________________________ OMEGA

  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // Writting to data files
  std::ofstream outputFile6("./OMEGA/omega.dat");
  for(i=0;i<Nx;i++){
    for(j=0;j<Ny;j++){
      for(k=0;k<Nz;k++){
	outputFile6 <<i<<" "<<j<<" "<<k<< " "<<w[0][i][j][k]<<" "<<w[1][i][j][k]<<" "<<w[2][i][j][k]<<" "<<w[3][i][j][k]<<" "<<w[4][i][j][k]<< " "<<w[5][i][j][k]<<" "<<w[6][i][j][k]<<" "<<w[7][i][j][k]<<std::endl;
      }
    }
  }
  outputFile6.close();
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



  // ________________________________________________ PHI

  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // Writting to data files
  std::ofstream outputFile16("./PHI/phi_xy.dat");
  for(i=0;i<Nx;i++){
    for(j=0;j<Ny;j++){
	outputFile16 <<i*dxyz[0]<<" "<<j*dxyz[1]<<" "<<phi[0][i][j][Nz/2]<<" "<<phi[1][i][j][Nz/2]<<" "<<phi[2][i][j][Nz/2]<<" "<<phi[3][i][j][Nz/2]<<" "<<phi[4][i][j][Nz/2]<< " "<<phi[5][i][j][Nz/2]<<" "<<phi[6][i][j][Nz/2]<<" "<<phi[7][i][j][Nz/2]<<std::endl;
    }
  }
  outputFile16.close();
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // Writting to data files
  std::ofstream outputFile26("./PHI/phi_xz.dat");
  for(i=0;i<Nx;i++){
    for(k=0;k<Nz;k++){
      outputFile26 <<i*dxyz[0]<<" "<<k*dxyz[2]<<" "<<phi[0][i][Ny/2][k]<<" "<<phi[1][i][Ny/2][k]<<" "<<phi[2][i][Ny/2][k]<<" "<<phi[3][i][Ny/2][k]<<" "<<phi[4][i][Ny/2][k]<< " "<<phi[5][i][Ny/2][k]<<" "<<phi[6][i][Ny/2][k]<<" "<<phi[7][i][Ny/2][k]<<std::endl;
    }
  }
  outputFile26.close();
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // Writting to data files
  std::ofstream outputFile36("./PHI/phi_yz.dat");
  for(j=0;j<Ny;j++){
    for(k=0;k<Nz;k++){
	outputFile36 <<j*dxyz[1]<<" "<<k*dxyz[2]<<" "<<phi[0][Nx/2][j][k]<<" "<<phi[1][Nx/2][j][k]<<" "<<phi[2][Nx/2][j][k]<<" "<<phi[3][Nx/2][j][k]<<" "<<phi[4][Nx/2][j][k]<< " "<<phi[5][Nx/2][j][k]<<" "<<phi[6][Nx/2][j][k]<<" "<<phi[7][Nx/2][j][k]<<std::endl;
    }
  }
  outputFile36.close();
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


 //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // Writting to data files
  std::ofstream outputFile46("./PHI/phi_x.dat");
  for(i=0;i<Nx;i++){
    outputFile46 <<i*dxyz[0]<<" "<<phi[0][i][Ny/2][Nz/2]<<" "<<phi[1][i][Ny/2][Nz/2]<<" "<<phi[2][i][Ny/2][Nz/2]<<" "<<phi[3][i][Ny/2][Nz/2]<<" "<<phi[4][i][Ny/2][Nz/2]<< " "<<phi[5][i][Ny/2][Nz/2]<<" "<<phi[6][i][Ny/2][Nz/2]<<" "<<phi[7][i][Ny/2][Nz/2]<<" "<<PHI_0[i][Ny/2][Nz/2]<<std::endl;
  }
  outputFile46.close();
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // Writting to data files
  std::ofstream outputFile47("./PHI/phi_y.dat");
  for(j=0;j<Ny;j++){
    outputFile47 <<j*dxyz[1]<<" "<<phi[0][Nx/2][j][Nz/2]<<" "<<phi[1][Nx/2][j][Nz/2]<<" "<<phi[2][Nx/2][j][Nz/2]<<" "<<phi[3][Nx/2][j][Nz/2]<<" "<<phi[4][Nx/2][j][Nz/2]<< " "<<phi[5][Nx/2][j][Nz/2]<<" "<<phi[6][Nx/2][j][Nz/2]<<" "<<phi[7][Nx/2][j][Nz/2]<<" "<<PHI_0[Nx/2][j][Nz/2]<<std::endl;
  }
  outputFile47.close();
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // Writting to data files
  std::ofstream outputFile48("./PHI/phi_z.dat");
  for(k=0;k<Nz;k++){
    outputFile48 <<k*dxyz[2]<<" "<<phi[0][Nx/2][Ny/2][k]<<" "<<phi[1][Nx/2][Ny/2][k]<<" "<<phi[2][Nx/2][Ny/2][k]<<" "<<phi[3][Nx/2][Ny/2][k]<<" "<<phi[4][Nx/2][Ny/2][k]<< " "<<phi[5][Nx/2][Ny/2][k]<<" "<<phi[6][Nx/2][Ny/2][k]<<" "<<phi[7][Nx/2][Ny/2][k]<<" "<<PHI_0[Nx/2][Ny/2][k]<<std::endl;
  }
  outputFile48.close();
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      
 
};



