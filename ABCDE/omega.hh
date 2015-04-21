void omega(double ****w){

  int i,j,k;
  int ii,jj,kk;
 

  if(Iomega==0){ // Read
    
    std::ifstream infile;
    if(LAM==1){
      if(HOR==1){infile.open("./OMEGA/READ/Omega_LAM_20_20_20_HOR.read");}
      if(VER==1){infile.open("./OMEGA/READ/Omega_LAM_20_20_20_VER.read");}
    }else if(HEX==1){
      if(HOR==1){infile.open("./OMEGA/READ/Omega_HEX_20_20_20_HOR.read");}
      if(VER==1){infile.open("./OMEGA/READ/Omega_HEX_20_20_20_VER.read");}
    }
    
    for(i=0;i<Nx;i++){
      for(j=0;j<Ny;j++){
	for(k=0;k<Nz;k++){
	  infile >> ii >> jj >> kk >> w[0][i][j][k] >> w[1][i][j][k] >> w[2][i][j][k] >> w[3][i][j][k] >> w[4][i][j][k] >> w[5][i][j][k] >> w[6][i][j][k]>> w[7][i][j][k]; 
	}
      }
    }
  
    infile.close();

  }else if(Iomega==1){
    // LAM ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if(LAM==1){
      if(HOR==1){
	for(i=0;i<Nx;i++){
	  for(j=0;j<Ny;j++){
	    for(k=0;k<Nz;k++){
	      w[0][i][j][k]=-1.0*cos(Numb_of_Periods*2.0*Pi*k/Nz); //A1
	      w[1][i][j][k]=-1.0*cos(Numb_of_Periods*2.0*Pi*k/Nz); //A2
	      w[2][i][j][k]=-1.0*cos(Numb_of_Periods*2.0*Pi*k/Nz); //A3
	      w[3][i][j][k]=-1.0*cos(Numb_of_Periods*2.0*Pi*k/Nz); //A4
	    }
	  }
	}
      } else if(VER==1){
	for(i=0;i<Nx;i++){
	  for(j=0;j<Ny;j++){
	    for(k=0;k<Nz;k++){
	      w[0][i][j][k]=-1.0*cos(4.0*Pi*i/Nx); //A1
	      w[1][i][j][k]=-1.0*cos(4.0*Pi*i/Nx); //A2
	      w[2][i][j][k]=-1.0*cos(4.0*Pi*i/Nx); //A3
	      w[3][i][j][k]=-1.0*cos(4.0*Pi*i/Nx); //A4
	    }
	  }
	}
      }
    }
    // HEX ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if(HEX==1){
      if(HOR==1){
	for(i=0;i<Nx;i++){
	  for(j=0;j<Ny;j++){
	    for(k=0;k<Nz;k++){
	      w[0][i][j][k]=-1.0*cos(2.0*Pi*i/Nx)*cos(Numb_of_Periods*2.0*Pi*k/Nz); //A1
	      w[1][i][j][k]=-1.0*cos(2.0*Pi*i/Nx)*cos(Numb_of_Periods*2.0*Pi*k/Nz); //A2
	      w[2][i][j][k]=-1.0*cos(2.0*Pi*i/Nx)*cos(Numb_of_Periods*2.0*Pi*k/Nz); //A3
	      w[3][i][j][k]=-1.0*cos(2.0*Pi*i/Nx)*cos(Numb_of_Periods*2.0*Pi*k/Nz); //A4
	    }
	  }
	}
      } else if(VER==1){
	for(i=0;i<Nx;i++){
	  for(j=0;j<Ny;j++){
	    for(k=0;k<Nz;k++){
	      w[0][i][j][k]=-1.0*cos(2.0*Pi*i/Nx)*cos(2.0*Pi*j/Ny); //A1
	      w[1][i][j][k]=-1.0*cos(2.0*Pi*i/Nx)*cos(2.0*Pi*j/Ny); //A2
	      w[2][i][j][k]=-1.0*cos(2.0*Pi*i/Nx)*cos(2.0*Pi*j/Ny); //A3
	      w[3][i][j][k]=-1.0*cos(2.0*Pi*i/Nx)*cos(2.0*Pi*j/Ny); //A4
	    }
	  }
	}
      }
    }

  }else if(Iomega==2){ // Random
    
    for(i=0;i<Nx;i++){
      for(j=0;j<Ny;j++){
	for(k=0;k<Nz;k++){

	  w[0][i][j][k]=-5.0*(drand48()-0.50); //A1
	  w[1][i][j][k]=-5.0*(drand48()-0.50); //A2
	  w[2][i][j][k]=-5.0*(drand48()-0.50); //A3
	  w[3][i][j][k]=-5.0*(drand48()-0.50); //A4
	  w[4][i][j][k]=-5.0*(drand48()-0.50); //B1
	  w[5][i][j][k]=-5.0*(drand48()-0.50); //B2
	  w[6][i][j][k]=-5.0*(drand48()-0.50); //B3
	  w[7][i][j][k]=-5.0*(drand48()-0.50); //B4

	}
      }
    }
  


  }
 
};



