void omega(double ****w){

  int i,j,k;
  int ii,jj,kk;
 

  if(Iomega==0){ // Read
    
    std::ifstream infile;
    infile.open("./Hex_Horz_32_32_32.read");
    
    for(i=0;i<Nx;i++){
      for(j=0;j<Ny;j++){
	for(k=0;k<Nz;k++){
	  infile >> ii >> jj >> kk >> w[0][i][j][k] >> w[1][i][j][k] >> w[2][i][j][k] >> w[3][i][j][k] >> w[4][i][j][k] >> w[5][i][j][k] >> w[6][i][j][k]>> w[7][i][j][k]; 
	}
      }
    }
  
    infile.close();

  }else if(Iomega==1){

    if(LAM==1){
      for(i=0;i<Nx;i++){
	for(j=0;j<Ny;j++){
	  for(k=0;k<Nz;k++){
	    w[0][i][j][k]=-1.0*cos(2.0*Pi*i/Nx); //A1
	    w[1][i][j][k]=-1.0*cos(2.0*Pi*i/Nx); //A2
	    w[2][i][j][k]=-1.0*cos(2.0*Pi*i/Nx); //A3
	    w[3][i][j][k]=-1.0*cos(2.0*Pi*i/Nx); //A4
	  }
	}
      }
    }

    if(HEX==1){
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



