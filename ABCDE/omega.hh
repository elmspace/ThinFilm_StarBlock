void omega(double ****w){

  int i,j,k;
  int ii,jj,kk;
 
  if(Round==1){
    if(Iomega==0){ // Read
      
      std::ifstream infile;
      if(Bulk_Calc==1){
	if(LAM==1){
	  infile.open("./OMEGA/READ_1_Period/Omega_LAM_20_20_20.read");
	}else if(HEX==1){
	  infile.open("./OMEGA/READ_1_Period/Omega_HEX_20_20_20.read");
	}else if(BCC==1){
	  infile.open("./OMEGA/READ_1_Period/Omega_BCC_20_20_20.read");
	}
      }else{
	if(LAM==1){
	  if(HOR==1){infile.open("./OMEGA/READ/Omega_LAM_40_1_40_HOR.read");}
	  if(VER==1){infile.open("./OMEGA/READ/Omega_LAM_40_1_40_VER.read");}
	  if(VER_2==1){infile.open("./OMEGA/READ/Omega_LAM_40_1_40_VER_2.read");}
	}else if(HEX==1){
	  if(HOR==1){infile.open("./OMEGA/READ/Omega_HEX_20_20_20_HOR.read");}
	  if(VER==1){infile.open("./OMEGA/READ/Omega_HEX_20_20_20_VER.read");}
	  if(VER_2==1){infile.open("./OMEGA/READ/Omega_HEX_20_20_20_VER_2.read");}
	}
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
	      for(k=3;k<(Nz-3);k++){
		w[0][i][j][k]=50.0*cos(Numb_of_Periods*2.0*Pi*k/Nz); //A1
		w[1][i][j][k]=w[0][i][j][k]; //A2
		w[2][i][j][k]=w[0][i][j][k]; //A3
		w[3][i][j][k]=w[0][i][j][k]; //A4
	      }
	    }
	  }
	} else if(VER==1){
	  for(i=0;i<Nx;i++){
	    for(j=0;j<Ny;j++){
	      for(k=2;k<(Nz-2);k++){
		w[0][i][j][k]=-10.0*cos(Numb_of_Periods*2.0*Pi*i/Nx); //A1
		w[1][i][j][k]=w[0][i][j][k]; //A2
		w[2][i][j][k]=w[0][i][j][k]; //A3
		w[3][i][j][k]=w[0][i][j][k]; //A4
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
	      for(k=2;k<(Nz-2);k++){
		w[0][i][j][k]=-10.0*cos(2.0*Pi*i/Nx)*cos(((Numb_of_Periods)*2.0*Pi*k/Nz)); //A1
		w[1][i][j][k]=w[0][i][j][k]; //A2
		w[2][i][j][k]=w[0][i][j][k]; //A3
		w[3][i][j][k]=w[0][i][j][k]; //A4
	      }
	    }
	  }
	} else if(VER==1){
	  for(i=0;i<Nx;i++){
	    for(j=0;j<Ny;j++){
	      for(k=2;k<(Nz-2);k++){
		w[0][i][j][k]=-10.0*cos(2.0*Pi*i/Nx)*cos(2.0*Pi*j/Ny); //A1
		w[1][i][j][k]=w[0][i][j][k]; //A2
		w[2][i][j][k]=w[0][i][j][k]; //A3
		w[3][i][j][k]=w[0][i][j][k]; //A4
	      }
	    }
	  }
	}
      }
       // BCC ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      if(BCC==1){

	i=0;
	j=0;
	k=0;
	w[0][i][j][k]=-10.0; //A1
	w[1][i][j][k]=-10.0; //A2
	w[2][i][j][k]=-10.0; //A3
	w[3][i][j][k]=-10.0; //A4

	i=Nx-1;
	j=0;
	k=0;
	w[0][i][j][k]=-10.0; //A1
	w[1][i][j][k]=-10.0; //A2
	w[2][i][j][k]=-10.0; //A3
	w[3][i][j][k]=-10.0; //A4

	i=0;
	j=Ny-1;
	k=0;
	w[0][i][j][k]=-10.0; //A1
	w[1][i][j][k]=-10.0; //A2
	w[2][i][j][k]=-10.0; //A3
	w[3][i][j][k]=-10.0; //A4

	i=Nx-1;
	j=Ny-1;
	k=0;
	w[0][i][j][k]=-1.0; //A1
	w[1][i][j][k]=-1.0; //A2
	w[2][i][j][k]=-1.0; //A3
	w[3][i][j][k]=-1.0; //A4

	i=0;
	j=0;
	k=Nz-1;
	w[0][i][j][k]=-1.0; //A1
	w[1][i][j][k]=-1.0; //A2
	w[2][i][j][k]=-1.0; //A3
	w[3][i][j][k]=-1.0; //A4

	i=Nx-1;
	j=0;
	k=Nz-1;
	w[0][i][j][k]=-1.0; //A1
	w[1][i][j][k]=-1.0; //A2
	w[2][i][j][k]=-1.0; //A3
	w[3][i][j][k]=-1.0; //A4

	i=0;
	j=Ny-1;
	k=Nz-1;
	w[0][i][j][k]=-1.0; //A1
	w[1][i][j][k]=-1.0; //A2
	w[2][i][j][k]=-1.0; //A3
	w[3][i][j][k]=-1.0; //A4

	i=Nx-1;
	j=Ny-1;
	k=Nz-1;
	w[0][i][j][k]=-1.0; //A1
	w[1][i][j][k]=-1.0; //A2
	w[2][i][j][k]=-1.0; //A3
	w[3][i][j][k]=-1.0; //A4

	i=Nx/2;
	j=Ny/2;
	k=Nz/2;
	w[0][i][j][k]=-1.0; //A1
	w[1][i][j][k]=-1.0; //A2
	w[2][i][j][k]=-1.0; //A3
	w[3][i][j][k]=-1.0; //A4

	
      }
      
    }else if(Iomega==2){ // Random
      
      for(i=0;i<Nx;i++){
	for(j=0;j<Ny;j++){
	  for(k=0;k<Nz;k++){
	    
	    w[0][i][j][k]=-1.1*(drand48()-0.50); //A1
	    w[1][i][j][k]=-1.1*(drand48()-0.50); //A2
	    w[2][i][j][k]=-1.1*(drand48()-0.50); //A3
	    w[3][i][j][k]=-1.1*(drand48()-0.50); //A4
	    w[4][i][j][k]=-1.1*(drand48()-0.50); //B1
	    w[5][i][j][k]=-1.1*(drand48()-0.50); //B2
	    w[6][i][j][k]=-1.1*(drand48()-0.50); //B3
	    w[7][i][j][k]=-1.1*(drand48()-0.50); //B4
	    
	  }
	}
      }
      
    }
  }else{
    // Will be using the previous omega files, since the code is looping.
  }
 
};



