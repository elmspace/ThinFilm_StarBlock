double ConcD(double ****phi){

  int         i,j,k;
  
  for(i=0;i<Nx;i++){
    for(j=0;j<Ny;j++){
      for(k=0;k<Nz;k++){
	
	if(forceD>0.5){
	  
	  phi[9][i][j][k]=Cube(i,j,k,2);
	  
	}else{

	  phi[9][i][j][k]=0.0;

	}

      }
    }
  }

  return 0.0;

};
