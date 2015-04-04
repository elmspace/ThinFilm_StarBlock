double ConcE(double ****phi){

  int         i,j,k;
 
  for(i=0;i<Nx;i++){
    for(j=0;j<Ny;j++){
      for(k=0;k<Nz;k++){
	
	if(forceE>0.5){
	  
	  phi[8][i][j][k]=Cube(i,j,k,3);

	}else{

	  phi[8][i][j][k]=0.0;

	}

      }
    }
  }
  
  return 0.0;

};
