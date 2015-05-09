void Incomp(double ***eta, double ****phi, double ***delphi){

  int     i,j,k;
  int     chain; 
  double  ptot;

  ptot=0.0;

  for(i=0;i<Nx;i++){
    for(j=0;j<Ny;j++){
      for(k=0;k<Nz;k++){
	
	ptot=0.0;
	delphi[i][j][k]=0.0;
    
	for(chain=0;chain<ChainType;chain++){
	  ptot+=phi[chain][i][j][k];
	}

	delphi[i][j][k]=1.0-ptot;
	eta[i][j][k]-=delphi[i][j][k];

      }
    }
  }
 


};
