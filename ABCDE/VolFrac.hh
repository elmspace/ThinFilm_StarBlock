double VolFracD(double *dxyz){

  double concenD;
  int    i,j,l;

  concenD=0.0;
  for(i=0;i<Nx;i++){
    for(j=0;j<Ny;j++){
      for(l=0;l<Nz;l++){

	concenD=concenD+Cube(i,j,l,2)*dxyz[0]*dxyz[1]*dxyz[2];

      }
    }
  }

  
  concenD=concenD/(Nx*Ny*Nz*dxyz[0]*dxyz[1]*dxyz[2]);

  return concenD;

};
double VolFracE(double *dxyz){

  double concenE;
  int    i,j,l;

  concenE=0.0;
  for(i=0;i<Nx;i++){
    for(j=0;j<Ny;j++){
      for(l=0;l<Nz;l++){
	
	concenE=concenE+Cube(i,j,l,3)*dxyz[0]*dxyz[1]*dxyz[2];

      }
    }
  }

  concenE=concenE/(Nx*Ny*Nz*dxyz[0]*dxyz[1]*dxyz[2]);


  return concenE;

};
