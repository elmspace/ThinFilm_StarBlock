double Set_h_function(double ****h, double *dxyz){

  int   i,j,k,ll;
  double Phi_Multi;

  Phi_Multi=0.0;
  // Setting up the surface field
  for(i=0;i<Nx;i++){
    for(j=0;j<Ny;j++){
      for(k=0;k<Nz;k++){
	for(ll=0;ll<ChainType;ll++){
	  
	  if(k*dxyz[2]<=epsilon){ // k=0 is the substrate surface
	    if(ll==9){
	      h[ll][i][j][k]=(cos(0.5*Pi*k*dxyz[2]/epsilon));
	    }else{
	      h[ll][i][j][k]=0.0;
	    }
	  }else if((k*dxyz[2]>=(FilmThickness-epsilon))){ // k=Nz-1 is the air interface
	    if(ll==8){
	      h[ll][i][j][k]=(cos(0.5*Pi*(FilmThickness-k*dxyz[2])/epsilon));
	    }else{
	      h[ll][i][j][k]=0.0;
	    }
	  }else{ // No surface interaction in bulk
	    h[ll][i][j][k]=0.0;
	  }
	  Phi_Multi+=h[ll][i][j][k]*dxyz[0]*dxyz[1]*dxyz[2];
	}
      }
    }
  }
  Phi_Multi/=Nx*Ny*Nz*dxyz[0]*dxyz[1]*dxyz[2];

  //std::cout<<Phi_Multi<<std::endl;
  //std::cin>>i;
  
  return (1.0-Phi_Multi);
  
};
