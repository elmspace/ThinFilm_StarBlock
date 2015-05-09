void set_PHI_0(double ***PHI_0,double *dxyz){

  int i,j,k;
  
  PHI_0_tot=0.0;
  for(i=0;i<Nx;i++){
    for(j=0;j<Ny;j++){
      for(k=0;k<Nz;k++){
	
	if((k*dxyz[2]>=0.0)&&(k*dxyz[2]<=epsilon)){// Bottom Surface
	  PHI_0[i][j][k]=0.5*(1.0-cos(Pi*k*dxyz[2]/epsilon));
	}else if((k*dxyz[2]>=(FilmThickness-epsilon))&&(k*dxyz[2]<=(FilmThickness+dxyz[2]))){// Top Surface
	  PHI_0[i][j][k]=0.5*(1.0-cos(Pi*(FilmThickness-k*dxyz[2])/epsilon));
	}else{// Middle part
	  PHI_0[i][j][k]=1.0;
	}

	PHI_0_tot+=PHI_0[i][j][k]*(dxyz[0]*dxyz[1]*dxyz[2]);
     
      } 
    }
  }
  PHI_0_tot/=((dxyz[0]*Nx)*(dxyz[1]*Ny)*(dxyz[2]*Nz));
  PHI_0_tot=(1.0-PHI_0_tot);
  
  
};



