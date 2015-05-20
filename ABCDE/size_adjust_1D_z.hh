double size_adjust_1D_z(double ****w, double ****phi, double ***eta, double *Ns, double ds, double ***k_vector, double *chi, double *dxyz, double **chiMatrix){

  int     k,l,ii,jj,ll,kk;
  double  delz;
  double  box, *dxyz_temp;
  double  *box_z,*box_fE;
  double  ****w_temp;

  w_temp=create_4d_double_array(ChainType,Nx,Ny,Nz,"w_temp");
  box_z=create_1d_double_array(2,"box_z");
  box_fE=create_1d_double_array(2,"box_fE");
  dxyz_temp=create_1d_double_array(3,"dxyz_temp");

  delz=0.2/Nz;
  // std::cout<<"++++++++++++++++++++++++++++++++++++++++++++++++++++"<<std::endl;
  l=0;
  for(k=-1;k<2;k++){
    
    if(k==0){
      
    }else{
 
      box_z[l]=k*delz;
      
      dxyz_temp[0]=dxyz[0];
      dxyz_temp[1]=dxyz[1];
      dxyz_temp[2]=dxyz[2];
      
      for(ii=0;ii<Nx;ii++){
	for(jj=0;jj<Ny;jj++){
	  for(kk=0;kk<Nz;kk++){
	    for(ll=0;ll<ChainType;ll++){
	      w_temp[ll][ii][jj][kk]=w[ll][ii][jj][kk];
	    }
	  }
	}
      }
      
      dxyz_temp[2]+=box_z[l];
      
      box_fE[l]=FreeEnergy_Box_Edition(w_temp,phi,eta,Ns,ds,k_vector,chi,dxyz_temp,chiMatrix); 
      l=l+1;
      
    }
    
  }
  //std::cout<<"++++++++++++++++++++++++++++++++++++++++++++++++++++"<<std::endl;

  if(box_fE[0]<box_fE[1]){
    dxyz[2]=dxyz[2]-delz;
  }else if(box_fE[1]<box_fE[0]){
    dxyz[2]=dxyz[2]+delz;
  }
    
  box=1;

  destroy_1d_double_array(dxyz_temp);
  destroy_1d_double_array(box_z);
  destroy_1d_double_array(box_fE);
  destroy_4d_double_array(w_temp);

  return box;
  
};

