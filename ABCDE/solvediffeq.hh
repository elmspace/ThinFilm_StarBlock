//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//                                Solving The Modified Diffusion Equation
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void solveModDiffEqn_FFT(double ****q, double ***w, double ***qint, double ds, int Ns, int sign, double ***k, double *dxyz){
  
  int            i,j,l,s,ss;  // some counters
  unsigned long  ijl; // This is used for the Fourier Transform
  double         ***wds, ***kds;

 

  wds=create_3d_double_array(Nx,Ny,Nz,"wds");
  kds=create_3d_double_array(Nx,Ny,Nz,"kds");

  for(i=0;i<Nx;i++){
    for(j=0;j<Ny;j++){
      for(l=0;l<Nz;l++){
	kds[i][j][l]=exp(-ds*k[i][j][l]);
	wds[i][j][l]=exp(-0.5*ds*w[i][j][l]);
      }
    }
  }
  
  if(sign==1){
   
    for(i=0;i<Nx;i++){
      for(j=0;j<Ny;j++){
	for(l=0;l<Nz;l++){
	  q[i][j][l][0]=qint[i][j][l];
	}
      }
    }
    
    for(s=0;s<((int)Ns);s++){               // Why not to Ns+1 ??
      
      for(i=0;i<Nx;i++){
	for(j=0;j<Ny;j++){
	  for(l=0;l<Nz;l++){
	    ss=l+Nz*(j+Ny*i);
	    input_q[ss]=q[i][j][l][s]*wds[i][j][l];
	  }
	}
      }
      fftw_execute(forward_plan);


      for(i=0;i<Nx;i++){
	for(j=0;j<Ny;j++){
	  for(l=0;l<Nz;l++){
	    ss=l+Nz*(j+Ny*i);
	    transformed_q[ss]*=kds[i][j][l];
	  }
	}
      }
      fftw_execute(inverse_plan);
      
      for(i=0;i<Nx;i++){
	for(j=0;j<Ny;j++){
	  for(l=0;l<Nz;l++){
	    ss=l+Nz*(j+Ny*i);
	    q[i][j][l][s+1]=((final_q[ss]*wds[i][j][l])/(8.0*Nx*Ny*Nz));
	  }
	}
      }
    }

  }else{

    for(i=0;i<Nx;i++){
      for(j=0;j<Ny;j++){
	for(l=0;l<Nz;l++){
	  q[i][j][l][0]=qint[i][j][l];
	}
      }
    }

    for(s=0;s<(Ns);s++){                   // why start from Ns and not 0 again?
      
      for(i=0;i<Nx;i++){
	for(j=0;j<Ny;j++){
	  for(l=0;l<Nz;l++){
	    ss=l+Nz*(j+Ny*i);
	    input_q[ss]=q[i][j][l][s]*wds[i][j][l];
	  }
	}
      }
      fftw_execute(forward_plan);

      for(i=0;i<Nx;i++){
	for(j=0;j<Ny;j++){
	  for(l=0;l<Nz;l++){
	    ss=l+Nz*(j+Ny*i);
	    transformed_q[ss]*=kds[i][j][l];
	  }
	}
      }
      fftw_execute(inverse_plan);
      
      for(i=0;i<Nx;i++){
	for(j=0;j<Ny;j++){
	  for(l=0;l<Nz;l++){
	    ss=l+Nz*(j+Ny*i);
	    q[i][j][l][s+1]=((final_q[ss]*wds[i][j][l])/(8.0*Nx*Ny*Nz));
	  }
	}
      }
    }
  }

    //std::cout<<Nx*Ny*Nz*dxyz[0]*dxyz[1]*dxyz[2]<<std::endl;

  destroy_3d_double_array(wds);
  destroy_3d_double_array(kds);

};
