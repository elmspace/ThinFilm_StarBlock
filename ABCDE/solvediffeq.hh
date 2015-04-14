//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//                                Solving The Modified Diffusion Equation
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void solveModDiffEqn_FFT(double ****q, double ***w, double ***qint, double ds, int Ns, int sign, double ***k, double *dxyz){
  
  int            i,j,l,s,ss;  // some counters
  unsigned long  ijl; // This is used for the Fourier Transform


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
    for(s=0;s<((int)Ns);s++){ 
      std::cout<<"1"<<std::endl;
      for(i=0;i<Nx;i++){
	for(j=0;j<Ny;j++){
	  for(l=0;l<Nz;l++){
	    ss=l+Nz*(j+Ny*i);
	    std::cout<<"s= "<<ss<<std::endl;
	    //input_q[ss]=q[i][j][l][s]*wds[i][j][l];
	    std::cout<<"inpt "<<input_q[ss]<<std::endl;
	    std::cout<<"q "<<q[i][j][l][s]<<std::endl;
	    std::cout<<"w "<<wds[i][j][l]<<std::endl;
	  }
	}
      }
      std::cout<<"1 after"<<std::endl;
      fftw_execute(forward_plan);
      std::cout<<"2"<<std::endl;
      for(i=0;i<Nx;i++){
	for(j=0;j<Ny;j++){
	  for(l=0;l<Nz;l++){
	    ss=l+Nz*(j+Ny*i);
	    transformed_q[ss]*=kds[i][j][l];
	  }
	}
      }
      std::cout<<"3"<<std::endl;
      fftw_execute(inverse_plan);
      for(i=0;i<Nx;i++){
	for(j=0;j<Ny;j++){
	  for(l=0;l<Nz;l++){
	    ss=l+Nz*(j+Ny*i);
	    q[i][j][l][s+1]=((final_q[ss]*wds[i][j][l])/(8.0*Nx*Ny*Nz));
	  }
	}
      }
      std::cout<<"4"<<std::endl;
      //std::cout<<s<<std::endl;
    }
  }else{

    for(i=0;i<Nx;i++){
      for(j=0;j<Ny;j++){
	for(l=0;l<Nz;l++){
	  q[i][j][l][0]=qint[i][j][l];
	}
      }
    }

    for(s=0;s<(Ns);s++){
      
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


};
