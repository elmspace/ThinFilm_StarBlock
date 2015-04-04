
void calculateWavenumbers(double kx[Nx], double lengthX){
	
	int i,j,k;
	//Next step is to calculate the exponential part in wave vectors in the three dimensions of the system.
	for(k=0;k<Nx/2;k++){kx[k]=(2*Pi*k)/(lengthX);}
	for(k=Nx/2;k<Nx;k++){kx[k]=(2*Pi*(Nx-k))/(lengthX);}
	for(k=0;k<Nx;k++){kx[k]*=kx[k];}
	
};

void solveModDiffEqn_FFT(double **q,double *w,double *qInitial, double ds, int Ns, int sign, double kxxx[Nx]){
	
	int staller;
	int i,j,k,i_s; //Counters
	unsigned long ijk; //To convert to row-major form for the FFTW transforms
	double ds_complex;  //Not really needed as it will just end up being ds
	double *wds, *kxxxds; 
	
	wds=create_1d_double_array(Nx,"wds");
	kxxxds=create_1d_double_array(Nx,"kxxxds");
	
	//Here we are predefining the factors to multiply the FFT's with
	for(i=0;i<Nx;i++){  //Note that from http://www.statmech.org/rsh/diffeqn.pdf the dt is our ds! 
		kxxxds[i]=exp(-ds*kxxx[i]);  //This is the (21-22) *Predefine* http://www.statmech.org/rsh/diffeqn.pdf
		wds[i]=exp(-0.5*w[i]*ds);  //This is the (19) from http://www.statmech.org/rsh/diffeqn.pdf
	}
	
	if(sign==1){
		for(i=0;i<Nx;i++){
			q[i][0]=qInitial[i]; //This is initializing the s=0 step in all space q(s+ds)
		}
		
		for(i_s=0;i_s<(Ns+1);i_s++){ //For all of the segment values q(s+ds) ...
			//std::cout << i_s << "   " <<  q[1][1][1][i_s] << std::endl;
			//Multiply the q values by the aux field w
			for(i=0;i<Nx;i++){  
				input_q[i][0]=q[i][i_s]*wds[i];
				input_q[i][1]=0.0;  
			}
			
			fftw_execute(forward_plan);
			
			//Multiply the q transformed values by the kxyzds field
			for(i=0;i<Nx;i++){  
				
				transformed_q[i][0]*=kxxxds[i];
				transformed_q[i][1]*=kxxxds[i];
				
			}
			
			//Do the inverse transform back
			fftw_execute(inverse_plan);
			
			//Multiply the q values by the aux field w again and divide by system volume as well as set this equal to the next q step
			for(i=0;i<Nx;i++){  
				q[i][i_s+1]=((final_q[i][0]*wds[i])/(Nx)); 
				
			}
			
		}
		
		
	}else{
		
		//Starting from the other side of the chain
		for(i=0;i<Nx;i++){
			
			q[i][Ns]=qInitial[i]; //This is initializing the s=0 step in all space q(s+ds)
			
		}
		for(i_s=Ns;i_s>0;i_s--){ //For all of the segment values q(s+ds) ...
			
			//Multiply the q values by the aux field w
			for(i=0;i<Nx;i++){  
				
				input_q[i][0]=q[i][i_s]*wds[i];
				input_q[i][1]=0.0;
				
			}
			
			//Do the forward FFT on the input data
			fftw_execute(forward_plan);
			
			//Multiply the q transformed values by the kxyzds field
			for(i=0;i<Nx;i++){  
				
				transformed_q[i][0]*=kxxxds[i];
				transformed_q[i][1]*=kxxxds[i];
				
			}
			
			//Do the inverse transform back
			fftw_execute(inverse_plan);
			
			//Multiply the q values by the aux field w again and divide by system volume as well as set this equal to the next q step
			for(i=0;i<Nx;i++){  
				
				q[i][i_s-1]=((final_q[i][0]*wds[i])/(Nx));
				
			}
			
			
		}
	}
	
	
	destroy_1d_double_array(wds);
	destroy_1d_double_array(kxxxds);
	
	
};

