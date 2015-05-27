void Arms(int Numb_of_Arms, double *Ns, int N_each_Arm){


  if(Numb_of_Arms==1){ // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    Ns[1]=0;  // A2                                                                                                                      
    Ns[2]=0;  // A3                                                                                                             
    Ns[3]=0;  // A4                                                                                                   
    Ns[4]=N_each_Arm-Ns[0];  // B1                                                                            
    Ns[5]=0;  // B2               
    Ns[6]=0;  // B3              
    Ns[7]=0;  // B4 
    // Setting the size of the period for the system
    Lam_Period=3.6;
    Film_LZ_Lam=9.1;
    Hex_Period=3.9;
    Film_LZ_Hex=9.1; 
  }else if(Numb_of_Arms==2){ // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    Ns[1]=Ns[0];  // A2                                                                                                                      
    Ns[2]=0;  // A3                                                                                                                      
    Ns[3]=0;  // A4                                                                                                                      
    Ns[4]=N_each_Arm-Ns[0];  // B1                                                                                                                  
    Ns[5]=N_each_Arm-Ns[0];  // B2                                                                                                                  
    Ns[6]=0;  // B3                                                                                                                  
    Ns[7]=0;  // B4
    // Setting the size of the period for the system
    Lam_Period=2.7;
    Film_LZ_Lam=7.3;
    Hex_Period=3.0;
    Film_LZ_Hex=6.8;
  }else if(Numb_of_Arms==3){ // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    Ns[1]=Ns[0];  // A2                                                                                                                      
    Ns[2]=Ns[0];  // A3                                                                                                                      
    Ns[3]=0;  // A4                                                                                                                      
    Ns[4]=N_each_Arm-Ns[0];  // B1                                                                                                                  
    Ns[5]=N_each_Arm-Ns[0];  // B2                                                                                                                  
    Ns[6]=N_each_Arm-Ns[0];  // B3                                                                                                                  
    Ns[7]=0;  // B4
    // Setting the size of the period for the system
    Lam_Period=2.3;
    Film_LZ_Lam=6.0;
    Hex_Period=2.5;
    Film_LZ_Hex=5.95; 
  }else if(Numb_of_Arms==4){ // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    Ns[1]=Ns[0];  // A2                                                                                                                      
    Ns[2]=Ns[0];  // A3                                                                                                                      
    Ns[3]=Ns[0];  // A4                                                                                                                      
    Ns[4]=N_each_Arm-Ns[0];  // B1                                                                                                                  
    Ns[5]=N_each_Arm-Ns[0];  // B2                                                                                                                  
    Ns[6]=N_each_Arm-Ns[0];  // B3                                                                                                                 
    Ns[7]=N_each_Arm-Ns[0];  // B4
    // Setting the size of the period for the system
    Lam_Period=2.0;
    Film_LZ_Lam=5.4;
    Hex_Period=2.2;
    Film_LZ_Hex=5.3;
  }else{
    // Error
    std::cout<<"The number of arms you have selected is incorrect! (Arms.hh)"<<std::endl;
  }




}
