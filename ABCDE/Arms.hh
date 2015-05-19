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
    Lam_Period=3.45; //xAB 14
    Hex_Period=3.75; //xAB 14
    BCC_Period=4.60; //xAB 14

    //Lam_Period=3.83; //xAB 20
    //Hex_Period=3.78; //xAB 20
    //BCC_Period=5.00; //xAB 20 
  }else if(Numb_of_Arms==2){ // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    Ns[1]=Ns[0];  // A2                                                                                                                      
    Ns[2]=0;  // A3                                                                                                                      
    Ns[3]=0;  // A4                                                                                                                      
    Ns[4]=N_each_Arm-Ns[0];  // B1                                                                                                                  
    Ns[5]=N_each_Arm-Ns[0];  // B2                                                                                                                  
    Ns[6]=0;  // B3                                                                                                                  
    Ns[7]=0;  // B4
    // Setting the size of the period for the system
    Lam_Period=2.72; //xAB14
    Hex_Period=2.98; //xAB14
    BCC_Period=3.80; //xAB14

    //Lam_Period=2.97; //xAB20
    //Hex_Period=2.90; //xAB20
    //BCC_Period=4.00; //xAB20 
  }else if(Numb_of_Arms==3){ // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    Ns[1]=Ns[0];  // A2                                                                                                                      
    Ns[2]=Ns[0];  // A3                                                                                                                      
    Ns[3]=0;  // A4                                                                                                                      
    Ns[4]=N_each_Arm-Ns[0];  // B1                                                                                                                  
    Ns[5]=N_each_Arm-Ns[0];  // B2                                                                                                                  
    Ns[6]=N_each_Arm-Ns[0];  // B3                                                                                                                  
    Ns[7]=0;  // B4
    // Setting the size of the period for the system
    Lam_Period=2.32; //xAB14
    Hex_Period=2.49; //xAB14
    BCC_Period=3.30; //xAB14

    //Lam_Period=2.49; //xAB20
    //Hex_Period=2.43; //xAB20
    //BCC_Period=3.00; //xAB20 
  }else if(Numb_of_Arms==4){ // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    Ns[1]=Ns[0];  // A2                                                                                                                      
    Ns[2]=Ns[0];  // A3                                                                                                                      
    Ns[3]=Ns[0];  // A4                                                                                                                      
    Ns[4]=N_each_Arm-Ns[0];  // B1                                                                                                                  
    Ns[5]=N_each_Arm-Ns[0];  // B2                                                                                                                  
    Ns[6]=N_each_Arm-Ns[0];  // B3                                                                                                                 
    Ns[7]=N_each_Arm-Ns[0];  // B4
    // Setting the size of the period for the system
    Lam_Period=2.0; //xAB14
    Hex_Period=2.2; //xAB14
    BCC_Period=3.1; //xAB14

    //Lam_Period=2.18; //xAB20
    //Hex_Period=2.13; //xAB20
    //BCC_Period=3.00; //xAB20  
  }else{
    // Error
    std::cout<<"The number of arms you have selected is incorrect! (Arms.hh)"<<std::endl;
  }




}
