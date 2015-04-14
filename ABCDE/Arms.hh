void Arms(int Numb_of_Arms, double *Ns){


  if(Numb_of_Arms==1){
    Ns[1]=1;  // A2                                                                                                                      
    Ns[2]=1;  // A3                                                                                                             
    Ns[3]=1;  // A4                                                                                                   
    Ns[4]=100-Ns[0];  // B1                                                                            
    Ns[5]=1;  // B2               
    Ns[6]=1;  // B3              
    Ns[7]=1;  // B4 
  }else if(Numb_of_Arms==2){
    Ns[1]=Ns[0];  // A2                                                                                                                      
    Ns[2]=1;  // A3                                                                                                                      
    Ns[3]=1;  // A4                                                                                                                      
    Ns[4]=100-Ns[0];  // B1                                                                                                                  
    Ns[5]=100-Ns[0];  // B2                                                                                                                  
    Ns[6]=1;  // B3                                                                                                                  
    Ns[7]=1;  // B4
  }else if(Numb_of_Arms==3){
    Ns[1]=Ns[0];  // A2                                                                                                                      
    Ns[2]=Ns[0];  // A3                                                                                                                      
    Ns[3]=1;  // A4                                                                                                                      
    Ns[4]=100-Ns[0];  // B1                                                                                                                  
    Ns[5]=100-Ns[0];  // B2                                                                                                                  
    Ns[6]=100-Ns[0];  // B3                                                                                                                  
    Ns[7]=1;  // B4
  }else if(Numb_of_Arms==4){
    Ns[1]=Ns[0];  // A2                                                                                                                      
    Ns[2]=Ns[0];  // A3                                                                                                                      
    Ns[3]=Ns[0];  // A4                                                                                                                      
    Ns[4]=100-Ns[0];  // B1                                                                                                                  
    Ns[5]=100-Ns[0];  // B2                                                                                                                  
    Ns[6]=100-Ns[0];  // B3                                                                                                                 
    Ns[7]=100-Ns[0];  // B4
  }else{
    // Error
    std::cout<<"The number of arms you have selected is incorrect! (Arms.hh)"<<std::endl;
  }




}