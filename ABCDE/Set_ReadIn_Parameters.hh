int Set_ReadIn_Parameters(int numb_of_arg, char* arg_input[]){

  // We give it a pass variable, unless a problem comes up
  int  pass_or_fail=0;
  
  // There are going to be 4 input arguments and 1 file path (not input)
  // So numb_of_arg should be 5
  if(numb_of_arg==5){
    
    int num_of_arms= atoi(arg_input[1]);
    if((num_of_arms==1)||(num_of_arms==2)||(num_of_arms==3)||(num_of_arms==4)){
      Numb_of_Arms=num_of_arms;
    }else{
      std::cout<<"Something wrong with: Number of arms!"<<std::endl;
      pass_or_fail=1;
    }

    if(strcmp( arg_input[2], "Lam") == 0){
      LAM=1;
      HEX=0;
      BCC=0;
    }else if(strcmp( arg_input[2], "Hex") == 0){
      LAM=0;
      HEX=1;
      BCC=0;
    }else if(strcmp( arg_input[2], "BCC") == 0){
      LAM=0;
      HEX=0;
      BCC=1;
    }else{
      std::cout<<"Something wrong with: Phase!"<<std::endl;
      pass_or_fail=1;
    }
    
    if(strcmp( arg_input[3], "Ver") == 0){
      HOR=0;
      VER=1;
      VER_2=0;
    }else if(strcmp( arg_input[3], "Hor") == 0){
      HOR=1;
      VER=0;
      VER_2=0;
    }else if(strcmp( arg_input[3], "Ver_2") == 0){
      HOR=0;
      VER=0;
      VER_2=1;
    }else{
      std::cout<<"Something wrong with: Direction!"<<std::endl;
      pass_or_fail=1;
    }

    double chi_BAir = atof(arg_input[4]);
    if((chi_BAir>0.0)&&(chi_BAir<1.0)){
      xBAir=chi_BAir;
    }else{
      std::cout<<"Something wrong with: xBAir!"<<std::endl;
      pass_or_fail=1;
    }
    
  } else{
    
    // It has failed
    std::cout<<"Something wrong with: the number of arguments!"<<std::endl;
    pass_or_fail=1;
    
  }

  return pass_or_fail;
  
};



