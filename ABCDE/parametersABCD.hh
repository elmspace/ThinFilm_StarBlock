void parametersAB(double *chi,double *f,double &ds,double *Ns,double *dxyz,double **chiMatrix, double ****h){
  
  int i,j,k;
  int Ds;
  double surface=0.0; // turn on=1 or off=0 the surface interactions
  double xAB;
  double chi_HA=-100.0;
  double chi_HS=-100.0;
  
  pMultiAve=0.85;
  pAirAve=0.075;
  pSubAve=0.075;
  
  if(Bulk_Calc==1){
    Numb_of_Periods=1.0;
    surface=0.0;
  }else{
    Numb_of_Periods=3.0;
  }
  
  // 0 read
  // 1 make 
  // 2 random     
  Iomega=0;

  // Minimize with respect to box size (yes=1, No=0)
  box_min=1;
  if(Bulk_Calc==1){
    box_min_xy_relax=0;
    box_min_xyz_relax=1;
  }else{
    box_min_xy_relax=0; //flip these two when done
    box_min_xyz_relax=1;
  }

  Ns[8]=50;
  Ns[9]=50;
  // Degree of polymerization (Each arm of the star is 100)
  if(LAM==1){
    Ns[0]=50;  // A1
  }else if(HEX==1){
    Ns[0]=35;  // A1
  } else if(BCC==1){
    Ns[0]=35;  // A1
  }else{
    std::cout<<"You have not chosen a phase yet."<<std::endl;
    Ns[0]=50;  // A1 
  }
  Arms(Numb_of_Arms,Ns);
 
  // Total length
  Ds=Ns[0]+Ns[1]+Ns[2]+Ns[3]+Ns[4]+Ns[5]+Ns[6]+Ns[7];

  // Setting the generic chi parameters
  xAB=(0.14)*Ds;
  // Setting up the xAB
  chi[0]=xAB;    
  //++++++++++++++++++++++++++++++++++++++++++++++++
  
  h_AAir=surface*(0.08)*Ds;
  h_BAir=surface*(xBAir)*Ds; // This is a variable
  h_ASub=surface*(0.08)*Ds;
  h_BSub=surface*(0.06)*Ds;
  
  // setting the global values:
  global_xAB=xAB;
  global_xAAir=h_AAir;
  global_xASub=h_ASub;
  global_xBAir=h_BAir;
  global_xBSub=h_BSub;
  //__________________________

  if(Bulk_Calc==1){
    if(LAM==1){
      Lx=Lam_Period;
      Ly=Lam_Period;
      Lz=Lam_Period;
    }else if(HEX==1){
      Lx=Hex_Period;
      Ly=Hex_Period*sqrt(3.0);
      Lz=Hex_Period;
    }else if(BCC==1){
      Lx=BCC_Period;
      Ly=BCC_Period;
      Lz=BCC_Period;
    }else{
      std::cout<<"You have not chosen a phase yet."<<std::endl;
    }
  }else{
    if(LAM==1){
      Lx=2.0*Lam_Period;
      Ly=2.0*Lam_Period;
      Lz=Numb_of_Periods*Lam_Period;
    }else if(HEX==1){
      if(VER==1){
	Lx=Hex_Period;
	Ly=Hex_Period*sqrt(3.0);
	Lz=Numb_of_Periods*Hex_Period;
      }else{
	Lx=Hex_Period*sqrt(3.0);
	Ly=Hex_Period;
	Lz=Numb_of_Periods*Hex_Period;
      }
    }else{
      std::cout<<"You have not chosen a phase yet."<<std::endl;
    }
  }
 
  
  if(Round==1){
    dxyz[0]=Lx/Nx;
    dxyz[1]=Ly/Ny;
    dxyz[2]=Lz/Nz;
  }else{
    // Will be using the previous dxyz values. The code is looping.
  }

 
  f[0]=Ns[0]/Ds;  // fA1
  f[1]=Ns[1]/Ds;  // fA2
  f[2]=Ns[2]/Ds;  // fA3
  f[3]=Ns[3]/Ds;  // fA4
  f[4]=Ns[4]/Ds;  // fB1
  f[5]=Ns[5]/Ds;  // fB2
  f[6]=Ns[6]/Ds;  // fB3
  f[7]=Ns[7]/Ds;  // fB4

  kappa_HA=Ns[8]/Ds;
  kappa_HS=Ns[9]/Ds;


  // Average Concentrations
  pA1ave=pMultiAve*f[0]; //A1
  pA2ave=pMultiAve*f[1]; //A2
  pA3ave=pMultiAve*f[2]; //A3
  pA4ave=pMultiAve*f[3]; //A4

  pB1ave=pMultiAve*f[4]; //B1
  pB2ave=pMultiAve*f[5]; //B2
  pB3ave=pMultiAve*f[6]; //B3
  pB4ave=pMultiAve*f[7]; //B4


  ds=1.0/Ds;
  

  // Setting up the chi matrix in this case 10 by 10
  
  chiMatrix[0][0]=0.0;  
  chiMatrix[0][1]=0.0;
  chiMatrix[0][2]=0.0;
  chiMatrix[0][3]=0.0;
  chiMatrix[0][4]=chi[0];
  chiMatrix[0][5]=chi[0];
  chiMatrix[0][6]=chi[0];
  chiMatrix[0][7]=chi[0];
  chiMatrix[0][8]=h_AAir;
  chiMatrix[0][9]=h_ASub;
  
  chiMatrix[1][0]=0.0;  
  chiMatrix[1][1]=0.0;
  chiMatrix[1][2]=0.0;
  chiMatrix[1][3]=0.0;
  chiMatrix[1][4]=chi[0];
  chiMatrix[1][5]=chi[0];
  chiMatrix[1][6]=chi[0];
  chiMatrix[1][7]=chi[0];
  chiMatrix[1][8]=h_AAir;
  chiMatrix[1][9]=h_ASub;

  chiMatrix[2][0]=0.0;  
  chiMatrix[2][1]=0.0;
  chiMatrix[2][2]=0.0;
  chiMatrix[2][3]=0.0;
  chiMatrix[2][4]=chi[0];
  chiMatrix[2][5]=chi[0];
  chiMatrix[2][6]=chi[0];
  chiMatrix[2][7]=chi[0];
  chiMatrix[2][8]=h_AAir;
  chiMatrix[2][9]=h_ASub;

  chiMatrix[3][0]=0.0;  
  chiMatrix[3][1]=0.0;
  chiMatrix[3][2]=0.0;
  chiMatrix[3][3]=0.0;
  chiMatrix[3][4]=chi[0];
  chiMatrix[3][5]=chi[0];
  chiMatrix[3][6]=chi[0];
  chiMatrix[3][7]=chi[0];
  chiMatrix[3][8]=h_AAir;
  chiMatrix[3][9]=h_ASub;

  chiMatrix[4][0]=chi[0];  
  chiMatrix[4][1]=chi[0];
  chiMatrix[4][2]=chi[0];
  chiMatrix[4][3]=chi[0];
  chiMatrix[4][4]=0.0;
  chiMatrix[4][5]=0.0;
  chiMatrix[4][6]=0.0;
  chiMatrix[4][7]=0.0;
  chiMatrix[4][8]=h_BAir;
  chiMatrix[4][9]=h_BSub;

  chiMatrix[5][0]=chi[0];  
  chiMatrix[5][1]=chi[0];
  chiMatrix[5][2]=chi[0];
  chiMatrix[5][3]=chi[0];
  chiMatrix[5][4]=0.0;
  chiMatrix[5][5]=0.0;
  chiMatrix[5][6]=0.0;
  chiMatrix[5][7]=0.0;
  chiMatrix[5][8]=h_BAir;
  chiMatrix[5][9]=h_BSub;

  chiMatrix[6][0]=chi[0];  
  chiMatrix[6][1]=chi[0];
  chiMatrix[6][2]=chi[0];
  chiMatrix[6][3]=chi[0];
  chiMatrix[6][4]=0.0;
  chiMatrix[6][5]=0.0;
  chiMatrix[6][6]=0.0;
  chiMatrix[6][7]=0.0;
  chiMatrix[6][8]=h_BAir;
  chiMatrix[6][9]=h_BSub;

  chiMatrix[7][0]=chi[0];  
  chiMatrix[7][1]=chi[0];
  chiMatrix[7][2]=chi[0];
  chiMatrix[7][3]=chi[0];
  chiMatrix[7][4]=0.0;
  chiMatrix[7][5]=0.0;
  chiMatrix[7][6]=0.0;
  chiMatrix[7][7]=0.0;
  chiMatrix[7][8]=h_BAir;
  chiMatrix[7][9]=h_BSub;


  chiMatrix[8][0]=h_AAir;  
  chiMatrix[8][1]=h_AAir;
  chiMatrix[8][2]=h_AAir;
  chiMatrix[8][3]=h_AAir;
  chiMatrix[8][4]=h_BAir;
  chiMatrix[8][5]=h_BAir;
  chiMatrix[8][6]=h_BAir;
  chiMatrix[8][7]=h_BAir;
  chiMatrix[8][8]=0.0;
  chiMatrix[8][9]=0.0;

  chiMatrix[9][0]=h_ASub;  
  chiMatrix[9][1]=h_ASub;
  chiMatrix[9][2]=h_ASub;
  chiMatrix[9][3]=h_ASub;
  chiMatrix[9][4]=h_BSub;
  chiMatrix[9][5]=h_BSub;
  chiMatrix[9][6]=h_BSub;
  chiMatrix[9][7]=h_BSub;
  chiMatrix[9][8]=0.0;
  chiMatrix[9][9]=0.0;


 
  // Setting up the surface field
  for(i=0;i<Nx;i++){
    for(j=0;j<Ny;j++){
      for(k=0;k<Nz;k++){
	
	if((k==0)||(k==1)){ // k=0 is the substrate surface
	  h[0][i][j][k]=0.0;
	  h[1][i][j][k]=0.0;
	  h[2][i][j][k]=0.0;
	  h[3][i][j][k]=0.0;
	  h[4][i][j][k]=0.0;
	  h[5][i][j][k]=0.0;
	  h[6][i][j][k]=0.0;
	  h[7][i][j][k]=0.0;
	  h[8][i][j][k]=0.0;
	  h[9][i][j][k]=chi_HS;
	}else if((k==(Nz-1))||(k==(Nz-2))){ // k=Nz-1 is the air interface
	  h[0][i][j][k]=0.0;
	  h[1][i][j][k]=0.0;
	  h[2][i][j][k]=0.0;
	  h[3][i][j][k]=0.0;
	  h[4][i][j][k]=0.0;
	  h[5][i][j][k]=0.0;
	  h[6][i][j][k]=0.0;
	  h[7][i][j][k]=0.0;
	  h[8][i][j][k]=chi_HA;
	  h[9][i][j][k]=0.0;
	}else{ // No surface interaction in bulk
	  h[0][i][j][k]=0.0;
	  h[1][i][j][k]=0.0;
	  h[2][i][j][k]=0.0;
	  h[3][i][j][k]=0.0;
	  h[4][i][j][k]=0.0;
	  h[5][i][j][k]=0.0;
	  h[6][i][j][k]=0.0;
	  h[7][i][j][k]=0.0;
	  h[8][i][j][k]=0.0;
	  h[9][i][j][k]=0.0;
	}

      }
    }
  }
  // This is print out of the code parameters, You can un-comment it to check the vraiables.

  std::cout<<"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<std::endl;
  std::cout<<"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<std::endl;
  std::cout<<"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<std::endl;
  std::cout<<"++++++++++++++++++ Thin Film Project +++++++++++++++++++++++++"<<std::endl;
  std::cout<<"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<std::endl;
  std::cout<<"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<std::endl;
  std::cout<<" "<<std::endl;
  std::cout<<" "<<std::endl;
  std::cout<<"Number of arms= "<<Numb_of_Arms<<std::endl;
  std::cout<<" "<<std::endl;
  std::cout<<" "<<std::endl;
  if(LAM==1){std::cout<<"Phase= Lam"<<std::endl;}
  else if(HEX==1){std::cout<<"Phase= Hex"<<std::endl;}
  else{std::cout<<"Warning: No Phase was chosen."<<std::endl;}
  std::cout<<" "<<std::endl;
  std::cout<<" "<<std::endl;
  if(HOR==1){std::cout<<"Direction: Horizontal"<<std::endl;}
  if(VER==1){std::cout<<"Direction: Vertical"<<std::endl;}
  std::cout<<" "<<std::endl;
  std::cout<<" "<<std::endl;
  std::cout<<"xAB="<<xAB<<std::endl;
  std::cout<<" "<<std::endl;
  std::cout<<" "<<std::endl;
  std::cout<<"xAAir="<<h_AAir<<std::endl;
  std::cout<<"xBAir="<<h_BAir<<std::endl;
  std::cout<<"xAsub="<<h_ASub<<std::endl;
  std::cout<<"xBSub="<<h_BSub<<std::endl;
  std::cout<<" "<<std::endl;
  std::cout<<" "<<std::endl;
  std::cout<<"Lx="<<Lx<<std::endl;
  std::cout<<"Ly="<<Ly<<std::endl;  
  std::cout<<"Lz="<<Lz<<std::endl;
  std::cout<<" "<<std::endl;
  std::cout<<" "<<std::endl;
  std::cout<<"fA1="<<f[0]<<std::endl;
  std::cout<<"fA2="<<f[1]<<std::endl;
  std::cout<<"fA3="<<f[2]<<std::endl;
  std::cout<<"fA4="<<f[3]<<std::endl;
  std::cout<<"fB1="<<f[4]<<std::endl;
  std::cout<<"fB2="<<f[5]<<std::endl;
  std::cout<<"fB3="<<f[6]<<std::endl;
  std::cout<<"fB4="<<f[7]<<std::endl;
  std::cout<<"Total="<<f[0]+f[1]+f[2]+f[3]+f[4]+f[5]+f[6]+f[7]<<std::endl;
  std::cout<<" "<<std::endl;
  std::cout<<"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<std::endl;
  std::cout<<"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<std::endl;

};
