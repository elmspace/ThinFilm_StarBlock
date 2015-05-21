void parametersAB(double *chi,double *f,double &ds,double *Ns,double *dxyz,double **chiMatrix, double ****h){
  
  int i,j,k;
  int Ds;
  int N_each_Arm=100;
  double surface; // turn on=1 or off=0 the surface interactions
  double xAB;
  double chi_HA;
  double chi_HS;

  epsilon=0.5;

  if(Bulk_Calc==1){
    Numb_of_Periods=1.0;
    surface=0.0;
    chi_HA=0.0;
    chi_HS=0.0;
  }else{
    Numb_of_Periods=2.0;
    surface=1.0;
    chi_HA=-100.0*Numb_of_Arms;
    chi_HS=-100.0*Numb_of_Arms;
  }
  
  // 0 read
  // 1 make 
  // 2 random     
  Iomega=0;

  // Minimize with respect to box size (yes=1, No=0)
  box_min=1;
  if(Bulk_Calc==1){
    box_min_xy_relax=0;
    box_min_z_relax=0;
    box_min_xyz_relax=1;
  }else{
    box_min_xy_relax=1;
    box_min_z_relax=0;
    box_min_xyz_relax=0; 
  }

  Ns[8]=Numb_of_Arms*N_each_Arm/10.0;
  Ns[9]=Numb_of_Arms*N_each_Arm/10.0;
  // Degree of polymerization (Each arm of the star is 100)
  if(LAM==1){
    Ns[0]=0.5*N_each_Arm;  // A1
  }else if(HEX==1){
    Ns[0]=0.35*N_each_Arm;  // A1
  } else if(BCC==1){
    Ns[0]=0.35*N_each_Arm;  // A1
  }
  Arms(Numb_of_Arms,Ns,N_each_Arm);
 
  // Total length
  Ds=Ns[0]+Ns[1]+Ns[2]+Ns[3]+Ns[4]+Ns[5]+Ns[6]+Ns[7];

  // Setting the generic chi parameters
  xAB=14.0/N_each_Arm;
  xAB*=Ds;
  chi[0]=xAB;    
  //++++++++++++++++++++++++++++++++++++++++++++++++
  
  h_AAir=surface*(0.0)*Ds;
  h_BAir=surface*(xBAir)*Ds; // This is a variable
  h_ASub=surface*(0.0)*Ds;
  h_BSub=surface*(0.0)*Ds;
  
  // setting the global values:
  global_xAB=xAB;
  global_xAAir=h_AAir;
  global_xASub=h_ASub;
  global_xBAir=h_BAir;
  global_xBSub=h_BSub;
  //__________________________

  if(Bulk_Calc==1){
    if(LAM==1){
      Lx=Lam_Period*Numb_of_Periods;
      Ly=Lam_Period*Numb_of_Periods;
      Lz=Lam_Period*Numb_of_Periods;
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
      Lx=Numb_of_Periods*Lam_Period;
      Ly=Numb_of_Periods*Lam_Period;
      Lz=Film_LZ_Lam;
    }else if(HEX==1){
      if((VER==1)||(VER_2==1)){
	Lx=Hex_Period;
	Ly=Hex_Period*sqrt(3.0);
	Lz=Numb_of_Periods*(Hex_Period+(0.0*(Numb_of_Periods*Hex_Period/Nz)));
      }else{
	Lx=Hex_Period*sqrt(3.0);
	Ly=Hex_Period*sqrt(3.0);
	Lz=Numb_of_Periods*(Hex_Period+(0.0*(Numb_of_Periods*Hex_Period/Nz)));
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
  FilmThickness=Lz-dxyz[2];

  if(Bulk_Calc==1){
    pMultiAve=1.0;
    pAirAve=0.0;
    pSubAve=0.0;
  }else{
    pMultiAve=Set_h_function(h,dxyz);
    pAirAve=0.5*(1.0-pMultiAve);
    pSubAve=pAirAve;
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

	h[8][i][j][k]*=chi_HA;
	h[9][i][j][k]*=chi_HS;
		
      }
    }
  }
  // This is print out of the code parameters, You can un-comment it to check the vraiables.
  /*
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
  std::cout<<" "<<std::endl;
  std::cout<<" "<<std::endl;
  std::cout<<"NA1="<<Ns[0]<<std::endl;
  std::cout<<"NA2="<<Ns[1]<<std::endl;
  std::cout<<"NA3="<<Ns[2]<<std::endl;
  std::cout<<"NA4="<<Ns[3]<<std::endl;
  std::cout<<"NB1="<<Ns[4]<<std::endl;
  std::cout<<"NB2="<<Ns[5]<<std::endl;
  std::cout<<"NB3="<<Ns[6]<<std::endl;
  std::cout<<"NB4="<<Ns[7]<<std::endl;
  std::cout<<"NHA="<<Ns[8]<<std::endl;
  std::cout<<"NHS="<<Ns[9]<<std::endl;
  std::cout<<"Total="<<Ns[0]+Ns[1]+Ns[2]+Ns[3]+Ns[4]+Ns[5]+Ns[6]+Ns[7]<<std::endl;
  std::cout<<"kappa_HA="<<kappa_HA<<std::endl;
  std::cout<<"kappa_HS="<<kappa_HS<<std::endl;
  std::cout<<" "<<std::endl;
  std::cout<<"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<std::endl;
  std::cout<<"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<std::endl;
  */
};
