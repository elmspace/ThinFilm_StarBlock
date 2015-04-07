void parametersAB(double *chi,double *f,double &ds,double *Ns,double *dxyz,double **chiMatrix){
  
  int Ds;
  double xAB;
  double zz;

  zz=3.0;

  // Minimize with respect to box size (yes=1.0, No=-1.0)
  box_min=1.0;

  // Setting the generic chi parameters

  xAB=(0.14)*400.0;
  xAD=(0.08)*400.0; // Substrate
  xAE=(0.08)*400.0; // Air
  xBD=(0.06)*400.0; // Substrate
  xBE=(0.08)*400.0; // Air

  xED=0.0; 

  // Degree of polymerization
  Ns[0]=35;  // A1
  Ns[1]=Ns[0];  // A2
  Ns[2]=Ns[0];  // A3
  Ns[3]=Ns[0];  // A4
  Ns[4]=100-Ns[0];  // B1
  Ns[5]=100-Ns[0];  // B2
  Ns[6]=100-Ns[0];  // B3
  Ns[7]=100-Ns[0];  // B4
  // E and D dont have lengths

  // Total length
  Ds=Ns[0]+Ns[1]+Ns[2]+Ns[3]+Ns[4]+Ns[5]+Ns[6]+Ns[7];
  
  //set to 1 to read from restart file
  //set to >0 <1 to set manually
  //set to negative to randomize
  Iomega=1.1;
  
 
  
  dxyz[0]=(2.0*2.2)/Nx;
  dxyz[1]=(2.0*2.2)/Ny;
  dxyz[2]=(zz*2.2)/22;
 
  f[0]=Ns[0]/Ds;  // fA1
  f[1]=Ns[1]/Ds;  // fA2
  f[2]=Ns[2]/Ds;  // fA3
  f[3]=Ns[3]/Ds;  // fA4
  f[4]=Ns[4]/Ds;  // fB1
  f[5]=Ns[5]/Ds;  // fB2
  f[6]=Ns[6]/Ds;  // fB3
  f[7]=Ns[7]/Ds;  // fB4

  // Setting up the individual chi values
  chi[0]=0.0;  
  chi[1]=0.0;  
  
  

  //++++++++++++++++++++++++++++++++++++++++++++++++

  // Average Concentrations
  pA1ave=p1ave*f[0]; //A1
  pA2ave=p1ave*f[1]; //A2
  pA3ave=p1ave*f[2]; //A3
  pA4ave=p1ave*f[3]; //A4

  pB1ave=p1ave*f[4]; //B1
  pB2ave=p1ave*f[5]; //B2
  pB3ave=p1ave*f[6]; //B3
  pB4ave=p1ave*f[7]; //B4

  pEave=p3ave;       //E
  pDave=p2ave;       //D


  ds=1.0/Ds;
  


  // Setting up the chi matrix in this case 8 by 8
  //1
  chiMatrix[0][0]=0.0;  
  chiMatrix[0][1]=chi[0];
  chiMatrix[0][2]=chi[1];
  chiMatrix[0][3]=chi[2];
  chiMatrix[0][4]=chi[3];
  chiMatrix[0][5]=chi[4];
  chiMatrix[0][6]=chi[5];
  chiMatrix[0][7]=chi[6];
  chiMatrix[0][8]=chi[7];
  chiMatrix[0][9]=chi[8];
  //2
  chiMatrix[1][0]=chi[0];  
  chiMatrix[1][1]=0.0;
  chiMatrix[1][2]=chi[9];
  chiMatrix[1][3]=chi[10];
  chiMatrix[1][4]=chi[11];
  chiMatrix[1][5]=chi[12];
  chiMatrix[1][6]=chi[13];
  chiMatrix[1][7]=chi[14];
  chiMatrix[1][8]=chi[15];
  chiMatrix[1][9]=chi[16];
  //3
  chiMatrix[2][0]=chi[1];  
  chiMatrix[2][1]=chi[9];
  chiMatrix[2][2]=0.0;
  chiMatrix[2][3]=chi[17];
  chiMatrix[2][4]=chi[18];
  chiMatrix[2][5]=chi[19];
  chiMatrix[2][6]=chi[20];
  chiMatrix[2][7]=chi[21];
  chiMatrix[2][8]=chi[22];
  chiMatrix[2][9]=chi[23];
  //4
  chiMatrix[3][0]=chi[2];  
  chiMatrix[3][1]=chi[10];
  chiMatrix[3][2]=chi[17];
  chiMatrix[3][3]=0.0;
  chiMatrix[3][4]=chi[24];
  chiMatrix[3][5]=chi[25];
  chiMatrix[3][6]=chi[26];
  chiMatrix[3][7]=chi[27];
  chiMatrix[3][8]=chi[28];
  chiMatrix[3][9]=chi[29];
  //5
  chiMatrix[4][0]=chi[3];  
  chiMatrix[4][1]=chi[11];
  chiMatrix[4][2]=chi[18];
  chiMatrix[4][3]=chi[24];
  chiMatrix[4][4]=0.0;
  chiMatrix[4][5]=chi[30];
  chiMatrix[4][6]=chi[31];
  chiMatrix[4][7]=chi[32];
  chiMatrix[4][8]=chi[33];
  chiMatrix[4][9]=chi[34];
  //6
  chiMatrix[5][0]=chi[4];  
  chiMatrix[5][1]=chi[12];
  chiMatrix[5][2]=chi[19];
  chiMatrix[5][3]=chi[25];
  chiMatrix[5][4]=chi[30];
  chiMatrix[5][5]=0.0;
  chiMatrix[5][6]=chi[35];
  chiMatrix[5][7]=chi[36];
  chiMatrix[5][8]=chi[37];
  chiMatrix[5][9]=chi[38];
  //7
  chiMatrix[6][0]=chi[5];  
  chiMatrix[6][1]=chi[13];
  chiMatrix[6][2]=chi[20];
  chiMatrix[6][3]=chi[26];
  chiMatrix[6][4]=chi[31];
  chiMatrix[6][5]=chi[35];
  chiMatrix[6][6]=0.0;
  chiMatrix[6][7]=chi[39];
  chiMatrix[6][8]=chi[40];
  chiMatrix[6][9]=chi[41];
  //8
  chiMatrix[7][0]=chi[6];  
  chiMatrix[7][1]=chi[14];
  chiMatrix[7][2]=chi[21];
  chiMatrix[7][3]=chi[27];
  chiMatrix[7][4]=chi[32];
  chiMatrix[7][5]=chi[36];
  chiMatrix[7][6]=chi[39];
  chiMatrix[7][7]=0.0;
  chiMatrix[7][8]=chi[42];
  chiMatrix[7][9]=chi[43];
  //9
  chiMatrix[8][0]=chi[7];  
  chiMatrix[8][1]=chi[15];
  chiMatrix[8][2]=chi[22];
  chiMatrix[8][3]=chi[28];
  chiMatrix[8][4]=chi[33];
  chiMatrix[8][5]=chi[37];
  chiMatrix[8][6]=chi[40];
  chiMatrix[8][7]=chi[42];
  chiMatrix[8][8]=0.0;
  chiMatrix[8][9]=chi[44];
  //10
  chiMatrix[9][0]=chi[8];  
  chiMatrix[9][1]=chi[16];
  chiMatrix[9][2]=chi[23];
  chiMatrix[9][3]=chi[29];
  chiMatrix[9][4]=chi[34];
  chiMatrix[9][5]=chi[38];
  chiMatrix[9][6]=chi[41];
  chiMatrix[9][7]=chi[43];
  chiMatrix[9][8]=chi[44];
  chiMatrix[9][9]=0.0;


};
