double homogenousfE(double **chiMatrix){

  int     i,j;
  double  *avphi;
  double  fE_homo;

  avphi=create_1d_double_array(ChainType,"avphi");

  avphi[0]=pA1ave; // A1 average
  avphi[1]=pA2ave; // A2 average
  avphi[2]=pA3ave; // A3 average
  avphi[3]=pA4ave; // A4 average

  avphi[4]=pB1ave; // B1 average
  avphi[5]=pB2ave; // B2 average
  avphi[6]=pB3ave; // B3 average
  avphi[7]=pB4ave; // B4 average

  fE_homo=0.0;

  for(i=0;i<ChainType;i++){
    for(j=i;j<ChainType;j++){

      fE_homo+=avphi[i]*avphi[j]*chiMatrix[i][j];

    }}
  
  return fE_homo;

  destroy_1d_double_array(avphi);

};
