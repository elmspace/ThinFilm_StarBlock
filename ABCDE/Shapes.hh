double Cube(int x, int y,int z,int kind){

   double conc_at_xyz;

  if(kind==2){

    conc_at_xyz=(normlD*(1.0+tanh((radDx-abs(x-sigmaDx))/(widthD))))*(normlD*(1.0+tanh((radDy-abs(y-sigmaDy))/(widthD))))*(normlD*(1.0+tanh((radDz-abs(z-sigmaDz))/(widthD))));

  }else if(kind==3){

    conc_at_xyz=(normlE*(1.0+tanh((radEx-abs(x-sigmaEx))/(widthE))))*(normlE*(1.0+tanh((radEy-abs(y-sigmaEy))/(widthE))))*(normlE*(1.0+tanh((radEz-abs(z-sigmaEz))/(widthE))));

  }else{

    std::cout<<"something bad happened in function Cube!";

  }
  
  return conc_at_xyz;


};

