void Mod1(double ****w, double ****phi, double ***eta, double ***PHI_0, double *Ns, double ds, double ***k_vector, double *chi, double *dxyz, double **chiMatrix, double ****h, double *f){

  double del_xBAir=(0.01);
  double max_xBAir=0.2;

  int do_1_calc=1; // if =1 this will just do one calculation, if =0 it will break and wont continue
  Bulk_Calc=0; // Bulk_Calc=1 will do a 1period system in bulk

  
  // Cleaning the .dat file
  std::ofstream outputFile37("./RESULTS/MOD1.dat");
  outputFile37 << std::endl;
  outputFile37.close();
  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  Round=1;
  do{
    
    parametersAB(chi,f,ds,Ns,dxyz,chiMatrix,h);
    set_PHI_0(PHI_0,dxyz);
    omega(w);
    std::cout<<"Round="<<Round<<"   xBAir="<<global_xBAir<<"   delx="<<global_xBAir-global_xAAir<<std::endl;
    FreeEnergy(w,phi,eta,PHI_0,Ns,ds,k_vector,chi,dxyz,chiMatrix,h);
    
    std::ofstream outputFile37("./RESULTS/MOD1.dat" , ios::app);
    outputFile37 <<global_fE<<" "<<flobal_HomfE<<" "<<global_xAB<<" "<<global_xAAir<<" "<<global_xASub<<" "<<global_xBAir<<" "<<global_xBSub<<std::endl;
    outputFile37.close();

    if(do_1_calc==1){break;}
    xBAir+=del_xBAir;

    
    Round++;
  }while(xBAir<max_xBAir);
    
}


