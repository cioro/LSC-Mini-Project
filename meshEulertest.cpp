#include"Mesh.hpp"
#include"Euler.hpp"
#include<vector>
#include<fstream>
#include<iostream>
#include<assert.h>
#include<string>

//Reflective boundary function
Euler::W_state Refl_Left_Bound(Euler::W_state w){
 
     
  Euler::W_state w_result;
  
  w_result.rho =   w.rho;
  w_result.u   = - w.u;
  w_result.P   =   w.P;

  return w_result;
}

Euler::W_state Refl_Right_Bound(Euler::W_state w ){
 
 Euler::W_state w_result;  

  w_result.rho =   w.rho;
  w_result.u   = - w.u;
  w_result.P   =   w.P;

  return w_result;
}

//Transmissive boundary function
Euler::W_state Trans_Left_Bound(Euler::W_state w){
 
 Euler::W_state w_result;
   
  w_result.rho =  w.rho;
  w_result.u   =  w.u;
  w_result.P   =  w.P;

  return w_result;
}


Euler::W_state Trans_Right_Bound(Euler::W_state w){
 
  
  Euler::W_state w_result;
  
  w_result.rho =   w.rho;
  w_result.u   =   w.u;
  w_result.P   =   w.P;

  return w_result;
}
//Initial condition function

Euler::U_state initial_test1(double x){

  const double x_0 = 0.5;
  //Should do this assert(x_0 < x_max && x_0 > x_min); to avoid errors
  //but I would have to add to more args to the fcn.

  Euler e;//To allow use of CfromP and the creation of wl and ul. Not good design.
  
  if (x<=x_0){
    
    const Euler::W_state wL(1.0,0.0,1.0);
    Euler::U_state uL = e.CfromP(wL);
    return uL;
  }
  if (x>x_0){
    const Euler::W_state wR(0.125,0.0,0.1);
    Euler::U_state uR = e.CfromP(wR);
    return uR;
  }
  else{
    std::cout<<"Something went wrong inside inital function in main.cpp"<<std::endl;
    Euler::U_state u;
    return u;
  }
}

Euler::U_state initial_test2(double x){
  
  const double x_0 = 0.5;
  
  Euler e;

 if (x<=x_0){
    
    const Euler::W_state wL(1.0,-2.0,0.4);
    Euler::U_state uL = e.CfromP(wL);
    return uL;
  }
  if (x>x_0){
    const Euler::W_state wR(1.0,2.0,0.4);
    Euler::U_state uR = e.CfromP(wR);
    return uR;
  }
  else{
    std::cout<<"Something went wrong inside inital function in main.cpp"<<std::endl;
    Euler::U_state u;
    return u;
  }
}

Euler::U_state initial_test3(double x){
  
  const double x_0 = 0.5;
  
  Euler e;

 if (x<=x_0){
    
    const Euler::W_state wL(1.0,0.0,1000.0);
    Euler::U_state uL = e.CfromP(wL);
    return uL;
  }
  if (x>x_0){
    const Euler::W_state wR(1.0,0.0,0.01);
    Euler::U_state uR = e.CfromP(wR);
    return uR;
  }
  else{
    std::cout<<"Something went wrong inside inital function in main.cpp"<<std::endl;
    Euler::U_state u;
    return u;
  }

}


Euler::U_state initial_test4(double x){
  
  const double x_0 = 0.4;
  
  Euler e;

 if (x<=x_0){
    
    const Euler::W_state wL(5.99924,19.5975,460.894);
    Euler::U_state uL = e.CfromP(wL);
    return uL;
  }
  if (x>x_0){
    const Euler::W_state wR(5.99242,-6.19633,46.0950);
    Euler::U_state uR = e.CfromP(wR);
    return uR;
  }
  else{
    std::cout<<"Something went wrong inside inital function in main.cpp"<<std::endl;
    Euler::U_state u;
    return u;
  }

}


Euler::U_state initial_test5(double x){
  
  const double x_0 = 0.8;
  
  Euler e;

 if (x<=x_0){
    
    const Euler::W_state wL(1.0,-19.59745,1000.0);
    Euler::U_state uL = e.CfromP(wL);
    return uL;
  }
  if (x>x_0){
    const Euler::W_state wR(1.0,-19.59745,0.01);
    Euler::U_state uR = e.CfromP(wR);
    return uR;
  }
  else{
    std::cout<<"Something went wrong inside inital function in main.cpp"<<std::endl;
    Euler::U_state u;
    return u;
  }

}


//Exact riemann solver function

Euler::U_state Exact_solver(double x){

  Euler::U_state u_empty;

  return u_empty;
}


int main(){
  //Parameters of the problem
 
  
  double x_min = 0, x_max = 1.0; //domain length
  double cfl = 0.9;
 
 
  int ncells = 100;
  int nGhost = 2;

  double dt;
  // double dx;

  double T_max = 0.2;
    
  //Initialise mesh with reflective BC
  Mesh m(ncells, x_min, x_max, cfl, initial_test1, Trans_Left_Bound, Trans_Right_Bound, nGhost); 
  
  //Print mesh to file
  std::string file_init = "initial_u";
  m.save_u_state(file_init,Exact_solver);

    
  //Initialise flux vector
  std::vector<Euler::U_state> flux(m.ncells+1);
 
  //Main loop over domain
  for(double i = 0; i<T_max; i+=dt){
    
    m.applyBC();
    dt = m.Calculate_dt();
    std::cout <<"dt " << dt <<" The time is "<< i << "\n";
    // flux = HLLC(m);
    flux = WAF(m,dt);
    //for( std::vector<Euler::U_state>::iterator itflux = flux.begin(); itflux!=flux.end(); itflux++){
    //std::cout << (*itflux).rho << "\t" << (*itflux).momentum << "\t" << (*itflux).energy << "\n";
    // }

    //for(int j = 29; j < 33; j++){
    //    std::cout << flux[j].rho << "\t" << flux[j].momentum << "\t" << flux[j].energy << "\n";
    //   }

    std::cout << "\n "; 
    Mesh_update(m,flux,dt);
   
    std::string snapshot = "Snap";
    m.save_w_state(snapshot, Exact_solver); 
    
    m.time++;
  }
  std::string output = "Output";
  m.save_w_state(output,Exact_solver);
  

}
