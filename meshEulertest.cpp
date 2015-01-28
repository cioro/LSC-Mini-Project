#include"Mesh.hpp"
#include"Euler.hpp"
#include<vector>
#include<fstream>
#include<iostream>
#include<assert.h>
#include<string>

//Reflective boundary function
Euler::W_state Refl_Left_Bound(Mesh &m){
 
  Euler::W_state w_left_end = m.ptr_euler-> PfromC(m.data[m.nGhost]);
  
  Euler::W_state w_result;
  
  w_result.rho =   w_left_end.rho;
  w_result.u   = - w_left_end.u;
  w_result.P   =   w_left_end.P;

  return w_result;
}

Euler::W_state Refl_Right_Bound(Mesh &m){
 
  Euler::W_state w_right_end = m.ptr_euler-> PfromC(m.data[m.ncells]);
  
  Euler::W_state w_result;
  
  w_result.rho =   w_right_end.rho;
  w_result.u   = - w_right_end.u;
  w_result.P   =   w_right_end.P;

  return w_result;
}

//Transmissive boundary function
Euler::W_state Trans_Left_Bound(Mesh &m){
 
  Euler::W_state w_left_end = m.ptr_euler-> PfromC(m.data[0]);
  
  Euler::W_state w_result;
  
  w_result.rho =  w_left_end.rho;
  w_result.u   =  w_left_end.u;
  w_result.P   =  w_left_end.P;

  return w_result;
}


Euler::W_state Trans_Right_Bound(Mesh &m){
 
  Euler::W_state w_right_end = m.ptr_euler-> PfromC(m.data[m.ncells]);
  
  Euler::W_state w_result;
  
  w_result.rho =   w_right_end.rho;
  w_result.u   =   w_right_end.u;
  w_result.P   =   w_right_end.P;

  return w_result;
}
//Initial condition function

Euler::U_state initial(double x){

  const double x_0 = 0.5;
  //Should do this assert(x_0 < x_max && x_0 > x_min); to avoid errors
  //but I would have to add to more args to the fcn.

  if (x<=x_0){
    
    Euler::U_state uL(1.0,0.75,1.0);
    return uL;
  }
  if (x>x_0){
    Euler::U_state uR(0.125,0.0,0.1);
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
  //double x_0 = 0.5;//Location of discontinuity at time 0
 
  int ncells = 100;
  int nGhost = 1;

  double dt;
  double dx;

  double T_max = 0.2;
    
  //Initialise mesh with reflective BC
  Mesh m(ncells, x_min, x_max, cfl, initial, Refl_Left_Bound, Refl_Right_Bound, nGhost); 
  
  //Print mesh to file
  std::string file_init = "initial_u";
  m.save_u_state(file_init,Exact_solver);

  //std::cout << m.ptr_euler->gamma << std::endl;
  //Apply BC() function
  // m.applyBC();

  //Print nGhost cells to check reflective BC
  // m.data[0].print();
  // m.data[1].print();
  // m.data[m.ncells].print();
  // m.data[m.ncells+m.nGhost].print();
  //Initialise mesh with transmissive BC

  //Calculate dt
  //dt = m.Calculate_dt();

  
  //Print dt
  //std::cout << "This is the time step dt: " << dt << "\n";
  std::vector<Euler::U_state> flux(m.ncells+1);
 
  
  
  for(double i = m.time; i<T_max; i+=dt){
    
    m.applyBC();
    dt = m.Calculate_dt();
    std::cout <<"dt " << dt << "\n";
    flux = HLLC(m);
    for( std::vector<Euler::U_state>::iterator itflux = flux.begin(); itflux!=flux.end(); itflux++){
    std::cout << (*itflux).rho << "\t" << (*itflux).momentum << "\t" << (*itflux).energy << "\n";
    }
    // Mesh_update(m,flux,dt);
  

  }
  std::string output = "Output";
  m.save_w_state(output,Exact_solver);
}
