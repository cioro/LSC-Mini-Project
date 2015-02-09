//Mesh class:
//Initialises a 1d vector of primitive variables,
//Calculates the boundary conditions
//Calculates the adequate dt for each time step
//HLLC flux calculator
//WAF flux calculator

///// ------Note: need to include a way of changing BC outside of constructor -----

#ifndef MESH_H
#define MESH_H

#include <cstdio>
#include <vector>
#include <string>
#include "Euler.hpp"

class Mesh{
public:
  //Pointer to Euler class
  Euler * ptr_euler;
  
  //Parameters
  std::vector<double> axis;
  int ncells;

  double x_min;
  double x_max;
  double dx;
  double time;
  double cfl;
  
  //Boundary functions
  Euler::W_state (*boundary1)(Euler::W_state  w);
  Euler::W_state (*boundary2)(Euler::W_state  w);
  int nGhost;

   //Data
  std::vector<Euler::U_state> data;
  std::vector<Euler::U_state> data2;

  //Constructor
  Mesh();
  Mesh(int ncells, double x_min, double x_max,double cfl, Euler::U_state (*f)(double x),Euler::W_state (*b1)(Euler::W_state w), Euler::W_state (*b2)(Euler::W_state w), int nGhost);
  ~Mesh();
  //print to screen and print to file
  void print()const;
  void save_u_state(std::string name, Euler::U_state (*exact)(double x))const;
  void save_w_state(std::string name, Euler::U_state (*exact)(double x))const;
 
 //Apply BCs
  void applyBC();
  
  //Calculate dt 
  double Calculate_dt();

  void reset(Euler::U_state (*f)(double x));
  // //Mesh update
  //  //  void Mesh_update(Mesh &m, std::vector<Euler::U_state> &flux,double dt);
};
//HLLC method
std::vector<Euler::U_state> HLLC(Mesh &m);


//TVD WAF
std::vector<Euler::U_state> WAF(Mesh &m,double dt,std::string limiter);
std::vector<Euler::U_state> HLLC_U_state(Euler::U_state U_state_L, Euler::U_state U_state_R);
//Limiter fcns
double minmod(double r, double c);
double superbee(double r, double c);

int sign(double c);

 //Mesh update
void Mesh_update(Mesh &m, std::vector<Euler::U_state> &flux,double dt);


#endif
