#include"Euler.hpp"
#include<vector>
#include<cmath>
#include<iostream>

Euler::W_state::W_state() : rho(0), u(0), P(0) {}

Euler::W_state::W_state(double rho_arg, double u_arg, double P_arg) : rho(rho_arg), u(u_arg), P(P_arg) {}

void Euler::W_state::print()const {
  std::cout << " rho\t " << rho << "\t u \t " << u << " \t P \t" << P << "\n";

}


Euler::U_state::U_state() : rho(0), momentum(0), energy(0) {}

Euler::U_state::U_state(double rho_arg, double momentum_arg, double energy_arg) : \
 rho(rho_arg), momentum(momentum_arg), energy(energy_arg) {}

void Euler::U_state::print()const {
  std::cout << " rho \t " << rho << "\t momentum \t " << momentum \
 << " \t energy \t" << energy << "\n";

}


double Euler::a(){
return 0;
}

double Euler::a(const Euler::W_state& w){
  if(w.P < 0.0 || w.rho < 0.0){
    std::cout << "Error speed  of sound can't be computed. Negative pressure or density" << "\n";
    w.print();
      }
  double a_result = sqrt(gamma*(w.P)/(w.rho));
  return a_result;
}
  
double Euler::int_energy(const Euler::W_state& w){

  double e = w.P/((gamma-1)*w.rho);
  return e;
  
}

//The Formula come from Toro(ed.2009) p.89
Euler::U_state Euler::flux(Euler::U_state& u){
  double f1 = u.momentum;
  double f2 = 0.5*(3-gamma)*(u.momentum*u.momentum/u.rho)+(gamma-1)*u.energy;
  double f3 = gamma*(u.momentum/u.rho)*u.energy - 0.5*(gamma-1)*(u.momentum*u.momentum*u.momentum)/(u.rho*u.rho);

  return U_state(f1,f2,f3);
  

}

Euler::U_state Euler::CfromP(const Euler::W_state& w){
  double u1 = w.rho;
  double u2 = w.rho*w.u;
  double u3 = w.rho*(0.5*w.u*w.u +int_energy(w));
  return U_state(u1,u2,u3);

}


Euler::W_state Euler::PfromC(const Euler::U_state& u){
  double w1 = u.rho;
  double w2 = u.momentum/u.rho;
  double w3 = (gamma - 1)*(u.energy -0.5*(u.momentum*u.momentum/u.rho));
  return W_state(w1,w2,w3);

}

Euler::Euler(){
  U_state();
  W_state();
}

