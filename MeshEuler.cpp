#include<vector>
#include"Mesh.hpp"
#include<fstream>
#include<cmath>
#include<iostream> 
#include<string>
#include<sstream>
#include<algorithm>
#include"Euler.hpp"

//Initialises the parameters of the grid and fills the vector of primitives with the initial conditions
Mesh::Mesh(int ncells, double x_min, double x_max,double cfl, Euler::U_state (*f)(double x), Euler::W_state (*b1)(Mesh &m), Euler::W_state (*b2)(Mesh &m), int nGhost) : ncells(ncells), x_min(x_min),x_max(x_max),cfl(cfl), time(0), boundary1(b1), boundary2(b2), nGhost(nGhost)
 { 
   ptr_euler = new Euler();
   dx = (x_max-x_min)/(double)ncells;
   
   data.resize(ncells + 2*nGhost);
   data2.resize(ncells + 2*nGhost);
   axis.resize(ncells + 2*nGhost);
   
   std::vector<Euler::U_state>::iterator itdata = data.begin()+nGhost;
   std::vector<Euler::U_state>::iterator itdata2= data2.begin()+nGhost;
   std::vector<double>::iterator itaxis= axis.begin()+nGhost;
   
   for(int i = 0; i<ncells; i++){
     (*itaxis) = x_min + i*dx;
     (*itdata) = f(*itaxis);
     itaxis++;
     itdata++;
    
     }
  
 }
//Destructor
Mesh::~Mesh()
{
  delete ptr_euler;
  /* delete data;
  delete axis;
  delete data2;*/
}

//Prints vector of conserved variables to screen
void Mesh::print()const
{
  
  std::vector<Euler::U_state>::const_iterator itdata = data.begin()+nGhost;
  std::vector<Euler::U_state>::const_iterator itdata2= data2.begin()+nGhost;
   std::vector<double>::const_iterator itaxis= axis.begin()+nGhost;
  
   for (int i = nGhost; i<(ncells+nGhost); i++){
     std::cout << *itaxis << "\t";
     (*itdata).print();
     itaxis++;
     itdata++;
   }
  
}
//Print to a file the 1D vector of conserved variables and the exact solution
void Mesh::save_u_state(std::string filename, Euler::U_state (*exact)(double x))const
{
  
  std::stringstream ss;
  ss << filename <<"_" << time;
  std::string tmppath = ss.str();
  

  FILE * outfile = fopen(tmppath.c_str(),"w");
  std::vector<Euler::U_state>::const_iterator itdata= data.begin()+nGhost;
  std::vector<double>::const_iterator itaxis= axis.begin()+nGhost;
  for(int i=1; i<ncells+nGhost; i++)
    {
      fprintf(outfile, "%.4f \t %.4f \t %.4f \t %.4f  \n", *itaxis, (*itdata).rho,(*itdata).momentum,(*itdata).energy );
      itaxis++;
      itdata++;
  }
  fclose(outfile);
  
}

//Prints to a file the 1D vector of primitive variables and the exact solution 
void Mesh::save_w_state(std::string filename, Euler::U_state (*exact)(double x))const
{
  
  std::stringstream ss;
  ss << filename << time;
  std::string tmppath = ss.str();
  

  FILE * outfile = fopen(tmppath.c_str(),"w");
  std::vector<Euler::U_state>::const_iterator itdata= data.begin()+nGhost;
  std::vector<double>::const_iterator itaxis= axis.begin()+nGhost;
 
  for(int i=1; i<ncells+nGhost; i++)
    {
      
      Euler::W_state w_print_approx = ptr_euler->PfromC(*itdata);
      //Euler::W_state w_print_exact = ptr_euler->PfromC(*itaxis);
      
      fprintf(outfile, "%.4f \t %.4f \t %.4f \t %.4f \n", *itaxis, w_print_approx.rho,w_print_approx.u,w_print_approx.P );
      itaxis++;
      itdata++;
  }
  fclose(outfile);
  
}

//Implements the boundary conditions. The actual boundary condition function should be in the main file
void Mesh::applyBC(){
  
  Euler::W_state w_Left_End = ptr_euler -> PfromC(data[nGhost]);
  Euler::W_state w_BC_Left = boundary1(*this);
  data[0]= ptr_euler -> CfromP(w_BC_Left);

  Euler::W_state w_Right_End = ptr_euler -> PfromC(data[ncells]);
  Euler::W_state w_BC_Right = boundary2(*this);
  data[ncells+nGhost]= ptr_euler -> CfromP(w_BC_Right);
  
}

//Calculates the adquate size of the time step dt. See page 183 from Toro(ed.2009)
double Mesh::Calculate_dt(){
  std::vector<Euler::U_state>::iterator itdata = data.begin()+nGhost;
  double speed=0.0;
  double speedtemp=0.0;
    
  for(int i=nGhost; i<nGhost+ncells;i++){
    Euler::W_state w = ptr_euler->PfromC(*itdata);
    speedtemp = ptr_euler->a(w);
    if(fabs(speedtemp) > fabs(speed)){
      speed = speedtemp;
    }          
    itdata++;
  }

  double dt;
  dt=(cfl*dx)/speed;
  return dt;
}

//HLLC flux calculator (this is a free function not a member function of class Mesh)
//Check Toro(ed.2009) p.331 for summary of HLLC method
std::vector<Euler::U_state> HLLC(Mesh &m){
 
  //Total vector of fluxes
  std::vector<Euler::U_state> flux(m.ncells+1);
  
  double gamma = m.ptr_euler->gamma;
  //Variables used

  double P_star,P_L,P_R; //Presures
  double rho_star,rho_L,rho_R;//densities
  double u_L,u_R;//particle/gas speed in cell

  double a_L, a_R;// sound speed in cell
  double S_L, S_R,S_star;// wave speed in cell left wave, right wave, contact wave

  Euler::U_state U_star_L, U_star_R; // Star states of conserved var
  
  double P_pvrs;
  double rho_bar, a_bar;//average density and average sound speed of left and right cells

  Euler::W_state w_temp_left;
  Euler::W_state w_temp_right;

  double q_R,q_L;

  Euler::U_state U_state_L;
  Euler::U_state U_state_R;
  Euler::W_state W_L;
  Euler::W_state W_R;
  Euler::U_state U_state_L_star;
  Euler::U_state U_state_R_star;

  double star_coef_left; // The coeficient in eq. 10.73 from Toro(ed.2009); (S_k-u_k)/(S_k-u_star_k);
  double star_coef_right; // The coeficient in eq. 10.73 from Toro(ed.2009);
  
  std::vector<Euler::U_state>::iterator itflux = flux.begin();
  
  // m.data[m.nGhost].print();
  // m.data[0].print();

  //Loop over whole domain
  for(int i = m.nGhost-1; i < m.ncells+m.nGhost; i++){
    
    //Select U_state and initialise W_state

    U_state_L = m.data[i];
    U_state_R = m.data[i+1];

    W_L = m.ptr_euler->PfromC(U_state_L);
    W_R = m.ptr_euler->PfromC(U_state_R);


    //---------Pressure estimate-------------------------

    P_L = W_L.P;
    P_R = W_R.P;

    u_L = W_L.u;
    u_R = W_R.u;

    rho_L = W_L.rho;
    rho_R = W_R.rho;
   
    a_L =m.ptr_euler->a(W_L);
    a_R = m.ptr_euler->a(W_R);

    std::cout <<"Inside the HLLC function" << "\n";
        
    rho_bar = 0.5*(rho_L + rho_R);
    a_bar = 0.5*(a_L + a_R);

    P_pvrs = 0.5*(P_L + P_R)-0.5*(u_R-u_L)*rho_bar*a_bar;

    P_star = std::max(0.0,P_pvrs);
      
    //---------------------------------------------------
   
    //----------Wave speed estimate----------------------
    
    //Calculate q_R

    if(P_star <= P_R){
      q_R = 1.0;
    }
    else{
      q_R = sqrt(1 + ((gamma+1)/(2*gamma))*((P_star/P_R)-1));
    }
        
    //Calculate q_L
    if(P_star <= P_L){
      q_L = 1.0;
    }
    else{
      q_L = sqrt(1 + ((gamma+1)/(2*gamma))*((P_star/P_L)-1));
    }
 
    //Calculate S_R and S_L
    S_L = u_L - a_L*q_L;
    S_R = u_R + a_R*q_R;
    
    double numerator = P_R - P_L + rho_L*u_L*(S_L-u_L) - rho_R*u_R*(S_R-u_R); 
    double denominator = rho_L*(S_L-u_L)-rho_R*(S_R-u_R);

    S_star = numerator/denominator; 

    //---------------------------------------------------

    //HLLC flux

    if ( 0.0 <= S_L){
      
      Euler::U_state F_L;
      F_L = m.ptr_euler->flux(U_state_L);
      *itflux = F_L;
    }

    if( (S_L <= 0.0) && (S_star >= 0.0)){

      star_coef_left = rho_L*((S_L-u_L)/(S_L-S_star));
      
      //Calculate U_state_L_star
      U_state_L_star.rho = star_coef_left;
      U_state_L_star.momentum = star_coef_left*S_star;
      U_state_L_star.energy = star_coef_left*(U_state_L.energy/U_state_L.rho + (S_star - u_L)*(S_star + P_L/(rho_L*(S_L-u_L))));
      
      Euler::U_state F_L;
      F_L = m.ptr_euler->flux(U_state_L);
      /* Consider overloading the + operator to write this in one line */
      (*itflux).rho = F_L.rho + S_L*(U_state_L_star.rho -U_state_L.rho);
      (*itflux).momentum = F_L.momentum + S_L*(U_state_L_star.momentum -U_state_L.momentum);
      (*itflux).energy = F_L.energy + S_L*(U_state_L_star.energy -U_state_L.energy);
   
							
    }

    if( (S_star <= 0.0) && (S_R >= 0.0)){


      star_coef_right = rho_R*((S_R-u_R)/(S_R-S_star));
      
      //Calculate U_state_R_star
      U_state_R_star.rho = star_coef_right;
      U_state_R_star.momentum = star_coef_right*S_star;
      U_state_R_star.energy = star_coef_right*(U_state_R.energy/U_state_R.rho + (S_star - u_R)*(S_star + P_R/(rho_L*(S_R-u_R))));

      Euler::U_state F_R;    
      F_R = m.ptr_euler->flux(U_state_R);
      /* Consider overloading the + operator to write this in one line */
      (*itflux).rho = F_R.rho + S_R*(U_state_R_star.rho -U_state_R.rho);
      (*itflux).momentum = F_R.momentum + S_R*(U_state_R_star.momentum -U_state_R.momentum);
      (*itflux).energy = F_R.energy + S_R*(U_state_R_star.energy -U_state_R.energy);


    }

    if( 0.0 >= S_R){
      
      Euler::U_state F_R;
      F_R = m.ptr_euler->flux(U_state_R);
      *itflux = F_R;
      
    }
					       

    itflux++;
  }

  return flux;
}

void Mesh_update(Mesh &m, std::vector<Euler::U_state> &flux, double dt){
 
  double dt_dx = dt/m.dx;
  std::vector<Euler::U_state>::iterator itflux = flux.begin();
  
  for(int i = m.nGhost; i <m.ncells+ m.nGhost;i++){

    m.data[i].rho = m.data[i].rho - dt_dx*((*(itflux+1)).rho-((*itflux).rho));
    m.data[i].momentum = m.data[i].momentum - dt_dx*((*(itflux+1)).momentum - (*itflux).momentum);
    m.data[i].energy = m.data[i].energy - dt_dx*((*(itflux+1)).energy-((*itflux).energy));
    
    itflux++;


  }




}
