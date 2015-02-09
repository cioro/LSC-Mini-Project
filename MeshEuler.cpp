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
Mesh::Mesh(){};
Mesh::Mesh(int Ancells, double Ax_min, double Ax_max,double Acfl, Euler::U_state (*f)(double x), Euler::W_state (*b1)(Euler::W_state  w), Euler::W_state (*b2)(Euler::W_state w), int AnGhost) : ncells(Ancells), x_min(Ax_min),x_max(Ax_max),cfl(Acfl), time(0), boundary1(b1), boundary2(b2), nGhost(AnGhost)
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
void Mesh::reset(Euler::U_state (*f)(double x)){

  std::vector<Euler::U_state>::iterator itdata = data.begin()+nGhost;
   std::vector<Euler::U_state>::iterator itdata2= data2.begin()+nGhost;
   std::vector<double>::iterator itaxis= axis.begin()+nGhost;

  
   for(int i = 0; i<ncells; i++){
     (*itaxis) = x_min + i*dx;
     (*itdata) = f(*itaxis);
  
     itaxis++;
     itdata++;
    
     }

};
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
  std::string dir = "data/";
  std::stringstream ss;
  ss << dir << filename << time;
  std::string tmppath = ss.str();
  

  FILE * outfile = fopen(tmppath.c_str(),"w");
  std::vector<Euler::U_state>::const_iterator itdata= data.begin()+nGhost;
  std::vector<double>::const_iterator itaxis= axis.begin()+nGhost;
  for(int i=nGhost; i<ncells+nGhost; i++)
    {
      fprintf(outfile, "%.4f \t %.4f \t %.4f \t %.4f \t %.4f  \n", *itaxis, (*itdata).rho,(*itdata).momentum,(*itdata).energy,ptr_euler->int_energy(ptr_euler->PfromC((*itdata)) ));
      itaxis++;
      itdata++;
  }
  fclose(outfile);
  
}

//Prints to a file the 1D vector of primitive variables and the exact solution 
void Mesh::save_w_state(std::string filename, Euler::U_state (*exact)(double x))const
{
  std::string dir = "data/";
  std::stringstream ss;
  ss << dir << filename << time;
  std::string tmppath = ss.str();
  

  FILE * outfile = fopen(tmppath.c_str(),"w");
  std::vector<Euler::U_state>::const_iterator itdata= data.begin()+nGhost;
  std::vector<double>::const_iterator itaxis= axis.begin()+nGhost;
 
  for(int i=1; i<ncells+nGhost; i++)
    {
      
      Euler::W_state w_print_approx = ptr_euler->PfromC(*itdata);
      //Euler::W_state w_print_exact = ptr_euler->PfromC(*itaxis);
      
      fprintf(outfile, "%.4f \t %.4f \t %.4f \t %.4f \t %.4f \n", *itaxis, w_print_approx.rho,w_print_approx.u,w_print_approx.P, ptr_euler->int_energy(ptr_euler->PfromC((*itdata)) ) );
      itaxis++;
      itdata++;
  }
  fclose(outfile);
  
}

//Implements the boundary conditions. The actual boundary condition function should be in the main file
void Mesh::applyBC(){
  
  for(int j = 0; j < nGhost; j++){

  Euler::W_state w_Left_End = ptr_euler -> PfromC(data[j+nGhost]);
  Euler::W_state w_BC_Left = boundary1(w_Left_End);
  data[nGhost-1-j]= ptr_euler -> CfromP(w_BC_Left);

  }

  for(int i = 0; i < nGhost; i++ ){
    Euler::W_state w_Right_End = ptr_euler -> PfromC(data[(ncells+nGhost-1)-i]);
  Euler::W_state w_BC_Right = boundary2(w_Right_End);
  data[(nGhost + ncells)+i]= ptr_euler -> CfromP(w_BC_Right);
  }
}

//Calculates the adquate size of the time step dt. See page 183 from Toro(ed.2009)
double Mesh::Calculate_dt(){
  std::vector<Euler::U_state>::iterator itdata = data.begin()+nGhost;
  double speed=0.0;
  double speedtemp=0.0;
    
  for(int i=nGhost; i<nGhost+ncells;i++){
    Euler::W_state w = ptr_euler->PfromC(*itdata);
    speedtemp = ptr_euler->a(w) + fabs(w.u);
    if(fabs(speedtemp) > fabs(speed)){
      speed = speedtemp;
    }          
    itdata++;
  }

  double dt;
  //If time < 5 then dt = 0.2
  if(time < 10){
    double  cfl_init = 0.2;
    dt=(cfl_init*dx)/speed;
    
    std::cout << "Inside calculate dt function, cfl = " << cfl_init << "\n"; 

    return dt;

  }
  std::cout << "Inside calculate dt function, cfl = " << cfl << "\n";

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

    // std::cout <<"Inside the HLLC function" << "\n";
        
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
    S_L = u_L -a_L*q_L;
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
      U_state_R_star.energy = star_coef_right*(U_state_R.energy/U_state_R.rho + (S_star - u_R)*(S_star + P_R/(rho_R*(S_R-u_R))));

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

std::vector<Euler::U_state> WAF(Mesh &m, double dt, std::string limiter){
 
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

  //Courant number for each wave speed
  double c_L,c_star,c_R;
 

  //Ratio for each wave speed
  double r_L,r_star,r_R;
  double dq_l,dq_l_right_interface, dq_l_left_interface;
  double dq_star, dq_star_right_interface, dq_star_left_interface;
  double dq_r, dq_r_right_interface, dq_r_left_interface;

  //Limiter functions
  double minmod_l,minmod_star,minmod_r;
  double superbee_l,superbee_star,superbee_r;

  std::vector<Euler::U_state> left_interface;
  std::vector<Euler::U_state> right_interface;
  
  // m.data[m.nGhost].print();
  // m.data[0].print();
  std::cout <<"Inside WAF flux function " << "\n";
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

    // std::cout <<"Inside the HLLC function" << "\n";
        
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
    S_L = u_L -a_L*q_L;
    S_R = u_R + a_R*q_R;
    
    double numerator = P_R - P_L + rho_L*u_L*(S_L-u_L) - rho_R*u_R*(S_R-u_R); 
    double denominator = rho_L*(S_L-u_L)-rho_R*(S_R-u_R);

    S_star = numerator/denominator; 

    //---------------------------------------------------

    //---------------HLLC fluxes -------------------------------------------------//

    //Calculate F_L     
      Euler::U_state F_L;
      F_L = m.ptr_euler->flux(U_state_L);
       

     //Calculate F_L_star
      Euler::U_state F_L_star;
      star_coef_left = rho_L*((S_L-u_L)/(S_L-S_star));
      
      //---Calculate U_state_L_star
      U_state_L_star.rho = star_coef_left;
      U_state_L_star.momentum = star_coef_left*S_star;
      U_state_L_star.energy = star_coef_left*(U_state_L.energy/U_state_L.rho + (S_star - u_L)*(S_star + P_L/(rho_L*(S_L-u_L))));
      /* Consider overloading the + operator to write this in one line */
      F_L_star.rho = F_L.rho + S_L*(U_state_L_star.rho -U_state_L.rho);
      F_L_star.momentum = F_L.momentum + S_L*(U_state_L_star.momentum -U_state_L.momentum);
      F_L_star.energy = F_L.energy + S_L*(U_state_L_star.energy -U_state_L.energy);
   
      //Calculate F_R
      Euler::U_state F_R;
      F_R = m.ptr_euler->flux(U_state_R);
     					

      //Calculate F_R_star
      Euler::U_state F_R_star;
      star_coef_right = rho_R*((S_R-u_R)/(S_R-S_star));
      
      //Calculate U_state_R_star
      U_state_R_star.rho = star_coef_right;
      U_state_R_star.momentum = star_coef_right*S_star;
      U_state_R_star.energy = star_coef_right*(U_state_R.energy/U_state_R.rho + (S_star - u_R)*(S_star + P_R/(rho_R*(S_R-u_R))));


      /* Consider overloading the + operator to write this in one line */
      F_R_star.rho = F_R.rho + S_R*(U_state_R_star.rho -U_state_R.rho);
      F_R_star.momentum = F_R.momentum + S_R*(U_state_R_star.momentum -U_state_R.momentum);
      F_R_star.energy = F_R.energy + S_R*(U_state_R_star.energy -U_state_R.energy);

      //---------------END of HLLC Fluxes, F_L, F_L_star, F_R, F_R_star------------------//


      //Compute the courant number
      c_L = dt*S_L/m.dx;
      c_star = dt*S_star/m.dx;
      c_R = dt*S_R/m.dx;
      

      //Compute the ratios
      //-----compute dq, dq_up, dq_down for each wave at the interface. Up for 
      
      //dq for the interface we are trying to solve for i+1/2
      dq_l = U_state_L.rho-U_state_L_star.rho;
      dq_star = U_state_L_star.rho - U_state_R_star.rho;
      dq_r = U_state_R_star.rho - U_state_R.rho;

      left_interface = HLLC_U_state(m.data[i-1], m.data[i]);
      right_interface = HLLC_U_state(m.data[i+1],m.data[i+2]);

      dq_l_left_interface = left_interface[0].rho-left_interface[1].rho;
      dq_star_left_interface = left_interface[1].rho - left_interface[2].rho;
      dq_r_left_interface = left_interface[2].rho - left_interface[3].rho;
      
      dq_l_right_interface = right_interface[0].rho-right_interface[1].rho;
      dq_star_right_interface = right_interface[1].rho - right_interface[2].rho;
      dq_r_right_interface = right_interface[2].rho - right_interface[3].rho;
     
      /*
      std::cout << "This is the " << i << " element" <<"\n";

      std::cout << "The 4 states at the interface are: " << "\n";
      U_state_L.print();
      U_state_L_star.print();
      U_state_R_star.print();
      U_state_R.print();
      std::cout << "The jumps are dq_l,dq_star ,dq_r " << dq_l <<"\t" << dq_star << "\t" << dq_r <<"\n";
      std::cout << "The nearby states are " << "\n";
      std::cout << "The interface to the left of i+1/2 " << "\n";
      std::cout << "U_L : " << "\n";
      left_interface[0].print();
      std::cout << "U_L_star : " << "\n";
      left_interface[1].print();
      std::cout << "U_R_star : " << "\n";
      left_interface[2].print();
      std::cout << "U_R : " << "\n";
      left_interface[3].print();


      std::cout << "The interface to the right of i+1/2 " << "\n";
      std::cout << "U_L : " << "\n";
      right_interface[0].print();
      std::cout << "U_L_star : " << "\n";
      right_interface[1].print();
      std::cout << "U_R_star : " << "\n";
      right_interface[2].print();
      std::cout << "U_R : " << "\n";
      right_interface[3].print();
      
      std::cout << "\n";
      */
   
     
      
      if(dq_l > 0){
	if (c_L > 0){
	  r_L = dq_l_left_interface/dq_l;
	}
	else{
	  r_L = dq_l_right_interface/dq_l;
	}
      }
      else{
	minmod_l = 0.0;
	superbee_l = 0.0;
	r_L = 0.0;
      }

	if(dq_star > 0){

	  if (c_star > 0){
	    r_star = dq_star_left_interface/dq_star;
	  }
	  else{
	    r_star = dq_star_right_interface/dq_star;
	  }
	}
      else{
	minmod_star = 0.0;
	superbee_star = 0.0;
	r_star = 0.0;
      }

	  if(dq_r > 0){

	    if (c_R > 0){
	      r_R = dq_r_left_interface/dq_r;
	    }
	    else{
	      r_R = dq_r_right_interface/dq_r;
	    }
	  }
	  else{
	    minmod_r = 0.0;
	    superbee_r = 0.0;
	    r_R = 0.0;
	  } 
   
      //Compute limiter functions
      minmod_l = minmod(r_L,fabs(c_L));
      minmod_star = minmod(r_star, fabs(c_star));
      minmod_r = minmod(r_R,fabs(c_R));

      superbee_l = superbee(r_L, fabs(c_L));
      superbee_star = superbee(r_star, fabs(c_star));
      superbee_r = superbee(r_R, fabs(c_R));

       std::cout << "The ratios are r_L, r_star, r_R " << r_L<<"\t" << r_star << "\t" << r_R <<"\n";
       std::cout << "The superbee lim m_l m_star m_r  " << superbee_l <<"\t" <<superbee_star << "\t" << superbee_r <<"\n";  
      //Compute intercell flux at i+1/2 (*itflux).rho, (*itflux).momentum, (*itflux).energy
      //For 1D N= 3, total number of waves
      
       if(limiter == std::string("minmod")){

	 (*itflux).rho = 0.5*(F_L.rho + F_R.rho) - 0.5* (sign(c_L)*minmod_l*(F_L_star.rho-F_L.rho) + \
						     sign(c_star)*minmod_star*(F_R_star.rho-F_L_star.rho) + \
						     sign(c_R)*minmod_r*(F_R.rho-F_R_star.rho));

	 (*itflux).momentum = 0.5*(F_L.momentum + F_R.momentum) - 0.5* (sign(c_L)*minmod_l*(F_L_star.momentum-F_L.momentum) + \
						     sign(c_star)*minmod_star*(F_R_star.momentum-F_L_star.momentum) + \
						     sign(c_R)*minmod_r*(F_R.momentum-F_R_star.momentum));
						     
	 (*itflux).energy = 0.5*(F_L.energy + F_R.energy) - 0.5* (sign(c_L)*minmod_l*(F_L_star.energy-F_L.energy) + \
						     sign(c_star)*minmod_star*(F_R_star.energy-F_L_star.energy) + \
						     sign(c_R)*minmod_r*(F_R.energy-F_R_star.energy));

       } else  if(limiter == std::string("superbee")){
	 (*itflux).rho = 0.5*(F_L.rho + F_R.rho) - 0.5* (sign(c_L)*superbee_l*(F_L_star.rho-F_L.rho) + \
						     sign(c_star)*superbee_star*(F_R_star.rho-F_L_star.rho) + \
						     sign(c_R)*superbee_r*(F_R.rho-F_R_star.rho));

	 (*itflux).momentum = 0.5*(F_L.momentum + F_R.momentum) - 0.5* (sign(c_L)*superbee_l*(F_L_star.momentum-F_L.momentum) + \
						     sign(c_star)*superbee_star*(F_R_star.momentum-F_L_star.momentum) + \
						     sign(c_R)*superbee_r*(F_R.momentum-F_R_star.momentum));
						     
	 (*itflux).energy = 0.5*(F_L.energy + F_R.energy) - 0.5* (sign(c_L)*superbee_l*(F_L_star.energy-F_L.energy) + \
						     sign(c_star)*superbee_star*(F_R_star.energy-F_L_star.energy) + \
						     sign(c_R)*superbee_r*(F_R.energy-F_R_star.energy));

       }else{
	     std::cout << "no limiter speficied. Fail to compute WAF fluxes"<< "\n";
	     (*itflux).rho = 0.0;
	     (*itflux).momentum = 0.0;
	     (*itflux).energy = 0.0;
       }
      //For debugging use empty flux
      /*
      (*itflux).rho = 0.0;
      (*itflux).momentum = 0.0;
      (*itflux).energy = 0.0;
      */
      //---------------------//      
      itflux++;
  }
     
  return flux;
      
      
}

int sign(double c){
  
  if(c >= 0){
    return 1;
  }
  else{
    return -1;
  }

}
double minmod(double r, double c){
  if (r <= 0){
    return 1.0;
  }
  if( r >= 0 && r <= 1){
    return (1-(1-c)*r);
  }
  if (r > 1 ){
    return c;
  }
  else{
    std::cout  << "Error in minmod limiter calculation " << r << "\n";
    return 0;
  }
}

double superbee(double r, double c){

  if(r <= 0){
    return 1;
  }
  else{
    double top = (1-c)*2*r;
    return (1 - top/(1+r));
  }

}

std::vector<Euler::U_state> HLLC_U_state(Euler::U_state U_state_L, Euler::U_state U_state_R){

  std::vector<Euler::U_state> U_interface(4);//vector of size 4 with U_L,U_L_star,U_R_star,U_star;
  
  U_interface[0]=U_state_L;
  U_interface[3]=U_state_R;

  Euler e;
  const double gamma = e.gamma;//This is a massive fudge!!!!! 
  
  //Variables used
  double P_star,P_L,P_R; //Presures
  double rho_star,rho_L,rho_R;//densities
  double u_L,u_R;//particle/gas speed in cell

  double a_L, a_R;// sound speed in cell
  double S_L, S_R,S_star;// wave speed in cell left wave, right wave, contact wave
    
  double P_pvrs;
  double rho_bar, a_bar;//average density and average sound speed of left and right cells

  Euler::W_state w_temp_left;
  Euler::W_state w_temp_right;

  double q_R,q_L;

 
  Euler::W_state W_L;
  Euler::W_state W_R;
  Euler::U_state U_state_L_star;
  Euler::U_state U_state_R_star;

  double star_coef_left; // The coeficient in eq. 10.73 from Toro(ed.2009); (S_k-u_k)/(S_k-u_star_k);
  double star_coef_right; // The coeficient in eq. 10.73 from Toro(ed.2009);
  
  
  W_L = e.PfromC(U_state_L);
  W_R = e.PfromC(U_state_R);


    //---------Pressure estimate-------------------------

    P_L = W_L.P;
    P_R = W_R.P;

    u_L = W_L.u;
    u_R = W_R.u;

    rho_L = W_L.rho;
    rho_R = W_R.rho;
   
    a_L = e.a(W_L);
    a_R = e.a(W_R);

    // std::cout <<"Inside the HLLC function" << "\n";
        
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
    S_L = u_L -a_L*q_L;
    S_R = u_R + a_R*q_R;
    
    double numerator = P_R - P_L + rho_L*u_L*(S_L-u_L) - rho_R*u_R*(S_R-u_R); 
    double denominator = rho_L*(S_L-u_L)-rho_R*(S_R-u_R);

    S_star = numerator/denominator; 

    //---------------------------------------------------

    //---------------HLLC fluxes -------------------------------------------------//

    //Calculate F_L     
    // Euler::U_state F_L;
    // F_L = m.ptr_euler->flux(U_state_L);
       

     //Calculate F_L_star
     // Euler::U_state F_L_star;
      star_coef_left = rho_L*((S_L-u_L)/(S_L-S_star));
      
      //---Calculate U_state_L_star
      U_state_L_star.rho = star_coef_left;
      U_state_L_star.momentum = star_coef_left*S_star;
      U_state_L_star.energy = star_coef_left*(U_state_L.energy/U_state_L.rho + (S_star - u_L)*(S_star + P_L/(rho_L*(S_L-u_L))));
      U_interface[1]=U_state_L_star;

      // /* Consider overloading the + operator to write this in one line */
	// F_L_star.rho = F_L.rho + S_L*(U_state_L_star.rho -U_state_L.rho);
      //F_L_star.momentum = F_L.momentum + S_L*(U_state_L_star.momentum -U_state_L.momentum);
      //F_L_star.energy = F_L.energy + S_L*(U_state_L_star.energy -U_state_L.energy);
   
      //Calculate F_R
      // Euler::U_state F_R;
      // F_R = m.ptr_euler->flux(U_state_R);
     					

      //Calculate F_R_star
      // Euler::U_state F_R_star;
      star_coef_right = rho_R*((S_R-u_R)/(S_R-S_star));
      
      //Calculate U_state_R_star
      U_state_R_star.rho = star_coef_right;
      U_state_R_star.momentum = star_coef_right*S_star;
      U_state_R_star.energy = star_coef_right*(U_state_R.energy/U_state_R.rho + (S_star - u_R)*(S_star + P_R/(rho_R*(S_R-u_R))));
      U_interface[2]=U_state_R_star;

      /* Consider overloading the + operator to write this in one line */
      // F_R_star.rho = F_R.rho + S_R*(U_state_R_star.rho -U_state_R.rho);
      // F_R_star.momentum = F_R.momentum + S_R*(U_state_R_star.momentum -U_state_R.momentum);
      // F_R_star.energy = F_R.energy + S_R*(U_state_R_star.energy -U_state_R.energy);

      return U_interface;

}
