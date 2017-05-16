

// evolution.h
//  
// Made by Pablo Galaviz
// e-mail  <Pablo.Galaviz@me.com>
// 



//  This file is part of PostNewtonian3BP
//
//  PostNewtonian3BP is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  any later version.
//
//  PostNewtonian3BP is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with PostNewtonian3BP.  If not, see <http://www.gnu.org/licenses/>.
//



#ifndef EVOLUTION_H

#define EVOLUTION_H

#include<utils.h>

class evolution
{

  src::severity_logger< severity_level > lg;

  bool verbose; 

  gsl_odeiv_evolve * evolve;

  gsl_odeiv_control * control;
  
  gsl_odeiv_step * step;

  gsl_odeiv_system system;

  double time; 

  double final_time; 

  double dt; 

  valarray<double> y;
  
  double dx; 
  
  size_t space_dimension;
  
  size_t number_of_variables; 
   
  size_t number_of_particles; 
  
  valarray<double> dy; 

  valarray<double> par;

  valarray<double> position;
  valarray<double> momentum;
  valarray<double> spin;

  
  static int rhs(double t, const double * y, double * f, void * param);

  int rhsN(double t, const double * y, double * f, void * param);

  int rhs1PN(double t, const double * y, double * f, void * param);

  int rhs2PN(double t, const double * y, double * f, void * param);

  int rhs2_5PN(double t, const double * y, double * f, void * param);

  int rhsSOloPN(double t, const double * y, double * f, void * param);

  int rhsSSloPN(double t, const double * y, double * f, void * param);

  int rhsS2loPN(double t, const double * y, double * f, void * param);

  int rhsSOnloPN(double t, const double * y, double * f, void * param);

  static int jac(double t, const double y[], double * dfdy, double dfdt[], void * param);

  int jacN(double t, const double y[], double * dfdy, double dfdt[], void * param);
  int jac1PN(double t, const double y[], double * dfdy, double dfdt[], void * param);
  int jac2PN(double t, const double y[], double * dfdy, double dfdt[], void * param);
  int jac2_5PN(double t, const double y[], double * dfdy, double dfdt[], void * param);
  int jacSOloPN(double t, const double y[], double * dfdy, double dfdt[], void * param);
  int jacSSloPN(double t, const double y[], double * dfdy, double dfdt[], void * param);
  int jacS2loPN(double t, const double y[], double * dfdy, double dfdt[], void * param);
  int jacSOnloPN(double t, const double y[], double * dfdy, double dfdt[], void * param);

  static int rhsC(double t, const double * y, double * f, void * param);

  
  typedef double (evolution::*pt2rhs) (double t, const double * y, double * param, int i);
  typedef double (evolution::*pt2jac) (double t, const double * y, double * param, int i, int j);

  pt2rhs rhsN_F; 
  pt2rhs rhs1PN_F; 
  pt2rhs rhs2PN_F; 
  pt2rhs rhs2_5PN_F; 
  pt2rhs rhsSOloPN_F; 
  pt2rhs rhsSSloPN_F; 
  pt2rhs rhsS2loPN_F; 
  pt2rhs rhsSOnloPN_F; 


  pt2jac jacN_F; 
  pt2jac jac1PN_F; 
  pt2jac jac2PN_F; 
  pt2jac jac2_5PN_F; 
  pt2jac jacSOloPN_F; 
  pt2jac jacSSloPN_F; 
  pt2jac jacS2loPN_F; 
  pt2jac jacSOnloPN_F; 

  inline int RHSN(double t, const double * y, double * f, void * param) 
  {return ((this->rhsN) (t,y,f,param));}
  
  inline int RHS1PN(double t, const double * y, double * f, void * param) 
  {return ((this->rhs1PN) (t,y,f,param));}

  inline int RHS2PN(double t, const double * y, double * f, void * param) 
  {return ((this->rhs2PN) (t,y,f,param));}

  inline int RHS2_5PN(double t, const double * y, double * f, void * param) 
  {return ((this->rhs2_5PN) (t,y,f,param));}

  inline int RHSSOloPN(double t, const double * y, double * f, void * param) 
  {return ((this->rhsSOloPN) (t,y,f,param));}

  inline int RHSSSloPN(double t, const double * y, double * f, void * param) 
  {return ((this->rhsSSloPN) (t,y,f,param));}

  inline int RHSS2loPN(double t, const double * y, double * f, void * param) 
  {return ((this->rhsS2loPN) (t,y,f,param));}

  inline int RHSSOnloPN(double t, const double * y, double * f, void * param) 
  {return ((this->rhsSOnloPN) (t,y,f,param));}


  double rhs_zero(double t, const double y[], double *param, int i)
  {return (0.0);};

  double jac_zero(double t, const double y[], double *param, int i,int j)
  {return (0.0);};

  
  inline int JACN(double t, const double y[], double * dfdy, double dfdt[], void * param) 
  {return ((this->jacN) (t,y,dfdy,dfdt,param));}
  
  inline int JAC1PN(double t, const double y[], double * dfdy, double dfdt[], void * param) 
  {return ((this->jac1PN) (t,y,dfdy,dfdt,param));}

  inline int JAC2PN(double t, const double y[], double * dfdy, double dfdt[], void * param) 
  {return ((this->jac2PN) (t,y,dfdy,dfdt,param));}

  inline int JAC2_5PN(double t, const double y[], double * dfdy, double dfdt[], void * param)
  {return ((this->jac2_5PN) (t,y,dfdy,dfdt,param));}

  inline int JACSOloPN(double t, const double y[], double * dfdy, double dfdt[], void * param)
  {return ((this->jacSOloPN) (t,y,dfdy,dfdt,param));}

  inline int JACSSloPN(double t, const double y[], double * dfdy, double dfdt[], void * param)
  {return ((this->jacSSloPN) (t,y,dfdy,dfdt,param));}

  inline int JACS2loPN(double t, const double y[], double * dfdy, double dfdt[], void * param)
  {return ((this->jacS2loPN) (t,y,dfdy,dfdt,param));}


  inline int JACSOnloPN(double t, const double y[], double * dfdy, double dfdt[], void * param)
  {return ((this->jacSOnloPN) (t,y,dfdy,dfdt,param));}


  //Number of bodies np=2
  //Dimenssion dim=2
  //Newtonian equations
  double rhs_Nnp2d2_FNS(double t, const double * y, double * param,int i);
  double jac_Nnp2d2_FaNS(double t, const double * y, double * param,int i, int j);
  //1 post-Newtonian equations
  double rhs_1PNnp2d2_FNS(double t, const double * y, double * param,int i);
  double jac_1PNnp2d2_FaNS(double t, const double * y, double * param,int i, int j);
  //2 post-Newtonian equations
  double rhs_2PNnp2d2_FNS(double t, const double * y, double * param,int i);
  //2.5 post-Newtonian equations
  double rhs_2_5PNnp2d2_FNS(double t, const double * y, double * param,int i);


  //Dimenssion dim=3 // include spinning particles
  //Newtonian equations
  double rhs_Nnp2d3_FS(double t, const double * y, double * param,int i);
  double jac_Nnp2d3_FaS(double t, const double * y, double * param,int i, int j);
  //1 post-Newtonian equations
  double rhs_1PNnp2d3_FS(double t, const double * y, double * param,int i);
  double jac_1PNnp2d3_FaS(double t, const double * y, double * param,int i, int j);
  //2 post-Newtonian equations
  double rhs_2PNnp2d3_FS(double t, const double * y, double * param,int i);
  //2.5 post-Newtonian equations
  double rhs_2_5PNnp2d3_FS(double t, const double * y, double * param,int i);

  //Spin-Orbit leading order post-Newtonian equations
  double rhs_SOloPNnp2d3_FS(double t, const double * y, double * param,int i);
  double jac_SOloPNnp2d3_FaS(double t, const double * y, double * param,int i, int j);

  //Spin(a)-Spin(b) leading order post-Newtonian equations
  double rhs_SSloPNnp2d3_FS(double t, const double * y, double * param,int i);
  double jac_SSloPNnp2d3_FaS(double t, const double * y, double * param,int i, int j);

  //Spin(a)-Spin(a) leading order post-Newtonian equations
  double rhs_S2loPNnp2d3_FS(double t, const double * y, double * param,int i);
  double jac_S2loPNnp2d3_FaS(double t, const double * y, double * param,int i, int j);

  //Spin-Orbit next to the leading order post-Newtonian equations
  double rhs_SOnloPNnp2d3_FS(double t, const double * y, double * param,int i);


  //Number of bodies np=3
  //Dimenssion dim=2
  //Newtonian equations
  double rhs_Nnp3d2_FNS(double t, const double * y, double * param,int i);
  double jac_Nnp3d2_FaNS(double t, const double * y, double * param,int i, int j);
  //1 post-Newtonian equations
  double rhs_1PNnp3d2_FNS(double t, const double * y, double * param,int i);
  double jac_1PNnp3d2_FaNS(double t, const double * y, double * param,int i, int j);
  //2 post-Newtonian equations
  double rhs_2PNnp3d2_FNS(double t, const double * y, double * param,int i);
  //2.5 post-Newtonian equations
  double rhs_2_5PNnp3d2_FNS(double t, const double * y, double * param,int i);

  //Dimenssion dim=3 // include spinning particles
  //Newtonian equations
  double rhs_Nnp3d3_FS(double t, const double * y, double * param,int i);
  double jac_Nnp3d3_FaS(double t, const double * y, double * param,int i, int j);
  //1 post-Newtonian equations
  double rhs_1PNnp3d3_FS(double t, const double * y, double * param,int i);
  double jac_1PNnp3d3_FaS(double t, const double * y, double * param,int i, int j);
  //2 post-Newtonian equations
  double rhs_2PNnp3d3_FS(double t, const double * y, double * param,int i);
  //2.5 post-Newtonian equations
  double rhs_2_5PNnp3d3_FS(double t, const double * y, double * param,int i);

  //Spin-Orbit leading order post-Newtonian equations
  double rhs_SOloPNnp3d3_FS(double t, const double * y, double * param,int i);
  double jac_SOloPNnp3d3_FaS(double t, const double * y, double * param,int i, int j);

  //Spin(a)-Spin(b) leading order post-Newtonian equations
  double rhs_SSloPNnp3d3_FS(double t, const double * y, double * param,int i);
  double jac_SSloPNnp3d3_FaS(double t, const double * y, double * param,int i, int j);

  //Spin(a)-Spin(a) leading order post-Newtonian equations
  double rhs_S2loPNnp3d3_FS(double t, const double * y, double * param,int i);
  double jac_S2loPNnp3d3_FaS(double t, const double * y, double * param,int i, int j);


 //Spin-Orbit next to the leading order post-Newtonian equations
  double rhs_SOnloPNnp3d3_FS(double t, const double * y, double * param,int i);


  //Dimenssion dim=3 // non spinning particles
  //Newtonian equations
  double rhs_Nnp3d3_FNS(double t, const double * y, double * param,int i);
  double jac_Nnp3d3_FaNS(double t, const double * y, double * param,int i, int j);
  //1 post-Newtonian equations
  double rhs_1PNnp3d3_FNS(double t, const double * y, double * param,int i);
  double jac_1PNnp3d3_FaNS(double t, const double * y, double * param,int i, int j);
  //2 post-Newtonian equations
  double rhs_2PNnp3d3_FNS(double t, const double * y, double * param,int i);
  //2.5 post-Newtonian equations
  double rhs_2_5PNnp3d3_FNS(double t, const double * y, double * param,int i);


  //Numerical Jacobians
  double jac_numN_F(double t, const double y[], double par[], int i, int j);
  double jac_num1PN_F(double t, const double y[], double par[], int i, int j);  
  double jac_num2PN_F(double t, const double y[], double par[], int i, int j);
  double jac_num2_5PN_F(double t, const double y[], double par[], int i, int j);
  double jac_numSOloPN_F(double t, const double y[], double par[], int i, int j);
  double jac_numSSloPN_F(double t, const double y[], double par[], int i, int j);
  double jac_numS2loPN_F(double t, const double y[], double par[], int i, int j);
  double jac_numSOnloPN_F(double t, const double y[], double par[], int i, int j);
 
  
  void set_rhs(terms_t _pn_terms, bool jacA, bool _chaos_test, bool _spin);


  
 public:

  size_t NDvar; 

  
  evolution(bool _evolution_verbose,
	    string _ode_method,
	    bool _chaos_test,
	    bool _JacA,
	    double _initial_dt,
	    double _Jacobian_dx,
	    double _factor_chaos_test,
	    double _scaling_variable,
	    double _scaling_derivative,
	    double _eps_rel,
	    double _eps_abs,
	    double _final_time,
	    valarray<double> &_id_variables,
	    valarray<double> &_mass,
	    bool _spin, 
	    terms_t _pn_terms
	    );

 ~evolution(){};

 bool update(double &t);
 

 void close();

 valarray<double> get_position(){

   int i=0;
   for(int a=0; a<number_of_particles; a++)
     for(int axis=0; axis<space_dimension; axis++)
       {
	 position[i]=y[r_index(a,axis,space_dimension)];
	 i++;
       }
	 
   return position;
 }

  valarray<double> get_momentum(){
    int i=0;
    for(int a=0; a<number_of_particles; a++)
      for(int axis=0; axis<space_dimension; axis++)
	{
	  momentum[i]=y[p_index(a,axis,space_dimension)];
	 i++;
	}
    
   return momentum; }

  valarray<double> get_spin(){

    int i=0;
    
    if(spin.size()>0)
    for(int a=0; a<number_of_particles; a++)
      for(int axis=0; axis<space_dimension; axis++)
	{
	  spin[i]=y[s_index(a,axis,space_dimension)];
	 i++;
	}
    
    
   return spin; }

 
};


#endif
