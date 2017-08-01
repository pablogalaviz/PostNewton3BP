

// evolution.cc
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



#include "evolution.h"

void* pt2evolution;


evolution::evolution(
		     bool _evolution_verbose,
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
		     )
{

  pt2evolution =(void*) this;
  
  verbose = _evolution_verbose;

  iterations=0;  
  dt = _initial_dt;
  final_time = _final_time;
  dx = _Jacobian_dx;

  waves.resize(39);
  
  size_t number_of_particles = _mass.size();

  const gsl_odeiv_step_type * step_type;
  
  if( _ode_method.compare("rk2") == 0 ) 
    step_type = gsl_odeiv_step_rk2;
  else if( _ode_method.compare("rk4") == 0 ) 
    step_type = gsl_odeiv_step_rk4;
  else if( _ode_method.compare("rkck") == 0 ) 
    step_type = gsl_odeiv_step_rkck;
  else if( _ode_method.compare("rk8pd") == 0 ) 
    step_type = gsl_odeiv_step_rk8pd;
  else if( _ode_method.compare("rk2imp") == 0 ) 
    step_type = gsl_odeiv_step_rk2imp;
  else if( _ode_method.compare("rk4imp") == 0 ) 
    step_type = gsl_odeiv_step_rk4imp;
  else if( _ode_method.compare("bsimp") == 0 ) 
    step_type = gsl_odeiv_step_bsimp;
  else if( _ode_method.compare("gear1") == 0 ) 
    step_type = gsl_odeiv_step_gear1;
  else if( _ode_method.compare("gear2") == 0 ) 
    step_type = gsl_odeiv_step_gear2;
  else
    step_type = gsl_odeiv_step_rkf45;

  time =0; 
  
  if (number_of_particles < 2 || number_of_particles > 3)
    {
      
      BOOST_LOG_SEV(lg, error) << "evolution.particles  should be  2 or 3.";
      BOOST_LOG_SEV(lg, error) << "Current value : "<< number_of_particles;

      exit(Finalize(0));
    }

  bool plane_constrain;
    
  if(number_of_particles==2 && !_spin)
    plane_constrain=true;
  else
    if(plane_constrain && _spin)
      plane_constrain=false;

  number_of_variables = _id_variables.size();
    
  size_t ode_size = number_of_variables;

  NDvar = number_of_variables; 

  
 
  if(_chaos_test) {

    ode_size *= 2;

    double scale_abs [ode_size];
    
    for(int i = 0; i < number_of_variables; i++)
      scale_abs[i] = 1.0;

    for(int i = number_of_variables; i < ode_size; i++)
      scale_abs[i] = _factor_chaos_test;

    control = gsl_odeiv_control_scaled_new (
					    _eps_abs,
					    _eps_rel,
					    _scaling_variable,
					    _scaling_derivative,
					    scale_abs,
					    ode_size);

  }
  else
    control = gsl_odeiv_control_standard_new (
					    _eps_abs,
					    _eps_rel,
					    _scaling_variable,
					    _scaling_derivative);


  step = gsl_odeiv_step_alloc (step_type, ode_size);

  evolve = gsl_odeiv_evolve_alloc (ode_size);


  
  if(verbose)
    {
      BOOST_LOG_SEV(lg, info)  << "Evolution info:";
      BOOST_LOG_SEV(lg, info)  << "======================================";
      if(plane_constrain)
	BOOST_LOG_SEV(lg,info)  << "Evolution in 2D";
      if(_chaos_test)
	BOOST_LOG_SEV(lg,info)  << "Testing for chaos";
      BOOST_LOG_SEV(lg, info)  << "Method: "  << gsl_odeiv_step_name(step);
      BOOST_LOG_SEV(lg, info)  << "Control: "  << gsl_odeiv_control_name(control);
      BOOST_LOG_SEV(lg, info)  << "Number of particles: "  << number_of_particles;
      BOOST_LOG_SEV(lg, info)  << "Final time: "  << final_time;
      BOOST_LOG_SEV(lg, info)  << "======================================";
    }

  y.resize(ode_size);
  for(int i =0; i < _id_variables.size(); i++)
    y[i]=_id_variables[i];
  
  par.resize(number_of_particles);
  for(int i =0; i < _mass.size(); i++)
    par[i]=_mass[i];

  dy.resize(ode_size);
  ym1.resize(ode_size);
  ym2.resize(ode_size);
  ddy.resize(ode_size);
  
  system.dimension = ode_size;

  system.params = &par[0];
  
  
  dy.resize(ode_size);
        

  set_rhs(_pn_terms, _JacA,_chaos_test,_spin);

  
};

void evolution::set_rhs(terms_t _pn_terms, bool jacA,bool _chaos_test, bool _spin )
{

  number_of_particles = par.size();

  space_dimension = _spin ? number_of_variables/(3*number_of_particles) : number_of_variables/(2*number_of_particles);

  position.resize(number_of_particles*space_dimension+1);
  momentum.resize(number_of_particles*space_dimension+1);
  dxdt.resize(number_of_particles*space_dimension+1);
  dpdt.resize(number_of_particles*space_dimension+1);
  ddxdt2.resize(number_of_particles*space_dimension+1);
  ddpdt2.resize(number_of_particles*space_dimension+1);
  if(_spin){
    spin.resize(number_of_particles*space_dimension+1);
    dsdt.resize(number_of_particles*space_dimension+1);
    ddsdt2.resize(number_of_particles*space_dimension+1);
  }
  else{
    spin.resize(1);
    dsdt.resize(1);
    ddsdt2.resize(1);
  }
  rhs1PN_F = &evolution::rhs_zero;
  jac1PN_F = &evolution::jac_zero;

  rhs2PN_F = &evolution::rhs_zero;
  jac2PN_F = &evolution::jac_zero;

  rhs2_5PN_F = &evolution::rhs_zero;
  jac2_5PN_F = &evolution::jac_zero;

  rhsSOloPN_F = &evolution::rhs_zero;
  jacSOloPN_F = &evolution::jac_zero;

  rhsSSloPN_F = &evolution::rhs_zero;
  jacSSloPN_F = &evolution::jac_zero;
 

  rhsS2loPN_F = &evolution::rhs_zero;
  jacS2loPN_F = &evolution::jac_zero;
 
  rhsSOnloPN_F = &evolution::rhs_zero;
  jacSOnloPN_F = &evolution::jac_zero;

  
  
  //Newtonian
  jacN_F= &evolution::jac_numN_F;
  
  if(number_of_particles==2 && space_dimension == 2)
    {
      rhsN_F = &evolution::rhs_Nnp2d2_FNS;
      if(jacA)
	jacN_F =  &evolution::jac_Nnp2d2_FaNS;
    }
  else     
    if(number_of_particles==2 && space_dimension ==3)
      {
	rhsN_F = &evolution::rhs_Nnp2d3_FS;
	BOOST_LOG_SEV(lg, debug)  << "Set rhs_Nnp2d3_FS";
	
	if(jacA)
	  jacN_F =  &evolution::jac_Nnp2d3_FaS;
      }
    else     
      if(number_of_particles==3 && space_dimension ==2)
	{
	  rhsN_F = &evolution::rhs_Nnp3d2_FNS;
	  BOOST_LOG_SEV(lg, debug)  << "Set rhs_Nnp3d2_FNS";

	  if(jacA)
	    jacN_F =  &evolution::jac_Nnp3d2_FaNS;
	}
      else     
	if(number_of_particles==3 && space_dimension ==3 && _spin)
	  {
	    rhsN_F = &evolution::rhs_Nnp3d3_FS;
	    BOOST_LOG_SEV(lg, debug)  << "Set rhs_Nnp3d3_FS";
	    if(jacA)
	      jacN_F =  &evolution::jac_Nnp3d3_FaS;
	  }
	else  //	if(number_of_particles==3 && space_dimension ==3 )
	  {
	    rhsN_F = &evolution::rhs_Nnp3d3_FNS;
	    BOOST_LOG_SEV(lg, debug)  << "Set rhs_Nnp3d3_FNS";
	    if(jacA)
	      jacN_F =  &evolution::jac_Nnp3d3_FaNS;
	  }
     
#ifdef odePN1
  if(_pn_terms.pn1)
    {
      jac1PN_F = &evolution::jac_num1PN_F;

      if(number_of_particles==2 && space_dimension == 2)
	{
	  rhs1PN_F = &evolution::rhs_1PNnp2d2_FNS;
	  if(jacA)
	    jac1PN_F =  &evolution::jac_1PNnp2d2_FaNS;
	}
      else     
	if(number_of_particles==2 && space_dimension ==3)
	  {
	    rhs1PN_F = &evolution::rhs_1PNnp2d3_FS;
	    if(jacA)
	      jac1PN_F =  &evolution::jac_1PNnp2d3_FaS;
	  }
	else     
	  if(number_of_particles==3 && space_dimension ==2)
	    {
	      rhs1PN_F = &evolution::rhs_1PNnp3d2_FNS;
	      if(jacA)
		jac1PN_F =  &evolution::jac_1PNnp3d2_FaNS;
	    }
	  else     
	    if(number_of_particles==3 && space_dimension ==3 && _spin)
	      {
		rhs1PN_F = &evolution::rhs_1PNnp3d3_FS;
		if(jacA)
		  jac1PN_F =  &evolution::jac_1PNnp3d3_FaS;
	  }
	else  //	if(number_of_particles==3 && space_dimension ==3 )
	  {
	    rhs1PN_F = &evolution::rhs_1PNnp3d3_FNS;
	    if(jacA)
	      jac1PN_F =  &evolution::jac_1PNnp3d3_FaNS;
	  }
            
    }
#endif

#ifdef odePN2
  if(_pn_terms.pn2)
    {
      jac2PN_F = &evolution::jac_num2PN_F;

      if(number_of_particles==2 && space_dimension == 2)
	{
	  rhs2PN_F = &evolution::rhs_2PNnp2d2_FNS;
	  if(jacA)
	    jac2PN_F =  &evolution::jac_2PNnp2d2_FaNS;
	}
      else     
	if(number_of_particles==2 && space_dimension ==3)
	  {
	    rhs2PN_F = &evolution::rhs_2PNnp2d3_FS;
	    if(jacA)
	      jac2PN_F =  &evolution::jac_2PNnp2d3_FaS;
	  }
	else     
	  if(number_of_particles==3 && space_dimension ==2)
	    {
	      rhs2PN_F = &evolution::rhs_2PNnp3d2_FNS;
	      if(jacA)
		jac2PN_F =  &evolution::jac_2PNnp3d2_FaNS;
	    }
	  else     
	    if(number_of_particles==3 && space_dimension ==3 && _spin)
	      {
		rhs2PN_F = &evolution::rhs_2PNnp3d3_FS;
		if(jacA)
		  jac2PN_F =  &evolution::jac_2PNnp3d3_FaS;
	  }
	else  //	if(number_of_particles==3 && space_dimension ==3 )
	  {
	    rhs2PN_F = &evolution::rhs_2PNnp3d3_FNS;
	    if(jacA)
	      jac2PN_F =  &evolution::jac_2PNnp3d3_FaNS;
	  }
            
    }
#endif

#ifdef odePN2_2
  if(_terms.pn2_5)
    {
      jac2_5PN_F = &evolution::jac_num2_5PN_F;

      if(number_of_particles==2 && space_dimension == 2)
	{
	  rhs2_5PN_F = &evolution::rhs_2_5PNnp2d2_FNS;
	  if(jacA)
	    jac2_5PN_F =  &evolution::jac_2_5PNnp2d2_FaNS;
	}
      else     
	if(number_of_particles==2 && space_dimension ==3)
	  {
	    rhs2_5PN_F = &evolution::rhs_2_5PNnp2d3_FS;
	    if(jacA)
	      jac2_5PN_F =  &evolution::jac_2_5PNnp2d3_FaS;
	  }
	else     
	  if(number_of_particles==3 && space_dimension ==2)
	    {
	      rhs2_5PN_F = &evolution::rhs_2_5PNnp3d2_FNS;
	      if(jacA)
		jac2_5PN_F =  &evolution::jac_2_5PNnp3d2_FaNS;
	    }
	  else     
	    if(number_of_particles==3 && space_dimension ==3 && _spin)
	      {
		rhs2_5PN_F = &evolution::rhs_2_5PNnp3d3_FS;
		if(jacA)
		  jac2_5PN_F =  &evolution::jac_2_5PNnp3d3_FaS;
	  }
	else  //	if(number_of_particles==3 && space_dimension ==3 )
	  {
	    rhs2_5PN_F = &evolution::rhs_2_5PNnp3d3_FNS;
	    if(jacA)
	      jac2_5PN_F =  &evolution::jac_2_5PNnp3d3_FaNS;
	  }
            
    }
#endif

#ifdef odePNSlo
  if(_terms.pnSOlo)
    {
    jacSOloPN_F = &evolution::jac_numSOloPN_F;
    
    if(number_of_particles==2 )
      {
	rhsSOloPN_F =  &evolution::rhs_SOloPNnp2d3_FS;
	if(jacA)
	  jacSOloPN_F = &evolution::jac_SOloPNnp2d3_FS;
      }
    else
      {
	rhsSOloPN_F =  &evolution::rhs_SOloPNnp3d3_FS;
	if(jacA)
	  jacSOloPN_F = &evolution::jac_SOloPNnp3d3_FS;
      }
      
    }
  
  if(_terms.pnSSlo)
    {
      jacSSloPN_F = &evolution::jac_numSSloPN_F;

      if(number_of_particles==2 )
	{
	  rhsSSloPN_F =  &evolution::rhs_SSloPNnp2d3_FS;
	  if(jacA)
	    jacSSloPN_F = &evolution::jac_SSloPNnp2d3_FS;
	}
      else
	{
	  rhsSSloPN_F =  &evolution::rhs_SSloPNnp3d3_FS;
	  if(jacA)
	    jacSSloPN_F = &evolution::jac_SSloPNnp3d3_FS;
	}
     
    }
  if(_terms.pnS2lo)
    {
    jacS2loPN_F = &evolution::jac_numS2loPN_F;

    if(number_of_particles==2 )
      {
	rhsS2loPN_F =  &evolution::rhs_S2loPNnp2d3_FS;
	if(jacA)
	  jacS2loPN_F = &evolution::jac_S2loPNnp2d3_FS;
      }
    else
      {
	rhsS2loPN_F =  &evolution::rhs_S2loPNnp3d3_FS;
	if(jacA)
	  jacS2loPN_F = &evolution::jac_S2loPNnp3d3_FS;
      }

    
    }
    
#endif
  
#ifdef odePNSnlo
  if(_terms.pnSOnlo)
    {
    jacSOnloPN_F = &evolution::jac_numSOnloPN_F;

    if(number_of_particles==2 )
      {
	rhsSOnloPN_F =  &evolution::rhs_SOnloPNnp2d3_FS;
	if(jacA)
	  jacSOnloPN_F = &evolution::jac_SOnloPNnp2d3_FS;
      }
    else
      {
	rhsSOnloPN_F =  &evolution::rhs_SOnloPNnp3d3_FS;
	if(jacA)
	  jacSOnloPN_F = &evolution::jac_SOnloPNnp3d3_FS;
      }

    }
#endif
  
  
  if(_chaos_test)
    
    system.function = rhsC;

  else
    
    system.function = rhs;


  system.jacobian = jac;


  

}

  






bool evolution::update(double &t)
{



  ym2=ym1;
  ym1=y;
  dtm1=dt;
  
  int status = gsl_odeiv_evolve_apply (evolve, control, step,
				       &system,
				       &time, final_time,
				       &dt, &y[0]);

  dt = time-t;
  
  for(int kk=0; kk < number_of_variables; kk++)
    dy[kk]=evolve->dydt_out[kk];

  ddy = 2*(dtm1*y-(dt+dtm1)*ym1+dt*ym2)/(dt*dtm1*(dt+dtm1));

  if (status != GSL_SUCCESS)
    {
      BOOST_LOG_SEV(lg, error) << "Step failed at time step:"<< time;

      exit(Finalize(0));
    }

  
  t=time;
  iterations++;
  
  return time < final_time; 
  
}



void evolution::close()
{

  gsl_odeiv_evolve_free (evolve);
  gsl_odeiv_control_free (control);
  gsl_odeiv_step_free (step);

}


int evolution::rhs(double t, const double * y, double * f, void * param)
{

  
  evolution* mySelf = (evolution*) pt2evolution;


  // call member
  mySelf->RHSN(t,y,f,param);

  mySelf->RHS1PN(t,y,f,param);

  mySelf->RHS2PN(t,y,f,param);

  mySelf->RHS2_5PN(t,y,f,param);

  mySelf->RHSSOloPN(t,y,f,param);

  mySelf->RHSSSloPN(t,y,f,param);

  mySelf->RHSS2loPN(t,y,f,param);

  mySelf->RHSSOnloPN(t,y,f,param);

  return (GSL_SUCCESS);

}



int evolution::jac(double t, const double y[], double * dfdy, double dfdt[], void * param)
{
  
  evolution* mySelf = (evolution*) pt2evolution;


  // call member
  mySelf->JACN(t,y,dfdy,dfdt,param);
  
  mySelf->JAC1PN(t,y,dfdy,dfdt,param);

  mySelf->JAC2PN(t,y,dfdy,dfdt,param);

  mySelf->JAC2_5PN(t,y,dfdy,dfdt,param);

  mySelf->JACSOloPN(t,y,dfdy,dfdt,param);

  mySelf->JACSSloPN(t,y,dfdy,dfdt,param);

  mySelf->JACS2loPN(t,y,dfdy,dfdt,param);

  mySelf->JACSOnloPN(t,y,dfdy,dfdt,param);
  
  return (GSL_SUCCESS);
}



int evolution::rhsC(double t, const double * y, double * f, void * param)
{

  
  evolution* mySelf = (evolution*) pt2evolution;

  // call member
  mySelf->RHSN(t,y,f,param);

  mySelf->RHS1PN(t,y,f,param);

  mySelf->RHS2PN(t,y,f,param);

  mySelf->RHS2_5PN(t,y,f,param);

  mySelf->RHSSOloPN(t,y,f,param);

  mySelf->RHSSSloPN(t,y,f,param);

  mySelf->RHSS2loPN(t,y,f,param);

  mySelf->RHSSOnloPN(t,y,f,param);

  double dfdy[mySelf->number_of_variables*mySelf->number_of_variables];

  double dfdt[mySelf->number_of_variables];

  mySelf->JACN(t,y,dfdy,dfdt,param);

  mySelf->JAC1PN(t,y,dfdy,dfdt,param);

  mySelf->JAC2PN(t,y,dfdy,dfdt,param);

  mySelf->JAC2_5PN(t,y,dfdy,dfdt,param);

  mySelf->JACSOloPN(t,y,dfdy,dfdt,param);

  mySelf->JACSSloPN(t,y,dfdy,dfdt,param);

  mySelf->JACS2loPN(t,y,dfdy,dfdt,param);

  mySelf->JACSOnloPN(t,y,dfdy,dfdt,param);
  
  for(int i=0; i < mySelf->number_of_variables; i++) 
    {

      f[mySelf->number_of_variables+i] = 0.0;
      
      for(int j=0; j < mySelf->number_of_variables; j++) 
	
	f[mySelf->number_of_variables+i] += dfdy[i*mySelf->number_of_variables+j]*y[mySelf->number_of_variables+j];

    }

  return (GSL_SUCCESS);

}


int evolution::rhsN (double t, const double y[], double f[], void *param){


  double *par = static_cast<double *>(param);
  

  for(int i=0; i < number_of_variables; i++)

    f[i] = (this->*rhsN_F)(t,y,par,i);


  return (GSL_SUCCESS); }




int evolution::rhs1PN (double t, const double y[], double f[], void *param){


  double *par = static_cast<double *>(param);
  

  for(int i=0; i < number_of_variables; i++)

    f[i] += (this->*rhs1PN_F)(t,y,par,i);


  return (GSL_SUCCESS); }



int evolution::rhs2PN (double t, const double y[], double f[], void *param){


  double *par = static_cast<double *>(param);

  
  for(int i=0; i < number_of_variables; i++)

    f[i] += (this->*rhs2PN_F)(t,y,par,i);


  return (GSL_SUCCESS); }




int evolution::rhs2_5PN (double t, const double y[], double f[], void *param){


  double *par = static_cast<double *>(param);
  
  for(int i=0; i < number_of_variables; i++)

    f[i] += (this->*rhs2_5PN_F)(t,y,par,i);


  return (GSL_SUCCESS); }


int evolution::rhsSOloPN (double t, const double y[], double f[], void *param){


  double *par = static_cast<double *>(param);
  
  for(int i=0; i < number_of_variables; i++)

    f[i] += (this->*rhsSOloPN_F)(t,y,par,i);


  return (GSL_SUCCESS); }


int evolution::rhsSSloPN (double t, const double y[], double f[], void *param){


  double *par = static_cast<double *>(param);
  
  for(int i=0; i < number_of_variables; i++)

    f[i] += (this->*rhsSSloPN_F)(t,y,par,i);


  return (GSL_SUCCESS); }


int evolution::rhsS2loPN (double t, const double y[], double f[], void *param){


  double *par = static_cast<double *>(param);
  
  for(int i=0; i < number_of_variables; i++)

    f[i] += (this->*rhsS2loPN_F)(t,y,par,i);


  return (GSL_SUCCESS); }


int evolution::rhsSOnloPN (double t, const double y[], double f[], void *param){


  double *par = static_cast<double *>(param);
  
  for(int i=0; i < number_of_variables; i++)

    f[i] += (this->*rhsSOnloPN_F)(t,y,par,i);


  return (GSL_SUCCESS); }


int evolution::jacN (double t, const double y[], double * dfdy, double dfdt[], void * param){


  double *par = static_cast<double *>(param);
  
  for(int i=0; i < number_of_variables; i++)
    {
      
      dfdt[i] = 0;

      for(int j=0; j < number_of_variables; j++)
      
	dfdy[i*number_of_variables+j] = (this->*jacN_F)(t,y,par,i,j);



    }


  return (GSL_SUCCESS); }


int evolution::jac1PN (double t, const double y[], double * dfdy, double dfdt[], void * param){


  double *par = static_cast<double *>(param);
  

  for(int i=0; i < number_of_variables; i++)
    
    for(int j=0; j < number_of_variables; j++)

      dfdy[i*number_of_variables+j] += (this->*jac1PN_F)(t,y,par,i,j);


  return (GSL_SUCCESS); }



int evolution::jac2PN (double t, const double y[], double * dfdy, double dfdt[], void * param){


  double *par = static_cast<double *>(param);
  

  for(int i=0; i < number_of_variables; i++)
    
    for(int j=0; j < number_of_variables; j++)
      
      dfdy[i*number_of_variables+j] += (this->*jac2PN_F)(t,y,par,i,j);


  return (GSL_SUCCESS); }



int evolution::jac2_5PN (double t, const double y[], double * dfdy, double dfdt[], void * param){


  double *par = static_cast<double *>(param);
  
  
  for(int i=0; i < number_of_variables; i++)
    
    for(int j=0; j < number_of_variables; j++)

      dfdy[i*number_of_variables+j] += (this->*jac2_5PN_F)(t,y,par,i,j);


  return (GSL_SUCCESS); }



int evolution::jacSOloPN (double t, const double y[], double * dfdy, double dfdt[], void * param){


  double *par = static_cast<double *>(param);
  
  
  for(int i=0; i < number_of_variables; i++)
    
    for(int j=0; j < number_of_variables; j++)

      dfdy[i*number_of_variables+j] += (this->*jacSOloPN_F)(t,y,par,i,j);


  return (GSL_SUCCESS); }


int evolution::jacSSloPN (double t, const double y[], double * dfdy, double dfdt[], void * param){


  double *par = static_cast<double *>(param);
  
  
  for(int i=0; i < number_of_variables; i++)
    
    for(int j=0; j < number_of_variables; j++)

      dfdy[i*number_of_variables+j] += (this->*jacSSloPN_F)(t,y,par,i,j);


  return (GSL_SUCCESS); }



int evolution::jacS2loPN (double t, const double y[], double * dfdy, double dfdt[], void * param){


  double *par = static_cast<double *>(param);
  
  
  for(int i=0; i < number_of_variables; i++)
    
    for(int j=0; j < number_of_variables; j++)

      dfdy[i*number_of_variables+j] += (this->*jacS2loPN_F)(t,y,par,i,j);


  return (GSL_SUCCESS); }



int evolution::jacSOnloPN (double t, const double y[], double * dfdy, double dfdt[], void * param){


  double *par = static_cast<double *>(param);
  
  
  for(int i=0; i < number_of_variables; i++)
    
    for(int j=0; j < number_of_variables; j++)

      dfdy[i*number_of_variables+j] += (this->*jacSOnloPN_F)(t,y,par,i,j);


  return (GSL_SUCCESS); }



double evolution::jac_numN_F(double t, const double y[], double par[], int i, int j)
{


  double yp1[number_of_variables], yp2[number_of_variables];  

  double ym1[number_of_variables], ym2[number_of_variables];  

  double r0, round, trunc;

  double error;

  double Linf = 0;


  for(int jj=0; jj<j; jj++)
    {    
      yp1[jj] = y[jj];
      ym1[jj] = y[jj];
      yp2[jj] = y[jj];
      ym2[jj] = y[jj];
    }
  
  yp1[j] = y[j]+0.5*dx;
  ym1[j] = y[j]-0.5*dx;
  yp2[j] = y[j]+dx;
  ym2[j] = y[j]-dx;

  for(int jj=j+1; jj<number_of_variables; jj++)
    {    
      yp1[jj] = y[jj];
      ym1[jj] = y[jj];
      yp2[jj] = y[jj];
      ym2[jj] = y[jj];
    }

  double fp1 = (this->*rhsN_F)(t,yp1,par,i);
  double fm1 = (this->*rhsN_F)(t,ym1,par,i);  
  double fp2 = (this->*rhsN_F)(t,yp2,par,i);
  double fm2 = (this->*rhsN_F)(t,ym2,par,i);  


  central_deriv (fm2, fp2,
		 fm1, fp1, 
		 y[j], dx, 
		 &r0, &round, &trunc);

  error = round+trunc;

  if (round < trunc && (round > 0 && trunc > 0)){

    double r_opt, round_opt, trunc_opt, error_opt;
	
    double h_opt = dx * pow (round / (2.0 * trunc), 1.0 / 3.0);

    yp1[j] = y[j]+0.5*h_opt;
    ym1[j] = y[j]-0.5*h_opt;
    yp2[j] = y[j]+h_opt;
    ym2[j] = y[j]-h_opt;

    fp1 = (this->*rhsN_F)(t,yp1,par,i);
    fm1 = (this->*rhsN_F)(t,ym1,par,i);  
    fp2 = (this->*rhsN_F)(t,yp2,par,i);
    fm2 = (this->*rhsN_F)(t,ym2,par,i);  
    
    central_deriv (fm2, fp2,
		   fm1, fp1, 
		   y[j], h_opt, 
		   &r_opt, &round_opt, &trunc_opt);
  
    error_opt = round_opt + trunc_opt;

    if (error_opt < error && fabs (r_opt - r0) < 4.0 * error)
      {
	r0 = r_opt;
	error = error_opt;
      }

  }

  

  return(r0);

}






double evolution::jac_num1PN_F(double t, const double y[], double par[], int i, int j)
{


  double yp1[number_of_variables], yp2[number_of_variables];  

  double ym1[number_of_variables], ym2[number_of_variables];  

  double r0, round, trunc;

  double error;

  double Linf = 0;


  for(int jj=0; jj<j; jj++)
    {    
      yp1[jj] = y[jj];
      ym1[jj] = y[jj];
      yp2[jj] = y[jj];
      ym2[jj] = y[jj];
    }
  
  yp1[j] = y[j]+0.5*dx;
  ym1[j] = y[j]-0.5*dx;
  yp2[j] = y[j]+dx;
  ym2[j] = y[j]-dx;

  for(int jj=j+1; jj<number_of_variables; jj++)
    {    
      yp1[jj] = y[jj];
      ym1[jj] = y[jj];
      yp2[jj] = y[jj];
      ym2[jj] = y[jj];
    }

  double fp1 = (this->*rhs1PN_F)(t,yp1,par,i);
  double fm1 = (this->*rhs1PN_F)(t,ym1,par,i);  
  double fp2 = (this->*rhs1PN_F)(t,yp2,par,i);
  double fm2 = (this->*rhs1PN_F)(t,ym2,par,i);  


  central_deriv (fm2, fp2,
		 fm1, fp1, 
		 y[j], dx, 
		 &r0, &round, &trunc);

  error = round+trunc;

  if (round < trunc && (round > 0 && trunc > 0)){

    double r_opt, round_opt, trunc_opt, error_opt;
	
    double h_opt = dx * pow (round / (2.0 * trunc), 1.0 / 3.0);

    yp1[j] = y[j]+0.5*h_opt;
    ym1[j] = y[j]-0.5*h_opt;
    yp2[j] = y[j]+h_opt;
    ym2[j] = y[j]-h_opt;

    fp1 = (this->*rhs1PN_F)(t,yp1,par,i);
    fm1 = (this->*rhs1PN_F)(t,ym1,par,i);  
    fp2 = (this->*rhs1PN_F)(t,yp2,par,i);
    fm2 = (this->*rhs1PN_F)(t,ym2,par,i);  
    
    central_deriv (fm2, fp2,
		   fm1, fp1, 
		   y[j], h_opt, 
		   &r_opt, &round_opt, &trunc_opt);
  
    error_opt = round_opt + trunc_opt;

    if (error_opt < error && fabs (r_opt - r0) < 4.0 * error)
      {
	r0 = r_opt;
	error = error_opt;
      }

  }

  

  return(r0);

}





double evolution::jac_num2PN_F(double t, const double y[], double par[], int i, int j)
{


  double yp1[number_of_variables], yp2[number_of_variables];  

  double ym1[number_of_variables], ym2[number_of_variables];  

  double r0, round, trunc;

  double error;

  double Linf = 0;


  for(int jj=0; jj<j; jj++)
    {    
      yp1[jj] = y[jj];
      ym1[jj] = y[jj];
      yp2[jj] = y[jj];
      ym2[jj] = y[jj];
    }
  
  yp1[j] = y[j]+0.5*dx;
  ym1[j] = y[j]-0.5*dx;
  yp2[j] = y[j]+dx;
  ym2[j] = y[j]-dx;

  for(int jj=j+1; jj<number_of_variables; jj++)
    {    
      yp1[jj] = y[jj];
      ym1[jj] = y[jj];
      yp2[jj] = y[jj];
      ym2[jj] = y[jj];
    }

  double fp1 = (this->*rhs2PN_F)(t,yp1,par,i);
  double fm1 = (this->*rhs2PN_F)(t,ym1,par,i);  
  double fp2 = (this->*rhs2PN_F)(t,yp2,par,i);
  double fm2 = (this->*rhs2PN_F)(t,ym2,par,i);  


  central_deriv (fm2, fp2,
		 fm1, fp1, 
		 y[j], dx, 
		 &r0, &round, &trunc);

  error = round+trunc;

  if (round < trunc && (round > 0 && trunc > 0)){

    double r_opt, round_opt, trunc_opt, error_opt;
	
    double h_opt = dx * pow (round / (2.0 * trunc), 1.0 / 3.0);

    yp1[j] = y[j]+0.5*h_opt;
    ym1[j] = y[j]-0.5*h_opt;
    yp2[j] = y[j]+h_opt;
    ym2[j] = y[j]-h_opt;

    fp1 = (this->*rhs2PN_F)(t,yp1,par,i);
    fm1 = (this->*rhs2PN_F)(t,ym1,par,i);  
    fp2 = (this->*rhs2PN_F)(t,yp2,par,i);
    fm2 = (this->*rhs2PN_F)(t,ym2,par,i);  
    
    central_deriv (fm2, fp2,
		   fm1, fp1, 
		   y[j], h_opt, 
		   &r_opt, &round_opt, &trunc_opt);
  
    error_opt = round_opt + trunc_opt;

    if (error_opt < error && fabs (r_opt - r0) < 4.0 * error)
      {
	r0 = r_opt;
	error = error_opt;
      }

  }

  

  return(r0);

}





double evolution::jac_num2_5PN_F(double t, const double y[], double par[], int i, int j)
{


  double yp1[number_of_variables], yp2[number_of_variables];  

  double ym1[number_of_variables], ym2[number_of_variables];  

  double r0, round, trunc;

  double error;

  double Linf = 0;


  for(int jj=0; jj<j; jj++)
    {    
      yp1[jj] = y[jj];
      ym1[jj] = y[jj];
      yp2[jj] = y[jj];
      ym2[jj] = y[jj];
    }
  
  yp1[j] = y[j]+0.5*dx;
  ym1[j] = y[j]-0.5*dx;
  yp2[j] = y[j]+dx;
  ym2[j] = y[j]-dx;

  for(int jj=j+1; jj<number_of_variables; jj++)
    {    
      yp1[jj] = y[jj];
      ym1[jj] = y[jj];
      yp2[jj] = y[jj];
      ym2[jj] = y[jj];
    }

  double fp1 = (this->*rhs2_5PN_F)(t,yp1,par,i);
  double fm1 = (this->*rhs2_5PN_F)(t,ym1,par,i);  
  double fp2 = (this->*rhs2_5PN_F)(t,yp2,par,i);
  double fm2 = (this->*rhs2_5PN_F)(t,ym2,par,i);  


  central_deriv (fm2, fp2,
		 fm1, fp1, 
		 y[j], dx, 
		 &r0, &round, &trunc);

  error = round+trunc;

  if (round < trunc && (round > 0 && trunc > 0)){

    double r_opt, round_opt, trunc_opt, error_opt;
	
    double h_opt = dx * pow (round / (2.0 * trunc), 1.0 / 3.0);

    yp1[j] = y[j]+0.5*h_opt;
    ym1[j] = y[j]-0.5*h_opt;
    yp2[j] = y[j]+h_opt;
    ym2[j] = y[j]-h_opt;

    fp1 = (this->*rhs2_5PN_F)(t,yp1,par,i);
    fm1 = (this->*rhs2_5PN_F)(t,ym1,par,i);  
    fp2 = (this->*rhs2_5PN_F)(t,yp2,par,i);
    fm2 = (this->*rhs2_5PN_F)(t,ym2,par,i);  
    
    central_deriv (fm2, fp2,
		   fm1, fp1, 
		   y[j], h_opt, 
		   &r_opt, &round_opt, &trunc_opt);
  
    error_opt = round_opt + trunc_opt;

    if (error_opt < error && fabs (r_opt - r0) < 4.0 * error)
      {
	r0 = r_opt;
	error = error_opt;
      }

  }

  

  return(r0);

}





double evolution::jac_numSOloPN_F(double t, const double y[], double par[], int i, int j)
{


  double yp1[number_of_variables], yp2[number_of_variables];  

  double ym1[number_of_variables], ym2[number_of_variables];  

  double r0, round, trunc;

  double error;

  double Linf = 0;


  for(int jj=0; jj<j; jj++)
    {    
      yp1[jj] = y[jj];
      ym1[jj] = y[jj];
      yp2[jj] = y[jj];
      ym2[jj] = y[jj];
    }
  
  yp1[j] = y[j]+0.5*dx;
  ym1[j] = y[j]-0.5*dx;
  yp2[j] = y[j]+dx;
  ym2[j] = y[j]-dx;

  for(int jj=j+1; jj<number_of_variables; jj++)
    {    
      yp1[jj] = y[jj];
      ym1[jj] = y[jj];
      yp2[jj] = y[jj];
      ym2[jj] = y[jj];
    }

  double fp1 = (this->*rhsSOloPN_F)(t,yp1,par,i);
  double fm1 = (this->*rhsSOloPN_F)(t,ym1,par,i);  
  double fp2 = (this->*rhsSOloPN_F)(t,yp2,par,i);
  double fm2 = (this->*rhsSOloPN_F)(t,ym2,par,i);  


  central_deriv (fm2, fp2,
		 fm1, fp1, 
		 y[j], dx, 
		 &r0, &round, &trunc);

  error = round+trunc;

  if (round < trunc && (round > 0 && trunc > 0)){

    double r_opt, round_opt, trunc_opt, error_opt;
	
    double h_opt = dx * pow (round / (2.0 * trunc), 1.0 / 3.0);

    yp1[j] = y[j]+0.5*h_opt;
    ym1[j] = y[j]-0.5*h_opt;
    yp2[j] = y[j]+h_opt;
    ym2[j] = y[j]-h_opt;

    fp1 = (this->*rhsSOloPN_F)(t,yp1,par,i);
    fm1 = (this->*rhsSOloPN_F)(t,ym1,par,i);  
    fp2 = (this->*rhsSOloPN_F)(t,yp2,par,i);
    fm2 = (this->*rhsSOloPN_F)(t,ym2,par,i);  
    
    central_deriv (fm2, fp2,
		   fm1, fp1, 
		   y[j], h_opt, 
		   &r_opt, &round_opt, &trunc_opt);
  
    error_opt = round_opt + trunc_opt;

    if (error_opt < error && fabs (r_opt - r0) < 4.0 * error)
      {
	r0 = r_opt;
	error = error_opt;
      }

  }

  

  return(r0);

}







double evolution::jac_numSSloPN_F(double t, const double y[], double par[], int i, int j)
{


  double yp1[number_of_variables], yp2[number_of_variables];  

  double ym1[number_of_variables], ym2[number_of_variables];  

  double r0, round, trunc;

  double error;

  double Linf = 0;


  for(int jj=0; jj<j; jj++)
    {    
      yp1[jj] = y[jj];
      ym1[jj] = y[jj];
      yp2[jj] = y[jj];
      ym2[jj] = y[jj];
    }
  
  yp1[j] = y[j]+0.5*dx;
  ym1[j] = y[j]-0.5*dx;
  yp2[j] = y[j]+dx;
  ym2[j] = y[j]-dx;

  for(int jj=j+1; jj<number_of_variables; jj++)
    {    
      yp1[jj] = y[jj];
      ym1[jj] = y[jj];
      yp2[jj] = y[jj];
      ym2[jj] = y[jj];
    }

  double fp1 = (this->*rhsSSloPN_F)(t,yp1,par,i);
  double fm1 = (this->*rhsSSloPN_F)(t,ym1,par,i);  
  double fp2 = (this->*rhsSSloPN_F)(t,yp2,par,i);
  double fm2 = (this->*rhsSSloPN_F)(t,ym2,par,i);  


  central_deriv (fm2, fp2,
		 fm1, fp1, 
		 y[j], dx, 
		 &r0, &round, &trunc);

  error = round+trunc;

  if (round < trunc && (round > 0 && trunc > 0)){

    double r_opt, round_opt, trunc_opt, error_opt;
	
    double h_opt = dx * pow (round / (2.0 * trunc), 1.0 / 3.0);

    yp1[j] = y[j]+0.5*h_opt;
    ym1[j] = y[j]-0.5*h_opt;
    yp2[j] = y[j]+h_opt;
    ym2[j] = y[j]-h_opt;

    fp1 = (this->*rhsSSloPN_F)(t,yp1,par,i);
    fm1 = (this->*rhsSSloPN_F)(t,ym1,par,i);  
    fp2 = (this->*rhsSSloPN_F)(t,yp2,par,i);
    fm2 = (this->*rhsSSloPN_F)(t,ym2,par,i);  
    
    central_deriv (fm2, fp2,
		   fm1, fp1, 
		   y[j], h_opt, 
		   &r_opt, &round_opt, &trunc_opt);
  
    error_opt = round_opt + trunc_opt;

    if (error_opt < error && fabs (r_opt - r0) < 4.0 * error)
      {
	r0 = r_opt;
	error = error_opt;
      }

  }

  

  return(r0);

}





double evolution::jac_numS2loPN_F(double t, const double y[], double par[], int i, int j)
{


  double yp1[number_of_variables], yp2[number_of_variables];  

  double ym1[number_of_variables], ym2[number_of_variables];  

  double r0, round, trunc;

  double error;

  double Linf = 0;


  for(int jj=0; jj<j; jj++)
    {    
      yp1[jj] = y[jj];
      ym1[jj] = y[jj];
      yp2[jj] = y[jj];
      ym2[jj] = y[jj];
    }
  
  yp1[j] = y[j]+0.5*dx;
  ym1[j] = y[j]-0.5*dx;
  yp2[j] = y[j]+dx;
  ym2[j] = y[j]-dx;

  for(int jj=j+1; jj<number_of_variables; jj++)
    {    
      yp1[jj] = y[jj];
      ym1[jj] = y[jj];
      yp2[jj] = y[jj];
      ym2[jj] = y[jj];
    }

  double fp1 = (this->*rhsS2loPN_F)(t,yp1,par,i);
  double fm1 = (this->*rhsS2loPN_F)(t,ym1,par,i);  
  double fp2 = (this->*rhsS2loPN_F)(t,yp2,par,i);
  double fm2 = (this->*rhsS2loPN_F)(t,ym2,par,i);  


  central_deriv (fm2, fp2,
		 fm1, fp1, 
		 y[j], dx, 
		 &r0, &round, &trunc);

  error = round+trunc;

  if (round < trunc && (round > 0 && trunc > 0)){

    double r_opt, round_opt, trunc_opt, error_opt;
	
    double h_opt = dx * pow (round / (2.0 * trunc), 1.0 / 3.0);

    yp1[j] = y[j]+0.5*h_opt;
    ym1[j] = y[j]-0.5*h_opt;
    yp2[j] = y[j]+h_opt;
    ym2[j] = y[j]-h_opt;

    fp1 = (this->*rhsS2loPN_F)(t,yp1,par,i);
    fm1 = (this->*rhsS2loPN_F)(t,ym1,par,i);  
    fp2 = (this->*rhsS2loPN_F)(t,yp2,par,i);
    fm2 = (this->*rhsS2loPN_F)(t,ym2,par,i);  
    
    central_deriv (fm2, fp2,
		   fm1, fp1, 
		   y[j], h_opt, 
		   &r_opt, &round_opt, &trunc_opt);
  
    error_opt = round_opt + trunc_opt;

    if (error_opt < error && fabs (r_opt - r0) < 4.0 * error)
      {
	r0 = r_opt;
	error = error_opt;
      }

  }

  

  return(r0);

}





double evolution::jac_numSOnloPN_F(double t, const double y[], double par[], int i, int j)
{


  double yp1[number_of_variables], yp2[number_of_variables];  

  double ym1[number_of_variables], ym2[number_of_variables];  

  double r0, round, trunc;

  double error;

  double Linf = 0;


  for(int jj=0; jj<j; jj++)
    {    
      yp1[jj] = y[jj];
      ym1[jj] = y[jj];
      yp2[jj] = y[jj];
      ym2[jj] = y[jj];
    }
  
  yp1[j] = y[j]+0.5*dx;
  ym1[j] = y[j]-0.5*dx;
  yp2[j] = y[j]+dx;
  ym2[j] = y[j]-dx;

  for(int jj=j+1; jj<number_of_variables; jj++)
    {    
      yp1[jj] = y[jj];
      ym1[jj] = y[jj];
      yp2[jj] = y[jj];
      ym2[jj] = y[jj];
    }

  double fp1 = (this->*rhsSOnloPN_F)(t,yp1,par,i);
  double fm1 = (this->*rhsSOnloPN_F)(t,ym1,par,i);  
  double fp2 = (this->*rhsSOnloPN_F)(t,yp2,par,i);
  double fm2 = (this->*rhsSOnloPN_F)(t,ym2,par,i);  


  central_deriv (fm2, fp2,
		 fm1, fp1, 
		 y[j], dx, 
		 &r0, &round, &trunc);

  error = round+trunc;

  if (round < trunc && (round > 0 && trunc > 0)){

    double r_opt, round_opt, trunc_opt, error_opt;
	
    double h_opt = dx * pow (round / (2.0 * trunc), 1.0 / 3.0);

    yp1[j] = y[j]+0.5*h_opt;
    ym1[j] = y[j]-0.5*h_opt;
    yp2[j] = y[j]+h_opt;
    ym2[j] = y[j]-h_opt;

    fp1 = (this->*rhsSOnloPN_F)(t,yp1,par,i);
    fm1 = (this->*rhsSOnloPN_F)(t,ym1,par,i);  
    fp2 = (this->*rhsSOnloPN_F)(t,yp2,par,i);
    fm2 = (this->*rhsSOnloPN_F)(t,ym2,par,i);  
    
    central_deriv (fm2, fp2,
		   fm1, fp1, 
		   y[j], h_opt, 
		   &r_opt, &round_opt, &trunc_opt);
  
    error_opt = round_opt + trunc_opt;

    if (error_opt < error && fabs (r_opt - r0) < 4.0 * error)
      {
	r0 = r_opt;
	error = error_opt;
      }

  }

  

  return(r0);

}

valarray<double> evolution::get_position(){

   position[0]=time;
   int i=1;
   for(int a=0; a<number_of_particles; a++)
     for(int axis=0; axis<space_dimension; axis++)
       {
	 position[i]=y[r_index(a,axis,space_dimension)];
	 i++;
       }
	 
   return position;
 }

valarray<double> evolution::get_momentum(){

    momentum[0]=time;
    int i=1;
    for(int a=0; a<number_of_particles; a++)
      for(int axis=0; axis<space_dimension; axis++)
	{
	  momentum[i]=y[p_index(a,axis,space_dimension)];
	 i++;
	}
    
   return momentum; }

valarray<double> evolution::get_spin(){

    spin[0]=time;
    int i=1;
    
    if(spin.size()>1)
    for(int a=0; a<number_of_particles; a++)
      for(int axis=0; axis<space_dimension; axis++)
	{
	  spin[i]=y[s_index(a,axis,space_dimension)];
	 i++;
	}
    
    
   return spin; }


valarray<double> evolution::get_dxdt(){

   dxdt[0]=time;
   int i=1;
   for(int a=0; a<number_of_particles; a++)
     for(int axis=0; axis<space_dimension; axis++)
       {
	 dxdt[i]=dy[r_index(a,axis,space_dimension)];
	 i++;
       }
	 
   return dxdt;
 }

valarray<double> evolution::get_dpdt(){

    dpdt[0]=time;
    int i=1;
    for(int a=0; a<number_of_particles; a++)
      for(int axis=0; axis<space_dimension; axis++)
	{
	  dpdt[i]=dy[p_index(a,axis,space_dimension)];
	 i++;
	}
    
   return dpdt; }

valarray<double> evolution::get_dsdt(){

    dsdt[0]=time;
    int i=1;
    
    if(dsdt.size()>1)
    for(int a=0; a<number_of_particles; a++)
      for(int axis=0; axis<space_dimension; axis++)
	{
	  dsdt[i]=dy[s_index(a,axis,space_dimension)];
	 i++;
	}
    
    
   return dsdt; }

valarray<double> evolution::get_ddxdt2(){

   ddxdt2[0]=time;
   int i=1;
   for(int a=0; a<number_of_particles; a++)
     for(int axis=0; axis<space_dimension; axis++)
       {
	 ddxdt2[i]=ddy[r_index(a,axis,space_dimension)];
	 i++;
       }
	 
   return ddxdt2;
 }

valarray<double> evolution::get_ddpdt2(){

    ddpdt2[0]=time;
    int i=1;
    for(int a=0; a<number_of_particles; a++)
      for(int axis=0; axis<space_dimension; axis++)
	{
	  ddpdt2[i]=ddy[p_index(a,axis,space_dimension)];
	 i++;
	}
    
   return ddpdt2; }

valarray<double> evolution::get_ddsdt2(){

    ddsdt2[0]=time;
    int i=1;
    
    if(ddsdt2.size()>1)
    for(int a=0; a<number_of_particles; a++)
      for(int axis=0; axis<space_dimension; axis++)
	{
	  ddsdt2[i]=ddy[s_index(a,axis,space_dimension)];
	 i++;
	}
    
    
   return ddsdt2; }


valarray<double> evolution::get_waves(){

  double Qtt[3][3];
  double Qttt[3][3][3];
  double Ctt[3][3][3];

  if(time==0)
    return waves;
  
  comp_Qtt(&y[0],&dy[0],&ddy[0],Qtt);

  comp_Qttt(&y[0],&dy[0],&ddy[0],Qttt);

  comp_Ctt(&y[0],&dy[0],&ddy[0],Ctt);

  double h_sh[38]; //decomposition of h+, hx in spherical harmonics
	
  comp_waves(Qtt,Qttt,Ctt,&h_sh[0]);

  waves[0]=time;
  for(int i=0; i < 38; i++)
    waves[i+1]=h_sh[i];

  return waves;
    
}


  void evolution::comp_waves(const double Qtt[3][3],const double Qttt[3][3][3],const double Ctt[3][3][3], double *h_sh)
  {
  //Quadrupole real part from m=-2 to m=2;
  h_sh[0] = -1.5853309190424043*(Qtt[0][0] - Qtt[1][1]);

  h_sh[1] = -0.2642218198404007*(11.*Qtt[0][2] + Qtt[2][0]);
  
  h_sh[2] = 1.2944172750371328*(Qtt[0][0] + Qtt[1][1] - 2.*Qtt[2][2]);

  h_sh[3] = -h_sh[1];

  h_sh[4] = h_sh[0];


  //Quadrupole imaginary part from m=-2 to m=2;
  h_sh[5] = -1.5853309190424043*(Qtt[0][1] + Qtt[1][0]);

  h_sh[6] = -0.2642218198404007*(11.*Qtt[1][2] + Qtt[2][1]);

  h_sh[7] = 0.;

  h_sh[8] = h_sh[6];

  h_sh[9] = -h_sh[5]; 


  //Octupole real part from m=-3 to m=3;

  h_sh[10] = 0.18233037789862086*(Qttt[0][0][0] - Qttt[0][1][1] - Qttt[1][0][1] - Qttt[1][1][0]);

  h_sh[11] = 0.007443606507674209*(20.*Qttt[0][0][2] + 41.*Qttt[0][2][0] - 20.*Qttt[1][1][2]-41.*Qttt[1][2][1] - Qttt[2][0][0] + Qttt[2][1][1]);

  h_sh[12] = 0.023538750570302115*(-6.*Qttt[0][0][0] - 9.*Qttt[0][1][1] + 15.*Qttt[0][2][2] + 5.*Qttt[1][0][1] - 2.*Qttt[1][1][0] + Qttt[2][0][2] + 8.*Qttt[2][2][0]);

  h_sh[13] = -0.16308124773781663*(Qttt[0][0][2] + Qttt[0][2][0] + Qttt[1][1][2] + Qttt[1][2][1] + Qttt[2][0][0] + Qttt[2][1][1] - 2.*Qttt[2][2][2]);

  h_sh[14] = -h_sh[12]; //0.023538750570302115*(6.*Qttt[0][0][0] + 9.*Qttt[0][1][1] - 15.*Qttt[0][2][2] - 5.*Qttt[1][0][1] + 2.*Qttt[1][1][0] - 1.*Qttt[2][0][2] - 8.*Qttt[2][2][0]);

  h_sh[15] = h_sh[11];// 0.007443606507674209*(20.*Qttt[0][0][2] + 41.*Qttt[0][2][0] - 20.*Qttt[1][1][2] - 41.*Qttt[1][2][1] - 1.*Qttt[2][0][0] + Qttt[2][1][1]);

  h_sh[16] = -h_sh[10]; //0.18233037789862086*(-Qttt[0][0][0] + Qttt[0][1][1] + Qttt[1][0][1] + Qttt[1][1][0]);


  //Octupole imaginary part from m=-3 to m=3;

  h_sh[17] = 0.18233037789862086*(Qttt[0][0][1] + Qttt[0][1][0] + Qttt[1][0][0] - Qttt[1][1][1]);

  h_sh[18] = -0.007443606507674209* (-20.*Qttt[0][1][2] - 41.*Qttt[0][2][1] - 20.*Qttt[1][0][2] - 41.*Qttt[1][2][0] + Qttt[2][0][1] + Qttt[2][1][0]);

  h_sh[19] = -0.023538750570302115*(2.*Qttt[0][0][1] - 5.*Qttt[0][1][0] + 9.*Qttt[1][0][0] + 6.*Qttt[1][1][1] - 15.*Qttt[1][2][2] - 1.*Qttt[2][1][2] - 8.*Qttt[2][2][1]);

  h_sh[20] = 0.;

  h_sh[21] = h_sh[19];//-0.023538750570302115*(2.*Qttt[0][0][1] - 5.*Qttt[0][1][0] + 9.*Qttt[1][0][0] + 6.*Qttt[1][1][1] - 15.*Qttt[1][2][2] - 1.*Qttt[2][1][2] - 8.*Qttt[2][2][1]);

  h_sh[22] = -h_sh[18];//-0.007443606507674209*(20.*Qttt[0][1][2] + 41.*Qttt[0][2][1] + 20.*Qttt[1][0][2] + 41.*Qttt[1][2][0] - Qttt[2][0][1] - Qttt[2][1][0]);

  h_sh[23] = h_sh[17];// 0.18233037789862086*(Qttt[0][0][1] + Qttt[0][1][0] + Qttt[1][0][0] - Qttt[1][1][1])


  //Current quadrupole real part from m=-3 to m=3;

  h_sh[24] =  0.08807393994680024*(-8.*Ctt[0][0][2] + 7.*Ctt[0][2][0] + 8.*Ctt[1][1][2] - 
			7.*Ctt[1][2][1] + Ctt[2][0][0] - Ctt[2][1][1]);

  h_sh[25] = 0.17614787989360048*(3.*Ctt[0][1][1] - 3.*Ctt[0][2][2] + Ctt[1][0][1] - 
			4.*Ctt[1][1][0] - Ctt[2][0][2] + 4.*Ctt[2][2][0]);
  h_sh[26] = 0.;

  h_sh[27] = 0.17614787989360048*(3.*Ctt[0][1][1] - 3.*Ctt[0][2][2] + Ctt[1][0][1] - 4.*Ctt[1][1][0] - Ctt[2][0][2] + 4.*Ctt[2][2][0]);

  h_sh[28] = 0.08807393994680024*(8.*Ctt[0][0][2] - 7.*Ctt[0][2][0] - 8.*Ctt[1][1][2] + 7.*Ctt[1][2][1] - Ctt[2][0][0] + Ctt[2][1][1]);


  h_sh[29] = -0.08807393994680024*(8.*Ctt[0][1][2] - 7.*Ctt[0][2][1] + 8.*Ctt[1][0][2] - 
				   7.*Ctt[1][2][0] - Ctt[2][0][1] - Ctt[2][1][0]);
  
  h_sh[30] = -0.17614787989360048*(4.*Ctt[0][0][1] - 1.*Ctt[0][1][0] - 3.*Ctt[1][0][0] + 
				   3.*Ctt[1][2][2] + Ctt[2][1][2] - 4.*Ctt[2][2][1]);
  
  h_sh[31] = 0.4314724250123776*(Ctt[0][1][2] - 4.*Ctt[0][2][1] - Ctt[1][0][2] + 
				 4.*Ctt[1][2][0]);

  h_sh[32] = 0.17614787989360048*(4.*Ctt[0][0][1] - Ctt[0][1][0] - 3.*Ctt[1][0][0] + 3.*Ctt[1][2][2] + Ctt[2][1][2] - 4.*Ctt[2][2][1]);
  
  h_sh[33] = -0.08807393994680024*(8.*Ctt[0][1][2] - 7.*Ctt[0][2][1] + 8.*Ctt[1][0][2] - 7.*Ctt[1][2][0] - 1.*Ctt[2][0][1] - 1.*Ctt[2][1][0]);

  h_sh[34] = 0;

  h_sh[35] = h_sh[33];
  
  h_sh[36] = -h_sh[32];

  h_sh[37] = h_sh[31];

  h_sh[38] = -h_sh[30];

  }


void evolution::comp_Qtt(const double y[], const double dydt[], const double ddydt2[], double (&Qtt)[3][3])
{

  size_t num_par = par.size();

  
  double Mtt[3][3];

  for(int i=0; i< 3;i++)
    for(int j=0; j< 3;j++)
      {
	Mtt[i][j] = 0;
 	if(i<space_dimension && j<space_dimension)
	  for(int a=0; a < num_par; a++){
	    size_t I = r_index(a,i,space_dimension);
	    size_t J = r_index(a,j,space_dimension);
	    
	    Mtt[i][j] += par[a]*(ddydt2[I]*y[J]+2*dydt[I]*dydt[J]+ddydt2[J]*y[I]);
	  }

      }


  double Mtt_kk=0;
  for(int k=0; k< 3; k++)
    Mtt_kk += Mtt[k][k];
  
  for(int i=0; i< 3; i++)
    for(int j=0; j< 3; j++)
      Qtt[i][j] = i==j ? Mtt[i][j]-Mtt_kk/3.0 : Mtt[i][j];
  
}


void evolution::comp_Qttt(const double y[], const double dydt[], const double ddydt2[], double (&Qttt)[3][3][3])
{
  
  size_t num_par = par.size();

  double Mttt[3][3][3];


  for(int i=0; i< 3;i++)
    for(int j=0; j< 3;j++)
      for(int k=0; k< 3;k++)
	{
	  Mttt[i][j][k] = 0;

	  if(i<space_dimension && j<space_dimension && k<space_dimension)
	    for(int a=0; a < num_par; a++){
	      size_t I = r_index(a,i,space_dimension);
	      size_t J = r_index(a,j,space_dimension);
	      size_t K = r_index(a,k,space_dimension);
	      
	      size_t Ip = p_index(a,i,space_dimension);
	      size_t Jp = p_index(a,j,space_dimension);
	      size_t Kp = p_index(a,k,space_dimension);

	      Mttt[i][j][k] += ddydt2[Ip]*y[J]*y[K] + y[I]*ddydt2[Jp]*y[K] + y[I]*y[J]*ddydt2[Kp] + 3*par[a]*(ddydt2[I]*dydt[J]*y[K] + ddydt2[I]*dydt[K]*y[J] +ddydt2[J]*dydt[I]*y[K] + ddydt2[K]*dydt[I]*y[J] + ddydt2[K]*dydt[J]*y[I] + ddydt2[J]*dydt[K]*y[I] + 2*dydt[I]*dydt[J]*dydt[K]);


	    }

	}



  double Mttt_llk[3]={0,0,0};
  double Mttt_ljl[3]={0,0,0};
  double Mttt_ill[3]={0,0,0};

  for(int ijk=0; ijk< 3; ijk++)
    for(int l=0; l< 3; l++){
      Mttt_llk[ijk] += Mttt[l][l][ijk];
      Mttt_ljl[ijk] += Mttt[l][ijk][l];
      Mttt_ill[ijk] += Mttt[ijk][l][l];
    }


  for(int i=0; i< 3; i++)
    for(int j=0; j< 3; j++)
      for(int k=0; k< 3; k++){

	Qttt[i][j][k] = Mttt[i][j][k];


	if(i==j)
	  Qttt[i][j][k] -= Mttt_llk[k]/5.0;

	if(i==k)
	  Qttt[i][j][k] -= Mttt_ljl[j]/5.0;

	if(j==k)
	  Qttt[i][j][k] -= Mttt_ill[i]/5.0;


      }

}

void evolution::comp_Ctt(const double y[], const double dydt[], const double ddydt2[], double (&Ctt)[3][3][3])
{

  size_t num_par = par.size();

  double Ptt[3][3][3];


  for(int i=0; i< 3;i++)
    for(int j=0; j< 3;j++)
      for(int k=0; k< 3;k++)
	{
	  Ptt[i][j][k] = 0;

	  if(i<space_dimension && j<space_dimension && k<space_dimension)
	    for(int a=0; a < num_par; a++){
	      size_t I = r_index(a,i,space_dimension);
	      size_t J = r_index(a,j,space_dimension);
	      size_t K = r_index(a,k,space_dimension);
	      
	      size_t Ip = p_index(a,i,space_dimension);
	      size_t Jp = p_index(a,j,space_dimension);
	      size_t Kp = p_index(a,k,space_dimension);

	      Ptt[i][j][k] += ddydt2[Ip]*y[J]*y[K]+ddydt2[K]*y[Ip]*y[J] +  ddydt2[J]*y[K]*y[Ip]+ 2*(dydt[Ip]*dydt[J]*y[K] + dydt[K]*dydt[Ip]*y[J] + dydt[J]*dydt[K]*y[Ip]);


	    }

	}



  for(int i=0; i< 3; i++)
    for(int j=0; j< 3; j++)
      for(int k=0; k< 3; k++)

	Ctt[i][j][k] = Ptt[i][j][k]+Ptt[j][i][k]-2*Ptt[k][i][j];

     
}
