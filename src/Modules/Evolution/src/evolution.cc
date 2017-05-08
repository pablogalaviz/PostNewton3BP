

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


evolution::evolution(po::variables_map &vm){

  pt2evolution =(void*) this;
  
  verbose = vm["evolution.verbose"].as<bool>();


  
  if( vm["evolution.ode_method"].as<string>().compare("rk2") == 0 ) 
    step_type = gsl_odeiv_step_rk2;
  else if( vm["evolution.ode_method"].as<string>().compare("rk4") == 0 ) 
    step_type = gsl_odeiv_step_rk4;
  else if( vm["evolution.ode_method"].as<string>().compare("rkck") == 0 ) 
    step_type = gsl_odeiv_step_rkck;
  else if( vm["evolution.ode_method"].as<string>().compare("rk8pd") == 0 ) 
    step_type = gsl_odeiv_step_rk8pd;
  else if( vm["evolution.ode_method"].as<string>().compare("rk2imp") == 0 ) 
    step_type = gsl_odeiv_step_rk2imp;
  else if( vm["evolution.ode_method"].as<string>().compare("rk4imp") == 0 ) 
    step_type = gsl_odeiv_step_rk4imp;
  else if( vm["evolution.ode_method"].as<string>().compare("bsimp") == 0 ) 
    step_type = gsl_odeiv_step_bsimp;
  else if( vm["evolution.ode_method"].as<string>().compare("gear1") == 0 ) 
    step_type = gsl_odeiv_step_gear1;
  else if( vm["evolution.ode_method"].as<string>().compare("gear2") == 0 ) 
    step_type = gsl_odeiv_step_gear2;
  else
    step_type = gsl_odeiv_step_rkf45;

  time =0; 
  dt = vm["evolution.initial_dt"].as<double>();
  final_time = vm["evolution.final_time"].as<double>();
  dx = vm["evolution.jacobian_dx"].as<double>();

  number_of_particles  = vm["evolution.particles"].as<size_t>();

  
  if (number_of_particles < 2 || number_of_particles > 3)
    {
      
      BOOST_LOG_SEV(lg, error) << "evolution.particles  should be  2 or 3.";
      BOOST_LOG_SEV(lg, error) << "Current value : "<< number_of_particles;

      exit(Finalize(0));
    }

  bool plane_constrain =  vm["evolution.plane_constrain"].as<bool>();


  bool pnSOlo =false;
  bool pnSSlo=false;
  bool pnS2lo=false;
#ifdef odePNSlo
  pnSOlo =vm["terms.pnSOlo"].as<bool>();
  pnSSlo=vm["terms.pnSSlo"].as<bool>();
  pnS2lo=vm["terms.pnS2lo"].as<bool>();
#endif
  bool pnSOnlo=false;
  bool pnSSnlo=false;
#ifdef odePNSnlo
  pnSOnlo =vm["terms.pnSOnlo"].as<bool>();
  pnSSnlo=vm["terms.pnSSnlo"].as<bool>();
#endif

  bool spin =  pnSOlo||pnSSlo||pnS2lo||pnSOnlo||pnSSnlo;
    
  if(number_of_particles==2 && ! plane_constrain  && !spin)
    plane_constrain=true;
  else
    if(plane_constrain && spin)
      plane_constrain=false;

  number_of_variables = plane_constrain ?  number_of_particles*4 : number_of_particles*6;

  if(spin) 
    number_of_variables += plane_constrain ?  number_of_particles*2 : number_of_particles*3;

  ode_size = number_of_variables ;

  NDvar = number_of_variables; 

  
  bool chaos_test = vm["evolution.chaos_test"].as<bool>() ;
 
  if(chaos_test) {

    ode_size *= 2;

    double scale_abs [ode_size];
    
    double fac =  vm["evolution.factor_chaos_test"].as<double>();
    

    for(int i = 0; i < number_of_variables; i++)
      scale_abs[i] = 1.0;

    for(int i = number_of_variables; i < ode_size; i++)
      scale_abs[i] = fac;

    control = gsl_odeiv_control_scaled_new (
					    vm["evolution.epsilon_abs"].as<double>(),
					    vm["evolution.epsilon_rel"].as<double>(),
					    vm["evolution.scaling_variable"].as<double>(),
					    vm["evolution.scaling_derivative"].as<double>(),
					    scale_abs,
					    ode_size);

  }
  else
    control = gsl_odeiv_control_standard_new (
					      vm["evolution.epsilon_abs"].as<double>(),
					      vm["evolution.epsilon_rel"].as<double>(),
					      vm["evolution.scaling_variable"].as<double>(),
					      vm["evolution.scaling_derivative"].as<double>());



  step = gsl_odeiv_step_alloc (step_type, ode_size);

  evolve = gsl_odeiv_evolve_alloc (ode_size);


  if(verbose)
    {
      BOOST_LOG_SEV(lg, info)  << "Evolution info:";
      BOOST_LOG_SEV(lg, info)  << "======================================";
      if(plane_constrain)
	BOOST_LOG_SEV(lg,info)  << "Evolution in 2D";
      if(chaos_test)
	BOOST_LOG_SEV(lg,info)  << "Testing for chaos";
      BOOST_LOG_SEV(lg, info)  << "Method: "  << gsl_odeiv_step_name(step);
      BOOST_LOG_SEV(lg, info)  << "Control: "  << gsl_odeiv_control_name(control);
      BOOST_LOG_SEV(lg, info)  << "Number of particles: "  << number_of_particles;
      BOOST_LOG_SEV(lg, info)  << "Final time: "  << final_time;
      BOOST_LOG_SEV(lg, info)  << "======================================";
    }

  y = new double [ode_size];
  
  dy = new double [ode_size];
  
  par = new double [number_of_particles];
      
  set_rhs(vm);


  
};

void evolution::set_rhs(po::variables_map &vm)
{

  size_t space_dimension = 3; 

  
  if(number_of_particles == 2 || vm["evolution.plane_constrain"].as<bool>() )
    space_dimension = 2; 

  
  system.dimension = ode_size;

  system.params = &par[0];

  jacN_F = &evolution::jac_numN_F;


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
  

#ifdef odePN1
  if(vm["terms.pn1"].as<bool>())
    jac1PN_F = &evolution::jac_num1PN_F;
#endif

#ifdef odePN2
  if(vm["terms.pn2"].as<bool>())
    jac2PN_F = &evolution::jac_num2PN_F;
#endif
  
#ifdef odePN2_5
  if(vm["terms.pn2_5"].as<bool>())
    jac2_5PN_F = &evolution::jac_num2_5PN_F;
#endif

#ifdef odePNSlo
  if(vm["terms.pnSOlo"].as<bool>())
    jacSOloPN_F = &evolution::jac_numSOloPN_F;

  
  if(vm["terms.pnSSlo"].as<bool>())
    jacSSloPN_F = &evolution::jac_numSSloPN_F;


  if(vm["terms.pnS2lo"].as<bool>())
    jacS2loPN_F = &evolution::jac_numS2loPN_F;

#endif
  
#ifdef odePNSnlo
  if(vm["terms.pnSOnlo"].as<bool>())
    jacSOnloPN_F = &evolution::jac_numSOnloPN_F;
#endif
  

switch ( number_of_particles ) {

  case 2 : 

    switch ( space_dimension ) {

    case 2 : 

      rhsN_F = &evolution::rhs_Nnp2d2_FNS;

      if(vm["evolution.JacA"].as<bool>())
	jacN_F = &evolution::jac_Nnp2d2_FaNS;
  

#ifdef odePN1
      if(vm["terms.pn1"].as<bool>()){
	rhs1PN_F = &evolution::rhs_1PNnp2d2_FNS;
	if(JacA)
	  jac1PN_F = &evolution::jac_1PNnp2d2_FaNS;
      }
#endif
    
#ifdef odePN2
      if(vm["terms.pn2"].as<bool>())
	rhs2PN_F = &evolution::rhs_2PNnp2d2_FNS;
#endif
	  

#ifdef odePN2_5
      if(vm["terms.pn2_5"].as<bool>())
	rhs2_5PN_F = &evolution::rhs_2_5PNnp2d2_FNS;
#endif
	  

      break;

    case 3 : 

      rhsN_F = &evolution::rhs_Nnp2d3_FS;

      if(vm["evolution.JacA"].as<bool>())
	jacN_F = &evolution::jac_Nnp2d3_FaS;

#ifdef odePN1
      if(vm["terms.pn1"].as<bool>()){
	rhs1PN_F = &evolution::rhs_1PNnp2d3_FS;
	if(JacA)
	  jac1PN_F = &evolution::jac_1PNnp2d3_FaS;
      }
#endif
    
      
#ifdef odePN2
      if(vm["terms.pn2"].as<bool>())
	rhs2PN_F = &evolution::rhs_2PNnp2d3_FS;
#endif

#ifdef odePN2_5
      if(vm["terms.pn2_5"].as<bool>())
	rhs2_5PN_F = &evolution::rhs_2_5PNnp2d3_FS;
#endif
      
#ifdef odePNSlo
      if(vm["terms.pnSOlo"].as<bool>()){
	rhsSOloPN_F = &evolution::rhs_SOloPNnp2d3_FS;
	if(JacA)
	  jacSOloPN_F = &evolution::jac_SOloPNnp2d3_FaS;
      }

      if(vm["terms.pnSSlo"].as<bool>()){
	rhsSSloPN_F = &evolution::rhs_SSloPNnp2d3_FS;
      	if(JacA)
	  jacSSloPN_F = &evolution::jac_SSloPNnp2d3_FaS;
      }
      if(vm["terms.pnS2lo"].as<bool>()){
	rhsS2loPN_F = &evolution::rhs_S2loPNnp2d3_FS;
	if(JacA)
	  jacS2loPN_F = &evolution::jac_S2loPNnp2d3_FaS;
      }
#endif

      
#ifdef odePNSnlo
      if(vm["terms.pnSOnlo"].as<bool>())
	rhsSOnloPN_F = &evolution::rhs_SOnloPNnp2d3_FS;
#endif
    

      break;
      
    default : 
      {
	BOOST_LOG_SEV(lg, error) << "Number of dimensions must be 2 or 3: "<< space_dimension ;

	exit(Finalize(0));
      }
      
    }

    break;

  case 3 : 


    switch ( space_dimension ) {

    case 2 : 
        

      if(vm["evolution.JacA"].as<bool>())
	jacN_F = &evolution::jac_Nnp3d2_FaNS;
  
#ifdef odePN1
      if(vm["terms.pn1"].as<bool>()){
	rhs1PN_F = &evolution::rhs_1PNnp3d2_FNS;
	if(vm["terms.JacA"].as<bool>())
	  jac1PN_F = &evolution::jac_1PNnp3d2_FaNS;
      } 
#endif

#ifdef odePN2
      if(vm["terms.pn2"].as<bool>())
	rhs2PN_F = &evolution::rhs_2PNnp3d2_FNS;
#endif
	  
#ifdef odePN2_5
      if(vm["terms.pn2_5"].as<bool>())
	rhs2_5PN_F = &evolution::rhs_2_5PNnp3d2_FNS;
#endif
	  
    
      break;

    case 3 : 
    

      if(vm["terms.pnSOlo"].as<bool>() ||
	 vm["terms.pnSSlo"].as<bool>() ||
	 vm["terms.pnS2lo"].as<bool>() ||
	 vm["terms.pnSOnlo"].as<bool>()||
	 vm["terms.pnSSnlo"].as<bool>()) {


	rhsN_F = &evolution::rhs_Nnp3d3_FS;
	
	if(vm["terms.JacA"].as<bool>())
	  jacN_F = &evolution::jac_Nnp3d3_FaS;

#ifdef odePN2
	if(vm["terms.pn1"].as<bool>()){
	  rhs1PN_F = &evolution::rhs_1PNnp3d3_FS;
	  if(vm["evolution.JacA"].as<bool>())
	    jac1PN_F = &evolution::jac_1PNnp3d3_FaS;
	}
#endif
    
      
#ifdef odePN2
	if(vm["terms.pn2"].as<bool>())
	  rhs2PN_F = &evolution::rhs_2PNnp3d3_FS;
#endif
	
#ifdef odePN2_5
	/*

	  if(vm["terms.pn2_5"].as<bool>())
	  rhs2_5PN_F = &evolution::rhs_2_5PNnp3d3_F;
	*/
#endif
	
      }
      else {
	

	rhsN_F = &evolution::rhs_Nnp3d3_FNS;
	
	if(vm["terms.JacA"].as<bool>())
	  jacN_F = &evolution::jac_Nnp3d3_FaNS;

#ifdef odePN1
	if(vm["terms.pn1"].as<bool>()){
	  rhs1PN_F = &evolution::rhs_1PNnp3d3_FNS;
	  
	  if(vm["terms.JacA"].as<bool>())
	    jac1PN_F = &evolution::jac_1PNnp3d3_FaNS;
	  
	}
#endif

	
#ifdef odePN2
	if(vm["terms.pn2"].as<bool>())
	  rhs2PN_F = &evolution::rhs_2PNnp3d3_FNS;
#endif
	
#ifdef odePN2_5
	if(vm["terms.pn2_5"].as<bool>())
	  rhs2_5PN_F = &evolution::rhs_2_5PNnp3d3_FNS;
#endif

	
#ifdef odePNSlo      
	if(vm["terms.pnSOlo"].as<bool>()){

	  rhsSOloPN_F = &evolution::rhs_SOloPNnp3d3_FS;
	  if(vm["evolution.JacA"].as<bool>())
	    jacSOloPN_F = &evolution::jac_SOloPNnp3d3_FaS;
	}

	if(vm["terms.pnSSlo"].as<bool>()){
	  rhsSSloPN_F = &evolution::rhs_SSloPNnp3d3_FS;
	  if(vm["evolution.JacA"].as<bool>())
	    jacSSloPN_F = &evolution::jac_SSloPNnp3d3_FaS;
	}
	if(vm["terms.pnS2lo"].as<bool>()){
	  rhsS2loPN_F = &evolution::rhs_S2loPNnp3d3_FS;
	  if(vm["evolution.JacA"].as<bool>())
	    jacS2loPN_F = &evolution::jac_S2loPNnp3d3_FaS;
	}
#endif
	     

#ifdef odePNSnlo
	/*
	  if(vm["terms.pnSOnlo)
	  rhsSOnloPN_F = &evolution::rhs_SOnloPNnp3d3_FS;
	*/
#endif


	break;
    
	default : 
	  {    
	    BOOST_LOG_SEV(lg, error) << "Number of dimensions must be 2 or 3: "<< space_dimension ;

	    exit(Finalize(0));
	  }
      }

    }

    break;
    
  default : 
    {
      BOOST_LOG_SEV(lg, error) << "evolution.particles  should be  2 or 3.";
      BOOST_LOG_SEV(lg, error) << "Current value : "<< number_of_particles;
      exit(Finalize(0));
      
    }
  }


  
  if(vm["evolution.chaos_test"].as<bool>())
    
    system.function = rhsC;

  else
    
    system.function = rhs;


  system.jacobian = jac;


  

}

  



void evolution::init(valarray<double>  &_y, valarray<double>  &_par){

  if( _y.size() != ode_size || _dy.size() != ode_size || _par.size() != number_of_particles )
    {

      BOOST_LOG_SEV(lg, error) << "Wrong initialization! On variables size";

      exit(Finalize(0));

    }

  for(int i =0; i < ode_size; i++)
    {
      y[i]=_y[i];
      dy[i] = _dy[i];
    }
  

  
}



bool evolution::update(double &t)
{


  int status = gsl_odeiv_evolve_apply (evolve, control, step,
				       &system,
				       &time, final_time,
				       &dt, y);

  if (status != GSL_SUCCESS)
    {
      BOOST_LOG_SEV(lg, error) << "Step failed at time step:"<< time;

      exit(Finalize(0));
    }

  t=time; 

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
