

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

  number_of_particles  = vm["evolution.particles"].as<size_t>();
  
  if (Nparticles < 2 || number_of_particles > 3)
    {
      
      BOOST_LOG_SEV(lg, error) << "evolution.particles  should be  2 or 3.";
      BOOST_LOG_SEV(lg, error) << "Current value : "<< number_of_particles;

      exit(Finalize(0));
    }

  bool plane_constrain =  vm["evolution.plane_constrain"].as<bool>();
  
  if(plane_constrain || number_of_particles == 2) 
    number_of_variables = number_of_particles*4;
  else
    number_of_variables = 18;

  ode_size=number_of_variables*2;

  bool chaos_test = vm["evolution.chaos_test"].as<bool>() ;
  
  if(chaos_test) {
    
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

  system.dimension = NTvar;

  system.params = &par[0];

  jacN_F = &tbh_pn::jac_numN_F;


  rhs1PN_F = &tbh_pn::rhs_zero;
  jac1PN_F = &tbh_pn::jac_zero;

  rhs2PN_F = &tbh_pn::rhs_zero;
  jac2PN_F = &tbh_pn::jac_zero;

  rhs2_5PN_F = &tbh_pn::rhs_zero;
  jac2_5PN_F = &tbh_pn::jac_zero;

  rhsSOloPN_F = &tbh_pn::rhs_zero;
  jacSOloPN_F = &tbh_pn::jac_zero;

  rhsSSloPN_F = &tbh_pn::rhs_zero;
  jacSSloPN_F = &tbh_pn::jac_zero;
 

  rhsS2loPN_F = &tbh_pn::rhs_zero;
  jacS2loPN_F = &tbh_pn::jac_zero;
 
  rhsSOnloPN_F = &tbh_pn::rhs_zero;
  jacSOnloPN_F = &tbh_pn::jac_zero;
  

#ifdef odePN1
  if(pn1)
    jac1PN_F = &tbh_pn::jac_num1PN_F;
#endif

#ifdef odePN2
  if(pn2)
    jac2PN_F = &tbh_pn::jac_num2PN_F;
#endif
  
#ifdef odePN2_5
  if(pn2_5)
    jac2_5PN_F = &tbh_pn::jac_num2_5PN_F;
#endif

#ifdef odePNSlo
  if(pnSOlo)
    jacSOloPN_F = &tbh_pn::jac_numSOloPN_F;

  
  if(pnSSlo)
    jacSSloPN_F = &tbh_pn::jac_numSSloPN_F;


  if(pnS2lo)
    jacS2loPN_F = &tbh_pn::jac_numS2loPN_F;

#endif
  
#ifdef odePNSnlo
  if(pnSOnlo)
    jacSOnloPN_F = &tbh_pn::jac_numSOnloPN_F;
#endif
  

  HLOSO = HZERO;

  HLOSS = HZERO;

  HLOS2 = HZERO;

  HNLOSO = HZERO;

  HNLOSS = HZERO;


  switch ( num_par ) {

  case 2 : 

    switch ( space_dim ) {

    case 2 : 
   

      HN = HN_NP2D2NS;

      H1PN = H1PN_NP2D2NS;

      H2PN = H2PN_NP2D2NS;

      Eesc = Eesc_NP2D2NS;


      rhsN_F = &tbh_pn::rhs_Nnp2d2_FNS;

      if(JacA)
	jacN_F = &tbh_pn::jac_Nnp2d2_FaNS;
  

#ifdef odePN1
      if(pn1){
	rhs1PN_F = &tbh_pn::rhs_1PNnp2d2_FNS;
	if(JacA)
	  jac1PN_F = &tbh_pn::jac_1PNnp2d2_FaNS;
      }
#endif
    
#ifdef odePN2
      if(pn2)
	rhs2PN_F = &tbh_pn::rhs_2PNnp2d2_FNS;
#endif
	  

#ifdef odePN2_5
      if(pn2_5)
	rhs2_5PN_F = &tbh_pn::rhs_2_5PNnp2d2_FNS;
#endif
	  

      break;

    case 3 : 


      HN = HN_NP2D3S;

      H1PN = H1PN_NP2D3S;

      H2PN = H2PN_NP2D3S;

      HLOSO = HLOSO_NP2D3;

      HLOSS = HLOSS_NP2D3;

      HLOS2 = HLOS2_NP2D3;

      HNLOSO = HNLOSO_NP2D3;

      HNLOSS = HNLOSS_NP2D3;

      Eesc = Eesc_NP2D3S;


      rhsN_F = &tbh_pn::rhs_Nnp2d3_FS;

      if(JacA)
	jacN_F = &tbh_pn::jac_Nnp2d3_FaS;

#ifdef odePN1
      if(pn1){
	rhs1PN_F = &tbh_pn::rhs_1PNnp2d3_FS;
	if(JacA)
	  jac1PN_F = &tbh_pn::jac_1PNnp2d3_FaS;
      }
#endif
    
      
#ifdef odePN2
      if(pn2)
	rhs2PN_F = &tbh_pn::rhs_2PNnp2d3_FS;
#endif

#ifdef odePN2_5
      if(pn2_5)
	rhs2_5PN_F = &tbh_pn::rhs_2_5PNnp2d3_FS;
#endif
      
#ifdef odePNSlo
      if(pnSOlo){
	rhsSOloPN_F = &tbh_pn::rhs_SOloPNnp2d3_FS;
	if(JacA)
	  jacSOloPN_F = &tbh_pn::jac_SOloPNnp2d3_FaS;
      }

      if(pnSSlo){
	rhsSSloPN_F = &tbh_pn::rhs_SSloPNnp2d3_FS;
      	if(JacA)
	  jacSSloPN_F = &tbh_pn::jac_SSloPNnp2d3_FaS;
      }
      if(pnS2lo){
	rhsS2loPN_F = &tbh_pn::rhs_S2loPNnp2d3_FS;
	if(JacA)
	  jacS2loPN_F = &tbh_pn::jac_S2loPNnp2d3_FaS;
      }
#endif

      
#ifdef odePNSnlo
      if(pnSOnlo)
	rhsSOnloPN_F = &tbh_pn::rhs_SOnloPNnp2d3_FS;
#endif
    

      break;
      
    default : 
      
      error_exit("The number of dimension must be 2 or 3.");
      
    }

    break;

  case 3 : 


    switch ( space_dim ) {

    case 2 : 
        

      HN = HN_NP3D2NS;

      H1PN = H1PN_NP3D2NS;

      H2PN = H2PN_NP3D2NS;

      Eesc = Eesc_NP3D2NS;


      rhsN_F = &tbh_pn::rhs_Nnp3d2_FNS;

      if(JacA)
	jacN_F = &tbh_pn::jac_Nnp3d2_FaNS;
  
#ifdef odePN1
      if(pn1){
	rhs1PN_F = &tbh_pn::rhs_1PNnp3d2_FNS;
	if(JacA)
	  jac1PN_F = &tbh_pn::jac_1PNnp3d2_FaNS;
      } 
#endif

#ifdef odePN2
      if(pn2)
	rhs2PN_F = &tbh_pn::rhs_2PNnp3d2_FNS;
#endif
	  
#ifdef odePN2_5
      if(pn2_5)
	rhs2_5PN_F = &tbh_pn::rhs_2_5PNnp3d2_FNS;
#endif
	  

    
      break;

    case 3 : 
    

      if(pnSOlo||pnSSlo||pnS2lo||pnSOnlo||pnSSnlo) {

	HN = HN_NP3D3S;

	H1PN = H1PN_NP3D3S;

	H2PN = H2PN_NP3D3S;
	
	//HLOSO = HLOSO_NP3D3;

	//HLOSS = HLOSS_NP3D3;

	//HLOS2 = HLOS2_NP3D3;

	//HNLOSO = HNLOSO_NP3D3;

	//HNLOSS = HNLOSS_NP3D3;

	Eesc = Eesc_NP3D3S;

	rhsN_F = &tbh_pn::rhs_Nnp3d3_FS;
	
	if(JacA)
	  jacN_F = &tbh_pn::jac_Nnp3d3_FaS;

#ifdef odePN2
	if(pn1){
	  rhs1PN_F = &tbh_pn::rhs_1PNnp3d3_FS;
	  if(JacA)
	    jac1PN_F = &tbh_pn::jac_1PNnp3d3_FaS;
	}
#endif
    
      
#ifdef odePN2
	if(pn2)
	  rhs2PN_F = &tbh_pn::rhs_2PNnp3d3_FS;
#endif
	
#ifdef odePN2_5
	/*

	  if(pn2_5)
	  rhs2_5PN_F = &tbh_pn::rhs_2_5PNnp3d3_F;
	*/
#endif
	
      }
      else {
	
	HN = HN_NP3D3NS;

	H1PN = H1PN_NP3D3NS;

	H2PN = H2PN_NP3D3NS;
	
	Eesc = Eesc_NP3D3NS;


	rhsN_F = &tbh_pn::rhs_Nnp3d3_FNS;
	
	if(JacA)
	  jacN_F = &tbh_pn::jac_Nnp3d3_FaNS;

#ifdef odePN1
	if(pn1){
	  rhs1PN_F = &tbh_pn::rhs_1PNnp3d3_FNS;
	  
	  if(JacA)
	    jac1PN_F = &tbh_pn::jac_1PNnp3d3_FaNS;
	  
	}
#endif

	
#ifdef odePN2
	if(pn2)
	  rhs2PN_F = &tbh_pn::rhs_2PNnp3d3_FNS;
#endif
	
#ifdef odePN2_5
	if(pn2_5)
	  rhs2_5PN_F = &tbh_pn::rhs_2_5PNnp3d3_FNS;
#endif

	
#ifdef odePNSlo      
	      if(pnSOlo){

	      rhsSOloPN_F = &tbh_pn::rhs_SOloPNnp3d3_FS;
	      if(JacA)
	      jacSOloPN_F = &tbh_pn::jac_SOloPNnp3d3_FaS;
	      }

	      if(pnSSlo){
	      rhsSSloPN_F = &tbh_pn::rhs_SSloPNnp3d3_FS;
	      if(JacA)
	      jacSSloPN_F = &tbh_pn::jac_SSloPNnp3d3_FaS;
	      }
	      if(pnS2lo){
	      rhsS2loPN_F = &tbh_pn::rhs_S2loPNnp3d3_FS;
	      if(JacA)
	      jacS2loPN_F = &tbh_pn::jac_S2loPNnp3d3_FaS;
	      }
#endif
	     

#ifdef odePNSnlo
	      /*
	      if(pnSOnlo)
	      rhsSOnloPN_F = &tbh_pn::rhs_SOnloPNnp3d3_FS;
      */
#endif


      break;
    
    default : 
    
      error_exit("The number of dimension must be 2 or 3.");

      }

    }

    break;
    
  default : 
    
    error_exit("The number of particles must be 2 or 3.");

  }

  if(vm["evolution.chaos_test"])
    
    system.function = rhsC;

  else
    
    system.function = rhs;


  system.jacobian = jac;


  

}

  



void evolution::Init(valarray<double> _y, valarray<double> _dy, valarray<double> _par){

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



void close()
{

  gsl_odeiv_evolve_free (evolve);
  gsl_odeiv_control_free (control);
  gsl_odeiv_step_free (step);

}


int evolution::rhs(double t, const double * y, double * f, void * param)
{
  
  evolution* mySelf = (evolution*) pt2evolution;


  // call member
  mySelf->RHSN(time,y,f,param);

  return (GSL_SUCCESS);

}
