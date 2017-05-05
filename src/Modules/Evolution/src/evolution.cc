

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


evolution::evolution(po::variables_map &vm){

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

  dt = vm["evolution.initial_dt"].as<double>();
  final_time = vm["evolution.final_time"].as<double>();


};



bool evolution::update(double &t)
{

  time += dt; 


  t=time; 

  return time < final_time; 
  
}
