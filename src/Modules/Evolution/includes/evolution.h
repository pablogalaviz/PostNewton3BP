

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


  const gsl_odeiv_step_type * step_type;

  gsl_odeiv_step * step;

  gsl_odeiv_control * control;

  gsl_odeiv_evolve * evolve;

  const gsl_rng_type *rng_T;

  gsl_rng *rng_r;

  double *par;

  gsl_odeiv_system system;

  double time; 

  double final_time; 

  double dt; 
  
 public:

 evolution(po::variables_map &vm);

 ~evolution(){};

 bool update(double &t);
 
 

};


#endif
