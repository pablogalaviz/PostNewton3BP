

// initial_data.h
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



#ifndef INITIAL_DATA_H

#define INITIAL_DATA_H

#include<utils.h>


class initialData
{

  src::severity_logger< severity_level > lg;

  size_t simulation_size; 
  size_t number_of_particles; 
  
  vector<valarray<double>> y; 
  vector<valarray<double>> par; 

  void load_file(string filename);

  bool verbose; 
  
 public:

 initialData(po::variables_map &vm);

 ~initialData(){};


 inline int sim_size(){ return simulation_size;}

 inline valarray<double> get_y(int i){return y[i]; }
 
 inline valarray<double> get_par(int i){return par[i]; }

};


#endif
