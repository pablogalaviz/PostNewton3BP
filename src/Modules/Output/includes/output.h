

// output.h
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



#ifndef OUTPUT_H

#define OUTPUT_H

#include<utils.h>
#include<hdf5.h>
#include<evolution.h>
#include<initial_data.h> 

class output
{

  src::severity_logger< logging::trivial::severity_level > lg;

  clock_t tStart;

  bool verbose;
  
  bool debug;
  
  double delta_time;
  
  string file_name;

  string monitor_file_name;
  
  hid_t file_id;

  int iteration;
  int iteration_set;

  double next_output;

  void save(double t,valarray<double> pos, valarray<double> mom,valarray<double> spin);

  
 public:

 output(string output_directory,bool _verbose, bool _debug, double _delta_time);

 ~output(){};

 
 void init(initialData &id, size_t index);
  
  void update(double t, int index,evolution &evo,  bool force=false);

 

};


#endif
