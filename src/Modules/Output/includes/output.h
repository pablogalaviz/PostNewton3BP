

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
#include<analysis.h>
#include<initial_data.h> 

#define NDIMS 2

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
  vector<size_t> iteration_set;

  double next_output;

  void save(valarray<double> pos,string group_name, size_t index);

  void create_dataset(string name, hid_t group_id, size_t ncols,int np, int dim, valarray<double> &mass);

  void saveField(analysis &an, string name, double *data, hid_t group_id,double time=0);

  
 public:

 output(string output_directory,bool _verbose, bool _debug, double _delta_time);

 ~output(){};

 
 void init(initialData &id, size_t index, analysis &an);
  
 void update(double t, int index,evolution &evo,  analysis &ana, bool force=false);

  inline void close(){hsize_t status = H5Fclose(file_id);}
  inline bool write_output(double t, bool force = false){  return t >= next_output || t==0 || force; }

};


#endif
