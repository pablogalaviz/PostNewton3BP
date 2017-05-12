

// output.cc
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



#include "output.h"


output::output(string output_directory,bool _verbose, bool _debug, double _delta_time){


  tStart = clock();

  verbose=_verbose;
  
  debug=_debug;
  
  delta_time=_delta_time;

  next_output = 0;

  
  
  BOOST_LOG_SEV(lg, logging::trivial::info) << "----------- output setup -------------";

  file_name = output_directory+"/output.h5";

  hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_id, MPI::COMM_WORLD, MPI::INFO_NULL);  
  file_id = H5Fcreate(file_name.c_str(),H5F_ACC_TRUNC,H5P_DEFAULT,plist_id);
  H5Pclose(plist_id);

  monitor_file_name = output_directory+"/monitor.gp";

  if(debug)
    BOOST_LOG_SEV(lg, logging::trivial::warning) << "debugging output activated";

  if(verbose)
    BOOST_LOG_SEV(lg, logging::trivial::warning) << "verbose output activated";

  BOOST_LOG_SEV(lg, logging::trivial::info) << "output delta set to: "<< delta_time;


}





void output::init(initialData &id, size_t i)
{


  BOOST_LOG_SEV(lg, logging::trivial::info) << "------------ output init -------------";

  int mpi_rank = MPI::COMM_WORLD.Get_rank ();
  int mpi_size = MPI::COMM_WORLD.Get_size ();

  stringstream ss_gpr;
  ss_gpr << "/"<< i;  
  hid_t group_id = H5Gcreate(file_id, ss_gpr.str().c_str(), H5P_DEFAULT,  H5P_DEFAULT, H5P_DEFAULT);
  H5Gclose(group_id);

  
  stringstream ss_position;
  ss_position << "/"<< i <<"/position";  
   group_id = H5Gcreate(file_id, ss_position.str().c_str(), H5P_DEFAULT,  H5P_DEFAULT, H5P_DEFAULT);
  H5Gclose(group_id);

  

  
  stringstream ss_momentum;
  ss_momentum << "/"<< i <<"/momentum";  
  group_id = H5Gcreate(file_id, ss_momentum.str().c_str(), H5P_DEFAULT,  H5P_DEFAULT, H5P_DEFAULT);
  H5Gclose(group_id);

  if(id.has_spin())
    {
      stringstream ss_spin;
      ss_spin << "/"<< i <<"/spin";  

      group_id = H5Gcreate(file_id, ss_spin.str().c_str(), H5P_DEFAULT,  H5P_DEFAULT, H5P_DEFAULT);
      H5Gclose(group_id);
    }

  
  MPI::COMM_WORLD.Barrier ();

}


void output::update(double t, int index,evolution &evo , bool force)
{

  iteration++;

  if(t >= next_output || t==0 || force){
    if(verbose)
      {
	BOOST_LOG_SEV(lg, logging::trivial::info) << "output at time t = " << t <<" - iteration "<< iteration;

      }

    save(t,evo.get_position(),evo.get_momentum(),evo.get_spin());
    
    next_output += delta_time;

    
  }
  
}

void output::save(double t,valarray<double> pos, valarray<double> mom,valarray<double> spin)
{

  cout << valarray2str(pos) << endl;

  

}
