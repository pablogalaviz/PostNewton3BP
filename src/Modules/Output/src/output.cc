

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

  iteration=0;
  
  
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

  /* Close the file. */
  herr_t status = H5Fclose(file_id);

  
}





void output::init(initialData &id, size_t i)
{


  BOOST_LOG_SEV(lg, logging::trivial::info) << "------------ output init -------------";

  int mpi_rank = MPI::COMM_WORLD.Get_rank ();
  int mpi_size = MPI::COMM_WORLD.Get_size ();

  stringstream ss_grp;
  ss_grp << "/"<< i;

  BOOST_LOG_SEV(lg, logging::trivial::debug) << "Group name: "<< i << " - " << ss_grp.str();

  hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_id, MPI::COMM_WORLD, MPI::INFO_NULL);  
  file_id = H5Fopen(file_name.c_str(),H5F_ACC_RDWR,plist_id);
  H5Pclose(plist_id);
  
  hid_t group_id = H5Gcreate(file_id, ss_grp.str().c_str(), H5P_DEFAULT,  H5P_DEFAULT, H5P_DEFAULT);

  size_t ncols=id.get_number_of_particles()*id.get_dimension()+1;

  
  create_dataset("position", group_id, ncols);
  create_dataset("momentum", group_id, ncols);
  create_dataset("spin", group_id, ncols);
  create_dataset("waves", group_id, 39);
  
 
  H5Gclose(group_id);
  
  /* Close the file. */
  herr_t status = H5Fflush(file_id,H5F_SCOPE_LOCAL);
  
  
  MPI::COMM_WORLD.Barrier ();

  iteration_set.push_back(0);
  
}


void output::create_dataset(string name, hid_t group_id, size_t ncols)
{

  
  hsize_t dimsf[NDIMS] = {0,ncols};
  
  hsize_t chunk_dims[NDIMS] = {1,ncols};

  hsize_t maxdims[NDIMS] = {H5S_UNLIMITED, ncols};
  
  
  hid_t filespace = H5Screate_simple(NDIMS, dimsf, maxdims); 

  
  hid_t plist_id = H5Pcreate(H5P_DATASET_CREATE);  
  H5Pset_chunk(plist_id, NDIMS, chunk_dims);

  hid_t dset_id = H5Dcreate(group_id, name.c_str(), H5T_NATIVE_DOUBLE, filespace,H5P_DEFAULT, plist_id, H5P_DEFAULT);

  hsize_t dims = 1;

  hid_t dataspace_id = H5Screate_simple(1, &dims, NULL);

  herr_t attribute_id = H5Acreate (dset_id, "columns",  H5T_IEEE_F64BE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);

  double zz=0.0;
  H5Awrite(attribute_id, H5T_NATIVE_DOUBLE, &zz);
  
  H5Aclose(attribute_id);
  
  
  H5Pclose(plist_id);
  H5Sclose(filespace);    
  H5Dclose(dset_id);

  
}


void output::update(double t, int index,evolution &evo , bool force)
{

  iteration++;

  if(t >= next_output || t==0 || force){
    if(verbose)
      {
	BOOST_LOG_SEV(lg, logging::trivial::info) << "output at time t = " << t <<" - iteration "<< iteration;

      }

    iteration_set[index]+=1;

    save(evo.get_position(),"position",index);
    save(evo.get_momentum(),"momentum",index);
    save(evo.get_spin(),"spin",index);
    save(evo.get_waves(),"waves",index);

    
    next_output += delta_time;

    
  }
  
}

void output::save(valarray<double> data, string group_name , size_t index)
{

  
  stringstream ss_grp;
  ss_grp << "/"<< index;  


  herr_t status;  
  
  hid_t group_id = H5Gopen2(file_id, ss_grp.str().c_str(), H5P_DEFAULT);
  

  hid_t dset_id = H5Dopen (group_id, group_name.c_str(), H5P_DEFAULT);
  size_t ncols=data.size();

  hsize_t      size[NDIMS]={iteration_set[index],ncols};
  
  status = H5Dset_extent (dset_id, size);

  /* Select a hyperslab in extended portion of dataset  */
  hid_t filespace = H5Dget_space (dset_id);

  hsize_t	offset[NDIMS] = {iteration_set[index]-1,0};
  hsize_t	count [NDIMS] = {1,ncols};	         
  
  status = H5Sselect_hyperslab (filespace, H5S_SELECT_SET, offset, NULL,
				count, NULL);  

  /* Define memory space */
  hid_t memspace = H5Screate_simple (NDIMS, count, NULL); 

  hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);

  BOOST_LOG_SEV(lg, logging::trivial::debug) << "data_size: "<< data.size();

  status = H5Dwrite(dset_id, //check
		    H5T_NATIVE_DOUBLE,
		    memspace, //check
		    filespace, //check
		    plist_id, //check
		    &data[0]); //check

  BOOST_LOG_SEV(lg, logging::trivial::debug) << "done";
  
  H5Dclose(dset_id);
  H5Gclose(group_id);
  
  
}
