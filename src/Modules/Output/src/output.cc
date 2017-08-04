

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





void output::init(initialData &id, size_t i,analysis &an)
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
  size_t np=id.get_number_of_particles();
  size_t dim=id.get_dimension();
  
  create_dataset("position", group_id, ncols,np,dim);
  create_dataset("momentum", group_id, ncols,np,dim);
  create_dataset("spin", group_id, ncols,np,dim);
  create_dataset("waves", group_id, 39,np,dim);

  create_dataset("dxdt", group_id, ncols,np,dim);
  create_dataset("dpdt", group_id, ncols,np,dim);
  create_dataset("dsdt", group_id, ncols,np,dim);
  
  create_dataset("ddxdt2", group_id, ncols,np,dim);
  create_dataset("ddpdt2", group_id, ncols,np,dim);
  create_dataset("ddsdt2", group_id, ncols,np,dim);

  
  H5Gclose(group_id);


  group_id = H5Gcreate(file_id,"coordinates", H5P_DEFAULT,  H5P_DEFAULT, H5P_DEFAULT);

  array3D  xx = an.getMeshX();
  saveField(an, "x", xx.data(), group_id);
  array3D  yy = an.getMeshY();
  saveField(an, "y", yy.data(), group_id);
  array3D  zz = an.getMeshZ();
  saveField(an, "z", zz.data(), group_id);

  H5Gclose(group_id);

  ss_grp << "/fields";
  group_id = H5Gcreate(file_id,ss_grp.str().c_str(), H5P_DEFAULT,  H5P_DEFAULT, H5P_DEFAULT);
  H5Gclose(group_id);

  
  /* Close the file. */
  herr_t status = H5Fflush(file_id,H5F_SCOPE_LOCAL);
  
  
  MPI::COMM_WORLD.Barrier ();

  iteration_set.push_back(0);
  
}


void output::create_dataset(string name, hid_t group_id, size_t ncols,int np,int dim)
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

  BOOST_LOG_SEV(lg, logging::trivial::debug) << "np: "<< np;
   
  herr_t attribute_id = H5Acreate (dset_id, "np",  H5T_NATIVE_INT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(attribute_id, H5T_NATIVE_INT, &np);
  H5Aclose(attribute_id);

  BOOST_LOG_SEV(lg, logging::trivial::debug) << "dim: "<< dim;

  attribute_id = H5Acreate (dset_id, "dim",  H5T_NATIVE_INT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(attribute_id, H5T_NATIVE_INT, &dim);
  H5Aclose(attribute_id);
  
  
  H5Sclose(dataspace_id);
  
  H5Pclose(plist_id);
  H5Sclose(filespace);    
  H5Dclose(dset_id);

  
}


void output::update(double t, int index,evolution &evo, analysis &an , bool force)
{


  if(write_output(t,force)){
    if(verbose)
      {
	BOOST_LOG_SEV(lg, logging::trivial::info) << "output at time t = " << t <<" - iteration "<< iteration;

      }


    save(evo.get_position(),"position",index);
    save(evo.get_momentum(),"momentum",index);
    save(evo.get_spin(),"spin",index);
    save(evo.get_waves(),"waves",index);

    save(evo.get_dxdt(),"dxdt",index);
    save(evo.get_dpdt(),"dpdt",index);
    save(evo.get_dsdt(),"dsdt",index);

    save(evo.get_ddxdt2(),"ddxdt2",index);
    save(evo.get_ddpdt2(),"ddpdt2",index);
    save(evo.get_ddsdt2(),"ddsdt2",index);


    stringstream ss_grp;
    ss_grp << "/"<< index << "/fields/"<< iteration_set[index];
    
    hid_t  group_id = H5Gcreate(file_id,ss_grp.str().c_str(), H5P_DEFAULT,  H5P_DEFAULT, H5P_DEFAULT);

    
    for(size_t i=0; i < DIMENSION; i++)
      for(size_t j=0; j < DIMENSION; j++)
	{
	  array3D  mm = an.getMetric(i,j);
	  stringstream metric_name;
	  metric_name << "metric_"<<i<<"_"<<j; 
	  saveField(an, metric_name.str(), mm.data(), group_id);
	}

    H5Gclose(group_id);
    
    next_output += delta_time;
    iteration++;
    iteration_set[index]+=1;
 
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

  hsize_t      size[NDIMS]={iteration_set[index]+1,ncols};
  
  status = H5Dset_extent (dset_id, size);

  /* Select a hyperslab in extended portion of dataset  */
  hid_t filespace = H5Dget_space (dset_id);

  hsize_t	offset[NDIMS] = {iteration_set[index],0};
  hsize_t	count [NDIMS] = {1,ncols};	         
  
  status = H5Sselect_hyperslab (filespace, H5S_SELECT_SET, offset, NULL,
				count, NULL);  

  /* Define memory space */
  hid_t memspace = H5Screate_simple (NDIMS, count, NULL); 

  hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);


  status = H5Dwrite(dset_id, //check
		    H5T_NATIVE_DOUBLE,
		    memspace, //check
		    filespace, //check
		    plist_id, //check
		    &data[0]); //check

  
  H5Dclose(dset_id);
  H5Gclose(group_id);
  
  
}



void output::saveField(analysis &mg, string name, double *data, hid_t group_id,double time)
{  

  BOOST_LOG_SEV(lg, logging::trivial::debug) << "Save field: "<< name;
  
  hsize_t dimsf[DIMENSION];
  dimsf[X] = (hsize_t)mg.getN(X);
  dimsf[Y] = (hsize_t)mg.getN(Y);
  dimsf[Z] = (hsize_t)mg.getN(Z);  

  hid_t dataspace_id = H5Screate_simple(DIMENSION, dimsf, NULL);  

  hid_t dset_id = H5Dcreate2(group_id, name.c_str(), H5T_NATIVE_DOUBLE, dataspace_id,
			    H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  herr_t  status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,  data);


  status = H5Sclose(dataspace_id);

   /* Close the first dataset. */
   status = H5Dclose(dset_id);

  
  
  
}
