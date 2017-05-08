

// utils.cc
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



#include "utils.h"



string PN3BP_Logo(){

  stringstream result;


  result << endl
	 << "--------------------------------------------------------------------------" << endl
	 <<  "_____          _   _   _               _              ____  ____  _____  " << endl
	 << "|  __ \\        | | | \\ | |             | |            |___ \\|  _ \\|  __ \\ " << endl
	 << "| |__) |__  ___| |_|  \\| | _____      _| |_ ___  _ __   __) | |_) | |__) |" << endl
	 << "|  ___/ _ \\/ __| __| . ` |/ _ \\ \\ /\\ / / __/ _ \\| '_ \\ |__ <|  _ <|  ___/ " << endl
	 << "| |  | (_) \\__ \\ |_| |\\  |  __/\\ \\V  \\V /| || (_) | | | |___) | |_) | |     " << endl
	 << "|_|   \\___/|___/\\__|_| \\_|\\___| \\_/\\_/  \\__\\___/|_| |_|____/|____/|_|     " << endl
	 << "--------------------------------------------------------------------------" << endl
	 << endl
	 << endl;

  result << "                   =======  PN3BP 16.01  =======" << endl
         << endl
         << "                   Author: Pablo Galaviz    " << endl
         << endl
         << "                   Email:  Pablo.Galaviz@me.com  " << endl
         << endl
         << "                   ==============================" << endl
         << endl;



  return result.str();



}



int Finalize(int e){

  MPI::COMM_WORLD.Barrier ();

  MPI::Finalize ();

  return e;

}


void setupLog(string dir_name, bool silent, bool debug){

  int mpi_rank = MPI::COMM_WORLD.Get_rank ();

  string log_file = "/output"+to_string(mpi_rank)+".log";


  logging::add_file_log
    (
     keywords::file_name = dir_name+log_file,                                        /*< file name pattern >*/
     keywords::rotation_size = 100 * 1024 * 1024,                                   /*< rotate files every 100 MiB... >*/
     keywords::format =  (expr::stream
			  << expr::format_date_time< boost::posix_time::ptime >("TimeStamp", "%Y-%m-%d %H:%M:%S")
			  << " PN3BP [cpu "<< mpi_rank << " " << logging::trivial::severity
			  << "] | " << expr::smessage),                            /*< log record format >*/
     boost::log::keywords::auto_flush = true
     );

  if(!silent && mpi_rank == MPI_ROOT_NODE)
    logging::add_console_log(cout, keywords::format = (expr::stream
						       << expr::format_date_time< boost::posix_time::ptime >("TimeStamp", "%Y-%m-%d %H:%M:%S")
						       << " PN3BP [cpu "<< mpi_rank << " " << logging::trivial::severity
						       << "] | " << expr::smessage)  );

  if(debug)
    logging::core::get()->set_filter
      (
       logging::trivial::severity >= logging::trivial::debug
       );
  else
    logging::core::get()->set_filter
      (
       logging::trivial::severity >= logging::trivial::info
       );


  logging::add_common_attributes();

  BOOST_LOG_TRIVIAL(info) << PN3BP_Logo();
  BOOST_LOG_TRIVIAL(info) << "Project directory: "<< dir_name;


}


void makeOutputDirectory(string dir_name){

  if(MPI::COMM_WORLD.Get_rank () == MPI_ROOT_NODE)
    {

      string rmdir = "rm -rf " + dir_name + "_prev";


      int sys_out = system(rmdir.c_str());



      string mvdir = "mv -f " + dir_name + " " + dir_name + "_prev";


      sys_out = system(mvdir.c_str());


      string mkdir = "mkdir " + dir_name;


      sys_out = system(mkdir.c_str());
    }

  MPI::COMM_WORLD.Barrier ();


}




void central_deriv (double fm1, double fp1, 
		    double fmh, double fph, 
		    double x,
		    double h, double *result, 
		    double *abserr_round, double *abserr_trunc)
{
  /* Compute the derivative using the 5-point rule (x-h, x-h/2, x,
     x+h/2, x+h). Note that the central point is not used.  

     Compute the error using the difference between the 5-point and
     the 3-point rule (x-h,x,x+h). Again the central point is not
     used. */

  
  double r3 = 0.5 * (fp1 - fm1);
  double r5 = (4.0 / 3.0) * (fph - fmh) - (1.0 / 3.0) * r3;

  double e3 = (fabs (fp1) + fabs (fm1)) * GSL_DBL_EPSILON;
  double e5 = 2.0 * (fabs (fph) + fabs (fmh)) * GSL_DBL_EPSILON + e3;

  double dy = GSL_MAX (fabs (r3), fabs (r5)) * fabs (x) * GSL_DBL_EPSILON;

  /* The truncation error in the r5 approximation itself is O(h^4).
     However, for safety, we estimate the error from r5-r3, which is
     O(h^2).  By scaling h we will minimise this estimated error, not
     the actual truncation error in r5. */

  *result = r5 / h;
  *abserr_trunc = fabs ((r5 - r3) / h); /* Estimated truncation error O(h^2) */
  *abserr_round = fabs (e5 / h) + dy;   /* Rounding error (cancellations) */
}



