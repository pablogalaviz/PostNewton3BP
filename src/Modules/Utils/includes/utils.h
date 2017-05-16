

// utils.h
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



#ifndef UTILS_H

#define UTILS_H

#include <mpi.h>
#include <fstream>
#include <valarray>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_rng.h>

#include <boost/log/utility/setup/console.hpp>
#include <boost/log/support/date_time.hpp>
#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/sinks/text_file_backend.hpp>
#include <boost/log/utility/setup/file.hpp>
#include <boost/log/utility/setup/common_attributes.hpp>
#include <boost/log/sources/severity_logger.hpp>
#include <boost/log/sources/record_ostream.hpp>
#include <boost/program_options.hpp>


#define MPI_ROOT_NODE 0

#define X 0
#define Y 1
#define Z 2

#define P1 0
#define P2 1
#define P3 2

namespace po = boost::program_options;
namespace logging = boost::log;
namespace src = boost::log::sources;
namespace sinks = boost::log::sinks;
namespace keywords = boost::log::keywords;
namespace expr = boost::log::expressions;


using namespace logging::trivial;
using namespace std;

string PN3BP_Logo();

int Finalize(int e=0);

void setupLog(string dir_name, bool silent, bool debug);

void makeOutputDirectory(string dir_name);

void central_deriv (double fm1, double fp1, 
		    double fmh, double fph, 
		    double x,
		    double h, double *result, 
		    double *abserr_round, double *abserr_trunc);

double get_id_value(vector<double> &v, int i, double def=0);

size_t r_index(size_t a, size_t i, size_t space_dim);
size_t p_index(size_t a, size_t i, size_t space_dim);
size_t s_index(size_t a, size_t i, size_t space_dim);

struct terms_t {
  bool pn1 = false; 
  bool pn2 = false; 
  bool pn2_5 = false; 
  bool pnSOlo = false; 
  bool pnSSlo = false; 
  bool pnS2lo = false; 
  bool pnSOnlo = false; 
  bool pnSSnlo = false;
};

template<class T> 
string vector2str(vector<T> v) {

  stringstream ss;
  for(int i=0; i < v.size()-1; i++)
    ss << v[i] << ", "; 
  ss << v[v.size()-1];
  
  return ss.str();

}


template<class T> 
string valarray2str(valarray<T> v) {

  stringstream ss;
  for(int i=0; i < v.size()-1; i++)
    ss << v[i] << ", "; 
  ss << v[v.size()-1];
  
  return ss.str();

}


#endif

