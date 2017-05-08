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


#include <main.h>



int main(int ac, char*av[])

{


  
  clock_t tStart = clock();
  
  string output_directory;
  
  double final_time;
  bool output_verbose;
  bool debug;
  double delta_time;

  
  
  try {

    MPI::Init ();


    po::options_description genericOptions("PostNewtonian 3BP \nAllowed options");

    genericOptions.add_options()
      ("help,h", "Shows a help message")
      ("input_file", po::value<string>(), "Is a parameter file")
      ("silent,s", "Shows additional debug messages in log")
      ("debug,d", "Shows additional debug messages in log");

    po::options_description evolutionOptions("Evolution options");

   evolutionOptions.add_options()
      ("evolution.verbose", po::value<bool>()->default_value(false),"Shows information about the evolution progress ")
      ("evolution.ode_method", po::value<string>()->default_value("rk8pd"),"ODE method, one of:[rk2,rk4,rkck,rk8pd,rk2imp,rk4imp,bsimp,gear1,gear2]")
     ("evolution.plane_constrain", po::value<bool>()->default_value(false),"Constrain evolution to a plane")
     ("evolution.chaos_test", po::value<bool>()->default_value(false),"Calculate  Lyapunov chaos indicator")
     ("evolution.jacA", po::value<bool>()->default_value(false),"Uses analytic Jacobian (in chaos test)")
     ("evolution.particles", po::value<size_t>()->default_value(3),"Number of particles [2,3]")
     ("evolution.initial_dt", po::value<double>()->default_value(1e-6),"Initial time-step")
     ("evolution.jacobian_dx", po::value<double>()->default_value(1e-6),"Numeric Jacobian initial dx ")
     ("evolution.factor_chaos_test", po::value<double>()->default_value(100),"Scaling factor for chaos test")
     ("evolution.scaling_variable", po::value<double>()->default_value(1),"Scale of variables on adaptive step size control")
     ("evolution.scaling_derivative", po::value<double>()->default_value(1),"Scale of deivative on adaptive  step size control")
     ("evolution.epsilon_rel", po::value<double>()->default_value(1e-6),"Relative error tolerance")
     ("evolution.epsilon_abs", po::value<double>()->default_value(1e-6),"Absolute error tolerance")
     ("evolution.final_time", po::value<double>(&final_time)->default_value(1),"Final evolution time");
    
    po::options_description outputOptions("Output options");

    outputOptions.add_options()
      ("output.directory,o", po::value<string>(&output_directory)->default_value("output"), "set output directory")
      ("output.debug", po::value<bool>(&debug)->default_value(false),"output aditional fields (debugging)")
      ("output.verbose", po::value<bool>(&output_verbose)->default_value(false),"output additional information")
      ("output.delta_time", po::value<double>(&delta_time)->default_value(1),"time interval of output");

    po::options_description termsOptions("Post-Newtonian terms options");

    termsOptions.add_options()
#ifdef odePN1
      ("terms.pn1", po::value<bool>()->default_value(false),"Activate post-Newtonian terms of order 1  ")
#endif
#ifdef odePN2
      ("terms.pn2", po::value<bool>()->default_value(false),"Activate post-Newtonian terms of order 2  ")
#endif
#ifdef odePN2_5
      ("terms.pn2_5", po::value<bool>()->default_value(false),"Activate post-Newtonian terms of order 2.5  ")
#endif
#ifdef odePNSlo
      ("terms.pnSOlo", po::value<bool>()->default_value(false),"Activate spin-orbit post-Newtonian terms of leading order  ")
      ("terms.pnSSlo", po::value<bool>()->default_value(false),"Activate spin-spin post-Newtonian terms of leading order  ")
      ("terms.pnS2lo", po::value<bool>()->default_value(false),"Activate self-spin post-Newtonian terms of leading order  ")
#endif
#ifdef odePNSnlo      
      ("terms.pnSOnlo", po::value<bool>()->default_value(false),"Activate spin-orbit post-Newtonian terms of next-to leading order  ")
      ("terms.pnSSnlo", po::value<bool>()->default_value(false),"Activate spin-spin post-Newtonian terms of next-to leading order  ")
#endif
      ;
    
    po::positional_options_description p;
    p.add("input_file", -1);

    
    po::options_description cmdline_options;
    cmdline_options.add(genericOptions).add(evolutionOptions).add(outputOptions).add(termsOptions);

    po::options_description config_file_options;
    config_file_options.add(evolutionOptions).add(outputOptions).add(termsOptions);


    po::variables_map vm;
    po::store(po::command_line_parser(ac, av).options(cmdline_options).positional(p).run(), vm);
    po::notify(vm);


    
    if (vm.count("help") || vm.count("input_file") == 0) {

      if (MPI::COMM_WORLD.Get_rank () == MPI_ROOT_NODE)
	{
	cout << cmdline_options << "\n";
	  if(vm.count("input_file") == 0 && vm.count("help")==0)
	    cout  << endl <<  "MISSING PARAMETER FILE!" << endl << endl;
	}
      return Finalize();
    }

    if (vm.count("input_file")){

      string input_file = vm["input_file"].as<string>();

      ifstream ifs(input_file.c_str());
      if (!ifs)
	{
	  cout << "can not open config file: " << input_file << "\n";
	  return 0;
	}
      else
	{
	  store(parse_config_file(ifs, config_file_options), vm);
	  notify(vm);
	}

    }


    makeOutputDirectory(output_directory);

    string cmd = "echo '#!/bin/bash' > " + output_directory + "/command.sh";
    int sys_out = system(cmd.c_str());

    cmd = "echo cd `pwd` >> "  + output_directory + "/command.sh";
    sys_out = system(cmd.c_str());
    
    stringstream param;
    
    for (int i = 0; i < ac; i++)

      param << av[i] << " ";

    cmd = "echo " + param.str() + " >> " + output_directory + "/command.sh;";
    sys_out = system(cmd.c_str());

    
    cmd = "chmod +x " + output_directory + "/command.sh";
    sys_out = system(cmd.c_str());

    
    setupLog(output_directory,vm.count("silent") > 0,vm.count("debug") > 0 );

    using namespace logging::trivial;
    src::severity_logger< severity_level > lg;

    
    evolution evol(vm);

       
    output my_output(output_directory,output_verbose, debug, delta_time);


    double t=0; 
    bool proceed = true;

    do{
      
      proceed=evol.update(t);

      BOOST_LOG_SEV(lg, info) << setprecision(5) << "time: "<< t;

      my_output.update(t);
      
    }while(proceed);

    my_output.update(t,true);
    
    
  }
  catch(exception& e) {
    cerr << "error: " << e.what() << "\n";
    Finalize();
    return 1;
  }
  catch(...) {
    cerr << "Exception of unknown type!\n";
  }

  using namespace logging::trivial;
  src::severity_logger< severity_level > lg;

  
  double ttime = (double)(clock() - tStart)/CLOCKS_PER_SEC;

  int thour = int(ttime /3600.0);
  int tmin = int( (ttime - thour*3600 )/60.0 );
  double tsec = ttime - thour*3600 - tmin*60;

  
  BOOST_LOG_SEV(lg, info) << "All done. Total execution time: "<< thour << " h " << tmin << " m " << tsec << " s" ;


  return Finalize();


}
