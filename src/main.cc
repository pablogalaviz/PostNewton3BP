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
  
  
  
  try {

    MPI::Init ();


    po::options_description genericOptions("PostNewtonian 3BP \nAllowed options");

    genericOptions.add_options()
      ("help,h", "Shows a help message")
      ("input_file", po::value<string>(), "Is a parameter file")
      ("silent,s", "Shows additional debug messages in log")
      ("debug,d", "Shows additional debug messages in log");


    bool initial_data_verbose; 
    
    po::options_description idOptions("Initial data options");
    idOptions.add_options()
      ("initial_data.type", po::value<string>()->default_value(""),"Type of initial data [read_file,binary_single] ")
      ("initial_data.verbose", po::value<bool>(&initial_data_verbose)->default_value(false),"output additional information")
      ("initial_data.particles", po::value<size_t>()->default_value(3),"Number of particles [2,3]")
      ("initial_data.file", po::value<string>()->default_value(""),"Data file with position, momentum and spin for each body and each simulation")
      ("initial_data.m1", po::value< vector<double> >()->default_value( vector<double>(), "1" ),"Mass of body 1 ")
      ("initial_data.m2", po::value< vector<double> >()->default_value( vector<double>(), "1" ),"Mass of body 2 ")
      ("initial_data.m3", po::value< vector<double> >()->default_value( vector<double>(), "1" ),"Mass of body 3 ")
      ("initial_data.separation_inner_binary", po::value< vector<double> >()->default_value( vector<double>(), "130" ),"Initial separation of the inner binary")
      ("initial_data.separation_outer_binary", po::value< vector<double> >()->default_value( vector<double>(), "400" ),"Initial separation of the outer binary")
      ("initial_data.theta_inner_binary", po::value< vector<double> >()->default_value( vector<double>(), "0" ),"Theta angle of the inner radio vector")
      ("initial_data.phi_inner_binary", po::value< vector<double> >()->default_value( vector<double>(), "0" ),"Phi angle of the inner radio vector")
      ("initial_data.theta_outer_binary", po::value< vector<double> >()->default_value( vector<double>(), "0" ),"Theta angle of the outer radio vector")
      ("initial_data.phi_outer_binary", po::value< vector<double> >()->default_value( vector<double>(), "1.5707963268" ),"Phi angle of the outer radio vector")
      ("initial_data.excentricity_inner_binary", po::value< vector<double> >()->default_value( vector<double>(), "0" ),"Newtonian eccentricity of the inner binary ")
      ("initial_data.excentricity_outer_binary", po::value< vector<double> >()->default_value( vector<double>(), "0" ),"Newtonian eccentricity of the outer binary ")
      ("initial_data.spin_1", po::value< vector<double> >()->default_value( vector<double>(), "0" ),"Spin magnitude of body 1")
      ("initial_data.spin_2", po::value< vector<double> >()->default_value( vector<double>(), "0" ),"Spin magnitude of body 2")
      ("initial_data.spin_3", po::value< vector<double> >()->default_value( vector<double>(), "0" ),"Spin magnitude of body 3")
      ("initial_data.theta_spin_1", po::value< vector<double> >()->default_value( vector<double>(), "0" ),"Spin theta angle of body 1")
      ("initial_data.phi_spin_1", po::value< vector<double> >()->default_value( vector<double>(), "0" ),"Spin phi angle of body 1")
      ("initial_data.theta_spin_2", po::value< vector<double> >()->default_value( vector<double>(), "0" ),"Spin theta angle of body 2")
      ("initial_data.phi_spin_2", po::value< vector<double> >()->default_value( vector<double>(), "0" ),"Spin phi angle of body 2")
      ("initial_data.theta_spin_3", po::value< vector<double> >()->default_value( vector<double>(), "0" ),"Spin theta angle of body 3")
      ("initial_data.phi_spin_3", po::value< vector<double> >()->default_value( vector<double>(), "0" ),"Spin phi angle of body 3")
      ("initial_data.epsilon", po::value< vector<double> >()->default_value( vector<double>(), "0" ),"epsilon");

    double cx;
    double cy;
    double cz;
    double radius; 
    size_t resolution; 
    bool calculate_metric;
    
    po::options_description analysisOptions("Analysis options");
    analysisOptions.add_options()
      ("analysis.calculate_metric", po::value< bool >(&calculate_metric)->default_value(false),"Whether compute the metric.")
      ("analysis.center_x", po::value< double >(&cx)->default_value(0),"Coordinate X of analysis domain.")
      ("analysis.center_y", po::value< double >(&cy)->default_value(0),"Coordinate Y of analysis domain.")
      ("analysis.center_z", po::value< double >(&cz)->default_value(0),"Coordinate Z of analysis domain.")
      ("analysis.resolution", po::value< size_t>(&resolution)->default_value(64),"Resolution of analysis domain.")
      ("analysis.radius", po::value< double >(&radius)->default_value(1000),"Radius of analysis domain.");

    
    bool evolution_verbose;
    string ode_method;
    bool chaos_test;
    bool JacA;
    double initial_dt;
    double Jacobian_dx; 
    double factor_chaos_test;
    double scaling_variable;
    double scaling_derivative;
    double eps_rel; 
    double eps_abs; 
    double final_time;
    
    po::options_description evolutionOptions("Evolution options");
    evolutionOptions.add_options()
      ("evolution.verbose", po::value<bool>(&evolution_verbose)->default_value(false),"Shows information about the evolution progress ")
      ("evolution.ode_method", po::value<string>(&ode_method)->default_value("rk8pd"),"ODE method, one of:[rk2,rk4,rkck,rk8pd,rk2imp,rk4imp,bsimp,gear1,gear2]")
      ("evolution.chaos_test", po::value<bool>(&chaos_test)->default_value(false),"Calculate  Lyapunov chaos indicator")
      ("evolution.JacA", po::value<bool>(&JacA)->default_value(false),"Uses analytic Jacobian (in chaos test)")
      ("evolution.initial_dt", po::value<double>(&initial_dt)->default_value(1e-6),"Initial time-step")
      ("evolution.Jacobian_dx", po::value<double>(&Jacobian_dx)->default_value(1e-6),"Numeric Jacobian initial dx ")
      ("evolution.factor_chaos_test", po::value<double>(&factor_chaos_test)->default_value(100),"Scaling factor for chaos test")
      ("evolution.scaling_variable", po::value<double>(&scaling_variable)->default_value(1),"Scale of variables on adaptive step size control")
      ("evolution.scaling_derivative", po::value<double>(&scaling_derivative)->default_value(1),"Scale of deivative on adaptive  step size control")
      ("evolution.epsilon_rel", po::value<double>(&eps_rel)->default_value(1e-6),"Relative error tolerance")
      ("evolution.epsilon_abs", po::value<double>(&eps_abs)->default_value(1e-6),"Absolute error tolerance")
      ("evolution.final_time", po::value<double>(&final_time)->default_value(1),"Final evolution time");

    string output_directory;
    bool output_debug;
    bool output_verbose;
    double delta_time;
    
    po::options_description outputOptions("Output options");    
    outputOptions.add_options()
      ("output.directory,o", po::value<string>(&output_directory)->default_value("output"), "set output directory")
      ("output.debug", po::value<bool>(&output_debug)->default_value(false),"output aditional fields (debugging)")
      ("output.verbose", po::value<bool>(&output_verbose)->default_value(false),"output additional information")
      ("output.delta_time", po::value<double>(&delta_time)->default_value(1),"time interval of output");

    terms_t pn_terms;
    
    bool terms_verbose; 
    
    po::options_description termsOptions("Post-Newtonian terms options");
    termsOptions.add_options()
#ifdef odePN1
      ("terms.pn1", po::value<bool>(&terms.pn1)->default_value(false),"Activate post-Newtonian terms of order 1  ")
#endif
#ifdef odePN2
      ("terms.pn2", po::value<bool>(&pn2)->default_value(false),"Activate post-Newtonian terms of order 2  ")
#endif
#ifdef odePN2_5
      ("terms.pn2_5", po::value<bool>(&pn2_5)->default_value(false),"Activate post-Newtonian terms of order 2.5  ")
#endif
#ifdef odePNSlo
      ("terms.pnSOlo", po::value<bool>(&pnSOlo)->default_value(false),"Activate spin-orbit post-Newtonian terms of leading order  ")
      ("terms.pnSSlo", po::value<bool>(&pnSSlo)->default_value(false),"Activate spin(a)-spin(b) post-Newtonian terms of leading order  ")
      ("terms.pnS2lo", po::value<bool>(&pnS2lo)->default_value(false),"Activate spin(a)-spin(a) post-Newtonian terms of leading order  ")
#endif
#ifdef odePNSnlo      
      ("terms.pnSOnlo", po::value<bool>(&pnSOnlo)->default_value(false),"Activate spin-orbit post-Newtonian terms of next-to leading order  ")
      ("terms.pnSSnlo", po::value<bool>(&pnSSnlo)->default_value(false),"Activate spin(a)-spin(b) post-Newtonian terms of next-to leading order  ")
#endif
      ("terms.verbose", po::value<bool>(&terms_verbose)->default_value(false),"output additional information");

    
    
    po::positional_options_description p;
    p.add("input_file", -1);

    
    po::options_description cmdline_options;
    cmdline_options.add(genericOptions).add(evolutionOptions).add(outputOptions).add(termsOptions).add(idOptions).add(analysisOptions);

    po::options_description config_file_options;
    config_file_options.add(evolutionOptions).add(outputOptions).add(termsOptions).add(idOptions).add(analysisOptions);


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

    initialData id(vm);

    int n =  id.sim_size();

    output my_output(output_directory,output_verbose, output_debug, delta_time);    

    analysis my_analysis(calculate_metric,cx,cy,cz,radius, resolution);
    
    for(int i=0; i < n; i++)
      {


	valarray<double> _id_variables=id.get_variables(i);
	valarray<double> _mass = id.get_mass(i);
	
	evolution evol(evolution_verbose,
		       ode_method,
		       chaos_test,
		       JacA,
		       initial_dt,
		       Jacobian_dx,
		       factor_chaos_test,
		       scaling_variable,
		       scaling_derivative,
		       eps_rel,
		       eps_abs,
		       final_time,
		       _id_variables,
		       _mass,
		       id.has_spin(),
		       pn_terms
		       );
    

	my_output.init(id,i,my_analysis);
    
	double t=0; 
	bool proceed = true;
	do{

	  if(my_output.write_output(t))
	    my_analysis.update(evol);
	      
	  my_output.update(t,i,evol,my_analysis);
	  
	  proceed=evol.update(t);
	  if(evolution_verbose)
	    BOOST_LOG_SEV(lg, info) << setprecision(5) << "time: "<< t;
	  
      
	}while(proceed);

	my_analysis.update(evol);
	my_output.update(t,i,evol,my_analysis,true);
	my_output.close();
	
      }    
    
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
