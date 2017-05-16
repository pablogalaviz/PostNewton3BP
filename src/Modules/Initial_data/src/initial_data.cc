

// initial_data.cc
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



#include "initial_data.h"



initialData::initialData(po::variables_map &vm){

  number_of_particles = vm["initial_data.particles"].as<size_t>();
  
  string id_type= vm["initial_data.type"].as<string>();

  verbose= vm["initial_data.verbose"].as<bool>();
      
  
  if (id_type.compare("read-file") == 0 ){

    string input_file = vm["initial_data.file"].as<string>();
      load_file(input_file); 
  }


}



void initialData::load_file(string filename){


  ifstream ifs(filename.c_str());
  if (!ifs)
    {

      BOOST_LOG_SEV(lg, error) << "Can not open config file : "<< filename;
      exit(Finalize(0));
      
    }
  else
    {
      po::options_description varOptions("Variables option");

      vector<double> x1; 
      vector<double> y1; 
      vector<double> z1; 

      vector<double> px1; 
      vector<double> py1; 
      vector<double> pz1; 
      
      vector<double> sx1; 
      vector<double> sy1; 
      vector<double> sz1; 

      vector<double> x2; 
      vector<double> y2; 
      vector<double> z2; 

      vector<double> px2; 
      vector<double> py2; 
      vector<double> pz2;
      
      vector<double> sx2; 
      vector<double> sy2; 
      vector<double> sz2; 

      vector<double> x3; 
      vector<double> y3; 
      vector<double> z3; 

      vector<double> px3; 
      vector<double> py3; 
      vector<double> pz3;

      vector<double> sx3; 
      vector<double> sy3; 
      vector<double> sz3; 


      vector<double> m1; 
      vector<double> m2; 
      vector<double> m3;

      
      varOptions.add_options()
	("position.x1", po::value< vector<double> >(&x1),"Position x of body 1 ")
	("position.y1", po::value< vector<double> >(&y1),"Position y of body 1 ")
	("position.z1", po::value< vector<double> >(&z1),"Position z of body 1 ")
	("position.x2", po::value< vector<double> >(&x2),"Position x of body 2 ")
	("position.y2", po::value< vector<double> >(&y2),"Position y of body 2 ")
	("position.z2", po::value< vector<double> >(&z2),"Position z of body 2 ")
	("position.x3", po::value< vector<double> >(&x3),"Position x of body 3 ")
	("position.y3", po::value< vector<double> >(&y3),"Position y of body 3 ")
	("position.z3", po::value< vector<double> >(&z3),"Position z of body 3 ")
	("momentum.x1", po::value< vector<double> >(&px1),"Momentum x of body 1 ")
	("momentum.y1", po::value< vector<double> >(&py1),"Momentum y of body 1 ")
	("momentum.z1", po::value< vector<double> >(&pz1),"Momentum z of body 1 ")
	("momentum.x2", po::value< vector<double> >(&px2),"Momentum x of body 2 ")
	("momentum.y2", po::value< vector<double> >(&py2),"Momentum y of body 2 ")
	("momentum.z2", po::value< vector<double> >(&pz2),"Momentum z of body 2 ")
	("momentum.x3", po::value< vector<double> >(&px3),"Momentum x of body 3 ")
	("momentum.y3", po::value< vector<double> >(&py3),"Momentum y of body 3 ")
	("momentum.z3", po::value< vector<double> >(&pz3),"Momentum z of body 3 ")
	("spin.x1", po::value< vector<double> >(&sx1),"Spin x of body 1 ")
	("spin.y1", po::value< vector<double> >(&sy1),"Spin y of body 1 ")
	("spin.z1", po::value< vector<double> >(&sz1),"Spin z of body 1 ")
	("spin.x2", po::value< vector<double> >(&sx2),"Spin x of body 2 ")
	("spin.y2", po::value< vector<double> >(&sy2),"Spin y of body 2 ")
	("spin.z2", po::value< vector<double> >(&sz2),"Spin z of body 2 ")
	("spin.x3", po::value< vector<double> >(&sx3),"Spin x of body 3 ")
	("spin.y3", po::value< vector<double> >(&sy3),"Spin y of body 3 ")
	("spin.z3", po::value< vector<double> >(&sz3),"Spin z of body 3 ")
	("mass.m1", po::value< vector<double> >(&m1),"Mass of body 1 ")
	("mass.m2", po::value< vector<double> >(&m2),"Mass of body 2 ")
	("mass.m3", po::value< vector<double> >(&m3),"Mass of body 3 ")
	;
	

      po::variables_map vm;      
      store(parse_config_file(ifs, varOptions), vm);
      notify(vm);


      int pos1 = max(x1.size(),max(y1.size(),z1.size())); 
      int pos2 = max(x2.size(),max(y2.size(),z2.size())); 
      int pos3 = max(x3.size(),max(y3.size(),z3.size())); 
      
      int mom1 = max(px1.size(),max(py1.size(),pz1.size())); 
      int mom2 = max(px2.size(),max(py2.size(),pz2.size())); 
      int mom3 = max(px3.size(),max(py3.size(),pz3.size())); 

      int s1 = max(sx1.size(),max(sy1.size(),sz1.size())); 
      int s2 = max(sx2.size(),max(sy2.size(),sz2.size())); 
      int s3 = max(sx3.size(),max(sy3.size(),sz3.size())); 

      if( m1.size() ==0 || m2.size() == 0)
	{
	  
	  BOOST_LOG_SEV(lg, error) << "Missing mass parameter for m1 or m2: "<< filename;
	  exit(Finalize(0));

	}

      if( number_of_particles == 3 && m3.size() == 0)
	{
	  BOOST_LOG_SEV(lg, error) << "Missing mass parameter for m3 "<< filename;
	  exit(Finalize(0));
	}
	
      dimension=3;

      if(s1+s2+s3+z1.size()+z2.size()+z3.size()+pz1.size()+pz2.size()+pz3.size() == 0)
	dimension=2;

      int num_var = dimension*number_of_particles*2;

      spin=s1+s2+s3 > 0; 
      
      if(spin)
	num_var+=dimension;


      simulation_size = max(m1.size(),max(m2.size(),m3.size()));

      simulation_size =max(simulation_size,max(pos1,max(pos2,max(pos3,max(mom1,max(mom2,max(mom3,max(s1,max(s2,s3)))))))));


      stringstream ss; 

      if(verbose)
	{
	  BOOST_LOG_SEV(lg, info)  << "Initial data info:";
	  BOOST_LOG_SEV(lg, info)  << "======================================";
	  BOOST_LOG_SEV(lg, info)  << "Number of simulations: " << simulation_size;
	  BOOST_LOG_SEV(lg, info)  << "Dimensions: " << dimension;
	  BOOST_LOG_SEV(lg, info)  << "Number of particles: " << number_of_particles;
	  BOOST_LOG_SEV(lg, info)  << "Spin: " << spin;
	  BOOST_LOG_SEV(lg, info)  << "Total number of variables: " << num_var;
	  if(m1.size()>0)
	    BOOST_LOG_SEV(lg, info)  << "Particle 1 mass: " << vector2str<double>(m1);
	  if(x1.size()>0)
	    BOOST_LOG_SEV(lg, info)  << "Particle 1 position  x: " << vector2str<double>(x1);
	  if(y1.size()>0)
	    BOOST_LOG_SEV(lg, info)  << "Particle 1 position  y: " << vector2str<double>(y1);
	  if(z1.size()>0)
	    BOOST_LOG_SEV(lg, info)  << "Particle 1 position  z: " << vector2str<double>(z1);
	  if(px1.size()>0)
	    BOOST_LOG_SEV(lg, info)  << "Particle 1 momentum  x: " << vector2str<double>(px1);
	  if(py1.size()>0)
	    BOOST_LOG_SEV(lg, info)  << "Particle 1 momentum  y: " << vector2str<double>(py1);
	  if(pz1.size()>0)
	    BOOST_LOG_SEV(lg, info)  << "Particle 1 momentum  z: " << vector2str<double>(pz1);
	  if(sx1.size()>0)
	    BOOST_LOG_SEV(lg, info)  << "Particle 1 spin  x: " << vector2str<double>(sx1);
	  if(sy1.size()>0)
	    BOOST_LOG_SEV(lg, info)  << "Particle 1 spin  y: " << vector2str<double>(sy1);
	  if(sz1.size()>0)
	    BOOST_LOG_SEV(lg, info)  << "Particle 1 spin  z: " << vector2str<double>(sz1);
	  if(m2.size()>0)
	    BOOST_LOG_SEV(lg, info)  << "Particle 2 mass: " << vector2str<double>(m2);
	  if(x2.size()>0)
	    BOOST_LOG_SEV(lg, info)  << "Particle 2 position  x: " << vector2str<double>(x2);
	  if(y2.size()>0)
	    BOOST_LOG_SEV(lg, info)  << "Particle 2 position  y: " << vector2str<double>(y2);
	  if(z2.size()>0)
	    BOOST_LOG_SEV(lg, info)  << "Particle 2 position  z: " << vector2str<double>(z2);
	  if(px2.size()>0)
	    BOOST_LOG_SEV(lg, info)  << "Particle 2 momentum  x: " << vector2str<double>(px2);
	  if(py2.size()>0)
	    BOOST_LOG_SEV(lg, info)  << "Particle 2 momentum  y: " << vector2str<double>(py2);
	  if(pz2.size()>0)
	    BOOST_LOG_SEV(lg, info)  << "Particle 2 momentum  z: " << vector2str<double>(pz2);
	  if(sx2.size()>0)
	    BOOST_LOG_SEV(lg, info)  << "Particle 2 spin  x: " << vector2str<double>(sx2);
	  if(sy2.size()>0)
	    BOOST_LOG_SEV(lg, info)  << "Particle 2 spin  y: " << vector2str<double>(sy2);
	  if(sz2.size()>0)
	    BOOST_LOG_SEV(lg, info)  << "Particle 2 spin  z: " << vector2str<double>(sz2);
	  if(m3.size()>0)
	    BOOST_LOG_SEV(lg, info)  << "Particle 2 mass: " << vector2str<double>(m3);
	  if(x3.size()>0)
	    BOOST_LOG_SEV(lg, info)  << "Particle 3 position  x: " << vector2str<double>(x3);
	  if(y3.size()>0)
	    BOOST_LOG_SEV(lg, info)  << "Particle 3 position  y: " << vector2str<double>(y3);
	  if(z3.size()>0)
	    BOOST_LOG_SEV(lg, info)  << "Particle 3 position  z: " << vector2str<double>(z3);
	  if(px3.size()>0)
	    BOOST_LOG_SEV(lg, info)  << "Particle 3 momentum  x: " << vector2str<double>(px3);
	  if(py3.size()>0)
	    BOOST_LOG_SEV(lg, info)  << "Particle 3 momentum  y: " << vector2str<double>(py3);
	  if(pz3.size()>0)
	    BOOST_LOG_SEV(lg, info)  << "Particle 3 momentum  z: " << vector2str<double>(pz3);
	  if(sx3.size()>0)
	    BOOST_LOG_SEV(lg, info)  << "Particle 3 spin  x: " << vector2str<double>(sx3);
	  if(sy3.size()>0)
	    BOOST_LOG_SEV(lg, info)  << "Particle 3 spin  y: " << vector2str<double>(sy3);
	  if(sz3.size()>0)
	    BOOST_LOG_SEV(lg, info)  << "Particle 3 spin  z: " << vector2str<double>(sz3);
	  
	}


      
      for(int i=0; i < simulation_size; i++)
	{
	  valarray<double> _y(num_var);
	  valarray<double> _par(number_of_particles);

	  _y[r_index(P1,X,dimension)]=get_id_value(x1,i);
	  _y[r_index(P1,Y,dimension)]=get_id_value(y1,i);

	  _y[r_index(P2,X,dimension)]=get_id_value(x2,i);
	  _y[r_index(P2,Y,dimension)]=get_id_value(y2,i);
	  
	  _y[p_index(P1,X,dimension)]=get_id_value(px1,i);
	  _y[p_index(P1,Y,dimension)]=get_id_value(py1,i);

	  _y[p_index(P2,X,dimension)]=get_id_value(px2,i);
	  _y[p_index(P2,Y,dimension)]=get_id_value(py2,i);

	  _par[P1] = get_id_value(m1,i);
	  _par[P2] = get_id_value(m2,i);
	  

	  if(dimension == 3)
	    {

	      _y[r_index(P1,Z,dimension)]=get_id_value(z1,i);
	      _y[r_index(P2,Z,dimension)]=get_id_value(z2,i);
	  
	      _y[p_index(P1,Z,dimension)]=get_id_value(pz1,i);
	      _y[p_index(P2,Z,dimension)]=get_id_value(pz2,i);

	      if(spin)
		{
		  _y[s_index(P1,X,dimension)]=get_id_value(sx1,i);
		  _y[s_index(P1,Y,dimension)]=get_id_value(sy1,i);
		  _y[s_index(P1,Z,dimension)]=get_id_value(sz1,i);

		  _y[s_index(P2,X,dimension)]=get_id_value(sx2,i);
		  _y[s_index(P2,Y,dimension)]=get_id_value(sy2,i);
		  _y[s_index(P2,Z,dimension)]=get_id_value(sz2,i);

	      }

	    }
	  
	  if(number_of_particles ==3)
	    {

	      _y[r_index(P3,X,dimension)]=get_id_value(x3,i);
	      _y[r_index(P3,Y,dimension)]=get_id_value(y3,i);
	      	      
	      _y[p_index(P3,X,dimension)]=get_id_value(px3,i);
	      _y[p_index(P3,Y,dimension)]=get_id_value(py3,i);

	      _par[P3] = get_id_value(m3,i);
		  

	      if(dimension == 3)
		{
		  _y[r_index(P3,Z,dimension)]=get_id_value(z3,i);
		  _y[p_index(P3,Z,dimension)]=get_id_value(pz3,i);

		  if(spin)
		    {
		      _y[s_index(P3,X,dimension)]=get_id_value(sx3,i);
		      _y[s_index(P3,Y,dimension)]=get_id_value(sy3,i);
		      _y[s_index(P3,Z,dimension)]=get_id_value(sz3,i);
		    }


		}
	  
	    }

	  BOOST_LOG_SEV(lg,debug)  << "Simulation: " << i;
	  BOOST_LOG_SEV(lg,debug)  << "Variables:  "<< valarray2str<double>(_y);
	  BOOST_LOG_SEV(lg,debug)  << "Parameters: "<< valarray2str<double>(_par);

	  
	  y.push_back(_y);
	  par.push_back(_par);
	  
	  
	}
    }
}


