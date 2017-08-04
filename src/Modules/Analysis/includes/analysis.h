

// analysis.h
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



#ifndef ANALYSIS_H

#define ANALYSIS_H

#include<utils.h>
#include<evolution.h>

class analysis
{

  bool w_metric; 
  vector<size_t> N; 

  vector<array3D> mesh;
  
  valarray<double> metric_data;
  vector< vector< array3Dref* > > metric; 

 public:

  analysis(bool _metric, double _cx, double _cy, double _cz, double _radius, vector<size_t> &resolution );
 ~analysis(){};

 void update(evolution &ev);

 inline  array3Dref & getMeshX() { return mesh[X];}
 inline  array3Dref & getMeshY() { return mesh[Y];}
 inline  array3Dref & getMeshZ() { return mesh[Z];}

 inline vector<size_t> getN(){ return N;}
 inline size_t getN(size_t x){ return N[x];}

 inline  array3Dref & getMetric(size_t i, size_t j ) { return *metric[i][j];}

 
};


#endif
