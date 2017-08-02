

// analysis.cc
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



#include "analysis.h"


analysis::analysis(bool _metric, double _cx, double _cy, double _cz, double _radius, size_t _resolution){
  w_metric = _metric;
  N=_resolution;

  double dxyz = 2*_radius/(N-1.0);
  
  array3D m_x(boost::extents[N][N][N]);
  array3D m_y(boost::extents[N][N][N]);
  array3D m_z(boost::extents[N][N][N]);

  for(index3D i=0; i<N; i++)
    for(index3D j=0; j<N; j++)
      for(index3D k=0; k<N; k++)
	{
	  m_x[i][j][k] = _cx-_radius+i*dxyz;
	  m_y[i][j][k] = _cy-_radius+j*dxyz;
	  m_z[i][j][k] = _cz-_radius+k*dxyz;
	}

    mesh.push_back(m_x);
    mesh.push_back(m_y);
    mesh.push_back(m_z);

    size_t N3= N*N*N;
    metric_data.resize(DIMENSION*DIMENSION*N3);

    size_t indx=0;
    for(int i=0; i<DIMENSION; i++)
      {
	vector< array3Dref* > new_row;      
	for(int j=0; j<DIMENSION; j++)
	  {
	    new_row.push_back(new array3Dref(&metric_data[indx],boost::extents[N][N][N]));
	    indx+=N3;
	  }
	metric.push_back(new_row);
      }
    
    
  };


void analysis::update(evolution &ev)
{

  valarray<double> x = ev.get_position();
  valarray<double> p = ev.get_momentum();
  valarray<double> m = ev.get_mass();

  size_t np=ev.get_number_of_particles();
  size_t dim=ev.get_space_dim();

  for(index3D i=0; i<N; i++)
    for(index3D j=0; j<N; j++)
      for(index3D k=0; k<N; k++)
	{
	  
	  double xx = mesh[X][i][j][k];
	  double yy = mesh[Y][i][j][k];
	  double zz = mesh[Z][i][j][k];
	  double phi2=0;
	  double phi4=0;
	  for(int a=0; a < np; a++)
	    {

	      double xa = x[a*dim+X+1];
	      double ya = x[a*dim+Y+1];
	      double za = dim==DIMENSION ? x[a*dim+Z+1] : 0;

	      double pxa = p[a*dim+X+1];
	      double pya = p[a*dim+Y+1];
	      double pza = dim==DIMENSION ? p[a*dim+Z+1] : 0;

	      double pa2=pxa*pxa+pya*pya+pza*pza; 
	      
	      double ra = sqrt(gsl_sf_pow_int(xx-xa,2)+gsl_sf_pow_int(yy-ya,2)+gsl_sf_pow_int(zz-za,2));
	      phi2 += m[a]/ra;

	      phi4 += pa2/(m[a]*ra);
	      for(int b=0; b < np; b++)
		if(a!=b)
		  {
		    double xb = x[b*dim+X+1];
		    double yb = x[b*dim+Y+1];
		    double zb = dim==DIMENSION ? x[b*dim+Z+1] : 0;
		    double rab = sqrt(gsl_sf_pow_int(xb-xa,2)+gsl_sf_pow_int(yb-ya,2)+gsl_sf_pow_int(zb-za,2));
		    phi4 -= m[a]*m[b]/(ra*rab);
		  }

	    }
	  phi2 *= 4;
	  phi4 *= 2; 
	  
	  double psi_PN4=gsl_sf_pow_int(1+phi2+phi4,4);


	  for(int I=0; I<DIMENSION; I++)
	    {
	      (*metric[I][I])[i][j][k] = psi_PN4; 
	    }
	  
	}  
  

}
