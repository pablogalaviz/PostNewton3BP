

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

  valarray<double> position = ev.get_position();
  valarray<double> momentum = ev.get_momentum();
  valarray<double> m = ev.get_mass();

  size_t np=ev.get_number_of_particles();
  size_t dim=ev.get_space_dim();

  double x[np][DIMENSION];
  double p[np][DIMENSION];
  double rv[np][DIMENSION];
  double n[np][DIMENSION];
  double rrv[np][np][DIMENSION];
  double nnv[np][np][DIMENSION];

  for(int a=0; a<np; a++)
    for(int d=0; d<DIMENSION; d++)
    {
      x[a][d]=0;
      p[a][d]=0;
      rv[a][d]=0;
      n[a][d]=0;
      for(int b=0; b<np; b++)
	{
	  rrv[a][b][d]=0;
	  nnv[a][b][d]=0;
	}
    }

  
  double p2[np];
  double r[np];
  double n_p[np];
  double rr[np][np];

  for( int a =0; a < np; a++)
    {
      p2[a]=0;
    for(int d=0; d < dim; d++)
      {
	x[a][d]=position[a*dim+d+1];
	p[a][d]=momentum[a*dim+d+1];
	p2[a]+=p[a][d];
      }

    for( int b =0; b < np; b++)
      {
	rr[a][b] =0;
	for(int d=0; d < dim; d++)
	  {
	  rrv[a][b][d]+=x[b][d]-x[a][d];
	  rr[a][b]+=gsl_sf_pow_int(rrv[a][b][d],2);
	  }
	rr[a][b] = sqrt(rr[a][b]);
      }

    for( int b =0; b < np; b++)
      for(int d=0; d < dim; d++)
	nnv[a][b][d]=rrv[a][b][d]/rr[a][b];
    
    }
    
  for(index3D i=0; i<N; i++)
    for(index3D j=0; j<N; j++)
      for(index3D k=0; k<N; k++)
	{
	  
	  double phi2=0;
	  double phi4=0;
	  double hTT4_t1=0;
	  for(int a=0; a < np; a++)
	    {	      
	      r[a]=0;
	      for(int d=0; d < dim; d++)
		{
		  rv[a][d]=mesh[d][i][j][k]-x[a][d];
		  r[a] += gsl_sf_pow_int(rv[a][d],2);
		}
	      r[a] = sqrt(r[a]);

	      n_p[a]=0;
	      for(int d=0; d < dim; d++)
		{
		  n[a][d] = rv[a][d]/r[a];
		  n_p[a] += n[a][d]*p[a][d]; 
		}
	      	      
	      phi2 += m[a]/r[a];
	      phi4 += p2[a]/(m[a]*r[a]);
	      hTT4_t1+=(p2[a]-5*n_p[a]*n_p[a])/(m[a]*r[a]);
	      
	      for(int b=0; b < np; b++)
		if(a!=b)
		  phi4 -= m[a]*m[b]/(r[a]*rr[a][b]);
		  
	    }
	  phi2 *= 4;
	  phi4 *= 2; 
	  hTT4_t1*=0.25;
	  
	  double psi_PN4=gsl_sf_pow_int(1+phi2+phi4,4);

	  for(int I=0; I<DIMENSION; I++)
	    for(int J=0; J<DIMENSION; J++)
	      {
		
		(*metric[I][J])[i][j][k]=0;
		if(I==J)
		  (*metric[I][I])[i][j][k] = psi_PN4+hTT4_t1;

		for(int a=0; a< np; a++)
		  {
		  (*metric[I][J])[i][j][k]+=2*p[a][I]*p[a][J] + (3*n_p[a]*n_p[a]-5*p2[a])*n[a][I]*n[a][J];
		  (*metric[I][J])[i][j][k]+=12*n_p[a]*(n[a][I]*p[a][J]+n[a][J]*p[a][I]);

		  for(int b=0; b<np; b++)
		    if(a!=b)
		    {
		      double sab=r[a]+r[b]+rr[a][b];
		      (*metric[I][J])[i][j][k] += m[a]*m[b]*0.125*(-32*(1./rr[a][b]+1./sab)*nnv[a][b][I]*nnv[a][b][J]/sab + 2*((r[a]+r[b])/gsl_sf_pow_int(rr[a][b],3)+12./(sab*sab) )*n[a][I]*n[a][J] );
		    }
		  }
	      }
	  
	}  
  

}
