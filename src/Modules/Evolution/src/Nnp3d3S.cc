#include "evolution.h"



double evolution::rhs_Nnp3d3_FS(double t, const double y[], double par[3], int i)
{
  
 
  switch ( i ) {

  case 0 : 
    return(y[3]/par[0]);
    break;

  case 1 : 
    return(y[4]/par[0]);
    break;

  case 2 : 
    return(y[5]/par[0]);
    break;

  case 3 : 
    return((-((par[0]*par[1]*(y[0] - y[9]))/pow(pow(y[0] - y[9],2) + pow(y[1] - y[10],2) + pow(y[2] - y[11],2),1.5)) + (par[0]*par[1]*(-y[0] + y[9]))/pow(pow(-y[0] + y[9],2) + pow(-y[1] + y[10],2) + pow(-y[2] + y[11],2),1.5) - (par[0]*par[2]*(y[0] - y[18]))/pow(pow(y[0] - y[18],2) + pow(y[1] - y[19],2) + pow(y[2] - y[20],2),1.5) + (par[0]*par[2]*(-y[0] + y[18]))/pow(pow(-y[0] + y[18],2) + pow(-y[1] + y[19],2) + pow(-y[2] + y[20],2),1.5))/2.);
    break;

  case 4 : 
    return((-((par[0]*par[1]*(y[1] - y[10]))/pow(pow(y[0] - y[9],2) + pow(y[1] - y[10],2) + pow(y[2] - y[11],2),1.5)) + (par[0]*par[1]*(-y[1] + y[10]))/pow(pow(-y[0] + y[9],2) + pow(-y[1] + y[10],2) + pow(-y[2] + y[11],2),1.5) - (par[0]*par[2]*(y[1] - y[19]))/pow(pow(y[0] - y[18],2) + pow(y[1] - y[19],2) + pow(y[2] - y[20],2),1.5) + (par[0]*par[2]*(-y[1] + y[19]))/pow(pow(-y[0] + y[18],2) + pow(-y[1] + y[19],2) + pow(-y[2] + y[20],2),1.5))/2.);
    break;

  case 5 : 
    return((-((par[0]*par[1]*(y[2] - y[11]))/pow(pow(y[0] - y[9],2) + pow(y[1] - y[10],2) + pow(y[2] - y[11],2),1.5)) + (par[0]*par[1]*(-y[2] + y[11]))/pow(pow(-y[0] + y[9],2) + pow(-y[1] + y[10],2) + pow(-y[2] + y[11],2),1.5) - (par[0]*par[2]*(y[2] - y[20]))/pow(pow(y[0] - y[18],2) + pow(y[1] - y[19],2) + pow(y[2] - y[20],2),1.5) + (par[0]*par[2]*(-y[2] + y[20]))/pow(pow(-y[0] + y[18],2) + pow(-y[1] + y[19],2) + pow(-y[2] + y[20],2),1.5))/2.);
    break;

  case 6 : 
    return(0);
    break;

  case 7 : 
    return(0);
    break;

  case 8 : 
    return(0);
    break;

  case 9 : 
    return(y[12]/par[1]);
    break;

  case 10 : 
    return(y[13]/par[1]);
    break;

  case 11 : 
    return(y[14]/par[1]);
    break;

  case 12 : 
    return(((par[0]*par[1]*(y[0] - y[9]))/pow(pow(y[0] - y[9],2) + pow(y[1] - y[10],2) + pow(y[2] - y[11],2),1.5) - (par[0]*par[1]*(-y[0] + y[9]))/pow(pow(-y[0] + y[9],2) + pow(-y[1] + y[10],2) + pow(-y[2] + y[11],2),1.5) - (par[1]*par[2]*(y[9] - y[18]))/pow(pow(y[9] - y[18],2) + pow(y[10] - y[19],2) + pow(y[11] - y[20],2),1.5) + (par[1]*par[2]*(-y[9] + y[18]))/pow(pow(-y[9] + y[18],2) + pow(-y[10] + y[19],2) + pow(-y[11] + y[20],2),1.5))/2.);
    break;

  case 13 : 
    return(((par[0]*par[1]*(y[1] - y[10]))/pow(pow(y[0] - y[9],2) + pow(y[1] - y[10],2) + pow(y[2] - y[11],2),1.5) - (par[0]*par[1]*(-y[1] + y[10]))/pow(pow(-y[0] + y[9],2) + pow(-y[1] + y[10],2) + pow(-y[2] + y[11],2),1.5) - (par[1]*par[2]*(y[10] - y[19]))/pow(pow(y[9] - y[18],2) + pow(y[10] - y[19],2) + pow(y[11] - y[20],2),1.5) + (par[1]*par[2]*(-y[10] + y[19]))/pow(pow(-y[9] + y[18],2) + pow(-y[10] + y[19],2) + pow(-y[11] + y[20],2),1.5))/2.);
    break;

  case 14 : 
    return(((par[0]*par[1]*(y[2] - y[11]))/pow(pow(y[0] - y[9],2) + pow(y[1] - y[10],2) + pow(y[2] - y[11],2),1.5) - (par[0]*par[1]*(-y[2] + y[11]))/pow(pow(-y[0] + y[9],2) + pow(-y[1] + y[10],2) + pow(-y[2] + y[11],2),1.5) - (par[1]*par[2]*(y[11] - y[20]))/pow(pow(y[9] - y[18],2) + pow(y[10] - y[19],2) + pow(y[11] - y[20],2),1.5) + (par[1]*par[2]*(-y[11] + y[20]))/pow(pow(-y[9] + y[18],2) + pow(-y[10] + y[19],2) + pow(-y[11] + y[20],2),1.5))/2.);
    break;

  case 15 : 
    return(0);
    break;

  case 16 : 
    return(0);
    break;

  case 17 : 
    return(0);
    break;

  case 18 : 
    return(y[21]/par[2]);
    break;

  case 19 : 
    return(y[22]/par[2]);
    break;

  case 20 : 
    return(y[23]/par[2]);
    break;

  case 21 : 
    return(((par[0]*par[2]*(y[0] - y[18]))/pow(pow(y[0] - y[18],2) + pow(y[1] - y[19],2) + pow(y[2] - y[20],2),1.5) + (par[1]*par[2]*(y[9] - y[18]))/pow(pow(y[9] - y[18],2) + pow(y[10] - y[19],2) + pow(y[11] - y[20],2),1.5) - (par[0]*par[2]*(-y[0] + y[18]))/pow(pow(-y[0] + y[18],2) + pow(-y[1] + y[19],2) + pow(-y[2] + y[20],2),1.5) - (par[1]*par[2]*(-y[9] + y[18]))/pow(pow(-y[9] + y[18],2) + pow(-y[10] + y[19],2) + pow(-y[11] + y[20],2),1.5))/2.);
    break;

  case 22 : 
    return(((par[0]*par[2]*(y[1] - y[19]))/pow(pow(y[0] - y[18],2) + pow(y[1] - y[19],2) + pow(y[2] - y[20],2),1.5) + (par[1]*par[2]*(y[10] - y[19]))/pow(pow(y[9] - y[18],2) + pow(y[10] - y[19],2) + pow(y[11] - y[20],2),1.5) - (par[0]*par[2]*(-y[1] + y[19]))/pow(pow(-y[0] + y[18],2) + pow(-y[1] + y[19],2) + pow(-y[2] + y[20],2),1.5) - (par[1]*par[2]*(-y[10] + y[19]))/pow(pow(-y[9] + y[18],2) + pow(-y[10] + y[19],2) + pow(-y[11] + y[20],2),1.5))/2.);
    break;

  case 23 : 
    return(((par[0]*par[2]*(y[2] - y[20]))/pow(pow(y[0] - y[18],2) + pow(y[1] - y[19],2) + pow(y[2] - y[20],2),1.5) + (par[1]*par[2]*(y[11] - y[20]))/pow(pow(y[9] - y[18],2) + pow(y[10] - y[19],2) + pow(y[11] - y[20],2),1.5) - (par[0]*par[2]*(-y[2] + y[20]))/pow(pow(-y[0] + y[18],2) + pow(-y[1] + y[19],2) + pow(-y[2] + y[20],2),1.5) - (par[1]*par[2]*(-y[11] + y[20]))/pow(pow(-y[9] + y[18],2) + pow(-y[10] + y[19],2) + pow(-y[11] + y[20],2),1.5))/2.);
    break;

  case 24 : 
    return(0);
    break;

  case 25 : 
    return(0);
    break;

  case 26 : 
    return(0);
    break;
}

 return 0; 
}


double evolution::jac_Nnp3d3_FaS(double t, const double y[], double par[3], int i, int j)
{

 switch ( j*NDvar+i ) {

  case 0 : 
    return(0);
    break;

  case 1 : 
    return(0);
    break;

  case 2 : 
    return(0);
    break;

  case 3 : 
    return(((3*par[0]*par[1]*pow(y[0] - y[9],2))/pow(pow(y[0] - y[9],2) + pow(y[1] - y[10],2) + pow(y[2] - y[11],2),2.5) - (par[0]*par[1])/pow(pow(y[0] - y[9],2) + pow(y[1] - y[10],2) + pow(y[2] - y[11],2),1.5) + (3*par[0]*par[1]*pow(-y[0] + y[9],2))/pow(pow(-y[0] + y[9],2) + pow(-y[1] + y[10],2) + pow(-y[2] + y[11],2),2.5) - (par[0]*par[1])/pow(pow(-y[0] + y[9],2) + pow(-y[1] + y[10],2) + pow(-y[2] + y[11],2),1.5) + (3*par[0]*par[2]*pow(y[0] - y[18],2))/pow(pow(y[0] - y[18],2) + pow(y[1] - y[19],2) + pow(y[2] - y[20],2),2.5) - (par[0]*par[2])/pow(pow(y[0] - y[18],2) + pow(y[1] - y[19],2) + pow(y[2] - y[20],2),1.5) + (3*par[0]*par[2]*pow(-y[0] + y[18],2))/pow(pow(-y[0] + y[18],2) + pow(-y[1] + y[19],2) + pow(-y[2] + y[20],2),2.5) - (par[0]*par[2])/pow(pow(-y[0] + y[18],2) + pow(-y[1] + y[19],2) + pow(-y[2] + y[20],2),1.5))/2.);
    break;

  case 4 : 
    return(((3*par[0]*par[1]*(y[0] - y[9])*(y[1] - y[10]))/pow(pow(y[0] - y[9],2) + pow(y[1] - y[10],2) + pow(y[2] - y[11],2),2.5) + (3*par[0]*par[1]*(-y[0] + y[9])*(-y[1] + y[10]))/pow(pow(-y[0] + y[9],2) + pow(-y[1] + y[10],2) + pow(-y[2] + y[11],2),2.5) + (3*par[0]*par[2]*(y[0] - y[18])*(y[1] - y[19]))/pow(pow(y[0] - y[18],2) + pow(y[1] - y[19],2) + pow(y[2] - y[20],2),2.5) + (3*par[0]*par[2]*(-y[0] + y[18])*(-y[1] + y[19]))/pow(pow(-y[0] + y[18],2) + pow(-y[1] + y[19],2) + pow(-y[2] + y[20],2),2.5))/2.);
    break;

  case 5 : 
    return(((3*par[0]*par[1]*(y[0] - y[9])*(y[2] - y[11]))/pow(pow(y[0] - y[9],2) + pow(y[1] - y[10],2) + pow(y[2] - y[11],2),2.5) + (3*par[0]*par[1]*(-y[0] + y[9])*(-y[2] + y[11]))/pow(pow(-y[0] + y[9],2) + pow(-y[1] + y[10],2) + pow(-y[2] + y[11],2),2.5) + (3*par[0]*par[2]*(y[0] - y[18])*(y[2] - y[20]))/pow(pow(y[0] - y[18],2) + pow(y[1] - y[19],2) + pow(y[2] - y[20],2),2.5) + (3*par[0]*par[2]*(-y[0] + y[18])*(-y[2] + y[20]))/pow(pow(-y[0] + y[18],2) + pow(-y[1] + y[19],2) + pow(-y[2] + y[20],2),2.5))/2.);
    break;

  case 6 : 
    return(0);
    break;

  case 7 : 
    return(0);
    break;

  case 8 : 
    return(0);
    break;

  case 9 : 
    return(0);
    break;

  case 10 : 
    return(0);
    break;

  case 11 : 
    return(0);
    break;

  case 12 : 
    return(((-3*par[0]*par[1]*pow(y[0] - y[9],2))/pow(pow(y[0] - y[9],2) + pow(y[1] - y[10],2) + pow(y[2] - y[11],2),2.5) + (par[0]*par[1])/pow(pow(y[0] - y[9],2) + pow(y[1] - y[10],2) + pow(y[2] - y[11],2),1.5) - (3*par[0]*par[1]*pow(-y[0] + y[9],2))/pow(pow(-y[0] + y[9],2) + pow(-y[1] + y[10],2) + pow(-y[2] + y[11],2),2.5) + (par[0]*par[1])/pow(pow(-y[0] + y[9],2) + pow(-y[1] + y[10],2) + pow(-y[2] + y[11],2),1.5))/2.);
    break;

  case 13 : 
    return(((-3*par[0]*par[1]*(y[0] - y[9])*(y[1] - y[10]))/pow(pow(y[0] - y[9],2) + pow(y[1] - y[10],2) + pow(y[2] - y[11],2),2.5) - (3*par[0]*par[1]*(-y[0] + y[9])*(-y[1] + y[10]))/pow(pow(-y[0] + y[9],2) + pow(-y[1] + y[10],2) + pow(-y[2] + y[11],2),2.5))/2.);
    break;

  case 14 : 
    return(((-3*par[0]*par[1]*(y[0] - y[9])*(y[2] - y[11]))/pow(pow(y[0] - y[9],2) + pow(y[1] - y[10],2) + pow(y[2] - y[11],2),2.5) - (3*par[0]*par[1]*(-y[0] + y[9])*(-y[2] + y[11]))/pow(pow(-y[0] + y[9],2) + pow(-y[1] + y[10],2) + pow(-y[2] + y[11],2),2.5))/2.);
    break;

  case 15 : 
    return(0);
    break;

  case 16 : 
    return(0);
    break;

  case 17 : 
    return(0);
    break;

  case 18 : 
    return(0);
    break;

  case 19 : 
    return(0);
    break;

  case 20 : 
    return(0);
    break;

  case 21 : 
    return(((-3*par[0]*par[2]*pow(y[0] - y[18],2))/pow(pow(y[0] - y[18],2) + pow(y[1] - y[19],2) + pow(y[2] - y[20],2),2.5) + (par[0]*par[2])/pow(pow(y[0] - y[18],2) + pow(y[1] - y[19],2) + pow(y[2] - y[20],2),1.5) - (3*par[0]*par[2]*pow(-y[0] + y[18],2))/pow(pow(-y[0] + y[18],2) + pow(-y[1] + y[19],2) + pow(-y[2] + y[20],2),2.5) + (par[0]*par[2])/pow(pow(-y[0] + y[18],2) + pow(-y[1] + y[19],2) + pow(-y[2] + y[20],2),1.5))/2.);
    break;

  case 22 : 
    return(((-3*par[0]*par[2]*(y[0] - y[18])*(y[1] - y[19]))/pow(pow(y[0] - y[18],2) + pow(y[1] - y[19],2) + pow(y[2] - y[20],2),2.5) - (3*par[0]*par[2]*(-y[0] + y[18])*(-y[1] + y[19]))/pow(pow(-y[0] + y[18],2) + pow(-y[1] + y[19],2) + pow(-y[2] + y[20],2),2.5))/2.);
    break;

  case 23 : 
    return(((-3*par[0]*par[2]*(y[0] - y[18])*(y[2] - y[20]))/pow(pow(y[0] - y[18],2) + pow(y[1] - y[19],2) + pow(y[2] - y[20],2),2.5) - (3*par[0]*par[2]*(-y[0] + y[18])*(-y[2] + y[20]))/pow(pow(-y[0] + y[18],2) + pow(-y[1] + y[19],2) + pow(-y[2] + y[20],2),2.5))/2.);
    break;

  case 24 : 
    return(0);
    break;

  case 25 : 
    return(0);
    break;

  case 26 : 
    return(0);
    break;

  case 27 : 
    return(0);
    break;

  case 28 : 
    return(0);
    break;

  case 29 : 
    return(0);
    break;

  case 30 : 
    return(((3*par[0]*par[1]*(y[0] - y[9])*(y[1] - y[10]))/pow(pow(y[0] - y[9],2) + pow(y[1] - y[10],2) + pow(y[2] - y[11],2),2.5) + (3*par[0]*par[1]*(-y[0] + y[9])*(-y[1] + y[10]))/pow(pow(-y[0] + y[9],2) + pow(-y[1] + y[10],2) + pow(-y[2] + y[11],2),2.5) + (3*par[0]*par[2]*(y[0] - y[18])*(y[1] - y[19]))/pow(pow(y[0] - y[18],2) + pow(y[1] - y[19],2) + pow(y[2] - y[20],2),2.5) + (3*par[0]*par[2]*(-y[0] + y[18])*(-y[1] + y[19]))/pow(pow(-y[0] + y[18],2) + pow(-y[1] + y[19],2) + pow(-y[2] + y[20],2),2.5))/2.);
    break;

  case 31 : 
    return(((3*par[0]*par[1]*pow(y[1] - y[10],2))/pow(pow(y[0] - y[9],2) + pow(y[1] - y[10],2) + pow(y[2] - y[11],2),2.5) - (par[0]*par[1])/pow(pow(y[0] - y[9],2) + pow(y[1] - y[10],2) + pow(y[2] - y[11],2),1.5) + (3*par[0]*par[1]*pow(-y[1] + y[10],2))/pow(pow(-y[0] + y[9],2) + pow(-y[1] + y[10],2) + pow(-y[2] + y[11],2),2.5) - (par[0]*par[1])/pow(pow(-y[0] + y[9],2) + pow(-y[1] + y[10],2) + pow(-y[2] + y[11],2),1.5) + (3*par[0]*par[2]*pow(y[1] - y[19],2))/pow(pow(y[0] - y[18],2) + pow(y[1] - y[19],2) + pow(y[2] - y[20],2),2.5) - (par[0]*par[2])/pow(pow(y[0] - y[18],2) + pow(y[1] - y[19],2) + pow(y[2] - y[20],2),1.5) + (3*par[0]*par[2]*pow(-y[1] + y[19],2))/pow(pow(-y[0] + y[18],2) + pow(-y[1] + y[19],2) + pow(-y[2] + y[20],2),2.5) - (par[0]*par[2])/pow(pow(-y[0] + y[18],2) + pow(-y[1] + y[19],2) + pow(-y[2] + y[20],2),1.5))/2.);
    break;

  case 32 : 
    return(((3*par[0]*par[1]*(y[1] - y[10])*(y[2] - y[11]))/pow(pow(y[0] - y[9],2) + pow(y[1] - y[10],2) + pow(y[2] - y[11],2),2.5) + (3*par[0]*par[1]*(-y[1] + y[10])*(-y[2] + y[11]))/pow(pow(-y[0] + y[9],2) + pow(-y[1] + y[10],2) + pow(-y[2] + y[11],2),2.5) + (3*par[0]*par[2]*(y[1] - y[19])*(y[2] - y[20]))/pow(pow(y[0] - y[18],2) + pow(y[1] - y[19],2) + pow(y[2] - y[20],2),2.5) + (3*par[0]*par[2]*(-y[1] + y[19])*(-y[2] + y[20]))/pow(pow(-y[0] + y[18],2) + pow(-y[1] + y[19],2) + pow(-y[2] + y[20],2),2.5))/2.);
    break;

  case 33 : 
    return(0);
    break;

  case 34 : 
    return(0);
    break;

  case 35 : 
    return(0);
    break;

  case 36 : 
    return(0);
    break;

  case 37 : 
    return(0);
    break;

  case 38 : 
    return(0);
    break;

  case 39 : 
    return(((-3*par[0]*par[1]*(y[0] - y[9])*(y[1] - y[10]))/pow(pow(y[0] - y[9],2) + pow(y[1] - y[10],2) + pow(y[2] - y[11],2),2.5) - (3*par[0]*par[1]*(-y[0] + y[9])*(-y[1] + y[10]))/pow(pow(-y[0] + y[9],2) + pow(-y[1] + y[10],2) + pow(-y[2] + y[11],2),2.5))/2.);
    break;

  case 40 : 
    return(((-3*par[0]*par[1]*pow(y[1] - y[10],2))/pow(pow(y[0] - y[9],2) + pow(y[1] - y[10],2) + pow(y[2] - y[11],2),2.5) + (par[0]*par[1])/pow(pow(y[0] - y[9],2) + pow(y[1] - y[10],2) + pow(y[2] - y[11],2),1.5) - (3*par[0]*par[1]*pow(-y[1] + y[10],2))/pow(pow(-y[0] + y[9],2) + pow(-y[1] + y[10],2) + pow(-y[2] + y[11],2),2.5) + (par[0]*par[1])/pow(pow(-y[0] + y[9],2) + pow(-y[1] + y[10],2) + pow(-y[2] + y[11],2),1.5))/2.);
    break;

  case 41 : 
    return(((-3*par[0]*par[1]*(y[1] - y[10])*(y[2] - y[11]))/pow(pow(y[0] - y[9],2) + pow(y[1] - y[10],2) + pow(y[2] - y[11],2),2.5) - (3*par[0]*par[1]*(-y[1] + y[10])*(-y[2] + y[11]))/pow(pow(-y[0] + y[9],2) + pow(-y[1] + y[10],2) + pow(-y[2] + y[11],2),2.5))/2.);
    break;

  case 42 : 
    return(0);
    break;

  case 43 : 
    return(0);
    break;

  case 44 : 
    return(0);
    break;

  case 45 : 
    return(0);
    break;

  case 46 : 
    return(0);
    break;

  case 47 : 
    return(0);
    break;

  case 48 : 
    return(((-3*par[0]*par[2]*(y[0] - y[18])*(y[1] - y[19]))/pow(pow(y[0] - y[18],2) + pow(y[1] - y[19],2) + pow(y[2] - y[20],2),2.5) - (3*par[0]*par[2]*(-y[0] + y[18])*(-y[1] + y[19]))/pow(pow(-y[0] + y[18],2) + pow(-y[1] + y[19],2) + pow(-y[2] + y[20],2),2.5))/2.);
    break;

  case 49 : 
    return(((-3*par[0]*par[2]*pow(y[1] - y[19],2))/pow(pow(y[0] - y[18],2) + pow(y[1] - y[19],2) + pow(y[2] - y[20],2),2.5) + (par[0]*par[2])/pow(pow(y[0] - y[18],2) + pow(y[1] - y[19],2) + pow(y[2] - y[20],2),1.5) - (3*par[0]*par[2]*pow(-y[1] + y[19],2))/pow(pow(-y[0] + y[18],2) + pow(-y[1] + y[19],2) + pow(-y[2] + y[20],2),2.5) + (par[0]*par[2])/pow(pow(-y[0] + y[18],2) + pow(-y[1] + y[19],2) + pow(-y[2] + y[20],2),1.5))/2.);
    break;

  case 50 : 
    return(((-3*par[0]*par[2]*(y[1] - y[19])*(y[2] - y[20]))/pow(pow(y[0] - y[18],2) + pow(y[1] - y[19],2) + pow(y[2] - y[20],2),2.5) - (3*par[0]*par[2]*(-y[1] + y[19])*(-y[2] + y[20]))/pow(pow(-y[0] + y[18],2) + pow(-y[1] + y[19],2) + pow(-y[2] + y[20],2),2.5))/2.);
    break;

  case 51 : 
    return(0);
    break;

  case 52 : 
    return(0);
    break;

  case 53 : 
    return(0);
    break;

  case 54 : 
    return(0);
    break;

  case 55 : 
    return(0);
    break;

  case 56 : 
    return(0);
    break;

  case 57 : 
    return(((3*par[0]*par[1]*(y[0] - y[9])*(y[2] - y[11]))/pow(pow(y[0] - y[9],2) + pow(y[1] - y[10],2) + pow(y[2] - y[11],2),2.5) + (3*par[0]*par[1]*(-y[0] + y[9])*(-y[2] + y[11]))/pow(pow(-y[0] + y[9],2) + pow(-y[1] + y[10],2) + pow(-y[2] + y[11],2),2.5) + (3*par[0]*par[2]*(y[0] - y[18])*(y[2] - y[20]))/pow(pow(y[0] - y[18],2) + pow(y[1] - y[19],2) + pow(y[2] - y[20],2),2.5) + (3*par[0]*par[2]*(-y[0] + y[18])*(-y[2] + y[20]))/pow(pow(-y[0] + y[18],2) + pow(-y[1] + y[19],2) + pow(-y[2] + y[20],2),2.5))/2.);
    break;

  case 58 : 
    return(((3*par[0]*par[1]*(y[1] - y[10])*(y[2] - y[11]))/pow(pow(y[0] - y[9],2) + pow(y[1] - y[10],2) + pow(y[2] - y[11],2),2.5) + (3*par[0]*par[1]*(-y[1] + y[10])*(-y[2] + y[11]))/pow(pow(-y[0] + y[9],2) + pow(-y[1] + y[10],2) + pow(-y[2] + y[11],2),2.5) + (3*par[0]*par[2]*(y[1] - y[19])*(y[2] - y[20]))/pow(pow(y[0] - y[18],2) + pow(y[1] - y[19],2) + pow(y[2] - y[20],2),2.5) + (3*par[0]*par[2]*(-y[1] + y[19])*(-y[2] + y[20]))/pow(pow(-y[0] + y[18],2) + pow(-y[1] + y[19],2) + pow(-y[2] + y[20],2),2.5))/2.);
    break;

  case 59 : 
    return((-((par[0]*par[1])/pow(pow(y[0] - y[9],2) + pow(y[1] - y[10],2) + pow(y[2] - y[11],2),1.5)) + (3*par[0]*par[1]*pow(y[2] - y[11],2))/pow(pow(y[0] - y[9],2) + pow(y[1] - y[10],2) + pow(y[2] - y[11],2),2.5) + (3*par[0]*par[1]*pow(-y[2] + y[11],2))/pow(pow(-y[0] + y[9],2) + pow(-y[1] + y[10],2) + pow(-y[2] + y[11],2),2.5) - (par[0]*par[1])/pow(pow(-y[0] + y[9],2) + pow(-y[1] + y[10],2) + pow(-y[2] + y[11],2),1.5) - (par[0]*par[2])/pow(pow(y[0] - y[18],2) + pow(y[1] - y[19],2) + pow(y[2] - y[20],2),1.5) + (3*par[0]*par[2]*pow(y[2] - y[20],2))/pow(pow(y[0] - y[18],2) + pow(y[1] - y[19],2) + pow(y[2] - y[20],2),2.5) + (3*par[0]*par[2]*pow(-y[2] + y[20],2))/pow(pow(-y[0] + y[18],2) + pow(-y[1] + y[19],2) + pow(-y[2] + y[20],2),2.5) - (par[0]*par[2])/pow(pow(-y[0] + y[18],2) + pow(-y[1] + y[19],2) + pow(-y[2] + y[20],2),1.5))/2.);
    break;

  case 60 : 
    return(0);
    break;

  case 61 : 
    return(0);
    break;

  case 62 : 
    return(0);
    break;

  case 63 : 
    return(0);
    break;

  case 64 : 
    return(0);
    break;

  case 65 : 
    return(0);
    break;

  case 66 : 
    return(((-3*par[0]*par[1]*(y[0] - y[9])*(y[2] - y[11]))/pow(pow(y[0] - y[9],2) + pow(y[1] - y[10],2) + pow(y[2] - y[11],2),2.5) - (3*par[0]*par[1]*(-y[0] + y[9])*(-y[2] + y[11]))/pow(pow(-y[0] + y[9],2) + pow(-y[1] + y[10],2) + pow(-y[2] + y[11],2),2.5))/2.);
    break;

  case 67 : 
    return(((-3*par[0]*par[1]*(y[1] - y[10])*(y[2] - y[11]))/pow(pow(y[0] - y[9],2) + pow(y[1] - y[10],2) + pow(y[2] - y[11],2),2.5) - (3*par[0]*par[1]*(-y[1] + y[10])*(-y[2] + y[11]))/pow(pow(-y[0] + y[9],2) + pow(-y[1] + y[10],2) + pow(-y[2] + y[11],2),2.5))/2.);
    break;

  case 68 : 
    return(((par[0]*par[1])/pow(pow(y[0] - y[9],2) + pow(y[1] - y[10],2) + pow(y[2] - y[11],2),1.5) - (3*par[0]*par[1]*pow(y[2] - y[11],2))/pow(pow(y[0] - y[9],2) + pow(y[1] - y[10],2) + pow(y[2] - y[11],2),2.5) - (3*par[0]*par[1]*pow(-y[2] + y[11],2))/pow(pow(-y[0] + y[9],2) + pow(-y[1] + y[10],2) + pow(-y[2] + y[11],2),2.5) + (par[0]*par[1])/pow(pow(-y[0] + y[9],2) + pow(-y[1] + y[10],2) + pow(-y[2] + y[11],2),1.5))/2.);
    break;

  case 69 : 
    return(0);
    break;

  case 70 : 
    return(0);
    break;

  case 71 : 
    return(0);
    break;

  case 72 : 
    return(0);
    break;

  case 73 : 
    return(0);
    break;

  case 74 : 
    return(0);
    break;

  case 75 : 
    return(((-3*par[0]*par[2]*(y[0] - y[18])*(y[2] - y[20]))/pow(pow(y[0] - y[18],2) + pow(y[1] - y[19],2) + pow(y[2] - y[20],2),2.5) - (3*par[0]*par[2]*(-y[0] + y[18])*(-y[2] + y[20]))/pow(pow(-y[0] + y[18],2) + pow(-y[1] + y[19],2) + pow(-y[2] + y[20],2),2.5))/2.);
    break;

  case 76 : 
    return(((-3*par[0]*par[2]*(y[1] - y[19])*(y[2] - y[20]))/pow(pow(y[0] - y[18],2) + pow(y[1] - y[19],2) + pow(y[2] - y[20],2),2.5) - (3*par[0]*par[2]*(-y[1] + y[19])*(-y[2] + y[20]))/pow(pow(-y[0] + y[18],2) + pow(-y[1] + y[19],2) + pow(-y[2] + y[20],2),2.5))/2.);
    break;

  case 77 : 
    return(((par[0]*par[2])/pow(pow(y[0] - y[18],2) + pow(y[1] - y[19],2) + pow(y[2] - y[20],2),1.5) - (3*par[0]*par[2]*pow(y[2] - y[20],2))/pow(pow(y[0] - y[18],2) + pow(y[1] - y[19],2) + pow(y[2] - y[20],2),2.5) - (3*par[0]*par[2]*pow(-y[2] + y[20],2))/pow(pow(-y[0] + y[18],2) + pow(-y[1] + y[19],2) + pow(-y[2] + y[20],2),2.5) + (par[0]*par[2])/pow(pow(-y[0] + y[18],2) + pow(-y[1] + y[19],2) + pow(-y[2] + y[20],2),1.5))/2.);
    break;

  case 78 : 
    return(0);
    break;

  case 79 : 
    return(0);
    break;

  case 80 : 
    return(0);
    break;

  case 81 : 
    return(1/par[0]);
    break;

  case 82 : 
    return(0);
    break;

  case 83 : 
    return(0);
    break;

  case 84 : 
    return(0);
    break;

  case 85 : 
    return(0);
    break;

  case 86 : 
    return(0);
    break;

  case 87 : 
    return(0);
    break;

  case 88 : 
    return(0);
    break;

  case 89 : 
    return(0);
    break;

  case 90 : 
    return(0);
    break;

  case 91 : 
    return(0);
    break;

  case 92 : 
    return(0);
    break;

  case 93 : 
    return(0);
    break;

  case 94 : 
    return(0);
    break;

  case 95 : 
    return(0);
    break;

  case 96 : 
    return(0);
    break;

  case 97 : 
    return(0);
    break;

  case 98 : 
    return(0);
    break;

  case 99 : 
    return(0);
    break;

  case 100 : 
    return(0);
    break;

  case 101 : 
    return(0);
    break;

  case 102 : 
    return(0);
    break;

  case 103 : 
    return(0);
    break;

  case 104 : 
    return(0);
    break;

  case 105 : 
    return(0);
    break;

  case 106 : 
    return(0);
    break;

  case 107 : 
    return(0);
    break;

  case 108 : 
    return(0);
    break;

  case 109 : 
    return(1/par[0]);
    break;

  case 110 : 
    return(0);
    break;

  case 111 : 
    return(0);
    break;

  case 112 : 
    return(0);
    break;

  case 113 : 
    return(0);
    break;

  case 114 : 
    return(0);
    break;

  case 115 : 
    return(0);
    break;

  case 116 : 
    return(0);
    break;

  case 117 : 
    return(0);
    break;

  case 118 : 
    return(0);
    break;

  case 119 : 
    return(0);
    break;

  case 120 : 
    return(0);
    break;

  case 121 : 
    return(0);
    break;

  case 122 : 
    return(0);
    break;

  case 123 : 
    return(0);
    break;

  case 124 : 
    return(0);
    break;

  case 125 : 
    return(0);
    break;

  case 126 : 
    return(0);
    break;

  case 127 : 
    return(0);
    break;

  case 128 : 
    return(0);
    break;

  case 129 : 
    return(0);
    break;

  case 130 : 
    return(0);
    break;

  case 131 : 
    return(0);
    break;

  case 132 : 
    return(0);
    break;

  case 133 : 
    return(0);
    break;

  case 134 : 
    return(0);
    break;

  case 135 : 
    return(0);
    break;

  case 136 : 
    return(0);
    break;

  case 137 : 
    return(1/par[0]);
    break;

  case 138 : 
    return(0);
    break;

  case 139 : 
    return(0);
    break;

  case 140 : 
    return(0);
    break;

  case 141 : 
    return(0);
    break;

  case 142 : 
    return(0);
    break;

  case 143 : 
    return(0);
    break;

  case 144 : 
    return(0);
    break;

  case 145 : 
    return(0);
    break;

  case 146 : 
    return(0);
    break;

  case 147 : 
    return(0);
    break;

  case 148 : 
    return(0);
    break;

  case 149 : 
    return(0);
    break;

  case 150 : 
    return(0);
    break;

  case 151 : 
    return(0);
    break;

  case 152 : 
    return(0);
    break;

  case 153 : 
    return(0);
    break;

  case 154 : 
    return(0);
    break;

  case 155 : 
    return(0);
    break;

  case 156 : 
    return(0);
    break;

  case 157 : 
    return(0);
    break;

  case 158 : 
    return(0);
    break;

  case 159 : 
    return(0);
    break;

  case 160 : 
    return(0);
    break;

  case 161 : 
    return(0);
    break;

  case 162 : 
    return(0);
    break;

  case 163 : 
    return(0);
    break;

  case 164 : 
    return(0);
    break;

  case 165 : 
    return(0);
    break;

  case 166 : 
    return(0);
    break;

  case 167 : 
    return(0);
    break;

  case 168 : 
    return(0);
    break;

  case 169 : 
    return(0);
    break;

  case 170 : 
    return(0);
    break;

  case 171 : 
    return(0);
    break;

  case 172 : 
    return(0);
    break;

  case 173 : 
    return(0);
    break;

  case 174 : 
    return(0);
    break;

  case 175 : 
    return(0);
    break;

  case 176 : 
    return(0);
    break;

  case 177 : 
    return(0);
    break;

  case 178 : 
    return(0);
    break;

  case 179 : 
    return(0);
    break;

  case 180 : 
    return(0);
    break;

  case 181 : 
    return(0);
    break;

  case 182 : 
    return(0);
    break;

  case 183 : 
    return(0);
    break;

  case 184 : 
    return(0);
    break;

  case 185 : 
    return(0);
    break;

  case 186 : 
    return(0);
    break;

  case 187 : 
    return(0);
    break;

  case 188 : 
    return(0);
    break;

  case 189 : 
    return(0);
    break;

  case 190 : 
    return(0);
    break;

  case 191 : 
    return(0);
    break;

  case 192 : 
    return(0);
    break;

  case 193 : 
    return(0);
    break;

  case 194 : 
    return(0);
    break;

  case 195 : 
    return(0);
    break;

  case 196 : 
    return(0);
    break;

  case 197 : 
    return(0);
    break;

  case 198 : 
    return(0);
    break;

  case 199 : 
    return(0);
    break;

  case 200 : 
    return(0);
    break;

  case 201 : 
    return(0);
    break;

  case 202 : 
    return(0);
    break;

  case 203 : 
    return(0);
    break;

  case 204 : 
    return(0);
    break;

  case 205 : 
    return(0);
    break;

  case 206 : 
    return(0);
    break;

  case 207 : 
    return(0);
    break;

  case 208 : 
    return(0);
    break;

  case 209 : 
    return(0);
    break;

  case 210 : 
    return(0);
    break;

  case 211 : 
    return(0);
    break;

  case 212 : 
    return(0);
    break;

  case 213 : 
    return(0);
    break;

  case 214 : 
    return(0);
    break;

  case 215 : 
    return(0);
    break;

  case 216 : 
    return(0);
    break;

  case 217 : 
    return(0);
    break;

  case 218 : 
    return(0);
    break;

  case 219 : 
    return(0);
    break;

  case 220 : 
    return(0);
    break;

  case 221 : 
    return(0);
    break;

  case 222 : 
    return(0);
    break;

  case 223 : 
    return(0);
    break;

  case 224 : 
    return(0);
    break;

  case 225 : 
    return(0);
    break;

  case 226 : 
    return(0);
    break;

  case 227 : 
    return(0);
    break;

  case 228 : 
    return(0);
    break;

  case 229 : 
    return(0);
    break;

  case 230 : 
    return(0);
    break;

  case 231 : 
    return(0);
    break;

  case 232 : 
    return(0);
    break;

  case 233 : 
    return(0);
    break;

  case 234 : 
    return(0);
    break;

  case 235 : 
    return(0);
    break;

  case 236 : 
    return(0);
    break;

  case 237 : 
    return(0);
    break;

  case 238 : 
    return(0);
    break;

  case 239 : 
    return(0);
    break;

  case 240 : 
    return(0);
    break;

  case 241 : 
    return(0);
    break;

  case 242 : 
    return(0);
    break;

  case 243 : 
    return(0);
    break;

  case 244 : 
    return(0);
    break;

  case 245 : 
    return(0);
    break;

  case 246 : 
    return(((-3*par[0]*par[1]*pow(y[0] - y[9],2))/pow(pow(y[0] - y[9],2) + pow(y[1] - y[10],2) + pow(y[2] - y[11],2),2.5) + (par[0]*par[1])/pow(pow(y[0] - y[9],2) + pow(y[1] - y[10],2) + pow(y[2] - y[11],2),1.5) - (3*par[0]*par[1]*pow(-y[0] + y[9],2))/pow(pow(-y[0] + y[9],2) + pow(-y[1] + y[10],2) + pow(-y[2] + y[11],2),2.5) + (par[0]*par[1])/pow(pow(-y[0] + y[9],2) + pow(-y[1] + y[10],2) + pow(-y[2] + y[11],2),1.5))/2.);
    break;

  case 247 : 
    return(((-3*par[0]*par[1]*(y[0] - y[9])*(y[1] - y[10]))/pow(pow(y[0] - y[9],2) + pow(y[1] - y[10],2) + pow(y[2] - y[11],2),2.5) - (3*par[0]*par[1]*(-y[0] + y[9])*(-y[1] + y[10]))/pow(pow(-y[0] + y[9],2) + pow(-y[1] + y[10],2) + pow(-y[2] + y[11],2),2.5))/2.);
    break;

  case 248 : 
    return(((-3*par[0]*par[1]*(y[0] - y[9])*(y[2] - y[11]))/pow(pow(y[0] - y[9],2) + pow(y[1] - y[10],2) + pow(y[2] - y[11],2),2.5) - (3*par[0]*par[1]*(-y[0] + y[9])*(-y[2] + y[11]))/pow(pow(-y[0] + y[9],2) + pow(-y[1] + y[10],2) + pow(-y[2] + y[11],2),2.5))/2.);
    break;

  case 249 : 
    return(0);
    break;

  case 250 : 
    return(0);
    break;

  case 251 : 
    return(0);
    break;

  case 252 : 
    return(0);
    break;

  case 253 : 
    return(0);
    break;

  case 254 : 
    return(0);
    break;

  case 255 : 
    return(((3*par[0]*par[1]*pow(y[0] - y[9],2))/pow(pow(y[0] - y[9],2) + pow(y[1] - y[10],2) + pow(y[2] - y[11],2),2.5) - (par[0]*par[1])/pow(pow(y[0] - y[9],2) + pow(y[1] - y[10],2) + pow(y[2] - y[11],2),1.5) + (3*par[0]*par[1]*pow(-y[0] + y[9],2))/pow(pow(-y[0] + y[9],2) + pow(-y[1] + y[10],2) + pow(-y[2] + y[11],2),2.5) - (par[0]*par[1])/pow(pow(-y[0] + y[9],2) + pow(-y[1] + y[10],2) + pow(-y[2] + y[11],2),1.5) + (3*par[1]*par[2]*pow(y[9] - y[18],2))/pow(pow(y[9] - y[18],2) + pow(y[10] - y[19],2) + pow(y[11] - y[20],2),2.5) - (par[1]*par[2])/pow(pow(y[9] - y[18],2) + pow(y[10] - y[19],2) + pow(y[11] - y[20],2),1.5) + (3*par[1]*par[2]*pow(-y[9] + y[18],2))/pow(pow(-y[9] + y[18],2) + pow(-y[10] + y[19],2) + pow(-y[11] + y[20],2),2.5) - (par[1]*par[2])/pow(pow(-y[9] + y[18],2) + pow(-y[10] + y[19],2) + pow(-y[11] + y[20],2),1.5))/2.);
    break;

  case 256 : 
    return(((3*par[0]*par[1]*(y[0] - y[9])*(y[1] - y[10]))/pow(pow(y[0] - y[9],2) + pow(y[1] - y[10],2) + pow(y[2] - y[11],2),2.5) + (3*par[0]*par[1]*(-y[0] + y[9])*(-y[1] + y[10]))/pow(pow(-y[0] + y[9],2) + pow(-y[1] + y[10],2) + pow(-y[2] + y[11],2),2.5) + (3*par[1]*par[2]*(y[9] - y[18])*(y[10] - y[19]))/pow(pow(y[9] - y[18],2) + pow(y[10] - y[19],2) + pow(y[11] - y[20],2),2.5) + (3*par[1]*par[2]*(-y[9] + y[18])*(-y[10] + y[19]))/pow(pow(-y[9] + y[18],2) + pow(-y[10] + y[19],2) + pow(-y[11] + y[20],2),2.5))/2.);
    break;

  case 257 : 
    return(((3*par[0]*par[1]*(y[0] - y[9])*(y[2] - y[11]))/pow(pow(y[0] - y[9],2) + pow(y[1] - y[10],2) + pow(y[2] - y[11],2),2.5) + (3*par[0]*par[1]*(-y[0] + y[9])*(-y[2] + y[11]))/pow(pow(-y[0] + y[9],2) + pow(-y[1] + y[10],2) + pow(-y[2] + y[11],2),2.5) + (3*par[1]*par[2]*(y[9] - y[18])*(y[11] - y[20]))/pow(pow(y[9] - y[18],2) + pow(y[10] - y[19],2) + pow(y[11] - y[20],2),2.5) + (3*par[1]*par[2]*(-y[9] + y[18])*(-y[11] + y[20]))/pow(pow(-y[9] + y[18],2) + pow(-y[10] + y[19],2) + pow(-y[11] + y[20],2),2.5))/2.);
    break;

  case 258 : 
    return(0);
    break;

  case 259 : 
    return(0);
    break;

  case 260 : 
    return(0);
    break;

  case 261 : 
    return(0);
    break;

  case 262 : 
    return(0);
    break;

  case 263 : 
    return(0);
    break;

  case 264 : 
    return(((-3*par[1]*par[2]*pow(y[9] - y[18],2))/pow(pow(y[9] - y[18],2) + pow(y[10] - y[19],2) + pow(y[11] - y[20],2),2.5) + (par[1]*par[2])/pow(pow(y[9] - y[18],2) + pow(y[10] - y[19],2) + pow(y[11] - y[20],2),1.5) - (3*par[1]*par[2]*pow(-y[9] + y[18],2))/pow(pow(-y[9] + y[18],2) + pow(-y[10] + y[19],2) + pow(-y[11] + y[20],2),2.5) + (par[1]*par[2])/pow(pow(-y[9] + y[18],2) + pow(-y[10] + y[19],2) + pow(-y[11] + y[20],2),1.5))/2.);
    break;

  case 265 : 
    return(((-3*par[1]*par[2]*(y[9] - y[18])*(y[10] - y[19]))/pow(pow(y[9] - y[18],2) + pow(y[10] - y[19],2) + pow(y[11] - y[20],2),2.5) - (3*par[1]*par[2]*(-y[9] + y[18])*(-y[10] + y[19]))/pow(pow(-y[9] + y[18],2) + pow(-y[10] + y[19],2) + pow(-y[11] + y[20],2),2.5))/2.);
    break;

  case 266 : 
    return(((-3*par[1]*par[2]*(y[9] - y[18])*(y[11] - y[20]))/pow(pow(y[9] - y[18],2) + pow(y[10] - y[19],2) + pow(y[11] - y[20],2),2.5) - (3*par[1]*par[2]*(-y[9] + y[18])*(-y[11] + y[20]))/pow(pow(-y[9] + y[18],2) + pow(-y[10] + y[19],2) + pow(-y[11] + y[20],2),2.5))/2.);
    break;

  case 267 : 
    return(0);
    break;

  case 268 : 
    return(0);
    break;

  case 269 : 
    return(0);
    break;

  case 270 : 
    return(0);
    break;

  case 271 : 
    return(0);
    break;

  case 272 : 
    return(0);
    break;

  case 273 : 
    return(((-3*par[0]*par[1]*(y[0] - y[9])*(y[1] - y[10]))/pow(pow(y[0] - y[9],2) + pow(y[1] - y[10],2) + pow(y[2] - y[11],2),2.5) - (3*par[0]*par[1]*(-y[0] + y[9])*(-y[1] + y[10]))/pow(pow(-y[0] + y[9],2) + pow(-y[1] + y[10],2) + pow(-y[2] + y[11],2),2.5))/2.);
    break;

  case 274 : 
    return(((-3*par[0]*par[1]*pow(y[1] - y[10],2))/pow(pow(y[0] - y[9],2) + pow(y[1] - y[10],2) + pow(y[2] - y[11],2),2.5) + (par[0]*par[1])/pow(pow(y[0] - y[9],2) + pow(y[1] - y[10],2) + pow(y[2] - y[11],2),1.5) - (3*par[0]*par[1]*pow(-y[1] + y[10],2))/pow(pow(-y[0] + y[9],2) + pow(-y[1] + y[10],2) + pow(-y[2] + y[11],2),2.5) + (par[0]*par[1])/pow(pow(-y[0] + y[9],2) + pow(-y[1] + y[10],2) + pow(-y[2] + y[11],2),1.5))/2.);
    break;

  case 275 : 
    return(((-3*par[0]*par[1]*(y[1] - y[10])*(y[2] - y[11]))/pow(pow(y[0] - y[9],2) + pow(y[1] - y[10],2) + pow(y[2] - y[11],2),2.5) - (3*par[0]*par[1]*(-y[1] + y[10])*(-y[2] + y[11]))/pow(pow(-y[0] + y[9],2) + pow(-y[1] + y[10],2) + pow(-y[2] + y[11],2),2.5))/2.);
    break;

  case 276 : 
    return(0);
    break;

  case 277 : 
    return(0);
    break;

  case 278 : 
    return(0);
    break;

  case 279 : 
    return(0);
    break;

  case 280 : 
    return(0);
    break;

  case 281 : 
    return(0);
    break;

  case 282 : 
    return(((3*par[0]*par[1]*(y[0] - y[9])*(y[1] - y[10]))/pow(pow(y[0] - y[9],2) + pow(y[1] - y[10],2) + pow(y[2] - y[11],2),2.5) + (3*par[0]*par[1]*(-y[0] + y[9])*(-y[1] + y[10]))/pow(pow(-y[0] + y[9],2) + pow(-y[1] + y[10],2) + pow(-y[2] + y[11],2),2.5) + (3*par[1]*par[2]*(y[9] - y[18])*(y[10] - y[19]))/pow(pow(y[9] - y[18],2) + pow(y[10] - y[19],2) + pow(y[11] - y[20],2),2.5) + (3*par[1]*par[2]*(-y[9] + y[18])*(-y[10] + y[19]))/pow(pow(-y[9] + y[18],2) + pow(-y[10] + y[19],2) + pow(-y[11] + y[20],2),2.5))/2.);
    break;

  case 283 : 
    return(((3*par[0]*par[1]*pow(y[1] - y[10],2))/pow(pow(y[0] - y[9],2) + pow(y[1] - y[10],2) + pow(y[2] - y[11],2),2.5) - (par[0]*par[1])/pow(pow(y[0] - y[9],2) + pow(y[1] - y[10],2) + pow(y[2] - y[11],2),1.5) + (3*par[0]*par[1]*pow(-y[1] + y[10],2))/pow(pow(-y[0] + y[9],2) + pow(-y[1] + y[10],2) + pow(-y[2] + y[11],2),2.5) - (par[0]*par[1])/pow(pow(-y[0] + y[9],2) + pow(-y[1] + y[10],2) + pow(-y[2] + y[11],2),1.5) + (3*par[1]*par[2]*pow(y[10] - y[19],2))/pow(pow(y[9] - y[18],2) + pow(y[10] - y[19],2) + pow(y[11] - y[20],2),2.5) - (par[1]*par[2])/pow(pow(y[9] - y[18],2) + pow(y[10] - y[19],2) + pow(y[11] - y[20],2),1.5) + (3*par[1]*par[2]*pow(-y[10] + y[19],2))/pow(pow(-y[9] + y[18],2) + pow(-y[10] + y[19],2) + pow(-y[11] + y[20],2),2.5) - (par[1]*par[2])/pow(pow(-y[9] + y[18],2) + pow(-y[10] + y[19],2) + pow(-y[11] + y[20],2),1.5))/2.);
    break;

  case 284 : 
    return(((3*par[0]*par[1]*(y[1] - y[10])*(y[2] - y[11]))/pow(pow(y[0] - y[9],2) + pow(y[1] - y[10],2) + pow(y[2] - y[11],2),2.5) + (3*par[0]*par[1]*(-y[1] + y[10])*(-y[2] + y[11]))/pow(pow(-y[0] + y[9],2) + pow(-y[1] + y[10],2) + pow(-y[2] + y[11],2),2.5) + (3*par[1]*par[2]*(y[10] - y[19])*(y[11] - y[20]))/pow(pow(y[9] - y[18],2) + pow(y[10] - y[19],2) + pow(y[11] - y[20],2),2.5) + (3*par[1]*par[2]*(-y[10] + y[19])*(-y[11] + y[20]))/pow(pow(-y[9] + y[18],2) + pow(-y[10] + y[19],2) + pow(-y[11] + y[20],2),2.5))/2.);
    break;

  case 285 : 
    return(0);
    break;

  case 286 : 
    return(0);
    break;

  case 287 : 
    return(0);
    break;

  case 288 : 
    return(0);
    break;

  case 289 : 
    return(0);
    break;

  case 290 : 
    return(0);
    break;

  case 291 : 
    return(((-3*par[1]*par[2]*(y[9] - y[18])*(y[10] - y[19]))/pow(pow(y[9] - y[18],2) + pow(y[10] - y[19],2) + pow(y[11] - y[20],2),2.5) - (3*par[1]*par[2]*(-y[9] + y[18])*(-y[10] + y[19]))/pow(pow(-y[9] + y[18],2) + pow(-y[10] + y[19],2) + pow(-y[11] + y[20],2),2.5))/2.);
    break;

  case 292 : 
    return(((-3*par[1]*par[2]*pow(y[10] - y[19],2))/pow(pow(y[9] - y[18],2) + pow(y[10] - y[19],2) + pow(y[11] - y[20],2),2.5) + (par[1]*par[2])/pow(pow(y[9] - y[18],2) + pow(y[10] - y[19],2) + pow(y[11] - y[20],2),1.5) - (3*par[1]*par[2]*pow(-y[10] + y[19],2))/pow(pow(-y[9] + y[18],2) + pow(-y[10] + y[19],2) + pow(-y[11] + y[20],2),2.5) + (par[1]*par[2])/pow(pow(-y[9] + y[18],2) + pow(-y[10] + y[19],2) + pow(-y[11] + y[20],2),1.5))/2.);
    break;

  case 293 : 
    return(((-3*par[1]*par[2]*(y[10] - y[19])*(y[11] - y[20]))/pow(pow(y[9] - y[18],2) + pow(y[10] - y[19],2) + pow(y[11] - y[20],2),2.5) - (3*par[1]*par[2]*(-y[10] + y[19])*(-y[11] + y[20]))/pow(pow(-y[9] + y[18],2) + pow(-y[10] + y[19],2) + pow(-y[11] + y[20],2),2.5))/2.);
    break;

  case 294 : 
    return(0);
    break;

  case 295 : 
    return(0);
    break;

  case 296 : 
    return(0);
    break;

  case 297 : 
    return(0);
    break;

  case 298 : 
    return(0);
    break;

  case 299 : 
    return(0);
    break;

  case 300 : 
    return(((-3*par[0]*par[1]*(y[0] - y[9])*(y[2] - y[11]))/pow(pow(y[0] - y[9],2) + pow(y[1] - y[10],2) + pow(y[2] - y[11],2),2.5) - (3*par[0]*par[1]*(-y[0] + y[9])*(-y[2] + y[11]))/pow(pow(-y[0] + y[9],2) + pow(-y[1] + y[10],2) + pow(-y[2] + y[11],2),2.5))/2.);
    break;

  case 301 : 
    return(((-3*par[0]*par[1]*(y[1] - y[10])*(y[2] - y[11]))/pow(pow(y[0] - y[9],2) + pow(y[1] - y[10],2) + pow(y[2] - y[11],2),2.5) - (3*par[0]*par[1]*(-y[1] + y[10])*(-y[2] + y[11]))/pow(pow(-y[0] + y[9],2) + pow(-y[1] + y[10],2) + pow(-y[2] + y[11],2),2.5))/2.);
    break;

  case 302 : 
    return(((par[0]*par[1])/pow(pow(y[0] - y[9],2) + pow(y[1] - y[10],2) + pow(y[2] - y[11],2),1.5) - (3*par[0]*par[1]*pow(y[2] - y[11],2))/pow(pow(y[0] - y[9],2) + pow(y[1] - y[10],2) + pow(y[2] - y[11],2),2.5) - (3*par[0]*par[1]*pow(-y[2] + y[11],2))/pow(pow(-y[0] + y[9],2) + pow(-y[1] + y[10],2) + pow(-y[2] + y[11],2),2.5) + (par[0]*par[1])/pow(pow(-y[0] + y[9],2) + pow(-y[1] + y[10],2) + pow(-y[2] + y[11],2),1.5))/2.);
    break;

  case 303 : 
    return(0);
    break;

  case 304 : 
    return(0);
    break;

  case 305 : 
    return(0);
    break;

  case 306 : 
    return(0);
    break;

  case 307 : 
    return(0);
    break;

  case 308 : 
    return(0);
    break;

  case 309 : 
    return(((3*par[0]*par[1]*(y[0] - y[9])*(y[2] - y[11]))/pow(pow(y[0] - y[9],2) + pow(y[1] - y[10],2) + pow(y[2] - y[11],2),2.5) + (3*par[0]*par[1]*(-y[0] + y[9])*(-y[2] + y[11]))/pow(pow(-y[0] + y[9],2) + pow(-y[1] + y[10],2) + pow(-y[2] + y[11],2),2.5) + (3*par[1]*par[2]*(y[9] - y[18])*(y[11] - y[20]))/pow(pow(y[9] - y[18],2) + pow(y[10] - y[19],2) + pow(y[11] - y[20],2),2.5) + (3*par[1]*par[2]*(-y[9] + y[18])*(-y[11] + y[20]))/pow(pow(-y[9] + y[18],2) + pow(-y[10] + y[19],2) + pow(-y[11] + y[20],2),2.5))/2.);
    break;

  case 310 : 
    return(((3*par[0]*par[1]*(y[1] - y[10])*(y[2] - y[11]))/pow(pow(y[0] - y[9],2) + pow(y[1] - y[10],2) + pow(y[2] - y[11],2),2.5) + (3*par[0]*par[1]*(-y[1] + y[10])*(-y[2] + y[11]))/pow(pow(-y[0] + y[9],2) + pow(-y[1] + y[10],2) + pow(-y[2] + y[11],2),2.5) + (3*par[1]*par[2]*(y[10] - y[19])*(y[11] - y[20]))/pow(pow(y[9] - y[18],2) + pow(y[10] - y[19],2) + pow(y[11] - y[20],2),2.5) + (3*par[1]*par[2]*(-y[10] + y[19])*(-y[11] + y[20]))/pow(pow(-y[9] + y[18],2) + pow(-y[10] + y[19],2) + pow(-y[11] + y[20],2),2.5))/2.);
    break;

  case 311 : 
    return((-((par[0]*par[1])/pow(pow(y[0] - y[9],2) + pow(y[1] - y[10],2) + pow(y[2] - y[11],2),1.5)) + (3*par[0]*par[1]*pow(y[2] - y[11],2))/pow(pow(y[0] - y[9],2) + pow(y[1] - y[10],2) + pow(y[2] - y[11],2),2.5) + (3*par[0]*par[1]*pow(-y[2] + y[11],2))/pow(pow(-y[0] + y[9],2) + pow(-y[1] + y[10],2) + pow(-y[2] + y[11],2),2.5) - (par[0]*par[1])/pow(pow(-y[0] + y[9],2) + pow(-y[1] + y[10],2) + pow(-y[2] + y[11],2),1.5) - (par[1]*par[2])/pow(pow(y[9] - y[18],2) + pow(y[10] - y[19],2) + pow(y[11] - y[20],2),1.5) + (3*par[1]*par[2]*pow(y[11] - y[20],2))/pow(pow(y[9] - y[18],2) + pow(y[10] - y[19],2) + pow(y[11] - y[20],2),2.5) + (3*par[1]*par[2]*pow(-y[11] + y[20],2))/pow(pow(-y[9] + y[18],2) + pow(-y[10] + y[19],2) + pow(-y[11] + y[20],2),2.5) - (par[1]*par[2])/pow(pow(-y[9] + y[18],2) + pow(-y[10] + y[19],2) + pow(-y[11] + y[20],2),1.5))/2.);
    break;

  case 312 : 
    return(0);
    break;

  case 313 : 
    return(0);
    break;

  case 314 : 
    return(0);
    break;

  case 315 : 
    return(0);
    break;

  case 316 : 
    return(0);
    break;

  case 317 : 
    return(0);
    break;

  case 318 : 
    return(((-3*par[1]*par[2]*(y[9] - y[18])*(y[11] - y[20]))/pow(pow(y[9] - y[18],2) + pow(y[10] - y[19],2) + pow(y[11] - y[20],2),2.5) - (3*par[1]*par[2]*(-y[9] + y[18])*(-y[11] + y[20]))/pow(pow(-y[9] + y[18],2) + pow(-y[10] + y[19],2) + pow(-y[11] + y[20],2),2.5))/2.);
    break;

  case 319 : 
    return(((-3*par[1]*par[2]*(y[10] - y[19])*(y[11] - y[20]))/pow(pow(y[9] - y[18],2) + pow(y[10] - y[19],2) + pow(y[11] - y[20],2),2.5) - (3*par[1]*par[2]*(-y[10] + y[19])*(-y[11] + y[20]))/pow(pow(-y[9] + y[18],2) + pow(-y[10] + y[19],2) + pow(-y[11] + y[20],2),2.5))/2.);
    break;

  case 320 : 
    return(((par[1]*par[2])/pow(pow(y[9] - y[18],2) + pow(y[10] - y[19],2) + pow(y[11] - y[20],2),1.5) - (3*par[1]*par[2]*pow(y[11] - y[20],2))/pow(pow(y[9] - y[18],2) + pow(y[10] - y[19],2) + pow(y[11] - y[20],2),2.5) - (3*par[1]*par[2]*pow(-y[11] + y[20],2))/pow(pow(-y[9] + y[18],2) + pow(-y[10] + y[19],2) + pow(-y[11] + y[20],2),2.5) + (par[1]*par[2])/pow(pow(-y[9] + y[18],2) + pow(-y[10] + y[19],2) + pow(-y[11] + y[20],2),1.5))/2.);
    break;

  case 321 : 
    return(0);
    break;

  case 322 : 
    return(0);
    break;

  case 323 : 
    return(0);
    break;

  case 324 : 
    return(0);
    break;

  case 325 : 
    return(0);
    break;

  case 326 : 
    return(0);
    break;

  case 327 : 
    return(0);
    break;

  case 328 : 
    return(0);
    break;

  case 329 : 
    return(0);
    break;

  case 330 : 
    return(0);
    break;

  case 331 : 
    return(0);
    break;

  case 332 : 
    return(0);
    break;

  case 333 : 
    return(1/par[1]);
    break;

  case 334 : 
    return(0);
    break;

  case 335 : 
    return(0);
    break;

  case 336 : 
    return(0);
    break;

  case 337 : 
    return(0);
    break;

  case 338 : 
    return(0);
    break;

  case 339 : 
    return(0);
    break;

  case 340 : 
    return(0);
    break;

  case 341 : 
    return(0);
    break;

  case 342 : 
    return(0);
    break;

  case 343 : 
    return(0);
    break;

  case 344 : 
    return(0);
    break;

  case 345 : 
    return(0);
    break;

  case 346 : 
    return(0);
    break;

  case 347 : 
    return(0);
    break;

  case 348 : 
    return(0);
    break;

  case 349 : 
    return(0);
    break;

  case 350 : 
    return(0);
    break;

  case 351 : 
    return(0);
    break;

  case 352 : 
    return(0);
    break;

  case 353 : 
    return(0);
    break;

  case 354 : 
    return(0);
    break;

  case 355 : 
    return(0);
    break;

  case 356 : 
    return(0);
    break;

  case 357 : 
    return(0);
    break;

  case 358 : 
    return(0);
    break;

  case 359 : 
    return(0);
    break;

  case 360 : 
    return(0);
    break;

  case 361 : 
    return(1/par[1]);
    break;

  case 362 : 
    return(0);
    break;

  case 363 : 
    return(0);
    break;

  case 364 : 
    return(0);
    break;

  case 365 : 
    return(0);
    break;

  case 366 : 
    return(0);
    break;

  case 367 : 
    return(0);
    break;

  case 368 : 
    return(0);
    break;

  case 369 : 
    return(0);
    break;

  case 370 : 
    return(0);
    break;

  case 371 : 
    return(0);
    break;

  case 372 : 
    return(0);
    break;

  case 373 : 
    return(0);
    break;

  case 374 : 
    return(0);
    break;

  case 375 : 
    return(0);
    break;

  case 376 : 
    return(0);
    break;

  case 377 : 
    return(0);
    break;

  case 378 : 
    return(0);
    break;

  case 379 : 
    return(0);
    break;

  case 380 : 
    return(0);
    break;

  case 381 : 
    return(0);
    break;

  case 382 : 
    return(0);
    break;

  case 383 : 
    return(0);
    break;

  case 384 : 
    return(0);
    break;

  case 385 : 
    return(0);
    break;

  case 386 : 
    return(0);
    break;

  case 387 : 
    return(0);
    break;

  case 388 : 
    return(0);
    break;

  case 389 : 
    return(1/par[1]);
    break;

  case 390 : 
    return(0);
    break;

  case 391 : 
    return(0);
    break;

  case 392 : 
    return(0);
    break;

  case 393 : 
    return(0);
    break;

  case 394 : 
    return(0);
    break;

  case 395 : 
    return(0);
    break;

  case 396 : 
    return(0);
    break;

  case 397 : 
    return(0);
    break;

  case 398 : 
    return(0);
    break;

  case 399 : 
    return(0);
    break;

  case 400 : 
    return(0);
    break;

  case 401 : 
    return(0);
    break;

  case 402 : 
    return(0);
    break;

  case 403 : 
    return(0);
    break;

  case 404 : 
    return(0);
    break;

  case 405 : 
    return(0);
    break;

  case 406 : 
    return(0);
    break;

  case 407 : 
    return(0);
    break;

  case 408 : 
    return(0);
    break;

  case 409 : 
    return(0);
    break;

  case 410 : 
    return(0);
    break;

  case 411 : 
    return(0);
    break;

  case 412 : 
    return(0);
    break;

  case 413 : 
    return(0);
    break;

  case 414 : 
    return(0);
    break;

  case 415 : 
    return(0);
    break;

  case 416 : 
    return(0);
    break;

  case 417 : 
    return(0);
    break;

  case 418 : 
    return(0);
    break;

  case 419 : 
    return(0);
    break;

  case 420 : 
    return(0);
    break;

  case 421 : 
    return(0);
    break;

  case 422 : 
    return(0);
    break;

  case 423 : 
    return(0);
    break;

  case 424 : 
    return(0);
    break;

  case 425 : 
    return(0);
    break;

  case 426 : 
    return(0);
    break;

  case 427 : 
    return(0);
    break;

  case 428 : 
    return(0);
    break;

  case 429 : 
    return(0);
    break;

  case 430 : 
    return(0);
    break;

  case 431 : 
    return(0);
    break;

  case 432 : 
    return(0);
    break;

  case 433 : 
    return(0);
    break;

  case 434 : 
    return(0);
    break;

  case 435 : 
    return(0);
    break;

  case 436 : 
    return(0);
    break;

  case 437 : 
    return(0);
    break;

  case 438 : 
    return(0);
    break;

  case 439 : 
    return(0);
    break;

  case 440 : 
    return(0);
    break;

  case 441 : 
    return(0);
    break;

  case 442 : 
    return(0);
    break;

  case 443 : 
    return(0);
    break;

  case 444 : 
    return(0);
    break;

  case 445 : 
    return(0);
    break;

  case 446 : 
    return(0);
    break;

  case 447 : 
    return(0);
    break;

  case 448 : 
    return(0);
    break;

  case 449 : 
    return(0);
    break;

  case 450 : 
    return(0);
    break;

  case 451 : 
    return(0);
    break;

  case 452 : 
    return(0);
    break;

  case 453 : 
    return(0);
    break;

  case 454 : 
    return(0);
    break;

  case 455 : 
    return(0);
    break;

  case 456 : 
    return(0);
    break;

  case 457 : 
    return(0);
    break;

  case 458 : 
    return(0);
    break;

  case 459 : 
    return(0);
    break;

  case 460 : 
    return(0);
    break;

  case 461 : 
    return(0);
    break;

  case 462 : 
    return(0);
    break;

  case 463 : 
    return(0);
    break;

  case 464 : 
    return(0);
    break;

  case 465 : 
    return(0);
    break;

  case 466 : 
    return(0);
    break;

  case 467 : 
    return(0);
    break;

  case 468 : 
    return(0);
    break;

  case 469 : 
    return(0);
    break;

  case 470 : 
    return(0);
    break;

  case 471 : 
    return(0);
    break;

  case 472 : 
    return(0);
    break;

  case 473 : 
    return(0);
    break;

  case 474 : 
    return(0);
    break;

  case 475 : 
    return(0);
    break;

  case 476 : 
    return(0);
    break;

  case 477 : 
    return(0);
    break;

  case 478 : 
    return(0);
    break;

  case 479 : 
    return(0);
    break;

  case 480 : 
    return(0);
    break;

  case 481 : 
    return(0);
    break;

  case 482 : 
    return(0);
    break;

  case 483 : 
    return(0);
    break;

  case 484 : 
    return(0);
    break;

  case 485 : 
    return(0);
    break;

  case 486 : 
    return(0);
    break;

  case 487 : 
    return(0);
    break;

  case 488 : 
    return(0);
    break;

  case 489 : 
    return(((-3*par[0]*par[2]*pow(y[0] - y[18],2))/pow(pow(y[0] - y[18],2) + pow(y[1] - y[19],2) + pow(y[2] - y[20],2),2.5) + (par[0]*par[2])/pow(pow(y[0] - y[18],2) + pow(y[1] - y[19],2) + pow(y[2] - y[20],2),1.5) - (3*par[0]*par[2]*pow(-y[0] + y[18],2))/pow(pow(-y[0] + y[18],2) + pow(-y[1] + y[19],2) + pow(-y[2] + y[20],2),2.5) + (par[0]*par[2])/pow(pow(-y[0] + y[18],2) + pow(-y[1] + y[19],2) + pow(-y[2] + y[20],2),1.5))/2.);
    break;

  case 490 : 
    return(((-3*par[0]*par[2]*(y[0] - y[18])*(y[1] - y[19]))/pow(pow(y[0] - y[18],2) + pow(y[1] - y[19],2) + pow(y[2] - y[20],2),2.5) - (3*par[0]*par[2]*(-y[0] + y[18])*(-y[1] + y[19]))/pow(pow(-y[0] + y[18],2) + pow(-y[1] + y[19],2) + pow(-y[2] + y[20],2),2.5))/2.);
    break;

  case 491 : 
    return(((-3*par[0]*par[2]*(y[0] - y[18])*(y[2] - y[20]))/pow(pow(y[0] - y[18],2) + pow(y[1] - y[19],2) + pow(y[2] - y[20],2),2.5) - (3*par[0]*par[2]*(-y[0] + y[18])*(-y[2] + y[20]))/pow(pow(-y[0] + y[18],2) + pow(-y[1] + y[19],2) + pow(-y[2] + y[20],2),2.5))/2.);
    break;

  case 492 : 
    return(0);
    break;

  case 493 : 
    return(0);
    break;

  case 494 : 
    return(0);
    break;

  case 495 : 
    return(0);
    break;

  case 496 : 
    return(0);
    break;

  case 497 : 
    return(0);
    break;

  case 498 : 
    return(((-3*par[1]*par[2]*pow(y[9] - y[18],2))/pow(pow(y[9] - y[18],2) + pow(y[10] - y[19],2) + pow(y[11] - y[20],2),2.5) + (par[1]*par[2])/pow(pow(y[9] - y[18],2) + pow(y[10] - y[19],2) + pow(y[11] - y[20],2),1.5) - (3*par[1]*par[2]*pow(-y[9] + y[18],2))/pow(pow(-y[9] + y[18],2) + pow(-y[10] + y[19],2) + pow(-y[11] + y[20],2),2.5) + (par[1]*par[2])/pow(pow(-y[9] + y[18],2) + pow(-y[10] + y[19],2) + pow(-y[11] + y[20],2),1.5))/2.);
    break;

  case 499 : 
    return(((-3*par[1]*par[2]*(y[9] - y[18])*(y[10] - y[19]))/pow(pow(y[9] - y[18],2) + pow(y[10] - y[19],2) + pow(y[11] - y[20],2),2.5) - (3*par[1]*par[2]*(-y[9] + y[18])*(-y[10] + y[19]))/pow(pow(-y[9] + y[18],2) + pow(-y[10] + y[19],2) + pow(-y[11] + y[20],2),2.5))/2.);
    break;

  case 500 : 
    return(((-3*par[1]*par[2]*(y[9] - y[18])*(y[11] - y[20]))/pow(pow(y[9] - y[18],2) + pow(y[10] - y[19],2) + pow(y[11] - y[20],2),2.5) - (3*par[1]*par[2]*(-y[9] + y[18])*(-y[11] + y[20]))/pow(pow(-y[9] + y[18],2) + pow(-y[10] + y[19],2) + pow(-y[11] + y[20],2),2.5))/2.);
    break;

  case 501 : 
    return(0);
    break;

  case 502 : 
    return(0);
    break;

  case 503 : 
    return(0);
    break;

  case 504 : 
    return(0);
    break;

  case 505 : 
    return(0);
    break;

  case 506 : 
    return(0);
    break;

  case 507 : 
    return(((3*par[0]*par[2]*pow(y[0] - y[18],2))/pow(pow(y[0] - y[18],2) + pow(y[1] - y[19],2) + pow(y[2] - y[20],2),2.5) - (par[0]*par[2])/pow(pow(y[0] - y[18],2) + pow(y[1] - y[19],2) + pow(y[2] - y[20],2),1.5) + (3*par[1]*par[2]*pow(y[9] - y[18],2))/pow(pow(y[9] - y[18],2) + pow(y[10] - y[19],2) + pow(y[11] - y[20],2),2.5) - (par[1]*par[2])/pow(pow(y[9] - y[18],2) + pow(y[10] - y[19],2) + pow(y[11] - y[20],2),1.5) + (3*par[0]*par[2]*pow(-y[0] + y[18],2))/pow(pow(-y[0] + y[18],2) + pow(-y[1] + y[19],2) + pow(-y[2] + y[20],2),2.5) - (par[0]*par[2])/pow(pow(-y[0] + y[18],2) + pow(-y[1] + y[19],2) + pow(-y[2] + y[20],2),1.5) + (3*par[1]*par[2]*pow(-y[9] + y[18],2))/pow(pow(-y[9] + y[18],2) + pow(-y[10] + y[19],2) + pow(-y[11] + y[20],2),2.5) - (par[1]*par[2])/pow(pow(-y[9] + y[18],2) + pow(-y[10] + y[19],2) + pow(-y[11] + y[20],2),1.5))/2.);
    break;

  case 508 : 
    return(((3*par[0]*par[2]*(y[0] - y[18])*(y[1] - y[19]))/pow(pow(y[0] - y[18],2) + pow(y[1] - y[19],2) + pow(y[2] - y[20],2),2.5) + (3*par[1]*par[2]*(y[9] - y[18])*(y[10] - y[19]))/pow(pow(y[9] - y[18],2) + pow(y[10] - y[19],2) + pow(y[11] - y[20],2),2.5) + (3*par[0]*par[2]*(-y[0] + y[18])*(-y[1] + y[19]))/pow(pow(-y[0] + y[18],2) + pow(-y[1] + y[19],2) + pow(-y[2] + y[20],2),2.5) + (3*par[1]*par[2]*(-y[9] + y[18])*(-y[10] + y[19]))/pow(pow(-y[9] + y[18],2) + pow(-y[10] + y[19],2) + pow(-y[11] + y[20],2),2.5))/2.);
    break;

  case 509 : 
    return(((3*par[0]*par[2]*(y[0] - y[18])*(y[2] - y[20]))/pow(pow(y[0] - y[18],2) + pow(y[1] - y[19],2) + pow(y[2] - y[20],2),2.5) + (3*par[1]*par[2]*(y[9] - y[18])*(y[11] - y[20]))/pow(pow(y[9] - y[18],2) + pow(y[10] - y[19],2) + pow(y[11] - y[20],2),2.5) + (3*par[0]*par[2]*(-y[0] + y[18])*(-y[2] + y[20]))/pow(pow(-y[0] + y[18],2) + pow(-y[1] + y[19],2) + pow(-y[2] + y[20],2),2.5) + (3*par[1]*par[2]*(-y[9] + y[18])*(-y[11] + y[20]))/pow(pow(-y[9] + y[18],2) + pow(-y[10] + y[19],2) + pow(-y[11] + y[20],2),2.5))/2.);
    break;

  case 510 : 
    return(0);
    break;

  case 511 : 
    return(0);
    break;

  case 512 : 
    return(0);
    break;

  case 513 : 
    return(0);
    break;

  case 514 : 
    return(0);
    break;

  case 515 : 
    return(0);
    break;

  case 516 : 
    return(((-3*par[0]*par[2]*(y[0] - y[18])*(y[1] - y[19]))/pow(pow(y[0] - y[18],2) + pow(y[1] - y[19],2) + pow(y[2] - y[20],2),2.5) - (3*par[0]*par[2]*(-y[0] + y[18])*(-y[1] + y[19]))/pow(pow(-y[0] + y[18],2) + pow(-y[1] + y[19],2) + pow(-y[2] + y[20],2),2.5))/2.);
    break;

  case 517 : 
    return(((-3*par[0]*par[2]*pow(y[1] - y[19],2))/pow(pow(y[0] - y[18],2) + pow(y[1] - y[19],2) + pow(y[2] - y[20],2),2.5) + (par[0]*par[2])/pow(pow(y[0] - y[18],2) + pow(y[1] - y[19],2) + pow(y[2] - y[20],2),1.5) - (3*par[0]*par[2]*pow(-y[1] + y[19],2))/pow(pow(-y[0] + y[18],2) + pow(-y[1] + y[19],2) + pow(-y[2] + y[20],2),2.5) + (par[0]*par[2])/pow(pow(-y[0] + y[18],2) + pow(-y[1] + y[19],2) + pow(-y[2] + y[20],2),1.5))/2.);
    break;

  case 518 : 
    return(((-3*par[0]*par[2]*(y[1] - y[19])*(y[2] - y[20]))/pow(pow(y[0] - y[18],2) + pow(y[1] - y[19],2) + pow(y[2] - y[20],2),2.5) - (3*par[0]*par[2]*(-y[1] + y[19])*(-y[2] + y[20]))/pow(pow(-y[0] + y[18],2) + pow(-y[1] + y[19],2) + pow(-y[2] + y[20],2),2.5))/2.);
    break;

  case 519 : 
    return(0);
    break;

  case 520 : 
    return(0);
    break;

  case 521 : 
    return(0);
    break;

  case 522 : 
    return(0);
    break;

  case 523 : 
    return(0);
    break;

  case 524 : 
    return(0);
    break;

  case 525 : 
    return(((-3*par[1]*par[2]*(y[9] - y[18])*(y[10] - y[19]))/pow(pow(y[9] - y[18],2) + pow(y[10] - y[19],2) + pow(y[11] - y[20],2),2.5) - (3*par[1]*par[2]*(-y[9] + y[18])*(-y[10] + y[19]))/pow(pow(-y[9] + y[18],2) + pow(-y[10] + y[19],2) + pow(-y[11] + y[20],2),2.5))/2.);
    break;

  case 526 : 
    return(((-3*par[1]*par[2]*pow(y[10] - y[19],2))/pow(pow(y[9] - y[18],2) + pow(y[10] - y[19],2) + pow(y[11] - y[20],2),2.5) + (par[1]*par[2])/pow(pow(y[9] - y[18],2) + pow(y[10] - y[19],2) + pow(y[11] - y[20],2),1.5) - (3*par[1]*par[2]*pow(-y[10] + y[19],2))/pow(pow(-y[9] + y[18],2) + pow(-y[10] + y[19],2) + pow(-y[11] + y[20],2),2.5) + (par[1]*par[2])/pow(pow(-y[9] + y[18],2) + pow(-y[10] + y[19],2) + pow(-y[11] + y[20],2),1.5))/2.);
    break;

  case 527 : 
    return(((-3*par[1]*par[2]*(y[10] - y[19])*(y[11] - y[20]))/pow(pow(y[9] - y[18],2) + pow(y[10] - y[19],2) + pow(y[11] - y[20],2),2.5) - (3*par[1]*par[2]*(-y[10] + y[19])*(-y[11] + y[20]))/pow(pow(-y[9] + y[18],2) + pow(-y[10] + y[19],2) + pow(-y[11] + y[20],2),2.5))/2.);
    break;

  case 528 : 
    return(0);
    break;

  case 529 : 
    return(0);
    break;

  case 530 : 
    return(0);
    break;

  case 531 : 
    return(0);
    break;

  case 532 : 
    return(0);
    break;

  case 533 : 
    return(0);
    break;

  case 534 : 
    return(((3*par[0]*par[2]*(y[0] - y[18])*(y[1] - y[19]))/pow(pow(y[0] - y[18],2) + pow(y[1] - y[19],2) + pow(y[2] - y[20],2),2.5) + (3*par[1]*par[2]*(y[9] - y[18])*(y[10] - y[19]))/pow(pow(y[9] - y[18],2) + pow(y[10] - y[19],2) + pow(y[11] - y[20],2),2.5) + (3*par[0]*par[2]*(-y[0] + y[18])*(-y[1] + y[19]))/pow(pow(-y[0] + y[18],2) + pow(-y[1] + y[19],2) + pow(-y[2] + y[20],2),2.5) + (3*par[1]*par[2]*(-y[9] + y[18])*(-y[10] + y[19]))/pow(pow(-y[9] + y[18],2) + pow(-y[10] + y[19],2) + pow(-y[11] + y[20],2),2.5))/2.);
    break;

  case 535 : 
    return(((3*par[0]*par[2]*pow(y[1] - y[19],2))/pow(pow(y[0] - y[18],2) + pow(y[1] - y[19],2) + pow(y[2] - y[20],2),2.5) - (par[0]*par[2])/pow(pow(y[0] - y[18],2) + pow(y[1] - y[19],2) + pow(y[2] - y[20],2),1.5) + (3*par[1]*par[2]*pow(y[10] - y[19],2))/pow(pow(y[9] - y[18],2) + pow(y[10] - y[19],2) + pow(y[11] - y[20],2),2.5) - (par[1]*par[2])/pow(pow(y[9] - y[18],2) + pow(y[10] - y[19],2) + pow(y[11] - y[20],2),1.5) + (3*par[0]*par[2]*pow(-y[1] + y[19],2))/pow(pow(-y[0] + y[18],2) + pow(-y[1] + y[19],2) + pow(-y[2] + y[20],2),2.5) - (par[0]*par[2])/pow(pow(-y[0] + y[18],2) + pow(-y[1] + y[19],2) + pow(-y[2] + y[20],2),1.5) + (3*par[1]*par[2]*pow(-y[10] + y[19],2))/pow(pow(-y[9] + y[18],2) + pow(-y[10] + y[19],2) + pow(-y[11] + y[20],2),2.5) - (par[1]*par[2])/pow(pow(-y[9] + y[18],2) + pow(-y[10] + y[19],2) + pow(-y[11] + y[20],2),1.5))/2.);
    break;

  case 536 : 
    return(((3*par[0]*par[2]*(y[1] - y[19])*(y[2] - y[20]))/pow(pow(y[0] - y[18],2) + pow(y[1] - y[19],2) + pow(y[2] - y[20],2),2.5) + (3*par[1]*par[2]*(y[10] - y[19])*(y[11] - y[20]))/pow(pow(y[9] - y[18],2) + pow(y[10] - y[19],2) + pow(y[11] - y[20],2),2.5) + (3*par[0]*par[2]*(-y[1] + y[19])*(-y[2] + y[20]))/pow(pow(-y[0] + y[18],2) + pow(-y[1] + y[19],2) + pow(-y[2] + y[20],2),2.5) + (3*par[1]*par[2]*(-y[10] + y[19])*(-y[11] + y[20]))/pow(pow(-y[9] + y[18],2) + pow(-y[10] + y[19],2) + pow(-y[11] + y[20],2),2.5))/2.);
    break;

  case 537 : 
    return(0);
    break;

  case 538 : 
    return(0);
    break;

  case 539 : 
    return(0);
    break;

  case 540 : 
    return(0);
    break;

  case 541 : 
    return(0);
    break;

  case 542 : 
    return(0);
    break;

  case 543 : 
    return(((-3*par[0]*par[2]*(y[0] - y[18])*(y[2] - y[20]))/pow(pow(y[0] - y[18],2) + pow(y[1] - y[19],2) + pow(y[2] - y[20],2),2.5) - (3*par[0]*par[2]*(-y[0] + y[18])*(-y[2] + y[20]))/pow(pow(-y[0] + y[18],2) + pow(-y[1] + y[19],2) + pow(-y[2] + y[20],2),2.5))/2.);
    break;

  case 544 : 
    return(((-3*par[0]*par[2]*(y[1] - y[19])*(y[2] - y[20]))/pow(pow(y[0] - y[18],2) + pow(y[1] - y[19],2) + pow(y[2] - y[20],2),2.5) - (3*par[0]*par[2]*(-y[1] + y[19])*(-y[2] + y[20]))/pow(pow(-y[0] + y[18],2) + pow(-y[1] + y[19],2) + pow(-y[2] + y[20],2),2.5))/2.);
    break;

  case 545 : 
    return(((par[0]*par[2])/pow(pow(y[0] - y[18],2) + pow(y[1] - y[19],2) + pow(y[2] - y[20],2),1.5) - (3*par[0]*par[2]*pow(y[2] - y[20],2))/pow(pow(y[0] - y[18],2) + pow(y[1] - y[19],2) + pow(y[2] - y[20],2),2.5) - (3*par[0]*par[2]*pow(-y[2] + y[20],2))/pow(pow(-y[0] + y[18],2) + pow(-y[1] + y[19],2) + pow(-y[2] + y[20],2),2.5) + (par[0]*par[2])/pow(pow(-y[0] + y[18],2) + pow(-y[1] + y[19],2) + pow(-y[2] + y[20],2),1.5))/2.);
    break;

  case 546 : 
    return(0);
    break;

  case 547 : 
    return(0);
    break;

  case 548 : 
    return(0);
    break;

  case 549 : 
    return(0);
    break;

  case 550 : 
    return(0);
    break;

  case 551 : 
    return(0);
    break;

  case 552 : 
    return(((-3*par[1]*par[2]*(y[9] - y[18])*(y[11] - y[20]))/pow(pow(y[9] - y[18],2) + pow(y[10] - y[19],2) + pow(y[11] - y[20],2),2.5) - (3*par[1]*par[2]*(-y[9] + y[18])*(-y[11] + y[20]))/pow(pow(-y[9] + y[18],2) + pow(-y[10] + y[19],2) + pow(-y[11] + y[20],2),2.5))/2.);
    break;

  case 553 : 
    return(((-3*par[1]*par[2]*(y[10] - y[19])*(y[11] - y[20]))/pow(pow(y[9] - y[18],2) + pow(y[10] - y[19],2) + pow(y[11] - y[20],2),2.5) - (3*par[1]*par[2]*(-y[10] + y[19])*(-y[11] + y[20]))/pow(pow(-y[9] + y[18],2) + pow(-y[10] + y[19],2) + pow(-y[11] + y[20],2),2.5))/2.);
    break;

  case 554 : 
    return(((par[1]*par[2])/pow(pow(y[9] - y[18],2) + pow(y[10] - y[19],2) + pow(y[11] - y[20],2),1.5) - (3*par[1]*par[2]*pow(y[11] - y[20],2))/pow(pow(y[9] - y[18],2) + pow(y[10] - y[19],2) + pow(y[11] - y[20],2),2.5) - (3*par[1]*par[2]*pow(-y[11] + y[20],2))/pow(pow(-y[9] + y[18],2) + pow(-y[10] + y[19],2) + pow(-y[11] + y[20],2),2.5) + (par[1]*par[2])/pow(pow(-y[9] + y[18],2) + pow(-y[10] + y[19],2) + pow(-y[11] + y[20],2),1.5))/2.);
    break;

  case 555 : 
    return(0);
    break;

  case 556 : 
    return(0);
    break;

  case 557 : 
    return(0);
    break;

  case 558 : 
    return(0);
    break;

  case 559 : 
    return(0);
    break;

  case 560 : 
    return(0);
    break;

  case 561 : 
    return(((3*par[0]*par[2]*(y[0] - y[18])*(y[2] - y[20]))/pow(pow(y[0] - y[18],2) + pow(y[1] - y[19],2) + pow(y[2] - y[20],2),2.5) + (3*par[1]*par[2]*(y[9] - y[18])*(y[11] - y[20]))/pow(pow(y[9] - y[18],2) + pow(y[10] - y[19],2) + pow(y[11] - y[20],2),2.5) + (3*par[0]*par[2]*(-y[0] + y[18])*(-y[2] + y[20]))/pow(pow(-y[0] + y[18],2) + pow(-y[1] + y[19],2) + pow(-y[2] + y[20],2),2.5) + (3*par[1]*par[2]*(-y[9] + y[18])*(-y[11] + y[20]))/pow(pow(-y[9] + y[18],2) + pow(-y[10] + y[19],2) + pow(-y[11] + y[20],2),2.5))/2.);
    break;

  case 562 : 
    return(((3*par[0]*par[2]*(y[1] - y[19])*(y[2] - y[20]))/pow(pow(y[0] - y[18],2) + pow(y[1] - y[19],2) + pow(y[2] - y[20],2),2.5) + (3*par[1]*par[2]*(y[10] - y[19])*(y[11] - y[20]))/pow(pow(y[9] - y[18],2) + pow(y[10] - y[19],2) + pow(y[11] - y[20],2),2.5) + (3*par[0]*par[2]*(-y[1] + y[19])*(-y[2] + y[20]))/pow(pow(-y[0] + y[18],2) + pow(-y[1] + y[19],2) + pow(-y[2] + y[20],2),2.5) + (3*par[1]*par[2]*(-y[10] + y[19])*(-y[11] + y[20]))/pow(pow(-y[9] + y[18],2) + pow(-y[10] + y[19],2) + pow(-y[11] + y[20],2),2.5))/2.);
    break;

  case 563 : 
    return((-((par[0]*par[2])/pow(pow(y[0] - y[18],2) + pow(y[1] - y[19],2) + pow(y[2] - y[20],2),1.5)) - (par[1]*par[2])/pow(pow(y[9] - y[18],2) + pow(y[10] - y[19],2) + pow(y[11] - y[20],2),1.5) + (3*par[0]*par[2]*pow(y[2] - y[20],2))/pow(pow(y[0] - y[18],2) + pow(y[1] - y[19],2) + pow(y[2] - y[20],2),2.5) + (3*par[1]*par[2]*pow(y[11] - y[20],2))/pow(pow(y[9] - y[18],2) + pow(y[10] - y[19],2) + pow(y[11] - y[20],2),2.5) + (3*par[0]*par[2]*pow(-y[2] + y[20],2))/pow(pow(-y[0] + y[18],2) + pow(-y[1] + y[19],2) + pow(-y[2] + y[20],2),2.5) - (par[0]*par[2])/pow(pow(-y[0] + y[18],2) + pow(-y[1] + y[19],2) + pow(-y[2] + y[20],2),1.5) + (3*par[1]*par[2]*pow(-y[11] + y[20],2))/pow(pow(-y[9] + y[18],2) + pow(-y[10] + y[19],2) + pow(-y[11] + y[20],2),2.5) - (par[1]*par[2])/pow(pow(-y[9] + y[18],2) + pow(-y[10] + y[19],2) + pow(-y[11] + y[20],2),1.5))/2.);
    break;

  case 564 : 
    return(0);
    break;

  case 565 : 
    return(0);
    break;

  case 566 : 
    return(0);
    break;

  case 567 : 
    return(0);
    break;

  case 568 : 
    return(0);
    break;

  case 569 : 
    return(0);
    break;

  case 570 : 
    return(0);
    break;

  case 571 : 
    return(0);
    break;

  case 572 : 
    return(0);
    break;

  case 573 : 
    return(0);
    break;

  case 574 : 
    return(0);
    break;

  case 575 : 
    return(0);
    break;

  case 576 : 
    return(0);
    break;

  case 577 : 
    return(0);
    break;

  case 578 : 
    return(0);
    break;

  case 579 : 
    return(0);
    break;

  case 580 : 
    return(0);
    break;

  case 581 : 
    return(0);
    break;

  case 582 : 
    return(0);
    break;

  case 583 : 
    return(0);
    break;

  case 584 : 
    return(0);
    break;

  case 585 : 
    return(1/par[2]);
    break;

  case 586 : 
    return(0);
    break;

  case 587 : 
    return(0);
    break;

  case 588 : 
    return(0);
    break;

  case 589 : 
    return(0);
    break;

  case 590 : 
    return(0);
    break;

  case 591 : 
    return(0);
    break;

  case 592 : 
    return(0);
    break;

  case 593 : 
    return(0);
    break;

  case 594 : 
    return(0);
    break;

  case 595 : 
    return(0);
    break;

  case 596 : 
    return(0);
    break;

  case 597 : 
    return(0);
    break;

  case 598 : 
    return(0);
    break;

  case 599 : 
    return(0);
    break;

  case 600 : 
    return(0);
    break;

  case 601 : 
    return(0);
    break;

  case 602 : 
    return(0);
    break;

  case 603 : 
    return(0);
    break;

  case 604 : 
    return(0);
    break;

  case 605 : 
    return(0);
    break;

  case 606 : 
    return(0);
    break;

  case 607 : 
    return(0);
    break;

  case 608 : 
    return(0);
    break;

  case 609 : 
    return(0);
    break;

  case 610 : 
    return(0);
    break;

  case 611 : 
    return(0);
    break;

  case 612 : 
    return(0);
    break;

  case 613 : 
    return(1/par[2]);
    break;

  case 614 : 
    return(0);
    break;

  case 615 : 
    return(0);
    break;

  case 616 : 
    return(0);
    break;

  case 617 : 
    return(0);
    break;

  case 618 : 
    return(0);
    break;

  case 619 : 
    return(0);
    break;

  case 620 : 
    return(0);
    break;

  case 621 : 
    return(0);
    break;

  case 622 : 
    return(0);
    break;

  case 623 : 
    return(0);
    break;

  case 624 : 
    return(0);
    break;

  case 625 : 
    return(0);
    break;

  case 626 : 
    return(0);
    break;

  case 627 : 
    return(0);
    break;

  case 628 : 
    return(0);
    break;

  case 629 : 
    return(0);
    break;

  case 630 : 
    return(0);
    break;

  case 631 : 
    return(0);
    break;

  case 632 : 
    return(0);
    break;

  case 633 : 
    return(0);
    break;

  case 634 : 
    return(0);
    break;

  case 635 : 
    return(0);
    break;

  case 636 : 
    return(0);
    break;

  case 637 : 
    return(0);
    break;

  case 638 : 
    return(0);
    break;

  case 639 : 
    return(0);
    break;

  case 640 : 
    return(0);
    break;

  case 641 : 
    return(1/par[2]);
    break;

  case 642 : 
    return(0);
    break;

  case 643 : 
    return(0);
    break;

  case 644 : 
    return(0);
    break;

  case 645 : 
    return(0);
    break;

  case 646 : 
    return(0);
    break;

  case 647 : 
    return(0);
    break;

  case 648 : 
    return(0);
    break;

  case 649 : 
    return(0);
    break;

  case 650 : 
    return(0);
    break;

  case 651 : 
    return(0);
    break;

  case 652 : 
    return(0);
    break;

  case 653 : 
    return(0);
    break;

  case 654 : 
    return(0);
    break;

  case 655 : 
    return(0);
    break;

  case 656 : 
    return(0);
    break;

  case 657 : 
    return(0);
    break;

  case 658 : 
    return(0);
    break;

  case 659 : 
    return(0);
    break;

  case 660 : 
    return(0);
    break;

  case 661 : 
    return(0);
    break;

  case 662 : 
    return(0);
    break;

  case 663 : 
    return(0);
    break;

  case 664 : 
    return(0);
    break;

  case 665 : 
    return(0);
    break;

  case 666 : 
    return(0);
    break;

  case 667 : 
    return(0);
    break;

  case 668 : 
    return(0);
    break;

  case 669 : 
    return(0);
    break;

  case 670 : 
    return(0);
    break;

  case 671 : 
    return(0);
    break;

  case 672 : 
    return(0);
    break;

  case 673 : 
    return(0);
    break;

  case 674 : 
    return(0);
    break;

  case 675 : 
    return(0);
    break;

  case 676 : 
    return(0);
    break;

  case 677 : 
    return(0);
    break;

  case 678 : 
    return(0);
    break;

  case 679 : 
    return(0);
    break;

  case 680 : 
    return(0);
    break;

  case 681 : 
    return(0);
    break;

  case 682 : 
    return(0);
    break;

  case 683 : 
    return(0);
    break;

  case 684 : 
    return(0);
    break;

  case 685 : 
    return(0);
    break;

  case 686 : 
    return(0);
    break;

  case 687 : 
    return(0);
    break;

  case 688 : 
    return(0);
    break;

  case 689 : 
    return(0);
    break;

  case 690 : 
    return(0);
    break;

  case 691 : 
    return(0);
    break;

  case 692 : 
    return(0);
    break;

  case 693 : 
    return(0);
    break;

  case 694 : 
    return(0);
    break;

  case 695 : 
    return(0);
    break;

  case 696 : 
    return(0);
    break;

  case 697 : 
    return(0);
    break;

  case 698 : 
    return(0);
    break;

  case 699 : 
    return(0);
    break;

  case 700 : 
    return(0);
    break;

  case 701 : 
    return(0);
    break;

  case 702 : 
    return(0);
    break;

  case 703 : 
    return(0);
    break;

  case 704 : 
    return(0);
    break;

  case 705 : 
    return(0);
    break;

  case 706 : 
    return(0);
    break;

  case 707 : 
    return(0);
    break;

  case 708 : 
    return(0);
    break;

  case 709 : 
    return(0);
    break;

  case 710 : 
    return(0);
    break;

  case 711 : 
    return(0);
    break;

  case 712 : 
    return(0);
    break;

  case 713 : 
    return(0);
    break;

  case 714 : 
    return(0);
    break;

  case 715 : 
    return(0);
    break;

  case 716 : 
    return(0);
    break;

  case 717 : 
    return(0);
    break;

  case 718 : 
    return(0);
    break;

  case 719 : 
    return(0);
    break;

  case 720 : 
    return(0);
    break;

  case 721 : 
    return(0);
    break;

  case 722 : 
    return(0);
    break;

  case 723 : 
    return(0);
    break;

  case 724 : 
    return(0);
    break;

  case 725 : 
    return(0);
    break;

  case 726 : 
    return(0);
    break;

  case 727 : 
    return(0);
    break;

  case 728 : 
    return(0);
    break;
}

 return 0; 
}

