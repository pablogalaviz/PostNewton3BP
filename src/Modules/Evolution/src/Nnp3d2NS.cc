#include "evolution.h"



double evolution::rhs_Nnp3d2_FNS(double t, const double y[], double par[3], int i)
{
  
  switch ( i ) {

  case 0 : 
    return(y[2]/par[0]);
    break;

  case 1 : 
    return(y[3]/par[0]);
    break;

  case 2 : 
    return((-((par[0]*par[1]*(y[0] - y[4]))/pow(pow(y[0] - y[4],2) + pow(y[1] - y[5],2),1.5)) + (par[0]*par[1]*(-y[0] + y[4]))/pow(pow(-y[0] + y[4],2) + pow(-y[1] + y[5],2),1.5) - (par[0]*par[2]*(y[0] - y[8]))/pow(pow(y[0] - y[8],2) + pow(y[1] - y[9],2),1.5) + (par[0]*par[2]*(-y[0] + y[8]))/pow(pow(-y[0] + y[8],2) + pow(-y[1] + y[9],2),1.5))/2.);
    break;

  case 3 : 
    return((-((par[0]*par[1]*(y[1] - y[5]))/pow(pow(y[0] - y[4],2) + pow(y[1] - y[5],2),1.5)) + (par[0]*par[1]*(-y[1] + y[5]))/pow(pow(-y[0] + y[4],2) + pow(-y[1] + y[5],2),1.5) - (par[0]*par[2]*(y[1] - y[9]))/pow(pow(y[0] - y[8],2) + pow(y[1] - y[9],2),1.5) + (par[0]*par[2]*(-y[1] + y[9]))/pow(pow(-y[0] + y[8],2) + pow(-y[1] + y[9],2),1.5))/2.);
    break;

  case 4 : 
    return(y[6]/par[1]);
    break;

  case 5 : 
    return(y[7]/par[1]);
    break;

  case 6 : 
    return(((par[0]*par[1]*(y[0] - y[4]))/pow(pow(y[0] - y[4],2) + pow(y[1] - y[5],2),1.5) - (par[0]*par[1]*(-y[0] + y[4]))/pow(pow(-y[0] + y[4],2) + pow(-y[1] + y[5],2),1.5) - (par[1]*par[2]*(y[4] - y[8]))/pow(pow(y[4] - y[8],2) + pow(y[5] - y[9],2),1.5) + (par[1]*par[2]*(-y[4] + y[8]))/pow(pow(-y[4] + y[8],2) + pow(-y[5] + y[9],2),1.5))/2.);
    break;

  case 7 : 
    return(((par[0]*par[1]*(y[1] - y[5]))/pow(pow(y[0] - y[4],2) + pow(y[1] - y[5],2),1.5) - (par[0]*par[1]*(-y[1] + y[5]))/pow(pow(-y[0] + y[4],2) + pow(-y[1] + y[5],2),1.5) - (par[1]*par[2]*(y[5] - y[9]))/pow(pow(y[4] - y[8],2) + pow(y[5] - y[9],2),1.5) + (par[1]*par[2]*(-y[5] + y[9]))/pow(pow(-y[4] + y[8],2) + pow(-y[5] + y[9],2),1.5))/2.);
    break;

  case 8 : 
    return(y[10]/par[2]);
    break;

  case 9 : 
    return(y[11]/par[2]);
    break;

  case 10 : 
    return(((par[0]*par[2]*(y[0] - y[8]))/pow(pow(y[0] - y[8],2) + pow(y[1] - y[9],2),1.5) + (par[1]*par[2]*(y[4] - y[8]))/pow(pow(y[4] - y[8],2) + pow(y[5] - y[9],2),1.5) - (par[0]*par[2]*(-y[0] + y[8]))/pow(pow(-y[0] + y[8],2) + pow(-y[1] + y[9],2),1.5) - (par[1]*par[2]*(-y[4] + y[8]))/pow(pow(-y[4] + y[8],2) + pow(-y[5] + y[9],2),1.5))/2.);
    break;

  case 11 : 
    return(((par[0]*par[2]*(y[1] - y[9]))/pow(pow(y[0] - y[8],2) + pow(y[1] - y[9],2),1.5) + (par[1]*par[2]*(y[5] - y[9]))/pow(pow(y[4] - y[8],2) + pow(y[5] - y[9],2),1.5) - (par[0]*par[2]*(-y[1] + y[9]))/pow(pow(-y[0] + y[8],2) + pow(-y[1] + y[9],2),1.5) - (par[1]*par[2]*(-y[5] + y[9]))/pow(pow(-y[4] + y[8],2) + pow(-y[5] + y[9],2),1.5))/2.);
    break;
}

 return 0; 
}


double evolution::jac_Nnp3d2_FaNS(double t, const double y[], double par[3], int i, int j)
{

 switch ( j*NDvar+i ) {

  case 0 : 
    return(0);
    break;

  case 1 : 
    return(0);
    break;

  case 2 : 
    return(((3*par[0]*par[1]*pow(y[0] - y[4],2))/pow(pow(y[0] - y[4],2) + pow(y[1] - y[5],2),2.5) - (par[0]*par[1])/pow(pow(y[0] - y[4],2) + pow(y[1] - y[5],2),1.5) + (3*par[0]*par[1]*pow(-y[0] + y[4],2))/pow(pow(-y[0] + y[4],2) + pow(-y[1] + y[5],2),2.5) - (par[0]*par[1])/pow(pow(-y[0] + y[4],2) + pow(-y[1] + y[5],2),1.5) + (3*par[0]*par[2]*pow(y[0] - y[8],2))/pow(pow(y[0] - y[8],2) + pow(y[1] - y[9],2),2.5) - (par[0]*par[2])/pow(pow(y[0] - y[8],2) + pow(y[1] - y[9],2),1.5) + (3*par[0]*par[2]*pow(-y[0] + y[8],2))/pow(pow(-y[0] + y[8],2) + pow(-y[1] + y[9],2),2.5) - (par[0]*par[2])/pow(pow(-y[0] + y[8],2) + pow(-y[1] + y[9],2),1.5))/2.);
    break;

  case 3 : 
    return(((3*par[0]*par[1]*(y[0] - y[4])*(y[1] - y[5]))/pow(pow(y[0] - y[4],2) + pow(y[1] - y[5],2),2.5) + (3*par[0]*par[1]*(-y[0] + y[4])*(-y[1] + y[5]))/pow(pow(-y[0] + y[4],2) + pow(-y[1] + y[5],2),2.5) + (3*par[0]*par[2]*(y[0] - y[8])*(y[1] - y[9]))/pow(pow(y[0] - y[8],2) + pow(y[1] - y[9],2),2.5) + (3*par[0]*par[2]*(-y[0] + y[8])*(-y[1] + y[9]))/pow(pow(-y[0] + y[8],2) + pow(-y[1] + y[9],2),2.5))/2.);
    break;

  case 4 : 
    return(0);
    break;

  case 5 : 
    return(0);
    break;

  case 6 : 
    return(((-3*par[0]*par[1]*pow(y[0] - y[4],2))/pow(pow(y[0] - y[4],2) + pow(y[1] - y[5],2),2.5) + (par[0]*par[1])/pow(pow(y[0] - y[4],2) + pow(y[1] - y[5],2),1.5) - (3*par[0]*par[1]*pow(-y[0] + y[4],2))/pow(pow(-y[0] + y[4],2) + pow(-y[1] + y[5],2),2.5) + (par[0]*par[1])/pow(pow(-y[0] + y[4],2) + pow(-y[1] + y[5],2),1.5))/2.);
    break;

  case 7 : 
    return(((-3*par[0]*par[1]*(y[0] - y[4])*(y[1] - y[5]))/pow(pow(y[0] - y[4],2) + pow(y[1] - y[5],2),2.5) - (3*par[0]*par[1]*(-y[0] + y[4])*(-y[1] + y[5]))/pow(pow(-y[0] + y[4],2) + pow(-y[1] + y[5],2),2.5))/2.);
    break;

  case 8 : 
    return(0);
    break;

  case 9 : 
    return(0);
    break;

  case 10 : 
    return(((-3*par[0]*par[2]*pow(y[0] - y[8],2))/pow(pow(y[0] - y[8],2) + pow(y[1] - y[9],2),2.5) + (par[0]*par[2])/pow(pow(y[0] - y[8],2) + pow(y[1] - y[9],2),1.5) - (3*par[0]*par[2]*pow(-y[0] + y[8],2))/pow(pow(-y[0] + y[8],2) + pow(-y[1] + y[9],2),2.5) + (par[0]*par[2])/pow(pow(-y[0] + y[8],2) + pow(-y[1] + y[9],2),1.5))/2.);
    break;

  case 11 : 
    return(((-3*par[0]*par[2]*(y[0] - y[8])*(y[1] - y[9]))/pow(pow(y[0] - y[8],2) + pow(y[1] - y[9],2),2.5) - (3*par[0]*par[2]*(-y[0] + y[8])*(-y[1] + y[9]))/pow(pow(-y[0] + y[8],2) + pow(-y[1] + y[9],2),2.5))/2.);
    break;

  case 12 : 
    return(0);
    break;

  case 13 : 
    return(0);
    break;

  case 14 : 
    return(((3*par[0]*par[1]*(y[0] - y[4])*(y[1] - y[5]))/pow(pow(y[0] - y[4],2) + pow(y[1] - y[5],2),2.5) + (3*par[0]*par[1]*(-y[0] + y[4])*(-y[1] + y[5]))/pow(pow(-y[0] + y[4],2) + pow(-y[1] + y[5],2),2.5) + (3*par[0]*par[2]*(y[0] - y[8])*(y[1] - y[9]))/pow(pow(y[0] - y[8],2) + pow(y[1] - y[9],2),2.5) + (3*par[0]*par[2]*(-y[0] + y[8])*(-y[1] + y[9]))/pow(pow(-y[0] + y[8],2) + pow(-y[1] + y[9],2),2.5))/2.);
    break;

  case 15 : 
    return((-((par[0]*par[1])/pow(pow(y[0] - y[4],2) + pow(y[1] - y[5],2),1.5)) + (3*par[0]*par[1]*pow(y[1] - y[5],2))/pow(pow(y[0] - y[4],2) + pow(y[1] - y[5],2),2.5) + (3*par[0]*par[1]*pow(-y[1] + y[5],2))/pow(pow(-y[0] + y[4],2) + pow(-y[1] + y[5],2),2.5) - (par[0]*par[1])/pow(pow(-y[0] + y[4],2) + pow(-y[1] + y[5],2),1.5) - (par[0]*par[2])/pow(pow(y[0] - y[8],2) + pow(y[1] - y[9],2),1.5) + (3*par[0]*par[2]*pow(y[1] - y[9],2))/pow(pow(y[0] - y[8],2) + pow(y[1] - y[9],2),2.5) + (3*par[0]*par[2]*pow(-y[1] + y[9],2))/pow(pow(-y[0] + y[8],2) + pow(-y[1] + y[9],2),2.5) - (par[0]*par[2])/pow(pow(-y[0] + y[8],2) + pow(-y[1] + y[9],2),1.5))/2.);
    break;

  case 16 : 
    return(0);
    break;

  case 17 : 
    return(0);
    break;

  case 18 : 
    return(((-3*par[0]*par[1]*(y[0] - y[4])*(y[1] - y[5]))/pow(pow(y[0] - y[4],2) + pow(y[1] - y[5],2),2.5) - (3*par[0]*par[1]*(-y[0] + y[4])*(-y[1] + y[5]))/pow(pow(-y[0] + y[4],2) + pow(-y[1] + y[5],2),2.5))/2.);
    break;

  case 19 : 
    return(((par[0]*par[1])/pow(pow(y[0] - y[4],2) + pow(y[1] - y[5],2),1.5) - (3*par[0]*par[1]*pow(y[1] - y[5],2))/pow(pow(y[0] - y[4],2) + pow(y[1] - y[5],2),2.5) - (3*par[0]*par[1]*pow(-y[1] + y[5],2))/pow(pow(-y[0] + y[4],2) + pow(-y[1] + y[5],2),2.5) + (par[0]*par[1])/pow(pow(-y[0] + y[4],2) + pow(-y[1] + y[5],2),1.5))/2.);
    break;

  case 20 : 
    return(0);
    break;

  case 21 : 
    return(0);
    break;

  case 22 : 
    return(((-3*par[0]*par[2]*(y[0] - y[8])*(y[1] - y[9]))/pow(pow(y[0] - y[8],2) + pow(y[1] - y[9],2),2.5) - (3*par[0]*par[2]*(-y[0] + y[8])*(-y[1] + y[9]))/pow(pow(-y[0] + y[8],2) + pow(-y[1] + y[9],2),2.5))/2.);
    break;

  case 23 : 
    return(((par[0]*par[2])/pow(pow(y[0] - y[8],2) + pow(y[1] - y[9],2),1.5) - (3*par[0]*par[2]*pow(y[1] - y[9],2))/pow(pow(y[0] - y[8],2) + pow(y[1] - y[9],2),2.5) - (3*par[0]*par[2]*pow(-y[1] + y[9],2))/pow(pow(-y[0] + y[8],2) + pow(-y[1] + y[9],2),2.5) + (par[0]*par[2])/pow(pow(-y[0] + y[8],2) + pow(-y[1] + y[9],2),1.5))/2.);
    break;

  case 24 : 
    return(1/par[0]);
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
    return(0);
    break;

  case 31 : 
    return(0);
    break;

  case 32 : 
    return(0);
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
    return(1/par[0]);
    break;

  case 38 : 
    return(0);
    break;

  case 39 : 
    return(0);
    break;

  case 40 : 
    return(0);
    break;

  case 41 : 
    return(0);
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
    return(0);
    break;

  case 49 : 
    return(0);
    break;

  case 50 : 
    return(((-3*par[0]*par[1]*pow(y[0] - y[4],2))/pow(pow(y[0] - y[4],2) + pow(y[1] - y[5],2),2.5) + (par[0]*par[1])/pow(pow(y[0] - y[4],2) + pow(y[1] - y[5],2),1.5) - (3*par[0]*par[1]*pow(-y[0] + y[4],2))/pow(pow(-y[0] + y[4],2) + pow(-y[1] + y[5],2),2.5) + (par[0]*par[1])/pow(pow(-y[0] + y[4],2) + pow(-y[1] + y[5],2),1.5))/2.);
    break;

  case 51 : 
    return(((-3*par[0]*par[1]*(y[0] - y[4])*(y[1] - y[5]))/pow(pow(y[0] - y[4],2) + pow(y[1] - y[5],2),2.5) - (3*par[0]*par[1]*(-y[0] + y[4])*(-y[1] + y[5]))/pow(pow(-y[0] + y[4],2) + pow(-y[1] + y[5],2),2.5))/2.);
    break;

  case 52 : 
    return(0);
    break;

  case 53 : 
    return(0);
    break;

  case 54 : 
    return(((3*par[0]*par[1]*pow(y[0] - y[4],2))/pow(pow(y[0] - y[4],2) + pow(y[1] - y[5],2),2.5) - (par[0]*par[1])/pow(pow(y[0] - y[4],2) + pow(y[1] - y[5],2),1.5) + (3*par[0]*par[1]*pow(-y[0] + y[4],2))/pow(pow(-y[0] + y[4],2) + pow(-y[1] + y[5],2),2.5) - (par[0]*par[1])/pow(pow(-y[0] + y[4],2) + pow(-y[1] + y[5],2),1.5) + (3*par[1]*par[2]*pow(y[4] - y[8],2))/pow(pow(y[4] - y[8],2) + pow(y[5] - y[9],2),2.5) - (par[1]*par[2])/pow(pow(y[4] - y[8],2) + pow(y[5] - y[9],2),1.5) + (3*par[1]*par[2]*pow(-y[4] + y[8],2))/pow(pow(-y[4] + y[8],2) + pow(-y[5] + y[9],2),2.5) - (par[1]*par[2])/pow(pow(-y[4] + y[8],2) + pow(-y[5] + y[9],2),1.5))/2.);
    break;

  case 55 : 
    return(((3*par[0]*par[1]*(y[0] - y[4])*(y[1] - y[5]))/pow(pow(y[0] - y[4],2) + pow(y[1] - y[5],2),2.5) + (3*par[0]*par[1]*(-y[0] + y[4])*(-y[1] + y[5]))/pow(pow(-y[0] + y[4],2) + pow(-y[1] + y[5],2),2.5) + (3*par[1]*par[2]*(y[4] - y[8])*(y[5] - y[9]))/pow(pow(y[4] - y[8],2) + pow(y[5] - y[9],2),2.5) + (3*par[1]*par[2]*(-y[4] + y[8])*(-y[5] + y[9]))/pow(pow(-y[4] + y[8],2) + pow(-y[5] + y[9],2),2.5))/2.);
    break;

  case 56 : 
    return(0);
    break;

  case 57 : 
    return(0);
    break;

  case 58 : 
    return(((-3*par[1]*par[2]*pow(y[4] - y[8],2))/pow(pow(y[4] - y[8],2) + pow(y[5] - y[9],2),2.5) + (par[1]*par[2])/pow(pow(y[4] - y[8],2) + pow(y[5] - y[9],2),1.5) - (3*par[1]*par[2]*pow(-y[4] + y[8],2))/pow(pow(-y[4] + y[8],2) + pow(-y[5] + y[9],2),2.5) + (par[1]*par[2])/pow(pow(-y[4] + y[8],2) + pow(-y[5] + y[9],2),1.5))/2.);
    break;

  case 59 : 
    return(((-3*par[1]*par[2]*(y[4] - y[8])*(y[5] - y[9]))/pow(pow(y[4] - y[8],2) + pow(y[5] - y[9],2),2.5) - (3*par[1]*par[2]*(-y[4] + y[8])*(-y[5] + y[9]))/pow(pow(-y[4] + y[8],2) + pow(-y[5] + y[9],2),2.5))/2.);
    break;

  case 60 : 
    return(0);
    break;

  case 61 : 
    return(0);
    break;

  case 62 : 
    return(((-3*par[0]*par[1]*(y[0] - y[4])*(y[1] - y[5]))/pow(pow(y[0] - y[4],2) + pow(y[1] - y[5],2),2.5) - (3*par[0]*par[1]*(-y[0] + y[4])*(-y[1] + y[5]))/pow(pow(-y[0] + y[4],2) + pow(-y[1] + y[5],2),2.5))/2.);
    break;

  case 63 : 
    return(((par[0]*par[1])/pow(pow(y[0] - y[4],2) + pow(y[1] - y[5],2),1.5) - (3*par[0]*par[1]*pow(y[1] - y[5],2))/pow(pow(y[0] - y[4],2) + pow(y[1] - y[5],2),2.5) - (3*par[0]*par[1]*pow(-y[1] + y[5],2))/pow(pow(-y[0] + y[4],2) + pow(-y[1] + y[5],2),2.5) + (par[0]*par[1])/pow(pow(-y[0] + y[4],2) + pow(-y[1] + y[5],2),1.5))/2.);
    break;

  case 64 : 
    return(0);
    break;

  case 65 : 
    return(0);
    break;

  case 66 : 
    return(((3*par[0]*par[1]*(y[0] - y[4])*(y[1] - y[5]))/pow(pow(y[0] - y[4],2) + pow(y[1] - y[5],2),2.5) + (3*par[0]*par[1]*(-y[0] + y[4])*(-y[1] + y[5]))/pow(pow(-y[0] + y[4],2) + pow(-y[1] + y[5],2),2.5) + (3*par[1]*par[2]*(y[4] - y[8])*(y[5] - y[9]))/pow(pow(y[4] - y[8],2) + pow(y[5] - y[9],2),2.5) + (3*par[1]*par[2]*(-y[4] + y[8])*(-y[5] + y[9]))/pow(pow(-y[4] + y[8],2) + pow(-y[5] + y[9],2),2.5))/2.);
    break;

  case 67 : 
    return((-((par[0]*par[1])/pow(pow(y[0] - y[4],2) + pow(y[1] - y[5],2),1.5)) + (3*par[0]*par[1]*pow(y[1] - y[5],2))/pow(pow(y[0] - y[4],2) + pow(y[1] - y[5],2),2.5) + (3*par[0]*par[1]*pow(-y[1] + y[5],2))/pow(pow(-y[0] + y[4],2) + pow(-y[1] + y[5],2),2.5) - (par[0]*par[1])/pow(pow(-y[0] + y[4],2) + pow(-y[1] + y[5],2),1.5) - (par[1]*par[2])/pow(pow(y[4] - y[8],2) + pow(y[5] - y[9],2),1.5) + (3*par[1]*par[2]*pow(y[5] - y[9],2))/pow(pow(y[4] - y[8],2) + pow(y[5] - y[9],2),2.5) + (3*par[1]*par[2]*pow(-y[5] + y[9],2))/pow(pow(-y[4] + y[8],2) + pow(-y[5] + y[9],2),2.5) - (par[1]*par[2])/pow(pow(-y[4] + y[8],2) + pow(-y[5] + y[9],2),1.5))/2.);
    break;

  case 68 : 
    return(0);
    break;

  case 69 : 
    return(0);
    break;

  case 70 : 
    return(((-3*par[1]*par[2]*(y[4] - y[8])*(y[5] - y[9]))/pow(pow(y[4] - y[8],2) + pow(y[5] - y[9],2),2.5) - (3*par[1]*par[2]*(-y[4] + y[8])*(-y[5] + y[9]))/pow(pow(-y[4] + y[8],2) + pow(-y[5] + y[9],2),2.5))/2.);
    break;

  case 71 : 
    return(((par[1]*par[2])/pow(pow(y[4] - y[8],2) + pow(y[5] - y[9],2),1.5) - (3*par[1]*par[2]*pow(y[5] - y[9],2))/pow(pow(y[4] - y[8],2) + pow(y[5] - y[9],2),2.5) - (3*par[1]*par[2]*pow(-y[5] + y[9],2))/pow(pow(-y[4] + y[8],2) + pow(-y[5] + y[9],2),2.5) + (par[1]*par[2])/pow(pow(-y[4] + y[8],2) + pow(-y[5] + y[9],2),1.5))/2.);
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
    return(0);
    break;

  case 76 : 
    return(1/par[1]);
    break;

  case 77 : 
    return(0);
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
    return(0);
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
    return(1/par[1]);
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
    return(((-3*par[0]*par[2]*pow(y[0] - y[8],2))/pow(pow(y[0] - y[8],2) + pow(y[1] - y[9],2),2.5) + (par[0]*par[2])/pow(pow(y[0] - y[8],2) + pow(y[1] - y[9],2),1.5) - (3*par[0]*par[2]*pow(-y[0] + y[8],2))/pow(pow(-y[0] + y[8],2) + pow(-y[1] + y[9],2),2.5) + (par[0]*par[2])/pow(pow(-y[0] + y[8],2) + pow(-y[1] + y[9],2),1.5))/2.);
    break;

  case 99 : 
    return(((-3*par[0]*par[2]*(y[0] - y[8])*(y[1] - y[9]))/pow(pow(y[0] - y[8],2) + pow(y[1] - y[9],2),2.5) - (3*par[0]*par[2]*(-y[0] + y[8])*(-y[1] + y[9]))/pow(pow(-y[0] + y[8],2) + pow(-y[1] + y[9],2),2.5))/2.);
    break;

  case 100 : 
    return(0);
    break;

  case 101 : 
    return(0);
    break;

  case 102 : 
    return(((-3*par[1]*par[2]*pow(y[4] - y[8],2))/pow(pow(y[4] - y[8],2) + pow(y[5] - y[9],2),2.5) + (par[1]*par[2])/pow(pow(y[4] - y[8],2) + pow(y[5] - y[9],2),1.5) - (3*par[1]*par[2]*pow(-y[4] + y[8],2))/pow(pow(-y[4] + y[8],2) + pow(-y[5] + y[9],2),2.5) + (par[1]*par[2])/pow(pow(-y[4] + y[8],2) + pow(-y[5] + y[9],2),1.5))/2.);
    break;

  case 103 : 
    return(((-3*par[1]*par[2]*(y[4] - y[8])*(y[5] - y[9]))/pow(pow(y[4] - y[8],2) + pow(y[5] - y[9],2),2.5) - (3*par[1]*par[2]*(-y[4] + y[8])*(-y[5] + y[9]))/pow(pow(-y[4] + y[8],2) + pow(-y[5] + y[9],2),2.5))/2.);
    break;

  case 104 : 
    return(0);
    break;

  case 105 : 
    return(0);
    break;

  case 106 : 
    return(((3*par[0]*par[2]*pow(y[0] - y[8],2))/pow(pow(y[0] - y[8],2) + pow(y[1] - y[9],2),2.5) - (par[0]*par[2])/pow(pow(y[0] - y[8],2) + pow(y[1] - y[9],2),1.5) + (3*par[1]*par[2]*pow(y[4] - y[8],2))/pow(pow(y[4] - y[8],2) + pow(y[5] - y[9],2),2.5) - (par[1]*par[2])/pow(pow(y[4] - y[8],2) + pow(y[5] - y[9],2),1.5) + (3*par[0]*par[2]*pow(-y[0] + y[8],2))/pow(pow(-y[0] + y[8],2) + pow(-y[1] + y[9],2),2.5) - (par[0]*par[2])/pow(pow(-y[0] + y[8],2) + pow(-y[1] + y[9],2),1.5) + (3*par[1]*par[2]*pow(-y[4] + y[8],2))/pow(pow(-y[4] + y[8],2) + pow(-y[5] + y[9],2),2.5) - (par[1]*par[2])/pow(pow(-y[4] + y[8],2) + pow(-y[5] + y[9],2),1.5))/2.);
    break;

  case 107 : 
    return(((3*par[0]*par[2]*(y[0] - y[8])*(y[1] - y[9]))/pow(pow(y[0] - y[8],2) + pow(y[1] - y[9],2),2.5) + (3*par[1]*par[2]*(y[4] - y[8])*(y[5] - y[9]))/pow(pow(y[4] - y[8],2) + pow(y[5] - y[9],2),2.5) + (3*par[0]*par[2]*(-y[0] + y[8])*(-y[1] + y[9]))/pow(pow(-y[0] + y[8],2) + pow(-y[1] + y[9],2),2.5) + (3*par[1]*par[2]*(-y[4] + y[8])*(-y[5] + y[9]))/pow(pow(-y[4] + y[8],2) + pow(-y[5] + y[9],2),2.5))/2.);
    break;

  case 108 : 
    return(0);
    break;

  case 109 : 
    return(0);
    break;

  case 110 : 
    return(((-3*par[0]*par[2]*(y[0] - y[8])*(y[1] - y[9]))/pow(pow(y[0] - y[8],2) + pow(y[1] - y[9],2),2.5) - (3*par[0]*par[2]*(-y[0] + y[8])*(-y[1] + y[9]))/pow(pow(-y[0] + y[8],2) + pow(-y[1] + y[9],2),2.5))/2.);
    break;

  case 111 : 
    return(((par[0]*par[2])/pow(pow(y[0] - y[8],2) + pow(y[1] - y[9],2),1.5) - (3*par[0]*par[2]*pow(y[1] - y[9],2))/pow(pow(y[0] - y[8],2) + pow(y[1] - y[9],2),2.5) - (3*par[0]*par[2]*pow(-y[1] + y[9],2))/pow(pow(-y[0] + y[8],2) + pow(-y[1] + y[9],2),2.5) + (par[0]*par[2])/pow(pow(-y[0] + y[8],2) + pow(-y[1] + y[9],2),1.5))/2.);
    break;

  case 112 : 
    return(0);
    break;

  case 113 : 
    return(0);
    break;

  case 114 : 
    return(((-3*par[1]*par[2]*(y[4] - y[8])*(y[5] - y[9]))/pow(pow(y[4] - y[8],2) + pow(y[5] - y[9],2),2.5) - (3*par[1]*par[2]*(-y[4] + y[8])*(-y[5] + y[9]))/pow(pow(-y[4] + y[8],2) + pow(-y[5] + y[9],2),2.5))/2.);
    break;

  case 115 : 
    return(((par[1]*par[2])/pow(pow(y[4] - y[8],2) + pow(y[5] - y[9],2),1.5) - (3*par[1]*par[2]*pow(y[5] - y[9],2))/pow(pow(y[4] - y[8],2) + pow(y[5] - y[9],2),2.5) - (3*par[1]*par[2]*pow(-y[5] + y[9],2))/pow(pow(-y[4] + y[8],2) + pow(-y[5] + y[9],2),2.5) + (par[1]*par[2])/pow(pow(-y[4] + y[8],2) + pow(-y[5] + y[9],2),1.5))/2.);
    break;

  case 116 : 
    return(0);
    break;

  case 117 : 
    return(0);
    break;

  case 118 : 
    return(((3*par[0]*par[2]*(y[0] - y[8])*(y[1] - y[9]))/pow(pow(y[0] - y[8],2) + pow(y[1] - y[9],2),2.5) + (3*par[1]*par[2]*(y[4] - y[8])*(y[5] - y[9]))/pow(pow(y[4] - y[8],2) + pow(y[5] - y[9],2),2.5) + (3*par[0]*par[2]*(-y[0] + y[8])*(-y[1] + y[9]))/pow(pow(-y[0] + y[8],2) + pow(-y[1] + y[9],2),2.5) + (3*par[1]*par[2]*(-y[4] + y[8])*(-y[5] + y[9]))/pow(pow(-y[4] + y[8],2) + pow(-y[5] + y[9],2),2.5))/2.);
    break;

  case 119 : 
    return((-((par[0]*par[2])/pow(pow(y[0] - y[8],2) + pow(y[1] - y[9],2),1.5)) - (par[1]*par[2])/pow(pow(y[4] - y[8],2) + pow(y[5] - y[9],2),1.5) + (3*par[0]*par[2]*pow(y[1] - y[9],2))/pow(pow(y[0] - y[8],2) + pow(y[1] - y[9],2),2.5) + (3*par[1]*par[2]*pow(y[5] - y[9],2))/pow(pow(y[4] - y[8],2) + pow(y[5] - y[9],2),2.5) + (3*par[0]*par[2]*pow(-y[1] + y[9],2))/pow(pow(-y[0] + y[8],2) + pow(-y[1] + y[9],2),2.5) - (par[0]*par[2])/pow(pow(-y[0] + y[8],2) + pow(-y[1] + y[9],2),1.5) + (3*par[1]*par[2]*pow(-y[5] + y[9],2))/pow(pow(-y[4] + y[8],2) + pow(-y[5] + y[9],2),2.5) - (par[1]*par[2])/pow(pow(-y[4] + y[8],2) + pow(-y[5] + y[9],2),1.5))/2.);
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
    return(1/par[2]);
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
    return(0);
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
    return(1/par[2]);
    break;

  case 142 : 
    return(0);
    break;

  case 143 : 
    return(0);
    break;
}

 return 0; 
}

