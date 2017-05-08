#include "evolution.h"



double evolution::rhs_Nnp3d3_FNS(double t, const double y[], double par[3], int i)
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
    return((-((par[0]*par[1]*(y[0] - y[6]))/pow(pow(y[0] - y[6],2) + pow(y[1] - y[7],2) + pow(y[2] - y[8],2),1.5)) + (par[0]*par[1]*(-y[0] + y[6]))/pow(pow(-y[0] + y[6],2) + pow(-y[1] + y[7],2) + pow(-y[2] + y[8],2),1.5) - (par[0]*par[2]*(y[0] - y[12]))/pow(pow(y[0] - y[12],2) + pow(y[1] - y[13],2) + pow(y[2] - y[14],2),1.5) + (par[0]*par[2]*(-y[0] + y[12]))/pow(pow(-y[0] + y[12],2) + pow(-y[1] + y[13],2) + pow(-y[2] + y[14],2),1.5))/2.);
    break;

  case 4 : 
    return((-((par[0]*par[1]*(y[1] - y[7]))/pow(pow(y[0] - y[6],2) + pow(y[1] - y[7],2) + pow(y[2] - y[8],2),1.5)) + (par[0]*par[1]*(-y[1] + y[7]))/pow(pow(-y[0] + y[6],2) + pow(-y[1] + y[7],2) + pow(-y[2] + y[8],2),1.5) - (par[0]*par[2]*(y[1] - y[13]))/pow(pow(y[0] - y[12],2) + pow(y[1] - y[13],2) + pow(y[2] - y[14],2),1.5) + (par[0]*par[2]*(-y[1] + y[13]))/pow(pow(-y[0] + y[12],2) + pow(-y[1] + y[13],2) + pow(-y[2] + y[14],2),1.5))/2.);
    break;

  case 5 : 
    return((-((par[0]*par[1]*(y[2] - y[8]))/pow(pow(y[0] - y[6],2) + pow(y[1] - y[7],2) + pow(y[2] - y[8],2),1.5)) + (par[0]*par[1]*(-y[2] + y[8]))/pow(pow(-y[0] + y[6],2) + pow(-y[1] + y[7],2) + pow(-y[2] + y[8],2),1.5) - (par[0]*par[2]*(y[2] - y[14]))/pow(pow(y[0] - y[12],2) + pow(y[1] - y[13],2) + pow(y[2] - y[14],2),1.5) + (par[0]*par[2]*(-y[2] + y[14]))/pow(pow(-y[0] + y[12],2) + pow(-y[1] + y[13],2) + pow(-y[2] + y[14],2),1.5))/2.);
    break;

  case 6 : 
    return(y[9]/par[1]);
    break;

  case 7 : 
    return(y[10]/par[1]);
    break;

  case 8 : 
    return(y[11]/par[1]);
    break;

  case 9 : 
    return(((par[0]*par[1]*(y[0] - y[6]))/pow(pow(y[0] - y[6],2) + pow(y[1] - y[7],2) + pow(y[2] - y[8],2),1.5) - (par[0]*par[1]*(-y[0] + y[6]))/pow(pow(-y[0] + y[6],2) + pow(-y[1] + y[7],2) + pow(-y[2] + y[8],2),1.5) - (par[1]*par[2]*(y[6] - y[12]))/pow(pow(y[6] - y[12],2) + pow(y[7] - y[13],2) + pow(y[8] - y[14],2),1.5) + (par[1]*par[2]*(-y[6] + y[12]))/pow(pow(-y[6] + y[12],2) + pow(-y[7] + y[13],2) + pow(-y[8] + y[14],2),1.5))/2.);
    break;

  case 10 : 
    return(((par[0]*par[1]*(y[1] - y[7]))/pow(pow(y[0] - y[6],2) + pow(y[1] - y[7],2) + pow(y[2] - y[8],2),1.5) - (par[0]*par[1]*(-y[1] + y[7]))/pow(pow(-y[0] + y[6],2) + pow(-y[1] + y[7],2) + pow(-y[2] + y[8],2),1.5) - (par[1]*par[2]*(y[7] - y[13]))/pow(pow(y[6] - y[12],2) + pow(y[7] - y[13],2) + pow(y[8] - y[14],2),1.5) + (par[1]*par[2]*(-y[7] + y[13]))/pow(pow(-y[6] + y[12],2) + pow(-y[7] + y[13],2) + pow(-y[8] + y[14],2),1.5))/2.);
    break;

  case 11 : 
    return(((par[0]*par[1]*(y[2] - y[8]))/pow(pow(y[0] - y[6],2) + pow(y[1] - y[7],2) + pow(y[2] - y[8],2),1.5) - (par[0]*par[1]*(-y[2] + y[8]))/pow(pow(-y[0] + y[6],2) + pow(-y[1] + y[7],2) + pow(-y[2] + y[8],2),1.5) - (par[1]*par[2]*(y[8] - y[14]))/pow(pow(y[6] - y[12],2) + pow(y[7] - y[13],2) + pow(y[8] - y[14],2),1.5) + (par[1]*par[2]*(-y[8] + y[14]))/pow(pow(-y[6] + y[12],2) + pow(-y[7] + y[13],2) + pow(-y[8] + y[14],2),1.5))/2.);
    break;

  case 12 : 
    return(y[15]/par[2]);
    break;

  case 13 : 
    return(y[16]/par[2]);
    break;

  case 14 : 
    return(y[17]/par[2]);
    break;

  case 15 : 
    return(((par[0]*par[2]*(y[0] - y[12]))/pow(pow(y[0] - y[12],2) + pow(y[1] - y[13],2) + pow(y[2] - y[14],2),1.5) + (par[1]*par[2]*(y[6] - y[12]))/pow(pow(y[6] - y[12],2) + pow(y[7] - y[13],2) + pow(y[8] - y[14],2),1.5) - (par[0]*par[2]*(-y[0] + y[12]))/pow(pow(-y[0] + y[12],2) + pow(-y[1] + y[13],2) + pow(-y[2] + y[14],2),1.5) - (par[1]*par[2]*(-y[6] + y[12]))/pow(pow(-y[6] + y[12],2) + pow(-y[7] + y[13],2) + pow(-y[8] + y[14],2),1.5))/2.);
    break;

  case 16 : 
    return(((par[0]*par[2]*(y[1] - y[13]))/pow(pow(y[0] - y[12],2) + pow(y[1] - y[13],2) + pow(y[2] - y[14],2),1.5) + (par[1]*par[2]*(y[7] - y[13]))/pow(pow(y[6] - y[12],2) + pow(y[7] - y[13],2) + pow(y[8] - y[14],2),1.5) - (par[0]*par[2]*(-y[1] + y[13]))/pow(pow(-y[0] + y[12],2) + pow(-y[1] + y[13],2) + pow(-y[2] + y[14],2),1.5) - (par[1]*par[2]*(-y[7] + y[13]))/pow(pow(-y[6] + y[12],2) + pow(-y[7] + y[13],2) + pow(-y[8] + y[14],2),1.5))/2.);
    break;

  case 17 : 
    return(((par[0]*par[2]*(y[2] - y[14]))/pow(pow(y[0] - y[12],2) + pow(y[1] - y[13],2) + pow(y[2] - y[14],2),1.5) + (par[1]*par[2]*(y[8] - y[14]))/pow(pow(y[6] - y[12],2) + pow(y[7] - y[13],2) + pow(y[8] - y[14],2),1.5) - (par[0]*par[2]*(-y[2] + y[14]))/pow(pow(-y[0] + y[12],2) + pow(-y[1] + y[13],2) + pow(-y[2] + y[14],2),1.5) - (par[1]*par[2]*(-y[8] + y[14]))/pow(pow(-y[6] + y[12],2) + pow(-y[7] + y[13],2) + pow(-y[8] + y[14],2),1.5))/2.);
    break;
}

 return 0; 
}


double evolution::jac_Nnp3d3_FaNS(double t, const double y[], double par[3], int i, int j)
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
    return(((3*par[0]*par[1]*pow(y[0] - y[6],2))/pow(pow(y[0] - y[6],2) + pow(y[1] - y[7],2) + pow(y[2] - y[8],2),2.5) - (par[0]*par[1])/pow(pow(y[0] - y[6],2) + pow(y[1] - y[7],2) + pow(y[2] - y[8],2),1.5) + (3*par[0]*par[1]*pow(-y[0] + y[6],2))/pow(pow(-y[0] + y[6],2) + pow(-y[1] + y[7],2) + pow(-y[2] + y[8],2),2.5) - (par[0]*par[1])/pow(pow(-y[0] + y[6],2) + pow(-y[1] + y[7],2) + pow(-y[2] + y[8],2),1.5) + (3*par[0]*par[2]*pow(y[0] - y[12],2))/pow(pow(y[0] - y[12],2) + pow(y[1] - y[13],2) + pow(y[2] - y[14],2),2.5) - (par[0]*par[2])/pow(pow(y[0] - y[12],2) + pow(y[1] - y[13],2) + pow(y[2] - y[14],2),1.5) + (3*par[0]*par[2]*pow(-y[0] + y[12],2))/pow(pow(-y[0] + y[12],2) + pow(-y[1] + y[13],2) + pow(-y[2] + y[14],2),2.5) - (par[0]*par[2])/pow(pow(-y[0] + y[12],2) + pow(-y[1] + y[13],2) + pow(-y[2] + y[14],2),1.5))/2.);
    break;

  case 4 : 
    return(((3*par[0]*par[1]*(y[0] - y[6])*(y[1] - y[7]))/pow(pow(y[0] - y[6],2) + pow(y[1] - y[7],2) + pow(y[2] - y[8],2),2.5) + (3*par[0]*par[1]*(-y[0] + y[6])*(-y[1] + y[7]))/pow(pow(-y[0] + y[6],2) + pow(-y[1] + y[7],2) + pow(-y[2] + y[8],2),2.5) + (3*par[0]*par[2]*(y[0] - y[12])*(y[1] - y[13]))/pow(pow(y[0] - y[12],2) + pow(y[1] - y[13],2) + pow(y[2] - y[14],2),2.5) + (3*par[0]*par[2]*(-y[0] + y[12])*(-y[1] + y[13]))/pow(pow(-y[0] + y[12],2) + pow(-y[1] + y[13],2) + pow(-y[2] + y[14],2),2.5))/2.);
    break;

  case 5 : 
    return(((3*par[0]*par[1]*(y[0] - y[6])*(y[2] - y[8]))/pow(pow(y[0] - y[6],2) + pow(y[1] - y[7],2) + pow(y[2] - y[8],2),2.5) + (3*par[0]*par[1]*(-y[0] + y[6])*(-y[2] + y[8]))/pow(pow(-y[0] + y[6],2) + pow(-y[1] + y[7],2) + pow(-y[2] + y[8],2),2.5) + (3*par[0]*par[2]*(y[0] - y[12])*(y[2] - y[14]))/pow(pow(y[0] - y[12],2) + pow(y[1] - y[13],2) + pow(y[2] - y[14],2),2.5) + (3*par[0]*par[2]*(-y[0] + y[12])*(-y[2] + y[14]))/pow(pow(-y[0] + y[12],2) + pow(-y[1] + y[13],2) + pow(-y[2] + y[14],2),2.5))/2.);
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
    return(((-3*par[0]*par[1]*pow(y[0] - y[6],2))/pow(pow(y[0] - y[6],2) + pow(y[1] - y[7],2) + pow(y[2] - y[8],2),2.5) + (par[0]*par[1])/pow(pow(y[0] - y[6],2) + pow(y[1] - y[7],2) + pow(y[2] - y[8],2),1.5) - (3*par[0]*par[1]*pow(-y[0] + y[6],2))/pow(pow(-y[0] + y[6],2) + pow(-y[1] + y[7],2) + pow(-y[2] + y[8],2),2.5) + (par[0]*par[1])/pow(pow(-y[0] + y[6],2) + pow(-y[1] + y[7],2) + pow(-y[2] + y[8],2),1.5))/2.);
    break;

  case 10 : 
    return(((-3*par[0]*par[1]*(y[0] - y[6])*(y[1] - y[7]))/pow(pow(y[0] - y[6],2) + pow(y[1] - y[7],2) + pow(y[2] - y[8],2),2.5) - (3*par[0]*par[1]*(-y[0] + y[6])*(-y[1] + y[7]))/pow(pow(-y[0] + y[6],2) + pow(-y[1] + y[7],2) + pow(-y[2] + y[8],2),2.5))/2.);
    break;

  case 11 : 
    return(((-3*par[0]*par[1]*(y[0] - y[6])*(y[2] - y[8]))/pow(pow(y[0] - y[6],2) + pow(y[1] - y[7],2) + pow(y[2] - y[8],2),2.5) - (3*par[0]*par[1]*(-y[0] + y[6])*(-y[2] + y[8]))/pow(pow(-y[0] + y[6],2) + pow(-y[1] + y[7],2) + pow(-y[2] + y[8],2),2.5))/2.);
    break;

  case 12 : 
    return(0);
    break;

  case 13 : 
    return(0);
    break;

  case 14 : 
    return(0);
    break;

  case 15 : 
    return(((-3*par[0]*par[2]*pow(y[0] - y[12],2))/pow(pow(y[0] - y[12],2) + pow(y[1] - y[13],2) + pow(y[2] - y[14],2),2.5) + (par[0]*par[2])/pow(pow(y[0] - y[12],2) + pow(y[1] - y[13],2) + pow(y[2] - y[14],2),1.5) - (3*par[0]*par[2]*pow(-y[0] + y[12],2))/pow(pow(-y[0] + y[12],2) + pow(-y[1] + y[13],2) + pow(-y[2] + y[14],2),2.5) + (par[0]*par[2])/pow(pow(-y[0] + y[12],2) + pow(-y[1] + y[13],2) + pow(-y[2] + y[14],2),1.5))/2.);
    break;

  case 16 : 
    return(((-3*par[0]*par[2]*(y[0] - y[12])*(y[1] - y[13]))/pow(pow(y[0] - y[12],2) + pow(y[1] - y[13],2) + pow(y[2] - y[14],2),2.5) - (3*par[0]*par[2]*(-y[0] + y[12])*(-y[1] + y[13]))/pow(pow(-y[0] + y[12],2) + pow(-y[1] + y[13],2) + pow(-y[2] + y[14],2),2.5))/2.);
    break;

  case 17 : 
    return(((-3*par[0]*par[2]*(y[0] - y[12])*(y[2] - y[14]))/pow(pow(y[0] - y[12],2) + pow(y[1] - y[13],2) + pow(y[2] - y[14],2),2.5) - (3*par[0]*par[2]*(-y[0] + y[12])*(-y[2] + y[14]))/pow(pow(-y[0] + y[12],2) + pow(-y[1] + y[13],2) + pow(-y[2] + y[14],2),2.5))/2.);
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
    return(((3*par[0]*par[1]*(y[0] - y[6])*(y[1] - y[7]))/pow(pow(y[0] - y[6],2) + pow(y[1] - y[7],2) + pow(y[2] - y[8],2),2.5) + (3*par[0]*par[1]*(-y[0] + y[6])*(-y[1] + y[7]))/pow(pow(-y[0] + y[6],2) + pow(-y[1] + y[7],2) + pow(-y[2] + y[8],2),2.5) + (3*par[0]*par[2]*(y[0] - y[12])*(y[1] - y[13]))/pow(pow(y[0] - y[12],2) + pow(y[1] - y[13],2) + pow(y[2] - y[14],2),2.5) + (3*par[0]*par[2]*(-y[0] + y[12])*(-y[1] + y[13]))/pow(pow(-y[0] + y[12],2) + pow(-y[1] + y[13],2) + pow(-y[2] + y[14],2),2.5))/2.);
    break;

  case 22 : 
    return(((3*par[0]*par[1]*pow(y[1] - y[7],2))/pow(pow(y[0] - y[6],2) + pow(y[1] - y[7],2) + pow(y[2] - y[8],2),2.5) - (par[0]*par[1])/pow(pow(y[0] - y[6],2) + pow(y[1] - y[7],2) + pow(y[2] - y[8],2),1.5) + (3*par[0]*par[1]*pow(-y[1] + y[7],2))/pow(pow(-y[0] + y[6],2) + pow(-y[1] + y[7],2) + pow(-y[2] + y[8],2),2.5) - (par[0]*par[1])/pow(pow(-y[0] + y[6],2) + pow(-y[1] + y[7],2) + pow(-y[2] + y[8],2),1.5) + (3*par[0]*par[2]*pow(y[1] - y[13],2))/pow(pow(y[0] - y[12],2) + pow(y[1] - y[13],2) + pow(y[2] - y[14],2),2.5) - (par[0]*par[2])/pow(pow(y[0] - y[12],2) + pow(y[1] - y[13],2) + pow(y[2] - y[14],2),1.5) + (3*par[0]*par[2]*pow(-y[1] + y[13],2))/pow(pow(-y[0] + y[12],2) + pow(-y[1] + y[13],2) + pow(-y[2] + y[14],2),2.5) - (par[0]*par[2])/pow(pow(-y[0] + y[12],2) + pow(-y[1] + y[13],2) + pow(-y[2] + y[14],2),1.5))/2.);
    break;

  case 23 : 
    return(((3*par[0]*par[1]*(y[1] - y[7])*(y[2] - y[8]))/pow(pow(y[0] - y[6],2) + pow(y[1] - y[7],2) + pow(y[2] - y[8],2),2.5) + (3*par[0]*par[1]*(-y[1] + y[7])*(-y[2] + y[8]))/pow(pow(-y[0] + y[6],2) + pow(-y[1] + y[7],2) + pow(-y[2] + y[8],2),2.5) + (3*par[0]*par[2]*(y[1] - y[13])*(y[2] - y[14]))/pow(pow(y[0] - y[12],2) + pow(y[1] - y[13],2) + pow(y[2] - y[14],2),2.5) + (3*par[0]*par[2]*(-y[1] + y[13])*(-y[2] + y[14]))/pow(pow(-y[0] + y[12],2) + pow(-y[1] + y[13],2) + pow(-y[2] + y[14],2),2.5))/2.);
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
    return(((-3*par[0]*par[1]*(y[0] - y[6])*(y[1] - y[7]))/pow(pow(y[0] - y[6],2) + pow(y[1] - y[7],2) + pow(y[2] - y[8],2),2.5) - (3*par[0]*par[1]*(-y[0] + y[6])*(-y[1] + y[7]))/pow(pow(-y[0] + y[6],2) + pow(-y[1] + y[7],2) + pow(-y[2] + y[8],2),2.5))/2.);
    break;

  case 28 : 
    return(((-3*par[0]*par[1]*pow(y[1] - y[7],2))/pow(pow(y[0] - y[6],2) + pow(y[1] - y[7],2) + pow(y[2] - y[8],2),2.5) + (par[0]*par[1])/pow(pow(y[0] - y[6],2) + pow(y[1] - y[7],2) + pow(y[2] - y[8],2),1.5) - (3*par[0]*par[1]*pow(-y[1] + y[7],2))/pow(pow(-y[0] + y[6],2) + pow(-y[1] + y[7],2) + pow(-y[2] + y[8],2),2.5) + (par[0]*par[1])/pow(pow(-y[0] + y[6],2) + pow(-y[1] + y[7],2) + pow(-y[2] + y[8],2),1.5))/2.);
    break;

  case 29 : 
    return(((-3*par[0]*par[1]*(y[1] - y[7])*(y[2] - y[8]))/pow(pow(y[0] - y[6],2) + pow(y[1] - y[7],2) + pow(y[2] - y[8],2),2.5) - (3*par[0]*par[1]*(-y[1] + y[7])*(-y[2] + y[8]))/pow(pow(-y[0] + y[6],2) + pow(-y[1] + y[7],2) + pow(-y[2] + y[8],2),2.5))/2.);
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
    return(((-3*par[0]*par[2]*(y[0] - y[12])*(y[1] - y[13]))/pow(pow(y[0] - y[12],2) + pow(y[1] - y[13],2) + pow(y[2] - y[14],2),2.5) - (3*par[0]*par[2]*(-y[0] + y[12])*(-y[1] + y[13]))/pow(pow(-y[0] + y[12],2) + pow(-y[1] + y[13],2) + pow(-y[2] + y[14],2),2.5))/2.);
    break;

  case 34 : 
    return(((-3*par[0]*par[2]*pow(y[1] - y[13],2))/pow(pow(y[0] - y[12],2) + pow(y[1] - y[13],2) + pow(y[2] - y[14],2),2.5) + (par[0]*par[2])/pow(pow(y[0] - y[12],2) + pow(y[1] - y[13],2) + pow(y[2] - y[14],2),1.5) - (3*par[0]*par[2]*pow(-y[1] + y[13],2))/pow(pow(-y[0] + y[12],2) + pow(-y[1] + y[13],2) + pow(-y[2] + y[14],2),2.5) + (par[0]*par[2])/pow(pow(-y[0] + y[12],2) + pow(-y[1] + y[13],2) + pow(-y[2] + y[14],2),1.5))/2.);
    break;

  case 35 : 
    return(((-3*par[0]*par[2]*(y[1] - y[13])*(y[2] - y[14]))/pow(pow(y[0] - y[12],2) + pow(y[1] - y[13],2) + pow(y[2] - y[14],2),2.5) - (3*par[0]*par[2]*(-y[1] + y[13])*(-y[2] + y[14]))/pow(pow(-y[0] + y[12],2) + pow(-y[1] + y[13],2) + pow(-y[2] + y[14],2),2.5))/2.);
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
    return(((3*par[0]*par[1]*(y[0] - y[6])*(y[2] - y[8]))/pow(pow(y[0] - y[6],2) + pow(y[1] - y[7],2) + pow(y[2] - y[8],2),2.5) + (3*par[0]*par[1]*(-y[0] + y[6])*(-y[2] + y[8]))/pow(pow(-y[0] + y[6],2) + pow(-y[1] + y[7],2) + pow(-y[2] + y[8],2),2.5) + (3*par[0]*par[2]*(y[0] - y[12])*(y[2] - y[14]))/pow(pow(y[0] - y[12],2) + pow(y[1] - y[13],2) + pow(y[2] - y[14],2),2.5) + (3*par[0]*par[2]*(-y[0] + y[12])*(-y[2] + y[14]))/pow(pow(-y[0] + y[12],2) + pow(-y[1] + y[13],2) + pow(-y[2] + y[14],2),2.5))/2.);
    break;

  case 40 : 
    return(((3*par[0]*par[1]*(y[1] - y[7])*(y[2] - y[8]))/pow(pow(y[0] - y[6],2) + pow(y[1] - y[7],2) + pow(y[2] - y[8],2),2.5) + (3*par[0]*par[1]*(-y[1] + y[7])*(-y[2] + y[8]))/pow(pow(-y[0] + y[6],2) + pow(-y[1] + y[7],2) + pow(-y[2] + y[8],2),2.5) + (3*par[0]*par[2]*(y[1] - y[13])*(y[2] - y[14]))/pow(pow(y[0] - y[12],2) + pow(y[1] - y[13],2) + pow(y[2] - y[14],2),2.5) + (3*par[0]*par[2]*(-y[1] + y[13])*(-y[2] + y[14]))/pow(pow(-y[0] + y[12],2) + pow(-y[1] + y[13],2) + pow(-y[2] + y[14],2),2.5))/2.);
    break;

  case 41 : 
    return((-((par[0]*par[1])/pow(pow(y[0] - y[6],2) + pow(y[1] - y[7],2) + pow(y[2] - y[8],2),1.5)) + (3*par[0]*par[1]*pow(y[2] - y[8],2))/pow(pow(y[0] - y[6],2) + pow(y[1] - y[7],2) + pow(y[2] - y[8],2),2.5) + (3*par[0]*par[1]*pow(-y[2] + y[8],2))/pow(pow(-y[0] + y[6],2) + pow(-y[1] + y[7],2) + pow(-y[2] + y[8],2),2.5) - (par[0]*par[1])/pow(pow(-y[0] + y[6],2) + pow(-y[1] + y[7],2) + pow(-y[2] + y[8],2),1.5) - (par[0]*par[2])/pow(pow(y[0] - y[12],2) + pow(y[1] - y[13],2) + pow(y[2] - y[14],2),1.5) + (3*par[0]*par[2]*pow(y[2] - y[14],2))/pow(pow(y[0] - y[12],2) + pow(y[1] - y[13],2) + pow(y[2] - y[14],2),2.5) + (3*par[0]*par[2]*pow(-y[2] + y[14],2))/pow(pow(-y[0] + y[12],2) + pow(-y[1] + y[13],2) + pow(-y[2] + y[14],2),2.5) - (par[0]*par[2])/pow(pow(-y[0] + y[12],2) + pow(-y[1] + y[13],2) + pow(-y[2] + y[14],2),1.5))/2.);
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
    return(((-3*par[0]*par[1]*(y[0] - y[6])*(y[2] - y[8]))/pow(pow(y[0] - y[6],2) + pow(y[1] - y[7],2) + pow(y[2] - y[8],2),2.5) - (3*par[0]*par[1]*(-y[0] + y[6])*(-y[2] + y[8]))/pow(pow(-y[0] + y[6],2) + pow(-y[1] + y[7],2) + pow(-y[2] + y[8],2),2.5))/2.);
    break;

  case 46 : 
    return(((-3*par[0]*par[1]*(y[1] - y[7])*(y[2] - y[8]))/pow(pow(y[0] - y[6],2) + pow(y[1] - y[7],2) + pow(y[2] - y[8],2),2.5) - (3*par[0]*par[1]*(-y[1] + y[7])*(-y[2] + y[8]))/pow(pow(-y[0] + y[6],2) + pow(-y[1] + y[7],2) + pow(-y[2] + y[8],2),2.5))/2.);
    break;

  case 47 : 
    return(((par[0]*par[1])/pow(pow(y[0] - y[6],2) + pow(y[1] - y[7],2) + pow(y[2] - y[8],2),1.5) - (3*par[0]*par[1]*pow(y[2] - y[8],2))/pow(pow(y[0] - y[6],2) + pow(y[1] - y[7],2) + pow(y[2] - y[8],2),2.5) - (3*par[0]*par[1]*pow(-y[2] + y[8],2))/pow(pow(-y[0] + y[6],2) + pow(-y[1] + y[7],2) + pow(-y[2] + y[8],2),2.5) + (par[0]*par[1])/pow(pow(-y[0] + y[6],2) + pow(-y[1] + y[7],2) + pow(-y[2] + y[8],2),1.5))/2.);
    break;

  case 48 : 
    return(0);
    break;

  case 49 : 
    return(0);
    break;

  case 50 : 
    return(0);
    break;

  case 51 : 
    return(((-3*par[0]*par[2]*(y[0] - y[12])*(y[2] - y[14]))/pow(pow(y[0] - y[12],2) + pow(y[1] - y[13],2) + pow(y[2] - y[14],2),2.5) - (3*par[0]*par[2]*(-y[0] + y[12])*(-y[2] + y[14]))/pow(pow(-y[0] + y[12],2) + pow(-y[1] + y[13],2) + pow(-y[2] + y[14],2),2.5))/2.);
    break;

  case 52 : 
    return(((-3*par[0]*par[2]*(y[1] - y[13])*(y[2] - y[14]))/pow(pow(y[0] - y[12],2) + pow(y[1] - y[13],2) + pow(y[2] - y[14],2),2.5) - (3*par[0]*par[2]*(-y[1] + y[13])*(-y[2] + y[14]))/pow(pow(-y[0] + y[12],2) + pow(-y[1] + y[13],2) + pow(-y[2] + y[14],2),2.5))/2.);
    break;

  case 53 : 
    return(((par[0]*par[2])/pow(pow(y[0] - y[12],2) + pow(y[1] - y[13],2) + pow(y[2] - y[14],2),1.5) - (3*par[0]*par[2]*pow(y[2] - y[14],2))/pow(pow(y[0] - y[12],2) + pow(y[1] - y[13],2) + pow(y[2] - y[14],2),2.5) - (3*par[0]*par[2]*pow(-y[2] + y[14],2))/pow(pow(-y[0] + y[12],2) + pow(-y[1] + y[13],2) + pow(-y[2] + y[14],2),2.5) + (par[0]*par[2])/pow(pow(-y[0] + y[12],2) + pow(-y[1] + y[13],2) + pow(-y[2] + y[14],2),1.5))/2.);
    break;

  case 54 : 
    return(1/par[0]);
    break;

  case 55 : 
    return(0);
    break;

  case 56 : 
    return(0);
    break;

  case 57 : 
    return(0);
    break;

  case 58 : 
    return(0);
    break;

  case 59 : 
    return(0);
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
    return(0);
    break;

  case 67 : 
    return(0);
    break;

  case 68 : 
    return(0);
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
    return(1/par[0]);
    break;

  case 74 : 
    return(0);
    break;

  case 75 : 
    return(0);
    break;

  case 76 : 
    return(0);
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
    return(0);
    break;

  case 90 : 
    return(0);
    break;

  case 91 : 
    return(0);
    break;

  case 92 : 
    return(1/par[0]);
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
    return(0);
    break;

  case 110 : 
    return(0);
    break;

  case 111 : 
    return(((-3*par[0]*par[1]*pow(y[0] - y[6],2))/pow(pow(y[0] - y[6],2) + pow(y[1] - y[7],2) + pow(y[2] - y[8],2),2.5) + (par[0]*par[1])/pow(pow(y[0] - y[6],2) + pow(y[1] - y[7],2) + pow(y[2] - y[8],2),1.5) - (3*par[0]*par[1]*pow(-y[0] + y[6],2))/pow(pow(-y[0] + y[6],2) + pow(-y[1] + y[7],2) + pow(-y[2] + y[8],2),2.5) + (par[0]*par[1])/pow(pow(-y[0] + y[6],2) + pow(-y[1] + y[7],2) + pow(-y[2] + y[8],2),1.5))/2.);
    break;

  case 112 : 
    return(((-3*par[0]*par[1]*(y[0] - y[6])*(y[1] - y[7]))/pow(pow(y[0] - y[6],2) + pow(y[1] - y[7],2) + pow(y[2] - y[8],2),2.5) - (3*par[0]*par[1]*(-y[0] + y[6])*(-y[1] + y[7]))/pow(pow(-y[0] + y[6],2) + pow(-y[1] + y[7],2) + pow(-y[2] + y[8],2),2.5))/2.);
    break;

  case 113 : 
    return(((-3*par[0]*par[1]*(y[0] - y[6])*(y[2] - y[8]))/pow(pow(y[0] - y[6],2) + pow(y[1] - y[7],2) + pow(y[2] - y[8],2),2.5) - (3*par[0]*par[1]*(-y[0] + y[6])*(-y[2] + y[8]))/pow(pow(-y[0] + y[6],2) + pow(-y[1] + y[7],2) + pow(-y[2] + y[8],2),2.5))/2.);
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
    return(((3*par[0]*par[1]*pow(y[0] - y[6],2))/pow(pow(y[0] - y[6],2) + pow(y[1] - y[7],2) + pow(y[2] - y[8],2),2.5) - (par[0]*par[1])/pow(pow(y[0] - y[6],2) + pow(y[1] - y[7],2) + pow(y[2] - y[8],2),1.5) + (3*par[0]*par[1]*pow(-y[0] + y[6],2))/pow(pow(-y[0] + y[6],2) + pow(-y[1] + y[7],2) + pow(-y[2] + y[8],2),2.5) - (par[0]*par[1])/pow(pow(-y[0] + y[6],2) + pow(-y[1] + y[7],2) + pow(-y[2] + y[8],2),1.5) + (3*par[1]*par[2]*pow(y[6] - y[12],2))/pow(pow(y[6] - y[12],2) + pow(y[7] - y[13],2) + pow(y[8] - y[14],2),2.5) - (par[1]*par[2])/pow(pow(y[6] - y[12],2) + pow(y[7] - y[13],2) + pow(y[8] - y[14],2),1.5) + (3*par[1]*par[2]*pow(-y[6] + y[12],2))/pow(pow(-y[6] + y[12],2) + pow(-y[7] + y[13],2) + pow(-y[8] + y[14],2),2.5) - (par[1]*par[2])/pow(pow(-y[6] + y[12],2) + pow(-y[7] + y[13],2) + pow(-y[8] + y[14],2),1.5))/2.);
    break;

  case 118 : 
    return(((3*par[0]*par[1]*(y[0] - y[6])*(y[1] - y[7]))/pow(pow(y[0] - y[6],2) + pow(y[1] - y[7],2) + pow(y[2] - y[8],2),2.5) + (3*par[0]*par[1]*(-y[0] + y[6])*(-y[1] + y[7]))/pow(pow(-y[0] + y[6],2) + pow(-y[1] + y[7],2) + pow(-y[2] + y[8],2),2.5) + (3*par[1]*par[2]*(y[6] - y[12])*(y[7] - y[13]))/pow(pow(y[6] - y[12],2) + pow(y[7] - y[13],2) + pow(y[8] - y[14],2),2.5) + (3*par[1]*par[2]*(-y[6] + y[12])*(-y[7] + y[13]))/pow(pow(-y[6] + y[12],2) + pow(-y[7] + y[13],2) + pow(-y[8] + y[14],2),2.5))/2.);
    break;

  case 119 : 
    return(((3*par[0]*par[1]*(y[0] - y[6])*(y[2] - y[8]))/pow(pow(y[0] - y[6],2) + pow(y[1] - y[7],2) + pow(y[2] - y[8],2),2.5) + (3*par[0]*par[1]*(-y[0] + y[6])*(-y[2] + y[8]))/pow(pow(-y[0] + y[6],2) + pow(-y[1] + y[7],2) + pow(-y[2] + y[8],2),2.5) + (3*par[1]*par[2]*(y[6] - y[12])*(y[8] - y[14]))/pow(pow(y[6] - y[12],2) + pow(y[7] - y[13],2) + pow(y[8] - y[14],2),2.5) + (3*par[1]*par[2]*(-y[6] + y[12])*(-y[8] + y[14]))/pow(pow(-y[6] + y[12],2) + pow(-y[7] + y[13],2) + pow(-y[8] + y[14],2),2.5))/2.);
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
    return(((-3*par[1]*par[2]*pow(y[6] - y[12],2))/pow(pow(y[6] - y[12],2) + pow(y[7] - y[13],2) + pow(y[8] - y[14],2),2.5) + (par[1]*par[2])/pow(pow(y[6] - y[12],2) + pow(y[7] - y[13],2) + pow(y[8] - y[14],2),1.5) - (3*par[1]*par[2]*pow(-y[6] + y[12],2))/pow(pow(-y[6] + y[12],2) + pow(-y[7] + y[13],2) + pow(-y[8] + y[14],2),2.5) + (par[1]*par[2])/pow(pow(-y[6] + y[12],2) + pow(-y[7] + y[13],2) + pow(-y[8] + y[14],2),1.5))/2.);
    break;

  case 124 : 
    return(((-3*par[1]*par[2]*(y[6] - y[12])*(y[7] - y[13]))/pow(pow(y[6] - y[12],2) + pow(y[7] - y[13],2) + pow(y[8] - y[14],2),2.5) - (3*par[1]*par[2]*(-y[6] + y[12])*(-y[7] + y[13]))/pow(pow(-y[6] + y[12],2) + pow(-y[7] + y[13],2) + pow(-y[8] + y[14],2),2.5))/2.);
    break;

  case 125 : 
    return(((-3*par[1]*par[2]*(y[6] - y[12])*(y[8] - y[14]))/pow(pow(y[6] - y[12],2) + pow(y[7] - y[13],2) + pow(y[8] - y[14],2),2.5) - (3*par[1]*par[2]*(-y[6] + y[12])*(-y[8] + y[14]))/pow(pow(-y[6] + y[12],2) + pow(-y[7] + y[13],2) + pow(-y[8] + y[14],2),2.5))/2.);
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
    return(((-3*par[0]*par[1]*(y[0] - y[6])*(y[1] - y[7]))/pow(pow(y[0] - y[6],2) + pow(y[1] - y[7],2) + pow(y[2] - y[8],2),2.5) - (3*par[0]*par[1]*(-y[0] + y[6])*(-y[1] + y[7]))/pow(pow(-y[0] + y[6],2) + pow(-y[1] + y[7],2) + pow(-y[2] + y[8],2),2.5))/2.);
    break;

  case 130 : 
    return(((-3*par[0]*par[1]*pow(y[1] - y[7],2))/pow(pow(y[0] - y[6],2) + pow(y[1] - y[7],2) + pow(y[2] - y[8],2),2.5) + (par[0]*par[1])/pow(pow(y[0] - y[6],2) + pow(y[1] - y[7],2) + pow(y[2] - y[8],2),1.5) - (3*par[0]*par[1]*pow(-y[1] + y[7],2))/pow(pow(-y[0] + y[6],2) + pow(-y[1] + y[7],2) + pow(-y[2] + y[8],2),2.5) + (par[0]*par[1])/pow(pow(-y[0] + y[6],2) + pow(-y[1] + y[7],2) + pow(-y[2] + y[8],2),1.5))/2.);
    break;

  case 131 : 
    return(((-3*par[0]*par[1]*(y[1] - y[7])*(y[2] - y[8]))/pow(pow(y[0] - y[6],2) + pow(y[1] - y[7],2) + pow(y[2] - y[8],2),2.5) - (3*par[0]*par[1]*(-y[1] + y[7])*(-y[2] + y[8]))/pow(pow(-y[0] + y[6],2) + pow(-y[1] + y[7],2) + pow(-y[2] + y[8],2),2.5))/2.);
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
    return(((3*par[0]*par[1]*(y[0] - y[6])*(y[1] - y[7]))/pow(pow(y[0] - y[6],2) + pow(y[1] - y[7],2) + pow(y[2] - y[8],2),2.5) + (3*par[0]*par[1]*(-y[0] + y[6])*(-y[1] + y[7]))/pow(pow(-y[0] + y[6],2) + pow(-y[1] + y[7],2) + pow(-y[2] + y[8],2),2.5) + (3*par[1]*par[2]*(y[6] - y[12])*(y[7] - y[13]))/pow(pow(y[6] - y[12],2) + pow(y[7] - y[13],2) + pow(y[8] - y[14],2),2.5) + (3*par[1]*par[2]*(-y[6] + y[12])*(-y[7] + y[13]))/pow(pow(-y[6] + y[12],2) + pow(-y[7] + y[13],2) + pow(-y[8] + y[14],2),2.5))/2.);
    break;

  case 136 : 
    return(((3*par[0]*par[1]*pow(y[1] - y[7],2))/pow(pow(y[0] - y[6],2) + pow(y[1] - y[7],2) + pow(y[2] - y[8],2),2.5) - (par[0]*par[1])/pow(pow(y[0] - y[6],2) + pow(y[1] - y[7],2) + pow(y[2] - y[8],2),1.5) + (3*par[0]*par[1]*pow(-y[1] + y[7],2))/pow(pow(-y[0] + y[6],2) + pow(-y[1] + y[7],2) + pow(-y[2] + y[8],2),2.5) - (par[0]*par[1])/pow(pow(-y[0] + y[6],2) + pow(-y[1] + y[7],2) + pow(-y[2] + y[8],2),1.5) + (3*par[1]*par[2]*pow(y[7] - y[13],2))/pow(pow(y[6] - y[12],2) + pow(y[7] - y[13],2) + pow(y[8] - y[14],2),2.5) - (par[1]*par[2])/pow(pow(y[6] - y[12],2) + pow(y[7] - y[13],2) + pow(y[8] - y[14],2),1.5) + (3*par[1]*par[2]*pow(-y[7] + y[13],2))/pow(pow(-y[6] + y[12],2) + pow(-y[7] + y[13],2) + pow(-y[8] + y[14],2),2.5) - (par[1]*par[2])/pow(pow(-y[6] + y[12],2) + pow(-y[7] + y[13],2) + pow(-y[8] + y[14],2),1.5))/2.);
    break;

  case 137 : 
    return(((3*par[0]*par[1]*(y[1] - y[7])*(y[2] - y[8]))/pow(pow(y[0] - y[6],2) + pow(y[1] - y[7],2) + pow(y[2] - y[8],2),2.5) + (3*par[0]*par[1]*(-y[1] + y[7])*(-y[2] + y[8]))/pow(pow(-y[0] + y[6],2) + pow(-y[1] + y[7],2) + pow(-y[2] + y[8],2),2.5) + (3*par[1]*par[2]*(y[7] - y[13])*(y[8] - y[14]))/pow(pow(y[6] - y[12],2) + pow(y[7] - y[13],2) + pow(y[8] - y[14],2),2.5) + (3*par[1]*par[2]*(-y[7] + y[13])*(-y[8] + y[14]))/pow(pow(-y[6] + y[12],2) + pow(-y[7] + y[13],2) + pow(-y[8] + y[14],2),2.5))/2.);
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
    return(((-3*par[1]*par[2]*(y[6] - y[12])*(y[7] - y[13]))/pow(pow(y[6] - y[12],2) + pow(y[7] - y[13],2) + pow(y[8] - y[14],2),2.5) - (3*par[1]*par[2]*(-y[6] + y[12])*(-y[7] + y[13]))/pow(pow(-y[6] + y[12],2) + pow(-y[7] + y[13],2) + pow(-y[8] + y[14],2),2.5))/2.);
    break;

  case 142 : 
    return(((-3*par[1]*par[2]*pow(y[7] - y[13],2))/pow(pow(y[6] - y[12],2) + pow(y[7] - y[13],2) + pow(y[8] - y[14],2),2.5) + (par[1]*par[2])/pow(pow(y[6] - y[12],2) + pow(y[7] - y[13],2) + pow(y[8] - y[14],2),1.5) - (3*par[1]*par[2]*pow(-y[7] + y[13],2))/pow(pow(-y[6] + y[12],2) + pow(-y[7] + y[13],2) + pow(-y[8] + y[14],2),2.5) + (par[1]*par[2])/pow(pow(-y[6] + y[12],2) + pow(-y[7] + y[13],2) + pow(-y[8] + y[14],2),1.5))/2.);
    break;

  case 143 : 
    return(((-3*par[1]*par[2]*(y[7] - y[13])*(y[8] - y[14]))/pow(pow(y[6] - y[12],2) + pow(y[7] - y[13],2) + pow(y[8] - y[14],2),2.5) - (3*par[1]*par[2]*(-y[7] + y[13])*(-y[8] + y[14]))/pow(pow(-y[6] + y[12],2) + pow(-y[7] + y[13],2) + pow(-y[8] + y[14],2),2.5))/2.);
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
    return(((-3*par[0]*par[1]*(y[0] - y[6])*(y[2] - y[8]))/pow(pow(y[0] - y[6],2) + pow(y[1] - y[7],2) + pow(y[2] - y[8],2),2.5) - (3*par[0]*par[1]*(-y[0] + y[6])*(-y[2] + y[8]))/pow(pow(-y[0] + y[6],2) + pow(-y[1] + y[7],2) + pow(-y[2] + y[8],2),2.5))/2.);
    break;

  case 148 : 
    return(((-3*par[0]*par[1]*(y[1] - y[7])*(y[2] - y[8]))/pow(pow(y[0] - y[6],2) + pow(y[1] - y[7],2) + pow(y[2] - y[8],2),2.5) - (3*par[0]*par[1]*(-y[1] + y[7])*(-y[2] + y[8]))/pow(pow(-y[0] + y[6],2) + pow(-y[1] + y[7],2) + pow(-y[2] + y[8],2),2.5))/2.);
    break;

  case 149 : 
    return(((par[0]*par[1])/pow(pow(y[0] - y[6],2) + pow(y[1] - y[7],2) + pow(y[2] - y[8],2),1.5) - (3*par[0]*par[1]*pow(y[2] - y[8],2))/pow(pow(y[0] - y[6],2) + pow(y[1] - y[7],2) + pow(y[2] - y[8],2),2.5) - (3*par[0]*par[1]*pow(-y[2] + y[8],2))/pow(pow(-y[0] + y[6],2) + pow(-y[1] + y[7],2) + pow(-y[2] + y[8],2),2.5) + (par[0]*par[1])/pow(pow(-y[0] + y[6],2) + pow(-y[1] + y[7],2) + pow(-y[2] + y[8],2),1.5))/2.);
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
    return(((3*par[0]*par[1]*(y[0] - y[6])*(y[2] - y[8]))/pow(pow(y[0] - y[6],2) + pow(y[1] - y[7],2) + pow(y[2] - y[8],2),2.5) + (3*par[0]*par[1]*(-y[0] + y[6])*(-y[2] + y[8]))/pow(pow(-y[0] + y[6],2) + pow(-y[1] + y[7],2) + pow(-y[2] + y[8],2),2.5) + (3*par[1]*par[2]*(y[6] - y[12])*(y[8] - y[14]))/pow(pow(y[6] - y[12],2) + pow(y[7] - y[13],2) + pow(y[8] - y[14],2),2.5) + (3*par[1]*par[2]*(-y[6] + y[12])*(-y[8] + y[14]))/pow(pow(-y[6] + y[12],2) + pow(-y[7] + y[13],2) + pow(-y[8] + y[14],2),2.5))/2.);
    break;

  case 154 : 
    return(((3*par[0]*par[1]*(y[1] - y[7])*(y[2] - y[8]))/pow(pow(y[0] - y[6],2) + pow(y[1] - y[7],2) + pow(y[2] - y[8],2),2.5) + (3*par[0]*par[1]*(-y[1] + y[7])*(-y[2] + y[8]))/pow(pow(-y[0] + y[6],2) + pow(-y[1] + y[7],2) + pow(-y[2] + y[8],2),2.5) + (3*par[1]*par[2]*(y[7] - y[13])*(y[8] - y[14]))/pow(pow(y[6] - y[12],2) + pow(y[7] - y[13],2) + pow(y[8] - y[14],2),2.5) + (3*par[1]*par[2]*(-y[7] + y[13])*(-y[8] + y[14]))/pow(pow(-y[6] + y[12],2) + pow(-y[7] + y[13],2) + pow(-y[8] + y[14],2),2.5))/2.);
    break;

  case 155 : 
    return((-((par[0]*par[1])/pow(pow(y[0] - y[6],2) + pow(y[1] - y[7],2) + pow(y[2] - y[8],2),1.5)) + (3*par[0]*par[1]*pow(y[2] - y[8],2))/pow(pow(y[0] - y[6],2) + pow(y[1] - y[7],2) + pow(y[2] - y[8],2),2.5) + (3*par[0]*par[1]*pow(-y[2] + y[8],2))/pow(pow(-y[0] + y[6],2) + pow(-y[1] + y[7],2) + pow(-y[2] + y[8],2),2.5) - (par[0]*par[1])/pow(pow(-y[0] + y[6],2) + pow(-y[1] + y[7],2) + pow(-y[2] + y[8],2),1.5) - (par[1]*par[2])/pow(pow(y[6] - y[12],2) + pow(y[7] - y[13],2) + pow(y[8] - y[14],2),1.5) + (3*par[1]*par[2]*pow(y[8] - y[14],2))/pow(pow(y[6] - y[12],2) + pow(y[7] - y[13],2) + pow(y[8] - y[14],2),2.5) + (3*par[1]*par[2]*pow(-y[8] + y[14],2))/pow(pow(-y[6] + y[12],2) + pow(-y[7] + y[13],2) + pow(-y[8] + y[14],2),2.5) - (par[1]*par[2])/pow(pow(-y[6] + y[12],2) + pow(-y[7] + y[13],2) + pow(-y[8] + y[14],2),1.5))/2.);
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
    return(((-3*par[1]*par[2]*(y[6] - y[12])*(y[8] - y[14]))/pow(pow(y[6] - y[12],2) + pow(y[7] - y[13],2) + pow(y[8] - y[14],2),2.5) - (3*par[1]*par[2]*(-y[6] + y[12])*(-y[8] + y[14]))/pow(pow(-y[6] + y[12],2) + pow(-y[7] + y[13],2) + pow(-y[8] + y[14],2),2.5))/2.);
    break;

  case 160 : 
    return(((-3*par[1]*par[2]*(y[7] - y[13])*(y[8] - y[14]))/pow(pow(y[6] - y[12],2) + pow(y[7] - y[13],2) + pow(y[8] - y[14],2),2.5) - (3*par[1]*par[2]*(-y[7] + y[13])*(-y[8] + y[14]))/pow(pow(-y[6] + y[12],2) + pow(-y[7] + y[13],2) + pow(-y[8] + y[14],2),2.5))/2.);
    break;

  case 161 : 
    return(((par[1]*par[2])/pow(pow(y[6] - y[12],2) + pow(y[7] - y[13],2) + pow(y[8] - y[14],2),1.5) - (3*par[1]*par[2]*pow(y[8] - y[14],2))/pow(pow(y[6] - y[12],2) + pow(y[7] - y[13],2) + pow(y[8] - y[14],2),2.5) - (3*par[1]*par[2]*pow(-y[8] + y[14],2))/pow(pow(-y[6] + y[12],2) + pow(-y[7] + y[13],2) + pow(-y[8] + y[14],2),2.5) + (par[1]*par[2])/pow(pow(-y[6] + y[12],2) + pow(-y[7] + y[13],2) + pow(-y[8] + y[14],2),1.5))/2.);
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
    return(1/par[1]);
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
    return(1/par[1]);
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
    return(1/par[1]);
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
    return(((-3*par[0]*par[2]*pow(y[0] - y[12],2))/pow(pow(y[0] - y[12],2) + pow(y[1] - y[13],2) + pow(y[2] - y[14],2),2.5) + (par[0]*par[2])/pow(pow(y[0] - y[12],2) + pow(y[1] - y[13],2) + pow(y[2] - y[14],2),1.5) - (3*par[0]*par[2]*pow(-y[0] + y[12],2))/pow(pow(-y[0] + y[12],2) + pow(-y[1] + y[13],2) + pow(-y[2] + y[14],2),2.5) + (par[0]*par[2])/pow(pow(-y[0] + y[12],2) + pow(-y[1] + y[13],2) + pow(-y[2] + y[14],2),1.5))/2.);
    break;

  case 220 : 
    return(((-3*par[0]*par[2]*(y[0] - y[12])*(y[1] - y[13]))/pow(pow(y[0] - y[12],2) + pow(y[1] - y[13],2) + pow(y[2] - y[14],2),2.5) - (3*par[0]*par[2]*(-y[0] + y[12])*(-y[1] + y[13]))/pow(pow(-y[0] + y[12],2) + pow(-y[1] + y[13],2) + pow(-y[2] + y[14],2),2.5))/2.);
    break;

  case 221 : 
    return(((-3*par[0]*par[2]*(y[0] - y[12])*(y[2] - y[14]))/pow(pow(y[0] - y[12],2) + pow(y[1] - y[13],2) + pow(y[2] - y[14],2),2.5) - (3*par[0]*par[2]*(-y[0] + y[12])*(-y[2] + y[14]))/pow(pow(-y[0] + y[12],2) + pow(-y[1] + y[13],2) + pow(-y[2] + y[14],2),2.5))/2.);
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
    return(((-3*par[1]*par[2]*pow(y[6] - y[12],2))/pow(pow(y[6] - y[12],2) + pow(y[7] - y[13],2) + pow(y[8] - y[14],2),2.5) + (par[1]*par[2])/pow(pow(y[6] - y[12],2) + pow(y[7] - y[13],2) + pow(y[8] - y[14],2),1.5) - (3*par[1]*par[2]*pow(-y[6] + y[12],2))/pow(pow(-y[6] + y[12],2) + pow(-y[7] + y[13],2) + pow(-y[8] + y[14],2),2.5) + (par[1]*par[2])/pow(pow(-y[6] + y[12],2) + pow(-y[7] + y[13],2) + pow(-y[8] + y[14],2),1.5))/2.);
    break;

  case 226 : 
    return(((-3*par[1]*par[2]*(y[6] - y[12])*(y[7] - y[13]))/pow(pow(y[6] - y[12],2) + pow(y[7] - y[13],2) + pow(y[8] - y[14],2),2.5) - (3*par[1]*par[2]*(-y[6] + y[12])*(-y[7] + y[13]))/pow(pow(-y[6] + y[12],2) + pow(-y[7] + y[13],2) + pow(-y[8] + y[14],2),2.5))/2.);
    break;

  case 227 : 
    return(((-3*par[1]*par[2]*(y[6] - y[12])*(y[8] - y[14]))/pow(pow(y[6] - y[12],2) + pow(y[7] - y[13],2) + pow(y[8] - y[14],2),2.5) - (3*par[1]*par[2]*(-y[6] + y[12])*(-y[8] + y[14]))/pow(pow(-y[6] + y[12],2) + pow(-y[7] + y[13],2) + pow(-y[8] + y[14],2),2.5))/2.);
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
    return(((3*par[0]*par[2]*pow(y[0] - y[12],2))/pow(pow(y[0] - y[12],2) + pow(y[1] - y[13],2) + pow(y[2] - y[14],2),2.5) - (par[0]*par[2])/pow(pow(y[0] - y[12],2) + pow(y[1] - y[13],2) + pow(y[2] - y[14],2),1.5) + (3*par[1]*par[2]*pow(y[6] - y[12],2))/pow(pow(y[6] - y[12],2) + pow(y[7] - y[13],2) + pow(y[8] - y[14],2),2.5) - (par[1]*par[2])/pow(pow(y[6] - y[12],2) + pow(y[7] - y[13],2) + pow(y[8] - y[14],2),1.5) + (3*par[0]*par[2]*pow(-y[0] + y[12],2))/pow(pow(-y[0] + y[12],2) + pow(-y[1] + y[13],2) + pow(-y[2] + y[14],2),2.5) - (par[0]*par[2])/pow(pow(-y[0] + y[12],2) + pow(-y[1] + y[13],2) + pow(-y[2] + y[14],2),1.5) + (3*par[1]*par[2]*pow(-y[6] + y[12],2))/pow(pow(-y[6] + y[12],2) + pow(-y[7] + y[13],2) + pow(-y[8] + y[14],2),2.5) - (par[1]*par[2])/pow(pow(-y[6] + y[12],2) + pow(-y[7] + y[13],2) + pow(-y[8] + y[14],2),1.5))/2.);
    break;

  case 232 : 
    return(((3*par[0]*par[2]*(y[0] - y[12])*(y[1] - y[13]))/pow(pow(y[0] - y[12],2) + pow(y[1] - y[13],2) + pow(y[2] - y[14],2),2.5) + (3*par[1]*par[2]*(y[6] - y[12])*(y[7] - y[13]))/pow(pow(y[6] - y[12],2) + pow(y[7] - y[13],2) + pow(y[8] - y[14],2),2.5) + (3*par[0]*par[2]*(-y[0] + y[12])*(-y[1] + y[13]))/pow(pow(-y[0] + y[12],2) + pow(-y[1] + y[13],2) + pow(-y[2] + y[14],2),2.5) + (3*par[1]*par[2]*(-y[6] + y[12])*(-y[7] + y[13]))/pow(pow(-y[6] + y[12],2) + pow(-y[7] + y[13],2) + pow(-y[8] + y[14],2),2.5))/2.);
    break;

  case 233 : 
    return(((3*par[0]*par[2]*(y[0] - y[12])*(y[2] - y[14]))/pow(pow(y[0] - y[12],2) + pow(y[1] - y[13],2) + pow(y[2] - y[14],2),2.5) + (3*par[1]*par[2]*(y[6] - y[12])*(y[8] - y[14]))/pow(pow(y[6] - y[12],2) + pow(y[7] - y[13],2) + pow(y[8] - y[14],2),2.5) + (3*par[0]*par[2]*(-y[0] + y[12])*(-y[2] + y[14]))/pow(pow(-y[0] + y[12],2) + pow(-y[1] + y[13],2) + pow(-y[2] + y[14],2),2.5) + (3*par[1]*par[2]*(-y[6] + y[12])*(-y[8] + y[14]))/pow(pow(-y[6] + y[12],2) + pow(-y[7] + y[13],2) + pow(-y[8] + y[14],2),2.5))/2.);
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
    return(((-3*par[0]*par[2]*(y[0] - y[12])*(y[1] - y[13]))/pow(pow(y[0] - y[12],2) + pow(y[1] - y[13],2) + pow(y[2] - y[14],2),2.5) - (3*par[0]*par[2]*(-y[0] + y[12])*(-y[1] + y[13]))/pow(pow(-y[0] + y[12],2) + pow(-y[1] + y[13],2) + pow(-y[2] + y[14],2),2.5))/2.);
    break;

  case 238 : 
    return(((-3*par[0]*par[2]*pow(y[1] - y[13],2))/pow(pow(y[0] - y[12],2) + pow(y[1] - y[13],2) + pow(y[2] - y[14],2),2.5) + (par[0]*par[2])/pow(pow(y[0] - y[12],2) + pow(y[1] - y[13],2) + pow(y[2] - y[14],2),1.5) - (3*par[0]*par[2]*pow(-y[1] + y[13],2))/pow(pow(-y[0] + y[12],2) + pow(-y[1] + y[13],2) + pow(-y[2] + y[14],2),2.5) + (par[0]*par[2])/pow(pow(-y[0] + y[12],2) + pow(-y[1] + y[13],2) + pow(-y[2] + y[14],2),1.5))/2.);
    break;

  case 239 : 
    return(((-3*par[0]*par[2]*(y[1] - y[13])*(y[2] - y[14]))/pow(pow(y[0] - y[12],2) + pow(y[1] - y[13],2) + pow(y[2] - y[14],2),2.5) - (3*par[0]*par[2]*(-y[1] + y[13])*(-y[2] + y[14]))/pow(pow(-y[0] + y[12],2) + pow(-y[1] + y[13],2) + pow(-y[2] + y[14],2),2.5))/2.);
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
    return(((-3*par[1]*par[2]*(y[6] - y[12])*(y[7] - y[13]))/pow(pow(y[6] - y[12],2) + pow(y[7] - y[13],2) + pow(y[8] - y[14],2),2.5) - (3*par[1]*par[2]*(-y[6] + y[12])*(-y[7] + y[13]))/pow(pow(-y[6] + y[12],2) + pow(-y[7] + y[13],2) + pow(-y[8] + y[14],2),2.5))/2.);
    break;

  case 244 : 
    return(((-3*par[1]*par[2]*pow(y[7] - y[13],2))/pow(pow(y[6] - y[12],2) + pow(y[7] - y[13],2) + pow(y[8] - y[14],2),2.5) + (par[1]*par[2])/pow(pow(y[6] - y[12],2) + pow(y[7] - y[13],2) + pow(y[8] - y[14],2),1.5) - (3*par[1]*par[2]*pow(-y[7] + y[13],2))/pow(pow(-y[6] + y[12],2) + pow(-y[7] + y[13],2) + pow(-y[8] + y[14],2),2.5) + (par[1]*par[2])/pow(pow(-y[6] + y[12],2) + pow(-y[7] + y[13],2) + pow(-y[8] + y[14],2),1.5))/2.);
    break;

  case 245 : 
    return(((-3*par[1]*par[2]*(y[7] - y[13])*(y[8] - y[14]))/pow(pow(y[6] - y[12],2) + pow(y[7] - y[13],2) + pow(y[8] - y[14],2),2.5) - (3*par[1]*par[2]*(-y[7] + y[13])*(-y[8] + y[14]))/pow(pow(-y[6] + y[12],2) + pow(-y[7] + y[13],2) + pow(-y[8] + y[14],2),2.5))/2.);
    break;

  case 246 : 
    return(0);
    break;

  case 247 : 
    return(0);
    break;

  case 248 : 
    return(0);
    break;

  case 249 : 
    return(((3*par[0]*par[2]*(y[0] - y[12])*(y[1] - y[13]))/pow(pow(y[0] - y[12],2) + pow(y[1] - y[13],2) + pow(y[2] - y[14],2),2.5) + (3*par[1]*par[2]*(y[6] - y[12])*(y[7] - y[13]))/pow(pow(y[6] - y[12],2) + pow(y[7] - y[13],2) + pow(y[8] - y[14],2),2.5) + (3*par[0]*par[2]*(-y[0] + y[12])*(-y[1] + y[13]))/pow(pow(-y[0] + y[12],2) + pow(-y[1] + y[13],2) + pow(-y[2] + y[14],2),2.5) + (3*par[1]*par[2]*(-y[6] + y[12])*(-y[7] + y[13]))/pow(pow(-y[6] + y[12],2) + pow(-y[7] + y[13],2) + pow(-y[8] + y[14],2),2.5))/2.);
    break;

  case 250 : 
    return(((3*par[0]*par[2]*pow(y[1] - y[13],2))/pow(pow(y[0] - y[12],2) + pow(y[1] - y[13],2) + pow(y[2] - y[14],2),2.5) - (par[0]*par[2])/pow(pow(y[0] - y[12],2) + pow(y[1] - y[13],2) + pow(y[2] - y[14],2),1.5) + (3*par[1]*par[2]*pow(y[7] - y[13],2))/pow(pow(y[6] - y[12],2) + pow(y[7] - y[13],2) + pow(y[8] - y[14],2),2.5) - (par[1]*par[2])/pow(pow(y[6] - y[12],2) + pow(y[7] - y[13],2) + pow(y[8] - y[14],2),1.5) + (3*par[0]*par[2]*pow(-y[1] + y[13],2))/pow(pow(-y[0] + y[12],2) + pow(-y[1] + y[13],2) + pow(-y[2] + y[14],2),2.5) - (par[0]*par[2])/pow(pow(-y[0] + y[12],2) + pow(-y[1] + y[13],2) + pow(-y[2] + y[14],2),1.5) + (3*par[1]*par[2]*pow(-y[7] + y[13],2))/pow(pow(-y[6] + y[12],2) + pow(-y[7] + y[13],2) + pow(-y[8] + y[14],2),2.5) - (par[1]*par[2])/pow(pow(-y[6] + y[12],2) + pow(-y[7] + y[13],2) + pow(-y[8] + y[14],2),1.5))/2.);
    break;

  case 251 : 
    return(((3*par[0]*par[2]*(y[1] - y[13])*(y[2] - y[14]))/pow(pow(y[0] - y[12],2) + pow(y[1] - y[13],2) + pow(y[2] - y[14],2),2.5) + (3*par[1]*par[2]*(y[7] - y[13])*(y[8] - y[14]))/pow(pow(y[6] - y[12],2) + pow(y[7] - y[13],2) + pow(y[8] - y[14],2),2.5) + (3*par[0]*par[2]*(-y[1] + y[13])*(-y[2] + y[14]))/pow(pow(-y[0] + y[12],2) + pow(-y[1] + y[13],2) + pow(-y[2] + y[14],2),2.5) + (3*par[1]*par[2]*(-y[7] + y[13])*(-y[8] + y[14]))/pow(pow(-y[6] + y[12],2) + pow(-y[7] + y[13],2) + pow(-y[8] + y[14],2),2.5))/2.);
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
    return(((-3*par[0]*par[2]*(y[0] - y[12])*(y[2] - y[14]))/pow(pow(y[0] - y[12],2) + pow(y[1] - y[13],2) + pow(y[2] - y[14],2),2.5) - (3*par[0]*par[2]*(-y[0] + y[12])*(-y[2] + y[14]))/pow(pow(-y[0] + y[12],2) + pow(-y[1] + y[13],2) + pow(-y[2] + y[14],2),2.5))/2.);
    break;

  case 256 : 
    return(((-3*par[0]*par[2]*(y[1] - y[13])*(y[2] - y[14]))/pow(pow(y[0] - y[12],2) + pow(y[1] - y[13],2) + pow(y[2] - y[14],2),2.5) - (3*par[0]*par[2]*(-y[1] + y[13])*(-y[2] + y[14]))/pow(pow(-y[0] + y[12],2) + pow(-y[1] + y[13],2) + pow(-y[2] + y[14],2),2.5))/2.);
    break;

  case 257 : 
    return(((par[0]*par[2])/pow(pow(y[0] - y[12],2) + pow(y[1] - y[13],2) + pow(y[2] - y[14],2),1.5) - (3*par[0]*par[2]*pow(y[2] - y[14],2))/pow(pow(y[0] - y[12],2) + pow(y[1] - y[13],2) + pow(y[2] - y[14],2),2.5) - (3*par[0]*par[2]*pow(-y[2] + y[14],2))/pow(pow(-y[0] + y[12],2) + pow(-y[1] + y[13],2) + pow(-y[2] + y[14],2),2.5) + (par[0]*par[2])/pow(pow(-y[0] + y[12],2) + pow(-y[1] + y[13],2) + pow(-y[2] + y[14],2),1.5))/2.);
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
    return(((-3*par[1]*par[2]*(y[6] - y[12])*(y[8] - y[14]))/pow(pow(y[6] - y[12],2) + pow(y[7] - y[13],2) + pow(y[8] - y[14],2),2.5) - (3*par[1]*par[2]*(-y[6] + y[12])*(-y[8] + y[14]))/pow(pow(-y[6] + y[12],2) + pow(-y[7] + y[13],2) + pow(-y[8] + y[14],2),2.5))/2.);
    break;

  case 262 : 
    return(((-3*par[1]*par[2]*(y[7] - y[13])*(y[8] - y[14]))/pow(pow(y[6] - y[12],2) + pow(y[7] - y[13],2) + pow(y[8] - y[14],2),2.5) - (3*par[1]*par[2]*(-y[7] + y[13])*(-y[8] + y[14]))/pow(pow(-y[6] + y[12],2) + pow(-y[7] + y[13],2) + pow(-y[8] + y[14],2),2.5))/2.);
    break;

  case 263 : 
    return(((par[1]*par[2])/pow(pow(y[6] - y[12],2) + pow(y[7] - y[13],2) + pow(y[8] - y[14],2),1.5) - (3*par[1]*par[2]*pow(y[8] - y[14],2))/pow(pow(y[6] - y[12],2) + pow(y[7] - y[13],2) + pow(y[8] - y[14],2),2.5) - (3*par[1]*par[2]*pow(-y[8] + y[14],2))/pow(pow(-y[6] + y[12],2) + pow(-y[7] + y[13],2) + pow(-y[8] + y[14],2),2.5) + (par[1]*par[2])/pow(pow(-y[6] + y[12],2) + pow(-y[7] + y[13],2) + pow(-y[8] + y[14],2),1.5))/2.);
    break;

  case 264 : 
    return(0);
    break;

  case 265 : 
    return(0);
    break;

  case 266 : 
    return(0);
    break;

  case 267 : 
    return(((3*par[0]*par[2]*(y[0] - y[12])*(y[2] - y[14]))/pow(pow(y[0] - y[12],2) + pow(y[1] - y[13],2) + pow(y[2] - y[14],2),2.5) + (3*par[1]*par[2]*(y[6] - y[12])*(y[8] - y[14]))/pow(pow(y[6] - y[12],2) + pow(y[7] - y[13],2) + pow(y[8] - y[14],2),2.5) + (3*par[0]*par[2]*(-y[0] + y[12])*(-y[2] + y[14]))/pow(pow(-y[0] + y[12],2) + pow(-y[1] + y[13],2) + pow(-y[2] + y[14],2),2.5) + (3*par[1]*par[2]*(-y[6] + y[12])*(-y[8] + y[14]))/pow(pow(-y[6] + y[12],2) + pow(-y[7] + y[13],2) + pow(-y[8] + y[14],2),2.5))/2.);
    break;

  case 268 : 
    return(((3*par[0]*par[2]*(y[1] - y[13])*(y[2] - y[14]))/pow(pow(y[0] - y[12],2) + pow(y[1] - y[13],2) + pow(y[2] - y[14],2),2.5) + (3*par[1]*par[2]*(y[7] - y[13])*(y[8] - y[14]))/pow(pow(y[6] - y[12],2) + pow(y[7] - y[13],2) + pow(y[8] - y[14],2),2.5) + (3*par[0]*par[2]*(-y[1] + y[13])*(-y[2] + y[14]))/pow(pow(-y[0] + y[12],2) + pow(-y[1] + y[13],2) + pow(-y[2] + y[14],2),2.5) + (3*par[1]*par[2]*(-y[7] + y[13])*(-y[8] + y[14]))/pow(pow(-y[6] + y[12],2) + pow(-y[7] + y[13],2) + pow(-y[8] + y[14],2),2.5))/2.);
    break;

  case 269 : 
    return((-((par[0]*par[2])/pow(pow(y[0] - y[12],2) + pow(y[1] - y[13],2) + pow(y[2] - y[14],2),1.5)) - (par[1]*par[2])/pow(pow(y[6] - y[12],2) + pow(y[7] - y[13],2) + pow(y[8] - y[14],2),1.5) + (3*par[0]*par[2]*pow(y[2] - y[14],2))/pow(pow(y[0] - y[12],2) + pow(y[1] - y[13],2) + pow(y[2] - y[14],2),2.5) + (3*par[1]*par[2]*pow(y[8] - y[14],2))/pow(pow(y[6] - y[12],2) + pow(y[7] - y[13],2) + pow(y[8] - y[14],2),2.5) + (3*par[0]*par[2]*pow(-y[2] + y[14],2))/pow(pow(-y[0] + y[12],2) + pow(-y[1] + y[13],2) + pow(-y[2] + y[14],2),2.5) - (par[0]*par[2])/pow(pow(-y[0] + y[12],2) + pow(-y[1] + y[13],2) + pow(-y[2] + y[14],2),1.5) + (3*par[1]*par[2]*pow(-y[8] + y[14],2))/pow(pow(-y[6] + y[12],2) + pow(-y[7] + y[13],2) + pow(-y[8] + y[14],2),2.5) - (par[1]*par[2])/pow(pow(-y[6] + y[12],2) + pow(-y[7] + y[13],2) + pow(-y[8] + y[14],2),1.5))/2.);
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
    return(0);
    break;

  case 274 : 
    return(0);
    break;

  case 275 : 
    return(0);
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
    return(1/par[2]);
    break;

  case 283 : 
    return(0);
    break;

  case 284 : 
    return(0);
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
    return(0);
    break;

  case 292 : 
    return(0);
    break;

  case 293 : 
    return(0);
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
    return(0);
    break;

  case 301 : 
    return(1/par[2]);
    break;

  case 302 : 
    return(0);
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
    return(0);
    break;

  case 310 : 
    return(0);
    break;

  case 311 : 
    return(0);
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
    return(0);
    break;

  case 319 : 
    return(0);
    break;

  case 320 : 
    return(1/par[2]);
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
}

 return 0; 
}

