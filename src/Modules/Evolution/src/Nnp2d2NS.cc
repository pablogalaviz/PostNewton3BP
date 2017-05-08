#include "evolution.h"



double evolution::rhs_Nnp2d2_FNS(double t, const double y[], double par[2], int i)
{
  
 
  switch ( i ) {

  case 0 : 
    return(y[2]/par[0]);
    break;

  case 1 : 
    return(y[3]/par[0]);
    break;

  case 2 : 
    return((-((par[0]*par[1]*(y[0] - y[4]))/pow(pow(y[0] - y[4],2) + pow(y[1] - y[5],2),1.5)) + (par[0]*par[1]*(-y[0] + y[4]))/pow(pow(-y[0] + y[4],2) + pow(-y[1] + y[5],2),1.5))/2.);
    break;

  case 3 : 
    return((-((par[0]*par[1]*(y[1] - y[5]))/pow(pow(y[0] - y[4],2) + pow(y[1] - y[5],2),1.5)) + (par[0]*par[1]*(-y[1] + y[5]))/pow(pow(-y[0] + y[4],2) + pow(-y[1] + y[5],2),1.5))/2.);
    break;

  case 4 : 
    return(y[6]/par[1]);
    break;

  case 5 : 
    return(y[7]/par[1]);
    break;

  case 6 : 
    return(((par[0]*par[1]*(y[0] - y[4]))/pow(pow(y[0] - y[4],2) + pow(y[1] - y[5],2),1.5) - (par[0]*par[1]*(-y[0] + y[4]))/pow(pow(-y[0] + y[4],2) + pow(-y[1] + y[5],2),1.5))/2.);
    break;

  case 7 : 
    return(((par[0]*par[1]*(y[1] - y[5]))/pow(pow(y[0] - y[4],2) + pow(y[1] - y[5],2),1.5) - (par[0]*par[1]*(-y[1] + y[5]))/pow(pow(-y[0] + y[4],2) + pow(-y[1] + y[5],2),1.5))/2.);
    break;
}

 return 0; 
}


double evolution::jac_Nnp2d2_FaNS(double t, const double y[], double par[2], int i, int j)
{

 switch ( j*NDvar+i ) {

  case 0 : 
    return(0);
    break;

  case 1 : 
    return(0);
    break;

  case 2 : 
    return(((3*par[0]*par[1]*pow(y[0] - y[4],2))/pow(pow(y[0] - y[4],2) + pow(y[1] - y[5],2),2.5) - (par[0]*par[1])/pow(pow(y[0] - y[4],2) + pow(y[1] - y[5],2),1.5) + (3*par[0]*par[1]*pow(-y[0] + y[4],2))/pow(pow(-y[0] + y[4],2) + pow(-y[1] + y[5],2),2.5) - (par[0]*par[1])/pow(pow(-y[0] + y[4],2) + pow(-y[1] + y[5],2),1.5))/2.);
    break;

  case 3 : 
    return(((3*par[0]*par[1]*(y[0] - y[4])*(y[1] - y[5]))/pow(pow(y[0] - y[4],2) + pow(y[1] - y[5],2),2.5) + (3*par[0]*par[1]*(-y[0] + y[4])*(-y[1] + y[5]))/pow(pow(-y[0] + y[4],2) + pow(-y[1] + y[5],2),2.5))/2.);
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
    return(((3*par[0]*par[1]*(y[0] - y[4])*(y[1] - y[5]))/pow(pow(y[0] - y[4],2) + pow(y[1] - y[5],2),2.5) + (3*par[0]*par[1]*(-y[0] + y[4])*(-y[1] + y[5]))/pow(pow(-y[0] + y[4],2) + pow(-y[1] + y[5],2),2.5))/2.);
    break;

  case 11 : 
    return((-((par[0]*par[1])/pow(pow(y[0] - y[4],2) + pow(y[1] - y[5],2),1.5)) + (3*par[0]*par[1]*pow(y[1] - y[5],2))/pow(pow(y[0] - y[4],2) + pow(y[1] - y[5],2),2.5) + (3*par[0]*par[1]*pow(-y[1] + y[5],2))/pow(pow(-y[0] + y[4],2) + pow(-y[1] + y[5],2),2.5) - (par[0]*par[1])/pow(pow(-y[0] + y[4],2) + pow(-y[1] + y[5],2),1.5))/2.);
    break;

  case 12 : 
    return(0);
    break;

  case 13 : 
    return(0);
    break;

  case 14 : 
    return(((-3*par[0]*par[1]*(y[0] - y[4])*(y[1] - y[5]))/pow(pow(y[0] - y[4],2) + pow(y[1] - y[5],2),2.5) - (3*par[0]*par[1]*(-y[0] + y[4])*(-y[1] + y[5]))/pow(pow(-y[0] + y[4],2) + pow(-y[1] + y[5],2),2.5))/2.);
    break;

  case 15 : 
    return(((par[0]*par[1])/pow(pow(y[0] - y[4],2) + pow(y[1] - y[5],2),1.5) - (3*par[0]*par[1]*pow(y[1] - y[5],2))/pow(pow(y[0] - y[4],2) + pow(y[1] - y[5],2),2.5) - (3*par[0]*par[1]*pow(-y[1] + y[5],2))/pow(pow(-y[0] + y[4],2) + pow(-y[1] + y[5],2),2.5) + (par[0]*par[1])/pow(pow(-y[0] + y[4],2) + pow(-y[1] + y[5],2),1.5))/2.);
    break;

  case 16 : 
    return(1/par[0]);
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
    return(0);
    break;

  case 22 : 
    return(0);
    break;

  case 23 : 
    return(0);
    break;

  case 24 : 
    return(0);
    break;

  case 25 : 
    return(1/par[0]);
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
    return(((-3*par[0]*par[1]*pow(y[0] - y[4],2))/pow(pow(y[0] - y[4],2) + pow(y[1] - y[5],2),2.5) + (par[0]*par[1])/pow(pow(y[0] - y[4],2) + pow(y[1] - y[5],2),1.5) - (3*par[0]*par[1]*pow(-y[0] + y[4],2))/pow(pow(-y[0] + y[4],2) + pow(-y[1] + y[5],2),2.5) + (par[0]*par[1])/pow(pow(-y[0] + y[4],2) + pow(-y[1] + y[5],2),1.5))/2.);
    break;

  case 35 : 
    return(((-3*par[0]*par[1]*(y[0] - y[4])*(y[1] - y[5]))/pow(pow(y[0] - y[4],2) + pow(y[1] - y[5],2),2.5) - (3*par[0]*par[1]*(-y[0] + y[4])*(-y[1] + y[5]))/pow(pow(-y[0] + y[4],2) + pow(-y[1] + y[5],2),2.5))/2.);
    break;

  case 36 : 
    return(0);
    break;

  case 37 : 
    return(0);
    break;

  case 38 : 
    return(((3*par[0]*par[1]*pow(y[0] - y[4],2))/pow(pow(y[0] - y[4],2) + pow(y[1] - y[5],2),2.5) - (par[0]*par[1])/pow(pow(y[0] - y[4],2) + pow(y[1] - y[5],2),1.5) + (3*par[0]*par[1]*pow(-y[0] + y[4],2))/pow(pow(-y[0] + y[4],2) + pow(-y[1] + y[5],2),2.5) - (par[0]*par[1])/pow(pow(-y[0] + y[4],2) + pow(-y[1] + y[5],2),1.5))/2.);
    break;

  case 39 : 
    return(((3*par[0]*par[1]*(y[0] - y[4])*(y[1] - y[5]))/pow(pow(y[0] - y[4],2) + pow(y[1] - y[5],2),2.5) + (3*par[0]*par[1]*(-y[0] + y[4])*(-y[1] + y[5]))/pow(pow(-y[0] + y[4],2) + pow(-y[1] + y[5],2),2.5))/2.);
    break;

  case 40 : 
    return(0);
    break;

  case 41 : 
    return(0);
    break;

  case 42 : 
    return(((-3*par[0]*par[1]*(y[0] - y[4])*(y[1] - y[5]))/pow(pow(y[0] - y[4],2) + pow(y[1] - y[5],2),2.5) - (3*par[0]*par[1]*(-y[0] + y[4])*(-y[1] + y[5]))/pow(pow(-y[0] + y[4],2) + pow(-y[1] + y[5],2),2.5))/2.);
    break;

  case 43 : 
    return(((par[0]*par[1])/pow(pow(y[0] - y[4],2) + pow(y[1] - y[5],2),1.5) - (3*par[0]*par[1]*pow(y[1] - y[5],2))/pow(pow(y[0] - y[4],2) + pow(y[1] - y[5],2),2.5) - (3*par[0]*par[1]*pow(-y[1] + y[5],2))/pow(pow(-y[0] + y[4],2) + pow(-y[1] + y[5],2),2.5) + (par[0]*par[1])/pow(pow(-y[0] + y[4],2) + pow(-y[1] + y[5],2),1.5))/2.);
    break;

  case 44 : 
    return(0);
    break;

  case 45 : 
    return(0);
    break;

  case 46 : 
    return(((3*par[0]*par[1]*(y[0] - y[4])*(y[1] - y[5]))/pow(pow(y[0] - y[4],2) + pow(y[1] - y[5],2),2.5) + (3*par[0]*par[1]*(-y[0] + y[4])*(-y[1] + y[5]))/pow(pow(-y[0] + y[4],2) + pow(-y[1] + y[5],2),2.5))/2.);
    break;

  case 47 : 
    return((-((par[0]*par[1])/pow(pow(y[0] - y[4],2) + pow(y[1] - y[5],2),1.5)) + (3*par[0]*par[1]*pow(y[1] - y[5],2))/pow(pow(y[0] - y[4],2) + pow(y[1] - y[5],2),2.5) + (3*par[0]*par[1]*pow(-y[1] + y[5],2))/pow(pow(-y[0] + y[4],2) + pow(-y[1] + y[5],2),2.5) - (par[0]*par[1])/pow(pow(-y[0] + y[4],2) + pow(-y[1] + y[5],2),1.5))/2.);
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
    return(0);
    break;

  case 52 : 
    return(1/par[1]);
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
    return(1/par[1]);
    break;

  case 62 : 
    return(0);
    break;

  case 63 : 
    return(0);
    break;
}

 return 0; 
}

