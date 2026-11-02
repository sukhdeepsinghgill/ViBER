# include <stdio.h>
# include <math.h>
# include <string.h>
# include <stdlib.h>
# include <gsl/gsl_sf_bessel.h>
# include <gsl/gsl_sf_trig.h>


double Beam_sinc(double theta_x,double theta_y, double freq)  // gsl_sinc is same as in python
{
  double D=4.0, nubyc, c; // D at FWHM=61 arc min.
  double be, be_x, be_y, arg;
  c=3.*pow(10,8.); //m/s
  
  nubyc=freq/c; //m^-1
  
  D=D*nubyc;
  
  be_x = gsl_sf_sinc(D*theta_x);
  be_y = gsl_sf_sinc(D*theta_y);
  
  be_x = be_x*be_x ;
  be_y = be_y*be_y ;

  be = be_x*be_y ;

  //be=exp(-pos*pos/(9.8e-3*9.8e-3)); //main lobe matched with Gaussian for this theta0 value 
  return(be);
}




double Beam_bessel(double pos,double freq)
{
  double D=40.0, nubyc, c; // D at FWHM=61 arc min.
  double be,arg;
  c=3.*pow(10,8.); //m/s
  
  nubyc=freq/c; //m^-1
  
  D=D*nubyc;
  arg=M_PI*pos*D;
  if(pos==0.)
    {
      be=1.;
    }
  else
    {
      be=(2.)*(gsl_sf_bessel_J1(M_PI*pos*D)/(M_PI*pos*D));

      be=be*be;
    }
 
  //be=exp(-pos*pos/(9.8e-3*9.8e-3)); //main lobe matched with Gaussian for this theta0 value 
  return(be);
}
