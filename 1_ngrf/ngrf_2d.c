//Generate diffuse image at different freq. in FITS fromat
// uses input angular power spectrum of brightness temperature fluctuation 

#include<stdlib.h>
#include<stdio.h>
#include<fftw3.h>
#include<math.h>
#include<fitsio.h>
#include<unistd.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

double k_B=1.38e3 /* in Jy*/, c=3.0e8;


//return input angular power spectrum
double  P_T(double UU)// return input brightness temperature PS in K^2
{
  double ell,pu;  

  double Av=1, betav=-1.0; // amplitude and multipole exponent of APS

  ell=2.0*M_PI*UU;// calculate ell --. multipole


  //pu=Av*pow(ell/1000.0,betav)*exp(-1.*pow(UU/500.,2));
  //pu=Av*pow(ell/1000.0,betav);
  
  /////////////////////////////////////////////////////////////

  //pu=1.0*exp(-1.*pow(UU/500.,2)); // unit angular power spectrum in K^2
  pu=1.0 ;

  /////////////////////////////////////////////////////////////

  return(pu);
}


int main(int argc, char *argv[])
{



  double *out, *img;
  int NN,length, ii, jj, indx, ia, ja;
  double  LL,pixel,ff,nu0,nu; 
  unsigned long int seed;
  char  OPT, task[128];
  int Nchan, Rlzn;
  double chan0, deltanu, spindex;
  double thetax, thetay, be;

  FILE *fp;

  void Fill_Modes(int NN,double LL,unsigned long int seed,double (*func)(double),double *out);
  void Make_Image(int NN,double LL,double *out);
  void Write_Fits_Image(int NN,double LL,double *out,char *outfile);
  void Non_Gauss(int NN,double LL,double *out,double ff);
  void Get_Modes(int NN,double LL,double *out);
  void power_spec(int NN,double LL,double *out,int NB,double *usum,double *pksum,int *no);
  //void PB_Flux(int NN,double LL, double *out, double nu); 
  void  PB_Flux_MF(int NN,double LL, double *out, double *img, int Nchan, double nu0,double chan0, double spindex, double deltanu);
  void PB_Flux_MF_BAIC(int NN,double LL, double *out, double *img, int Nchan, double nu0,double chan0, double spindex,double deltanu);
  double Beam(double theta,double freq);
  double Beam_sinc(double theta_x,double theta_y, double freq);
  void Write_Fits_MF_Image(int NN,double LL, int Nchan, double nu0,double chan0, double deltanu, double *img,char *outfile);



  
  // reads name of parameter file and output FITS file at runtime
  if(argc!=13)
   {
      printf("Usage: %s <N no. of pixs.><L gridsize (arcmin)><seed int <0><non-Gassian parameter ff><spectral index (Intensity)><nu0><channel width><chan0><Number of channels><ouput FITS file SF><ouput FITS file MF>\n", argv[0]);
      return 1;
    }
    
  //system(task);
  sscanf(argv[1],"%d",&NN); // 2048
  sscanf(argv[2],"%lf",&LL); //0.515= 1.453e-3 radian
  sscanf(argv[3],"%ld",&seed); //7999
  sscanf(argv[4],"%lf",&ff); //0.5
  sscanf(argv[5],"%lf",&spindex); //0.8
  sscanf(argv[6],"%lf",&nu0); //150.00e06 Hz
  sscanf(argv[7],"%lf",&deltanu); //62.50e3 in Hz// channel width
  sscanf(argv[8],"%lf",&chan0); //1
  sscanf(argv[9],"%d",&Nchan); //can vary  8/20/64
  sscanf(argv[12],"%d",&Rlzn);
  printf("%d %lf %ld %lf %lf %lf %lf %lf %d \n",NN,LL,seed,ff,spindex,nu0,deltanu,chan0,Nchan);
   
 
  pixel=LL;// pixel resolution in arcmin, we have used this in HEADER
 
  LL=M_PI*LL/(180.*60.);// in rad, added
  
  //memory allocation for out and in  NN*(NN+2) array 

  out=(double*)calloc ((NN*(NN+2)),sizeof(double));
  img=(double*)calloc (NN*NN*Nchan,sizeof(double));
  // call the functions
  
  Fill_Modes(NN, LL, seed,(*P_T),out);
  
  Make_Image(NN,LL,out);


  Non_Gauss(NN,LL,out,ff);
   
  //Write_Fits_Image(NN,LL,out,argv[10]); //In temerature unit at central freq
 
  PB_Flux_MF(NN, LL, out, img, Nchan, nu0, chan0, spindex, deltanu);// multi-frequency//beam at N/2,N/2
  //PB_Flux_MF_BAIC(NN, LL, out, img, Nchan, nu0, chan0, spindex, deltanu); //beam at bottom left corner
 
  Write_Fits_MF_Image(NN, LL, Nchan, nu0,  chan0, deltanu, img, argv[11]); //In Intensity

 




  //----------- Multiply with beam -------------//

    /*  for(ii=0;ii<NN;++ii)
	for(jj=0;jj<NN;++jj)
	  {
	    indx =ii*(NN+2)+jj;
	    ia = (ii>NN/2) ? (NN-ii) : ii ;
	    ja = (jj>NN/2) ? (NN-jj) : jj ;
	    thetax= ia*LL;
	    thetay= ja*LL;

	    be = Beam_sinc(thetax,thetay,nu0);
	    out[indx]*=(be); 
	  }*/
	  
	/* for(ii=0;ii<1;++ii)
	  for(jj=0;jj<10;++jj){
		indx=ii*(NN+2)+jj;
		printf("%f\t\t",out[indx]); }*/
	  
	  // --- below is beam assuming N/2,N/2 is center ---//
     /*  for(ii=0;ii<NN;++ii)
	for(jj=0;jj<NN;++jj)
	  {
	    indx=ii*(NN+2)+jj;
	    thetax= (ii-NN/2)*LL;
	    thetay= (jj-NN/2)*LL;

	    be=Beam_sinc(thetax,thetay,nu0);
	    out[indx]*=be; // in function : Beam = beam*beam
	  }*/


  //----------- Following is for validation -------------//
  
  Get_Modes(NN,LL,out); 
  
       /*	for(ii=0;ii<1;++ii)
	  for(jj=0;jj<10;++jj){
		indx=ii*(NN+2)+jj;
		printf("%f\t\t",out[indx]); }*/
  
  // --- Writing in binary file to find bispectrum later //
  
    fftw_complex *comp_ro;
    comp_ro = (fftw_complex *)&(out[0]);
    FILE *pp;
    char fname[128];
    snprintf(fname,128*sizeof(char),"output_data/2D_ngrf_field_R%d",Rlzn);
    pp=fopen(fname,"w");
    
  /*  for(ii=0;ii<N1;++ii)
      for(jj=0;jj<N2+2;++jj){
      	index=ii*(NN+2)+jj;
	fwrite(&out[index],sizeof(double),1,pp);
       }*/
        
    for(ii=0;ii<NN;++ii)
      for(jj=0;jj<=NN/2;++jj){       //use if beam multiplied left bottom corner..(use it mostly)
	indx=ii*(NN/2+1)+jj;
	fwrite(&comp_ro[indx][0],sizeof(double),1,pp);
	fwrite(&comp_ro[indx][1],sizeof(double),1,pp);}
	
   /*  double re_part, im_part;
     for(ii=0;ii<NN;++ii)
      for(jj=0;jj<=NN/2;++jj){
	indx=ii*(NN/2+1)+jj;
	re_part = pow(-1.,ii+jj)*comp_ro[indx][0];
	im_part = pow(-1.,ii+jj)*(-1.)*comp_ro[indx][1];
	fwrite(&re_part,sizeof(double),1,pp); // if beam multiplied at center then it shifts image center at 0,0
	fwrite(&im_part,sizeof(double),1,pp);}*/

    fclose(pp);
    
    // -----------------------------------------------//
  
  
  //length=LL*NN; // angular dimension of one edge of image in radian
  //############### for binned power spectrum calculation #################### //
  int i,NB=15,*num; 
  double *binu,*binps;
 
  num=(int*)calloc(NB,sizeof(int)); // stores number of U modes in each bin
  binu=(double*)calloc(NB,sizeof(double));// stores sum of |U| in each bin
  binps=(double*)calloc(NB,sizeof(double));// PS in each bin
    
  power_spec(NN,LL,out,NB,binu,binps,num);// Calls the power spectrum function
  //power_spec(NN,LL,img,NB,binu,binps,num);// Calls the power spectrum function

  //Write power spectrum in data file
  //fp=fopen("psout.dat","w");
  snprintf(fname,128*sizeof(char),"output_data/psout_R%d.dat",Rlzn);
  fp=fopen(fname,"w");

  for(i=0;i<NB;++i)
    if(num[i]>0)
      {
	fprintf(fp,"%e %e %e %d \n",binu[i],binps[i],P_T(binu[i]),num[i]); // These are in temperature unit
      }
  fclose(fp);
  //############################################################
  
  free(binu);
  free(num);
  free(binps);
  
  printf("Realization = %d done.\n", Rlzn); 
}				  

