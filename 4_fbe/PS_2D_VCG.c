/*=========================================================================
       This code uses a convlved gridded visibility data in 2D for a 
       simulated non-Gaussian random field and computes the corresponding 
       power spectra and bispectra. Bispectrum computation uses FFT method.
       Authors : Abinash Kumar Shaw (IIT Kharagpur)
                 Prof. Somnath Bharadwaj (IIT Kharagpur)
                 Sk. Siyad Ali 
===========================================================================*/
#include<stdlib.h>
#include<stdio.h>
#include<fftw3.h>
#include<math.h>
#include<unistd.h>
#include<omp.h>
#include<fitsio.h>
#include<gsl/gsl_errno.h>
#include<gsl/gsl_spline.h>

/*---- Global variables ----*/

long N1,N2; // Dimension
int Nthreads, F_Task;
double Du, pi=M_PI; //Du is grid spacing in uv plane
double **ro, **ro_u; // for power spectrum
double **b_ro, **b_ro_u; // for bispectrum
fftw_plan p_ro, q_ro;


int nc=1,Nm,Nu;        
double fx,FWHM,dUx;

double *nua;                   
double **GVRRre,**GVRRim,**GVLLre,**GVLLim,**SCRR,**SCLL;
char GV[256];//char GV[128]="../3_vcg/GV_data";
/*---------------------------*/

int main(int argc, char *argv[])
{
  int Nbin, Rlzn;
  double sigma1, sigma2, ff1, LL, Area,ff2,Av=1.,beta=-2., Umax1, Umin1; 
  unsigned long int TotN=0;
  char  OPT;
  char file[128],file1[128],file2[128];
  FILE *fp,*fp1,*fp2;
  double t,T=omp_get_wtime();
  
  int i, j, l, NB, *num;
  double *binu,*binps,*logU,*logPU;
  
  int ii, jj;

  // void read_uvimage(char filename1[], char filename2[]);
 
  void power_spec_GV(int NB, double umax, double umin, double *usum, double *pksum, int *no);
  void power_spec_GVRL(int NB, double umax, double umin, double *usum, double *pksum, int *no);
  double **allocate_2d_double(int N1,int N2);
  
  
  
   sscanf(argv[1],"%d",&Rlzn);
  sscanf(argv[2],"%s",GV); 
  
  
  
  
  
  int initialize();
  void ReadGVnAC();
  
  /*----------------------------------------------------------*/
  initialize(); 
  printf("\nArrays initialized\n");
  //GV = "../3_vcg/UGV_data";
  ReadGVnAC();
  printf("\nRead gridded visibilities\n");
  printf("\nTotal no. of channels (nc) = %d",nc);
  printf("\nFWHM = %.2f arcmin.",FWHM);
  printf("\nGrid spacing, dU = %.2f wavelengths\n",dUx);
  printf("\nFreq. nu_0 = %e",nua[0]);
  printf("\nNm = %d",Nm);
  printf("\nNu = %d",Nu);
  printf("\nN1 = %ld",N1);
  
  Du = dUx;
  
  
  
  
  
 /* 
  int nn=0, index;
  int MM = 2*Nm;
  for(i=MM; i<Nu-MM;++i)
    for(j=MM;j<Nu-MM;++j)
      {
  	index=i*((Nu-2*MM)/2+1)+j;
  	if (GVRRre[index][nn]!=0)
  	 {printf("\n%f",GVRRre[index][nn]); }
  }
  */
  
  
  
  
  /*================ binned power spectrum computation ===============*/ 
  NB = 20;
  Umax1 = 500.;
  Umin1 =1;
  
  num=(int*)calloc(NB,sizeof(int));
  binu=(double*)calloc(NB,sizeof(double));
  binps=(double*)calloc(NB,sizeof(double));
  logU=(double*)calloc(NB,sizeof(double));
  logPU=(double*)calloc(NB,sizeof(double));
  
 // t=omp_get_wtime();
  
  power_spec_GVRL(NB, Umax1, Umin1, binu, binps, num);
  //power_spec_GV(NB, Umax1, Umin1, binu, binps, num);

  //Write power spectrum in data file
 // fp=fopen("apsout.dat","w");
    snprintf(file2, 512*sizeof(char),"PS_results/psout_%s",GV); // open file to write bispectrum
   fp=fopen(file2,"w");

  for(i=0;i<NB;++i)
   // if(num[i]>0)
      {
  	fprintf(fp,"%e\t%e\t%d\n",binu[i],binps[i],num[i]);
  	logU[i]=log10(binu[i]); // log is done for the interpolation used in B_Ana
  	logPU[i]=log10(binps[i]);
      }
  fclose(fp);
  fprintf(stdout,"\n Power spectrum time = %e sec\n",omp_get_wtime()-t);
  /*===================================================================*/
  

}
