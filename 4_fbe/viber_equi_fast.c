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
double **GVRRre,**GVRRim,**GVLLre,**GVLLim,**SCRR,**SCLL,**NTiG;
char GV[256];//="../3_vcg/output_data/GV_data_1";
/*---------------------------*/

int main(int argc, char *argv[])
{
  int Nbin, Rlzn, index;
  double sigma1, sigma2, ff1, LL, Area,ff2,Av=1.,beta=-1.; 
  unsigned long int TotN=0;
  char  OPT;
  char file[512],file1[512],file2[512];
  FILE *fp,*fp1,*fp2;
  double t,T=omp_get_wtime();
  
  int i, j, l, NB, *num;
  double *binu,*binps,*logU,*logPU;
  
  int ii, jj;
  double Bispec, Biana,Denom,Bierr, Num, tt, mu, PU1, PU2, PU3;
  double U1, U2, U3, dU, *ULow, *UHigh, *Uavg;
  double Umax, Umin, Umax1, Umin1;
  double **I1, **I2, **I3, **D1, **D2, **D3, **O1, **O2, **O3; ;

  
  // void read_uvimage(char filename1[], char filename2[]);
  void read_uvdata(char filename1[], char filename2[]);
  void power_spec_GV(int NB, double umax, double umin, double *usum, double *pksum, int *no);
  void fft_init(int Nbin, double *uL,double *uH);
  void fft_init_orig(int Nbin, double *uL,double *uH);
  void shapes(int nbin, double *uavg);
  double **allocate_2d_double(int N1,int N2);
  double P_T(double UU);
  
  int initialize_gridding();
  void ReadGVnAC_gridding();
  
  /*------ reading parameters from runtime arguments --------*/
 /*if(argc!=12)
   {
      printf("Usage: %s </full/path/for/data/FITSfile></full/path/for/UAPS/FITSfile><U1_max><U1_min><No. of U1 bins><# of threads><Bispec Flag (1/0)><fNG1(data)> ><fNG2(uaps)><sigma1(data)><sigma2(uaps)>\n", argv[0]);
      return 1;
    }
  
  sscanf(argv[3],"%lf",&Umax);
  sscanf(argv[4],"%lf",&Umin);
  sscanf(argv[5],"%d",&Nbin);
  sscanf(argv[6],"%d",&Nthreads);
  sscanf(argv[7],"%d",&F_Task);
  sscanf(argv[8],"%lf",&ff1);
  sscanf(argv[9],"%lf",&ff2);
  sscanf(argv[10],"%lf",&sigma1);
  sscanf(argv[11],"%lf",&sigma2);*/
  
  sscanf(argv[1],"%d",&Rlzn);
  sscanf(argv[2],"%s",GV); 

  
    
  
  //GV="../3_vcg/GV_data";
 // strcpy(GV, "../3_vcg/GV_data");
 
  Du = 1.0254384902642049;//((180.*60.)/(3.14159265359*2048.*1.719));//dUx; 
  N1 = 1176;
  N2 = 1176;
  
  Umax  	=  N1/3;//270;
  Umin 		=  1;
  Nbin   	=  25;//31//22;
  Nthreads	=  4;
  F_Task	=  1;
  
  /*--------- reading input file. ----------*/

  // LL= 2.50000e-01; //In arcmin , Origial Image resulation here 15"
  // LL*=(pi/(60.*180.)); // changing from aremin to radian # check this first
  // Area=(2048*LL)*(2048*LL); //Area Str for 2048*2048 
  
  
    /*----------------------------------------------------------*/
 
  
  initialize_gridding(); 
  printf("\nArrays initialized\n");
  ReadGVnAC_gridding();
  
  
  /*---------Data normalization-------------------*/
  
  
  /*  for(ii=0; ii<N1;ii++)
     for(jj=0; jj<N2/2+1;jj++)
      {
	index = ii*(N2/2+1)+jj; 
	for(i=0; i<nc; i++)
	  {
	    if (NTiG[index][i]>0)
	    {
	    GVRRre[index][i] /= NTiG[index][i];
	    GVRRim[index][i] /= NTiG[index][i];
	    }
	  }
      }
  
    for(ii=0; ii<N1;ii++)
     for(jj=0; jj<N2/2+1;jj++)
      {
	index = ii*(N2/2+1)+jj; 
	for(i=0; i<nc; i++)
	  {
	  if (NTiG[index][i]>0)
	   {
	    NTiG[index][i]=1;
	   }
	  }
      }
  */
  /*----------------------------------------------------------*/
  
  
  printf("\nRead gridded visibilities\n");
  printf("\nTotal no. of channels (nc) = %d",nc);
  printf("\nFWHM = %.2f arcmin.",FWHM);
  printf("\nGrid spacing, dU = %.2f wavelengths\n",dUx);
  printf("\nFreq. nu_0 = %e",nua[0]);
  printf("\nNm = %d",Nm);
  printf("\nNu = %d",Nu);
  printf("\nN1 = %ld",N1);
  
   
  
  fprintf(stdout,"N1=%ld  N2=%ld  dU=%e (in rad^{-1}) Area = %lf \n",N1,N2,Du,Area);
 
  fprintf(stdout,"Umax= %e\nUmin= %e\nU1 bin #= %d\nOPENMP Threads= %d\nF_Task= %d\n", Umax, Umin, Nbin, Nthreads, F_Task);

  /*------- runtime arguments reading done. ---------*/
  
  /*------------ Checking the limit of U1 ------------*/
  
   /*  if(Umax>Du*N1/3)
  {
   printf("Umax is out of the bound dU*(N1/3).\n");
   return 1; 
   } */

  /*----------- Dividing U1 plane for bispectra -----------*/
  ULow=(double*)calloc(Nbin,sizeof(double));
  UHigh=(double*)calloc(Nbin,sizeof(double));
  Uavg=(double*)calloc(Nbin,sizeof(double));
  
 /*------------- Dividing k1 space into linear bins -------------*/
 /* dU=(Umax-Umin)/(1.*Nbin); // goes max upto dU*(N/3) to avoid modes outside boxes
   fprintf(stdout,"Binsize= %e\n", dU);

  for(i=0;i<Nbin;++i)
  {
   ULow[i]=i*dU+Umin;
   UHigh[i]=(i+1)*dU+Umin;
   Uavg[i]=(2./3.)*(pow(UHigh[i],3.)-pow(ULow[i],3.))/(pow(UHigh[i],2.)-pow(ULow[i],2.)); // bin average value
   
   printf("Sl=%d \t ULow=%e\t UHigh=%e \t Uave=%e\n",i+1,ULow[i],UHigh[i], Uavg[i]);

  }
  */
  
  
  
   // Define intervals
    //int intervals[][4] = {{0, 50, 5}, {50, 100, 10}, {100, 200, 25}, {200, 350, 50}};
    
    //int intervals[][4] = {{0, 50, 5}, {50, 100, 10}, {100, 200, 20}, {200, 270, 35}};
    //int num_intervals = 4;
    
    int intervals[][5] = {{0, 50, 5}, {50, 100, 10}, {100, 200, 20}, {200, 305, 35}, {305, 393, 44}};
    int num_intervals = 5;
    
    //int intervals[][5] = {{0, 50, 5}, {50, 100, 10}, {100, 200, 20}, {200, 410, 35}, {410, 660, 50}};
    //int num_intervals = 5;

    // Index for ULow array
    index = 0;

    // Loop through intervals and fill ULow and UHigh
    for (i = 0; i < num_intervals; ++i) {
        int start = intervals[i][0];
        int end = intervals[i][1];
        int step = intervals[i][2];

        for (j = start; j < end; j += step) {
            if (index < Nbin) {
                ULow[index] = j;
                if (index > 0) {
                    UHigh[index - 1] = ULow[index];
                }
                ++index;
            }
        }
    }
 
  ULow[0] = 1.;
  UHigh[Nbin-1]=1.*Umax;
  
  
  for(i=0;i<Nbin;++i)
  { //k_perp_avg[i]=(k_perp_High[i]+k_perp_Low[i])/2.;

/*if (i<10)
   {k_perp_avg[i]=(k_perp_High[i]+k_perp_Low[i])/2.;}
else
   {k_perp_avg[i]=(2./3.)*(pow(k_perp_High[i],3.)-pow(k_perp_Low[i],3.))/(pow(k_perp_High[i],2.)-pow(k_perp_Low[i],2.)); }// bin average value */
   
   Uavg[i]=(2./3.)*(pow(UHigh[i],3.)-pow(ULow[i],3.))/(pow(UHigh[i],2.)-pow(ULow[i],2.));
   
   printf("\t%f\t%f\t%f\n",ULow[i],UHigh[i],Uavg[i]) ;
  }
  
  
  
  
  /*-------- varying dkperp in binning ----------*/
  /*int no_ring_set = 14;
  
  int currentNumber = 1;
  for ( i = 0; i < Nbin; i += no_ring_set) {
    for ( j = i; j < i + no_ring_set; j++) {
        ULow[j] = currentNumber;
	currentNumber += (i / no_ring_set) + 1;}}
  for(i=1;i<Nbin;++i)
  {
   // ULow[i] -= .1;  // otherwise mu=1, stretched bs=nan; could do +=.1 but then bs for k=0 mode is difficult to handle.so keep 0th ring from 0-.9
    UHigh[i-1]=ULow[i];
  }
  UHigh[Nbin-1]=1.*Umax;
  
  
  for(i=0;i<Nbin;++i)
  { //k_perp_avg[i]=(k_perp_High[i]+k_perp_Low[i])/2.;

//if (i<10)
 //  {k_perp_avg[i]=(k_perp_High[i]+k_perp_Low[i])/2.;}
//else
//   {k_perp_avg[i]=(2./3.)*(pow(k_perp_High[i],3.)-pow(k_perp_Low[i],3.))/(pow(k_perp_High[i],2.)-pow(k_perp_Low[i],2.)); }// bin average value 
   
   Uavg[i]=(2./3.)*(pow(UHigh[i],3.)-pow(ULow[i],3.))/(pow(UHigh[i],2.)-pow(ULow[i],2.));
   
   printf("\t%f\t%f\t%f\n",ULow[i],UHigh[i],Uavg[i]) ;
  }
  
*/
  
 /*---------------------- Binning done --------------------------*/
 
 if(F_Task==0)
 {
  fprintf(stdout,"\nWriting the various shapes of closed triangles of all possible sizes.\n");
  shapes(Nbin, Uavg);
 }
 else if(F_Task==1)
 {
  fprintf(stdout,"\nComputing the power spectrum and bispectrum of closed triangles of all possible shapes and sizes.\n");  
  
  /*------------ Initializing fft plan and openmp ------------*/
  fftw_init_threads();
  fftw_plan_with_nthreads(Nthreads);
  omp_set_num_threads(Nthreads);
  

  /*================ binned power spectrum computation ===============*/ 
  NB = 15;
  Umax1 = 500.;
  Umin1 =1.;
  
  num=(int*)calloc(NB,sizeof(int));
  binu=(double*)calloc(NB,sizeof(double));
  binps=(double*)calloc(NB,sizeof(double));
  logU=(double*)calloc(NB,sizeof(double));
  logPU=(double*)calloc(NB,sizeof(double));
  
  t=omp_get_wtime();
  
  power_spec_GV(NB, Umax1, Umin1, binu, binps, num);

  //Write power spectrum in data file

  
   snprintf(file2, 512*sizeof(char),"results/psout_%s",GV); // open file to write bispectrum
   fp2=fopen(file2,"w");

  for(i=0;i<NB;++i)
    if(num[i]>0)
      {
  	fprintf(fp2,"%e\t%e\t%d\n",binu[i],binps[i],num[i]);
  	logU[i]=log10(binu[i]); // log is done for the interpolation used in B_Ana
  	logPU[i]=log10(binps[i]);
      }
  fclose(fp2);
  fprintf(stdout,"Power spectrum time = %e sec\n",omp_get_wtime()-t);
  /*===================================================================*/
  
  /*=====================  Bispectrum computation ====================*/
  
  /*------- for power spectrum interpolation -------*/
  // gsl_interp_accel *acc= gsl_interp_accel_alloc();
  // gsl_spline *spline= gsl_spline_alloc (gsl_interp_cspline, NB);
  // gsl_spline_init (spline, logU, logPU, NB);
  
  t=omp_get_wtime();
  
  /*-------- Initializing arrays to store D(Un,x) and I(Un,x) ---------*/
  I1=allocate_2d_double(N1,N2+2); 
  D1=allocate_2d_double(N1,N2+2);

 fftw_complex *i1,*d1;
  i1=(fftw_complex*)&(I1[0][0]);
  d1=(fftw_complex*)&(D1[0][0]);
  
  fftw_plan p1, pd1;
  p1= fftw_plan_dft_c2r_2d(N1, N2, (fftw_complex*)&(I1[0][0]), &(I1[0][0]), FFTW_ESTIMATE);
  pd1= fftw_plan_dft_c2r_2d(N1, N2, (fftw_complex*)&(D1[0][0]), &(D1[0][0]), FFTW_ESTIMATE);
 double uu;
 int m, ia, ja;

  /*-------------------------------------------------------*/
  
 ////////// fft_init_orig(Nbin, ULow, UHigh); // perform FFTs to compute D and I for all rings
  
   snprintf(file, 512*sizeof(char),"results/bsout_%s",GV); // open file to write bispectrum
   fp1=fopen(file,"w");

  int index_1, nn=0;
  int ic, jc, iagv;
 for(m=0;m<Nbin;++m)
  {
  tt = 1.;
  mu = 0.5;
   U1=Du*Uavg[m];
   PU1=P_T(U1);
   //count=0;
   for(i=0;i<N1;++i)
    for(j=0;j<N2/3;++j)
      {
      
        ia=(i>N1/2) ? (N1-i) : i ;
        index=i*(N2/2+1) + j;
  	uu=sqrt(1.*(ia*ia+j*j));
  	
  	/*ic = i+N2/2;
        jc = j+N2/2;
        iagv=(ic>=N1) ? (ic-N1) : ic ;
        index_1 = iagv*N1+jc;*/
      
  	if(uu>=ULow[m] && uu<UHigh[m] && uu>0.)
  	{
  	 //i1[index][0]=1.; i1[index][1]=0.;
  	 i1[index][0]=NTiG[index][nn]; i1[index][1]=0.;
  	 //d1[index][0]=GVRRre[index_1][nn]; d1[index][1]=GVRRim[index_1][nn];// use if center at N/2,N/2
  	 d1[index][0]=GVRRre[index][nn]; d1[index][1]=GVRRim[index][nn];// use if center at 0,0
  	 //count++; // count modes within the bin
  	}
      }
      
      
    fftw_execute(p1); // perform FT
    fftw_execute(pd1); // perform FT
        
        
        Bispec=0.; 
        //Denom=0.; 
	Num=0.;
  	 for(ii=0;ii<N1;++ii)
    	   for(jj=0;jj<N2;++jj)
       	{
         	 Bispec+=D1[ii][jj]*D1[ii][jj]*D1[ii][jj]; // numerator (eq. 7 Shaw et al.)
         	 Num+=I1[ii][jj]*I1[ii][jj]*I1[ii][jj]; // denominator (eq. 10 Shaw et al.) 
        	}
        	
	  Bispec/= Num;//(Num/(6.*ff2/sigma2)); // FBE Bispectrum (eq. 24 Ali et al.)
	     
	  Num/=(N1*N2*1.); // No. of triangles
	  Biana=6.*(PU1*PU1);
	  //Biana= 2.*(PU1*PU2+PU2*PU3+PU3*PU1); // eq. 27 Ali et al.
	  //Biana/=(6.*ff1/sigma1);//testing for UABS=1

	 
	
	 //	Bierr=sqrt(Area*PU1*PU2*PU3+3.*Biana*Biana)/sqrt(Num); // cosmic varianc
	//	if(Bispec>0.0)

	  //fprintf(fp1,"%e\t%e\t%e\t%e\t%e\t%e\n",tt,mu,Bispec, Biana, ((fabs(Bispec)-Biana)*100)/Biana, Num); // writing output
	  fprintf(fp1,"%d\t%e\t%e\t%e\t%e\t%e\t%.0lf\n",m,tt,mu,Bispec,Biana, ((fabs(Bispec)-Biana)*100)/Biana,Num); 
     	TotN+=(long)round(Num);
     	
     for(i=0;i<N1;++i)
      for(j=0;j<N2+2;++j)
          {
            I1[i][j]=0.;
            D1[i][j]=0.;
          }
 
   //fclose(fp1);
  } fclose(fp1);
  
  

  fprintf(stdout,"Total triangles = %ld\n",TotN);
  fprintf(stdout,"Bispectra computation time = %e sec\n",omp_get_wtime()-t);
 
 //free(num);
// free(binps);
 //free(binu);

 
 
fprintf(stdout,"Total time taken = %d hr %d min %d sec\n",(int)((omp_get_wtime()-T)/3600), (int)((omp_get_wtime()-T)/60)%60, (int)(omp_get_wtime()-T)%60);
 }
 else
 {
   fprintf(stderr,"\nWrong \'F_Task\' value.\nF_Task= 0 --> Writes shapes and sizes of the triangles for given Nmax and Nbin.\nF_Task= 1 --> Computes the powerspectrum and bispectrum values.\n");
 }
}
