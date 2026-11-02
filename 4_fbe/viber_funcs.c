/*=========================================================================
       This supports in computing the power spectra and bispectra from
       a convolved gridded visiblity data of a 2D non-Gaussian field generated
       from simulations. Bispectrum computation uses FFT method.
       Authors : Abinash Kumar Shaw (IIT Kharagpur)
                 Prof. Somnath Bharadwaj (IIT Kharagpur) 
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

//------ Global variables ---------//

extern long N1,N2;
extern int Nthreads, F_Task;
extern double Du, pi;
extern double **ro, **ro_u; // for power spectrum
extern double **b_ro, **b_ro_u; // for bispectrum
extern fftw_plan p_ro, q_ro;


extern int nc,Nm,Nu;        
extern double fx,FWHM,dUx;
extern double *nua;                   
extern double **GVRRre,**GVRRim,**GVLLre,**GVLLim,**SCRR,**SCLL,**NTiG;

extern char GV[128];

/*-------------- Array Allocation functions -----------------*/
double **allocate_2d_double(int N1,int N2)
{
  int ii;
  long asize;
  double **phia, *phi;

  phia=(double **)calloc(N1,sizeof(double *));
  
  asize = N1*N2;

  if(!(phi = (double *) calloc(asize,sizeof(double))))
    {
      fprintf(stderr,"Error in allocate_2d_double.\n");
      exit(0);
    }

  for(ii=0;ii<N1;++ii)
	phia[ii]=phi+ N2*ii;

  return(phia);
}
/*-------------------------------------------------------*/


double P_T(double UU)// return input brightness temperature PS in K^2
{
  double ell,pu;  
  double Av=1.0, betav=-1.0; // amplitude and exponent of PS

  ell=2.0*M_PI*UU;// calculate ell  
  //pu=Av*pow(ell/1000.0,betav);  // use this if you are reading 2D_ngrf_field file
  pu=Av*pow(ell/1000.0,betav)*exp(-1.*pow(UU/500.,2));
  //pu=1.;   // use this if you are reading 2D_Ungrf_field file
  return(pu);
}



void fft_init(int Nbin, double *uL,double *uH)
{
  double uu;
  double **I, **D;
  int i,j,k,m,index,ia,ja,index_1;
  unsigned long int count;
  FILE *pp,*pf;
  char fname[128];
  int nn=0;
  //int index1;
  
  I=allocate_2d_double(N1,N2+2); // for denominator
  D=allocate_2d_double(N1,N2+2); // for numerator
  
  fftw_complex *i1,*d1, *comp_ro;
  
  i1=(fftw_complex*)&(I[0][0]);
  d1=(fftw_complex*)&(D[0][0]);
 
  /*---------------------------------------------------------*/
  
  fftw_plan p1, pd1,po1;
  p1= fftw_plan_dft_c2r_2d(N1, N2, (fftw_complex*)&(I[0][0]), &(I[0][0]), FFTW_ESTIMATE);
  pd1= fftw_plan_dft_c2r_2d(N1, N2, (fftw_complex*)&(D[0][0]), &(D[0][0]), FFTW_ESTIMATE);
 
 //comp_ro = (fftw_complex *)&(b_ro[0][0]); // storing data visibilities 
 
  
  pf=fopen("modes","w"); // contains modes within U bins
  for(m=0;m<Nbin;++m)
  {
   count=0;
   for(i=0;i<N1;++i)
      for(j=0;j<N2/3;++j)//for(j=(N2-1)/2;j<N2;++j) 
      {
	index=i*N1+j;
	
	
  	//ia=((N1-1)/2-i);
  	//ja=((N2-1)/2-j);
  	
        ia=(i>N1/2) ? (N1-i) : i ;
        index_1=i*(N2/2+1) + j;
  	
  	//index_1=ia*(N2/2-1) + ja;
  
  	uu = sqrt(1.*(ia*ia+j*j))*Du; 
  	//uu = sqrt(1.*(ia*ia+ja*ja))*Du; 

  	if(uu>=uL[m] && uu<uH[m] && uu>0.) // checking which grid is within m-th ring
  	{
  	 d1[index_1][0]=GVRRre[index][nn]; d1[index_1][1]=GVRRim[index][nn]; // filling NG data into array
	 i1[index_1][0]=1.; i1[index_1][1]=0.;
	
  	 count++; // count modes within the bin
  	}
      }
    fftw_execute(p1); // perform FT to obtain I
    fftw_execute(pd1); // perform FT to obtain D
 
    /*--------- Writing D and I for every rings in files ----------*/
    snprintf(fname,128*sizeof(char),"masked_%d",m);
    pp=fopen(fname,"w");
    
    for(i=0;i<N1;++i)
      for(j=0;j<N2;++j)
           fwrite(&D[i][j],sizeof(double),1,pp);
           
    for(i=0;i<N1;++i)
      for(j=0;j<N2;++j)
          fwrite(&I[i][j],sizeof(double),1,pp);

    fclose(pp);
    /*-----------------------------------------------------------*/
    
    /*------- Initializing for the next bin --------*/
    for(i=0;i<N1;++i)
      for(j=0;j<N2;++j)
          {
            I[i][j]=0.;
            D[i][j]=0.;
	  }
    /*----------------------------------------------*/
    fprintf(pf,"%d\t%ld\n",m,count);
   }
 fclose(pf);
 fftw_destroy_plan(p1);
 fftw_destroy_plan(pd1);
}








void fft_init_orig(int Nbin, double *uL,double *uH)
{
  double uu;
  double **I, **D;
  int i,j,k,m,index,index_1,ia,ja, nn=0;
  int ic, jc, iagv;
  unsigned long int count;
  FILE *pp,*pf;
  char fname[128];

  I=allocate_2d_double(N1,N2+2); 
  D=allocate_2d_double(N1,N2+2);
  
  fftw_complex *i1,*d1, *comp_ro;
  i1=(fftw_complex*)&(I[0][0]);
  d1=(fftw_complex*)&(D[0][0]);
  
  fftw_plan p1, pd1;
  p1= fftw_plan_dft_c2r_2d(N1, N2, (fftw_complex*)&(I[0][0]), &(I[0][0]), FFTW_ESTIMATE);
  pd1= fftw_plan_dft_c2r_2d(N1, N2, (fftw_complex*)&(D[0][0]), &(D[0][0]), FFTW_ESTIMATE);

  //comp_ro = (fftw_complex *)&(ro[0][0]);
 
  pf=fopen("modes","w"); // contains modes within k bins
  for(m=0;m<Nbin;++m)
  {
   count=0;
   for(i=0;i<N1;++i)
    for(j=0;j<N2/3;++j)    //(j=0;j<N2/2+1;++j)
      {
  	ia=(i>N1/2) ? (N1-i) : i ;
        index=i*(N2/2+1) + j;
  	uu=sqrt(1.*(ia*ia+j*j));
  	
  	/*ic = i+N2/2;
        jc = j+N2/2;
        iagv=(ic>=N1) ? (ic-N1) : ic ;
        index_1 = iagv*N1+jc;*/
  	
  	//igv = i+(N2/2+1);
  	//jgv = j+(N2/2+1);
  	//iagv=(igv>N1/2) ? (N1-igv) : igv ;
  	//index_1=iagv*N1+jgv;
  	
  	
  	if(uu>=uL[m] && uu<uH[m] && uu>0.)
  	{
  	 //i1[index][0]=1.; i1[index][1]=0.;
  	 i1[index][0]=NTiG[index][nn]; i1[index][1]=0.;
  	 //d1[index][0]=GVRRre[index_1][nn]; d1[index][1]=GVRRim[index_1][nn];
  	 d1[index][0]=GVRRre[index][nn]; d1[index][1]=GVRRim[index][nn];// use if center at 0,0
  	 count++; // count modes within the bin
  	}
      }
    fftw_execute(p1); // perform FT
    fftw_execute(pd1); // perform FT
    
    snprintf(fname,128*sizeof(char),"masked_%d",m);
    pp=fopen(fname,"w");
    
    for(i=0;i<N1;++i)
      for(j=0;j<N2;++j)
           fwrite(&D[i][j],sizeof(double),1,pp);
           
    for(i=0;i<N1;++i)
      for(j=0;j<N2;++j)
          fwrite(&I[i][j],sizeof(double),1,pp);
    
    fclose(pp);
    /*------- Initializing for the next bin --------*/
    for(i=0;i<N1;++i)
      for(j=0;j<N2+2;++j)
          {
            I[i][j]=0.;
            D[i][j]=0.;
          }
    /*----------------------------------------------*/
    fprintf(pf,"%d\t%ld\n",m,count);
   }
 fclose(pf);
 fftw_destroy_plan(p1);
 fftw_destroy_plan(pd1);
}




















/*-------------------------------------------------------*/

//=======================================================================//
// Function to calculate the power spectrum from Fourier modes
void power_spec(int NB, double umax, double umin, double *usum,double *pksum,int *no)
// NB -- number of bins in output PS, *usum -- array with  mean U of each bin, *psum array with -- binned PS for each bin, *no --array with  number of modes in each bin  
{
  
  fftw_complex *in, *in_u;
  
  in=(fftw_complex*)&b_ro[0][0];
  in_u=(fftw_complex*)&b_ro_u[0][0];
  
  // local variables 
  double binsiz,uu, numer, denom;
  int i,j,index,ia,tmp;

  for(i=0;i<NB;++i)
    {
      no[i]=0;
      usum[i]=0.0;
      pksum[i]=0.0;
    }
  
  binsiz=log10(umax/umin)/(1.*NB);

  for(i=0;i<N1;++i)
    for(j=0;j<N2/2+1;++j)
      {
  	index=i*(N2/2+1)+j;
  	ia=(i>N1/2) ? (N1-i) : i ;
  
  	uu = sqrt(1.*(ia*ia+j*j))*Du; //Du is grid spacing in uv plane

  	if((uu>=umin)&&(uu<umax))
  	  {
  	    denom = in_u[index][0]*in_u[index][0]+in_u[index][1]*in_u[index][1];
  	    if(denom >= 1.)
  	    {
  	      numer = in[index][0]*in[index][0]+in[index][1]*in[index][1];
  	      tmp=(int) floor(log10(uu/umin)/binsiz);
  	      no[tmp]++;
	      pksum[tmp]+=(numer/denom);
  	      usum[tmp]+=uu;
  	    }
  	  }
      }
      

  for(i=0;i<NB;++i)
    if(no[i] != 0)
      {
	usum[i]/=(1.*no[i]);
	pksum[i]/=(1.*no[i]);
      }
}
/*-------------------------------------------------------------------------*/




void power_spec_GVRL(int NB, double umax, double umin, double *usum,double *pksum,int *no)
// NB -- number of bins in output PS, *usum -- array with  mean U of each bin, *psum array with -- binned PS for each bin, *no --array with  number of modes in each bin  
{
  
 // fftw_complex *in, *in_u;
  
 // in=(fftw_complex*)&b_ro[0][0];
  //in_u=(fftw_complex*)&b_ro_u[0][0];
  
  // local variables 
  
  double binsiz,uu, numer, denom, re_part, im_part;
  int i,j,index,ia,ja,tmp;
  int nn=0;
  int MM = 2*Nm;
  
  int Nlim;
  
  Nlim=Nu-2*MM;
   
  printf("\n Nlim=%d", Nlim); 
  printf("\n hiii");
  for(i=0;i<NB;++i)
    {
      no[i]=0;
      usum[i]=0.0;
      pksum[i]=0.0;
    }
  
  binsiz=log10(umax/umin)/(1.*NB);
  printf("\n %f",binsiz);
  for(i=0;i<Nlim;++i)
    for(j=0;j<Nlim;++j)
      {
  	index=i*Nlim+j;//i*((Nu-2*MM)/2+1)+j;
  	//ia=(i>(Nlim-1)/2) ? ((Nlim-1)/2-i) : i ;
  	//ja=(j>(Nlim-1)/2) ? ((Nlim-1)/2-j) : j ;
  	
  	ia=((Nlim-1)/2-i);
  	ja=((Nlim-1)/2-j);
  
  	uu = sqrt(1.*(ia*ia+ja*ja))*Du; //Du is grid spacing in uv plane

  	if((uu>=umin)&&(uu<umax))
  	  {
  	     re_part = GVRRre[index][nn] + GVLLre[index][nn];
  	     im_part = GVRRim[index][nn] + GVLLim[index][nn];
  	     numer =  re_part*re_part+im_part*im_part;//in[index][0]*in[index][0]+in[index][1]*in[index][1];
  	     tmp=(int) floor(log10(uu/umin)/binsiz);
  	     no[tmp]++;
	     pksum[tmp]+=(numer);
	     //if (GVRRre[index][nn]!=0 ){
	     //printf("\n %d",tmp);}
	     //if (tmp<5){
	     //printf("\n %d",tmp);}
  	     usum[tmp]+=uu;
  	  }
      }
      

  for(i=0;i<NB;++i)
    if(no[i] != 0)
      {
	usum[i]/=(1.*no[i]);
	pksum[i]/=(1.*no[i]);
      }
}
/*-------------------------------------------------------------------------*/













void power_spec_GV(int NB, double umax, double umin, double *usum,double *pksum,int *no)
// NB -- number of bins in output PS, *usum -- array with  mean U of each bin, *psum array with -- binned PS for each bin, *no --array with  number of modes in each bin  
{
  
 // fftw_complex *in, *in_u;
  
 // in=(fftw_complex*)&b_ro[0][0];
  //in_u=(fftw_complex*)&b_ro_u[0][0];
  
  // local variables 
  
  double binsiz,uu, numer, denom;
  int i,j,index,ia,ja,tmp;
  int nn=0;
  int MM = 2*Nm;
  
  int Nlim;
  
  Nlim=Nu-2*MM;
   
  printf("\n Nlim=%d", Nlim); 
  printf("\n hiii");
  for(i=0;i<NB;++i)
    {
      no[i]=0;
      usum[i]=0.0;
      pksum[i]=0.0;
    }
  
  binsiz=log10(umax/umin)/(1.*NB);
  printf("\n %f",binsiz);
  for(i=0;i<Nlim;++i)
    for(j=0;j<Nlim;++j)
      {
  	index=i*Nlim+j;//i*((Nu-2*MM)/2+1)+j;
  	//ia=(i>(Nlim-1)/2) ? ((Nlim-1)/2-i) : i ;
  	//ja=(j>(Nlim-1)/2) ? ((Nlim-1)/2-j) : j ;
  	
  	ia=((Nlim-1)/2-i);
  	ja=((Nlim-1)/2-j);
  
  	uu = sqrt(1.*(ia*ia+ja*ja))*Du; //Du is grid spacing in uv plane

  	if((uu>=umin)&&(uu<umax))
  	  {
  	     numer =  GVRRre[index][nn]*GVRRre[index][nn]+GVRRim[index][nn]*GVRRim[index][nn];//in[index][0]*in[index][0]+in[index][1]*in[index][1];
  	     tmp=(int) floor(log10(uu/umin)/binsiz);
  	     no[tmp]++;
	     pksum[tmp]+=(numer);
	     //if (GVRRre[index][nn]!=0 ){
	     //printf("\n %d",tmp);}
	     //if (tmp<5){
	     //printf("\n %d",tmp);}
  	     usum[tmp]+=uu;
  	  }
      }
      

  for(i=0;i<NB;++i)
    if(no[i] != 0)
      {
	usum[i]/=(1.*no[i]);
	pksum[i]/=(1.*no[i]);
      }
}
/*-------------------------------------------------------------------------*/



void power_spec_GV_Ntig(int NB, double umax, double umin, double *usum,double *pksum,int *no)
// NB -- number of bins in output PS, *usum -- array with  mean U of each bin, *psum array with -- binned PS for each bin, *no --array with  number of modes in each bin  
{
  
 // fftw_complex *in, *in_u;
  
 // in=(fftw_complex*)&b_ro[0][0];
  //in_u=(fftw_complex*)&b_ro_u[0][0];
  
  // local variables 
  
  double binsiz,uu, numer, denom;
  int i,j,index,ia,ja,tmp;
  int nn=0;
  int MM = 2*Nm;
  
  int Nlim;
  
  Nlim=Nu-2*MM;
   
  printf("\n Nlim=%d", Nlim); 
  printf("\n hiii");
  for(i=0;i<NB;++i)
    {
      no[i]=0;
      usum[i]=0.0;
      pksum[i]=0.0;
    }
    
    
      /*for(i=0;i<N1;++i)
    for(j=0;j<N2/2+1;++j)
      {
  	index=i*(N2/2+1)+j;
  	ia=(i>N1/2) ? (N1-i) : i ;*/
    
    
  
  binsiz=log10(umax/umin)/(1.*NB);
  printf("\n %f",binsiz);
  for(i=0;i<N1;++i)
    for(j=0;j<N2/2+1;++j)
      {
  	index=i*(N2/2+1)+j;
  	
	ia=(i>N1/2) ? (N1-i) : i ;
  
  	uu = sqrt(1.*(ia*ia+j*j))*Du; //Du is grid spacing in uv plane

  	if((uu>=umin)&&(uu<umax))
  	  {
  	     numer =  GVRRre[index][nn]*GVRRre[index][nn]+GVRRim[index][nn]*GVRRim[index][nn];//in[index][0]*in[index][0]+in[index][1]*in[index][1];
  	     numer *= NTiG[index][nn];
  	     tmp=(int) floor(log10(uu/umin)/binsiz);
  	     no[tmp]+= NTiG[index][nn];
	     pksum[tmp]+=(numer);
	     //if (GVRRre[index][nn]!=0 ){
	     //printf("\n %d",tmp);}
	     //if (tmp<5){
	     //printf("\n %d",tmp);}
  	     usum[tmp]+=(uu*NTiG[index][nn]);
  	  }
      }
      

  for(i=0;i<NB;++i)
    if(no[i] != 0)
      {
	usum[i]/=(1.*no[i]);
	pksum[i]/=(1.*no[i]);
      }
}
/*-------------------------------------------------------------------------*/






/*-------- Computing shapes for various sizes of closed triangles ---------*/
void shapes(int nbin, double *uavg)
{
  int i,j,l; 
  double MU,TT, uu1;
  char file[128];
  FILE *fp1;
 /*----- Evaluating all possible unique shapes of closed triangles -----*/  
  for(i=0;i<nbin;++i)
  {
   uu1=uavg[i];
   
   snprintf(file, 128*sizeof(char),"shape_u%d_%.3lf",i,uu1);
   fp1=fopen(file,"w");
   
   for(j=0;j<=i;++j)
    {
     TT=uavg[j]/uavg[i];
     if(TT>=0.5 && TT<=1.)
     {
      for(l=0;l<=j;++l)
      {
       MU=0.5*(TT+1./TT-(uavg[l]/uavg[i])*(uavg[l]/uavg[j]));
       if(MU>=0.5/TT && MU<=1.) // this incorporates 2*mu*t>=1
          fprintf(fp1,"%e\t%e\n",TT,MU); // writing output
      }
     }
    }
   fclose(fp1);
  } 
}

/*--------------------------------------------------------------------------*/

/*****************************************************/
/* Print out cfitsio error messages and exit program */
/*****************************************************/
void printerror( int status)
{
    if (status)
    {
       fits_report_error(stderr, status); /* print error report */

       exit( status );    /* terminate the program, returning error status */
    }
    return;
}
/*--------------------------------------------------------------------------*/


/************************************************************************/
/*           Read a UV data in Multifrequency FITS format               */
/*         This module reads both the NG data and UAPS data             */
/************************************************************************/
void read_uvdata(char filename1[], char filename2[])
{
    fitsfile *fptr, *fptr1;       /* pointer to the FITS file, defined in fitsio.h */
    int status, nfound, anynull;
    long naxes[3], fpixel[3], lpixel[3], inc[3]={1,1,1}, nbuffer, npixels, ii, jj, index, index1;
    long Nu, Nv, Nz, Nvis;
    double *buffer1, *buffer2;
    double cdelt[3], crpix[3], crval[3];
    void *nullval;
    int ia, ja;
    
    fftw_complex *in, *in_u, *in1, *in_u1;
    
    if(access(filename1, F_OK)!=0)
    {
      printf("Input File %s does not exists\n",filename1);
      exit(0);
    }
  
    if(access(filename2, F_OK)!=0)
    {
      printf("Input File %s does not exists\n",filename2);
      exit(0);
    }
    
    /*------------- Reading NG data file --------------*/
    status = 0;

    if ( fits_open_file(&fptr, filename1, READONLY, &status) )
         printerror( status );

    /* read the NAXIS1, NAXIS2 and NAXIS3 keywords to get image size */
    if ( fits_read_keys_lng(fptr, "NAXIS", 1, 3, naxes, &nfound, &status) )
         printerror( status );
    
    /* read the CDELT1, CDELT2 and CDELT3 keywords to get data size */
    if ( fits_read_keys_dbl(fptr, "CDELT", 1, 3, cdelt, &nfound, &status) )
         printerror( status );
         
    if ( fits_read_keys_dbl(fptr, "CRVAL", 1, 3, crval, &nfound, &status) )
         printerror( status );
    
    if ( fits_read_keys_dbl(fptr, "CRPIX", 1, 3, crpix, &nfound, &status) )
         printerror( status );     

    Nu = naxes[0]; Nv = naxes[1]; Nz = naxes[2];          /* number of grids along the three sides */
    Nvis  = Nu*Nv;         /* number of visiblities in one channel */
    //nullval  = 0;                /* don't check for null values in the image */
    
    printf("dU=%lf naxes vis={%ld,%ld,%ld}\n",cdelt[0],naxes[0],naxes[1],naxes[2]);
    
    //memory allocation for rho[][] array
    N1 = Nu;
    N2 = Nu;
   
    Du = fabs(cdelt[0]); // grid size in Fourier space in the unit of lambda
    
    ro = allocate_2d_double(N1,N2+1); //   ro = allocate_2d_double(N1,N2+1);,N2 is odd
    ro_u = allocate_2d_double(N1,N2+1);
    
    b_ro = allocate_2d_double(N1,N2+1);
    b_ro_u = allocate_2d_double(N1,N2+1);
    
    in = (fftw_complex*)&ro[0][0];
    in_u = (fftw_complex*)&ro_u[0][0];
    
    in1 = (fftw_complex*)&b_ro[0][0];
    in_u1 = (fftw_complex*)&b_ro_u[0][0];
    	
    buffer1 = (double*)calloc(Nvis,sizeof(double));
    buffer2 = (double*)calloc(Nvis,sizeof(double));
    
    fpixel[0]=1; fpixel[1]=1; fpixel[2]=1;
    lpixel[0]= Nu; lpixel[1]= Nv; lpixel[2]=1;
    
    if( fits_read_subset(fptr,TDOUBLE,fpixel,lpixel,inc,nullval,buffer1,&anynull,&status) )
         printerror( status );
         
         
    fpixel[0]=1; fpixel[1]=1; fpixel[2]=2;
    lpixel[0]= Nu; lpixel[1]= Nv; lpixel[2]=2;
    
    if( fits_read_subset(fptr,TDOUBLE,fpixel,lpixel,inc,nullval,buffer2,&anynull,&status) )
         printerror( status );
         
         
    for(ii=0;ii<N1;ii++)
      for(jj=0;jj<N2/2+1;jj++)
        {
         index = jj*N1+ii;
         in[index][0] = buffer1[index];
         in[index][1] = buffer2[index];
        }
    
   /*------------------ Shifting the U-plane origin for FFTW ------------------*/
    for(ii=0;ii<N1;++ii)
    {
      ia = ii-N1/2;
      ja = (ia<0) ? (ia+N1) : ia;
      for(jj=0;jj<N2/2+1;++jj)
          {
           index = jj*N1 + ii;
           index1 = jj + ja*(N2/2+1);
           in1[index1][0]=in[index][0];
           in1[index1][1]=in[index][1];
          }
     }
  /*---------------------------------------------------------*/
        
    free(buffer1);
    free(buffer2);
    if ( fits_close_file(fptr, &status) )
         printerror( status );
   /*------------- Reading UAPS data file --------------*/
   
   status = 0;

    if ( fits_open_file(&fptr1, filename2, READONLY, &status) )
         printerror( status );

    
    buffer1 = (double*)calloc(Nvis,sizeof(double));
    buffer2 = (double*)calloc(Nvis,sizeof(double));
    
    fpixel[0]=1; fpixel[1]=1; fpixel[2]=1;
    lpixel[0]= Nu; lpixel[1]= Nv; lpixel[2]=1;
    
    if( fits_read_subset(fptr1,TDOUBLE,fpixel,lpixel,inc,nullval,buffer1,&anynull,&status) )
         printerror( status );
         
         
    fpixel[0]=1; fpixel[1]=1; fpixel[2]=2;
    lpixel[0]= Nu; lpixel[1]= Nv; lpixel[2]=2;
    
    if( fits_read_subset(fptr1,TDOUBLE,fpixel,lpixel,inc,nullval,buffer2,&anynull,&status) )
         printerror( status );
    
    
    for(ii=0;ii<N1;ii++)
      for(jj=0;jj<N2/2+1;jj++)
        {
         index = jj*N1+ii;
         in_u[index][0] = buffer1[index];
         in_u[index][1] = buffer2[index];
        }
        
    /*------------------ Shifting the U-plane origin for FFTW------------------*/
    for(ii=0;ii<N1;++ii)
    {
      ia = ii-N1/2;
      ja = (ia<0) ? (ia+N1) : ia;
      for(jj=0;jj<N2/2+1;++jj)
          {
           index = jj*N1 + ii;
	   index1 = jj + ja*(N2/2+1);
           in_u1[index1][0]=in_u[index][0];
           in_u1[index1][1]=in_u[index][1];
          }
     }
  /*---------------------------------------------------------*/
    free(buffer1);
    free(buffer2);
    
    if ( fits_close_file(fptr1, &status) )
         printerror( status );
}
/*--------------------------------------------------------------------------*/






int initialize()
{
  int ii,jj,nn,sf,Ngrid, Ng; 

  sf = 2; //sampling factor
 // dU= 1.; // sampling the effective PB twice.
  Nm=3*sf; //sqrt(ln(500)/ln(2) = 3)

  Ng = 696;//474;//change thise //(int) ceil(Umax/dU) + Nm; // padding for convolution. 
  Nu = 2*Ng +1;                        
  Ngrid = (Nu-4*Nm)*(Nu-4*Nm);
  printf("\nNumber of grids = %d",(Nu-4*Nm));

  nua = (double*)calloc(nc,sizeof(double));
  for(nn=0;nn<nc;nn++)
    {
      nua[nn]= 154.41e6;//(nu_chan0 + del_chan*(nn - (chan0-1)));
    }
                    
  
  GVRRre = (double**)calloc(Ngrid,sizeof(double*));  
  GVRRim = (double**)calloc(Ngrid,sizeof(double*));  
  GVLLre = (double**)calloc(Ngrid,sizeof(double*));  
  GVLLim = (double**)calloc(Ngrid,sizeof(double*));  
  for(ii=0;ii<Ngrid;ii++)
    {
      GVRRre[ii]=(double*)calloc(nc,sizeof(double)); 
      GVRRim[ii]=(double*)calloc(nc,sizeof(double)); 
      GVLLre[ii]=(double*)calloc(nc,sizeof(double)); 
      GVLLim[ii]=(double*)calloc(nc,sizeof(double)); 
    }

  SCRR = (double**)calloc(Ngrid,sizeof(double*));  
  SCLL = (double**)calloc(Ngrid,sizeof(double*));  
  for(ii=0;ii<Ngrid;ii++)
    {
      SCRR[ii]=(double*)calloc(nc,sizeof(double)); 
      SCLL[ii]=(double*)calloc(nc,sizeof(double)); 
    }

 // printf("\nTotal no. of channels (nc) = %d",nc);
 // printf("\nFWHM = %.2f arcmin.",FWHM);
 // printf("\nGrid spacing, dU = %.2f wavelengths",dUx);
}


int initialize_gridding()
{
  int ii,jj,nn,sf,Ngrid,NN; 

  //theta_0=0.6*FWHM*M_PI/(180.*60.);         

  //sf = 2; //sampling factor
  //dU=sqrt(log(2))/(M_PI*theta_eff*sf); // sampling the effective PB twice.
  //Nm=3*sf; //sqrt(ln(500)/ln(2) = 3)

  //Ng = NN;//(int) ceil(Umax/dU) ; // padding for convolution. 
  //Nu = NN;//2*Ng +1;  
  //NN=(int) N1;     
  //Ngrid = NN*(NN/2+1);                
  Ngrid = N1*(N2/2+1);//Nu*Nu;
  printf("\n Ngrid = %d",Ngrid);

  nua = (double*)calloc(nc,sizeof(double));
  for(nn=0;nn<nc;nn++)
    {
      nua[nn]=154.41e6;// (nu_chan0 + del_chan*(nn - (chan0-1)));
    }
                    
  
  GVRRre = (double**)calloc(Ngrid,sizeof(double*));  
  GVRRim = (double**)calloc(Ngrid,sizeof(double*));  
  GVLLre = (double**)calloc(Ngrid,sizeof(double*));  
  GVLLim = (double**)calloc(Ngrid,sizeof(double*));  
  NTiG   = (double**)calloc(Ngrid,sizeof(double*));  
  for(ii=0;ii<Ngrid;ii++)
    {
      GVRRre[ii]=(double*)calloc(nc,sizeof(double)); 
      GVRRim[ii]=(double*)calloc(nc,sizeof(double)); 
      GVLLre[ii]=(double*)calloc(nc,sizeof(double)); 
      GVLLim[ii]=(double*)calloc(nc,sizeof(double)); 
      NTiG[ii]=(double*)calloc(nc,sizeof(double));
    }

  /*SCRR = (double**)calloc(Ngrid,sizeof(double*));  
  SCLL = (double**)calloc(Ngrid,sizeof(double*));  
  for(ii=0;ii<Ngrid;ii++)
    {
      SCRR[ii]=(double*)calloc(nc,sizeof(double)); 
      SCLL[ii]=(double*)calloc(nc,sizeof(double)); 
    }*/

}








void ReadGVnAC()
{ 
  int ii,jj,nn, index;
  char file_GV[256];
  FILE *fp;
  
  int MM = 2*Nm, Nlim;
  
  Nlim=Nu-2*MM;
 // printf("\nSize of the output GV: %d*%d*%d.",Nu-2*MM,Nu-2*MM,nc);  



   snprintf(file_GV, 256*sizeof(char),"../3_vcg/output_data/%s",GV); // open file to write bispectrum
   fp=fopen(file_GV,"rb");
  //fp=fopen(GV,"rb"); 
  fread(&FWHM,sizeof(FWHM),1,fp);
  fread(&dUx,sizeof(dUx),1,fp);
  fread(&fx,sizeof(fx),1,fp);
  for(ii=0;ii<nc;ii++)
    fread(&nua[ii],sizeof(nua[ii]),1,fp);

  for(ii=0; ii<Nlim;ii++)
    for(jj=0; jj<Nlim;jj++)
      {
	index=ii*Nlim+jj;
	for(nn=0; nn<nc; nn++)
	  {
	    fread(&GVRRre[index][nn],sizeof(GVRRre[index][jj]),1,fp);
	    fread(&GVRRim[index][nn],sizeof(GVRRim[index][jj]),1,fp);
	   // printf("\ni= %d, j= %d",ii,jj);
	   // printf("\nRe= %e, Im= %e",GVRRre[index][nn],GVRRim[index][nn]);
	  }
      }
      
      
   
  for(ii=0; ii<Nlim;ii++)
    for(jj=0; jj<Nlim;jj++)
      {
	index=ii*Nlim+jj;
	for(nn=0; nn<nc; nn++)
	  {
	    fread(&GVLLre[index][nn],sizeof(GVLLre[index][jj]),1,fp);
	    fread(&GVLLim[index][nn],sizeof(GVLLim[index][jj]),1,fp);
	  }
      }

  //Self-correlation. (real) 
  for(ii=0; ii<Nlim;ii++)
    for(jj=0; jj<Nlim;jj++)
      {
	index=ii*Nlim+jj;
	for(nn=0; nn<nc; nn++)
	  {
	    fread(&SCRR[index][nn],sizeof(SCRR[index][jj]),1,fp);
	    fread(&SCLL[index][nn],sizeof(SCLL[index][jj]),1,fp);
	    //printf("\nRe= %e, Im= %e",SCRR[index][nn],SCRR[index][nn]);
	  }
      }
  fclose(fp); 
}











void ReadGVnAC_gridding()
{
  int ii,jj,nn, index;
  char file_GV[256];
  FILE *fp;
  
   snprintf(file_GV, 256*sizeof(char),"../3_vcg/output_data/%s",GV); // open file to write bispectrum
   fp=fopen(file_GV,"rb"); 
 
 // fwrite(&FWHM,sizeof(FWHM),1,fp);
 // fwrite(&dU,sizeof(dU),1,fp);
 // fwrite(&f,sizeof(f),1,fp);
 // for(ii=0;ii<nc;ii++)
   // fwrite(&nua[ii],sizeof(nua[ii]),1,fp);

  for(ii=0; ii<N1;ii++)
    for(jj=0; jj<N2/2+1;jj++)
      {
	index = ii*(N2/2+1)+jj; 
	for(nn=0; nn<nc; nn++)
	  {
	    fread(&GVRRre[index][nn],sizeof(GVRRre[index][nn]),1,fp);
	    fread(&GVRRim[index][nn],sizeof(GVRRim[index][nn]),1,fp);
	    //printf("\nRe= %e, Im= %e",GVRRre[index][nn],GVRRim[index][nn]);
	  }
      }
      
       printf("GV[0][0]=%f\n",GVRRre[0][0]);
       fclose(fp); 
       
       
       
    fp=fopen("../3_vcg/output_data/sampling","rb");    
   for(ii=0; ii<N1;ii++)
    for(jj=0; jj<N2/2+1;jj++)
      {
	index = ii*(N2/2+1)+jj; 
	for(nn=0; nn<nc; nn++)
	  {
	    fread(&NTiG[index][nn],sizeof(NTiG[index][nn]),1,fp);
	  }
      }
      
       fclose(fp); 
       
       
    
 /* for(ii=0; ii<NN;ii++)
    for(jj=0; jj<NN/2+1;jj++)
      {
	index = ii*(NN/2+1)+jj;
	for(nn=0; nn<nc; nn++)
	  {
	    fwrite(&GVLLre[index][nn],sizeof(GVLLre[index][jj]),1,fp);
	    fwrite(&GVLLim[index][nn],sizeof(GVLLim[index][jj]),1,fp);
	  }
      }

  //Self-correlation. (real) 
  for(ii=0; ii<NN;ii++)
    for(jj=0; jj<NN/2+1;jj++)
      {
	index = ii*(NN/2+1)+jj;
	for(nn=0; nn<nc; nn++)
	  {
	    fwrite(&SCRR[index][nn],sizeof(SCRR[index][jj]),1,fp);
	    fwrite(&SCLL[index][nn],sizeof(SCLL[index][jj]),1,fp);
	    //printf("\nRe= %e, Im= %e",SCRR[index][nn],SCRR[index][nn]);
	  }
      }*/
      
      
 // fclose(fp); 
}
