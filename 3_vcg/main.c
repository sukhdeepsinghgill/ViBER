# include <stdio.h>
# include <math.h>
# include <string.h>
# include <stdlib.h>
# include <unistd.h>
# include "tgefuncs.h"

char INFITS[128], input[128], GV[128],file_ntig[128];
int NN; 
double dU,theta_0,FWHM,Umax;


int main(int argc, char *argv[])
{ 
  int key_ntig;
  
  if(argc!=6) 
    {
      printf("Usage: %s <input FITS file> <input parameters file> <output GV file> <want sampling_fle?> <output sampling_file>\n", argv[0]);
      return 1;
    }
  sscanf(argv[1],"%s",INFITS);
  sscanf(argv[2],"%s",input); 
  sscanf(argv[3],"%s",GV); 
  sscanf(argv[4],"%d",&key_ntig); 
  sscanf(argv[5],"%s",file_ntig); 
  
  printf("\nkey_ntig=%d\n",key_ntig);
  
  read_inputs(input); 
 
  
  //initialize(); 
  initialize_gridding();

  //group_loop();
  group_loop_gridding();
  
  //WriteGVnAC();
  WriteGVnAC_gridding();
  
  if (key_ntig==1)
   {Write_NTiG();}
  
 
  
  
  printf("\nEnd.\n");
  return 0;   
} 
