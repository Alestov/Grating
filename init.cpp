#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <conio.h>
#include <istream>
#include <fstream> 
#include <ios>
#include <string.h>
#include <iostream>

#include "equation.h"
#include "init.h"
#include "outdata.h"
#include "init.h"
#include "psd.h"
#include "imsl.h"

using namespace std;

int sawtooth_period;


void SawToothInterpolationL ( int NSlope1, int NSlope2, int NSlopeD1, int NSlopeD2)
{
  int i, j, k;
  // sawtooth parameters 
  //number of sawtooths
  int saw_num;
  //angles of sawtooths
  double saw_angle1,    // lest sawtooth angle
         saw_angle2;    // right sawtooth angle
  // length of the sawtooth basis
  double saw_length; 
  //number of points for the left sawtooth slope
  int saw_slope;
  //number of points per one sawtooth
  int saw_points; 
  //multipling factor = saw_points  / saw_length(nm)  
  double saw_factor;
  long int NData;
  double XData[6], FData[6];
  //long int NumberDelta = 37; //NSlope1 + NSlope2 + NSlopeD1 + NSlopeD2; //NSlope1 + NSlope2 + NSlopeD1 + NSlopeD2;
  int NumberDelta = NSlope1 + NSlope2 + NSlopeD1 + NSlopeD2; //NSlope1 + NSlope2 + NSlopeD1 + NSlopeD2;
  long int NSize;
  //double XVec[38], NValue[38];
  double XVec[NumberDelta+1], NValue[NumberDelta+1];
  double CurrhOld[10000];
  FILE * hO;
  int NCurrS;
  
  
  for(i = 0; i < N; i++)
    CurrhOld[i] = hOld[i];
    
  
  
   hO = fopen("profileSpline","w+");   
     //printf("%f", CurrhOld[225]);
     //getch();
    //calculation of the spline
    N = N + 2 * (NSlopeD1 + NSlopeD2);
    
	NData = 6;
	XData[0] = double(N_one_period - NSlope1);
    XData[1] = double(N_one_period + NSlopeD1);
    XData[2] = double(N_one_period + NSlopeD1 + 2);
    XData[3] = double(N_one_period + NSlopeD1 + 4);
    XData[4] = double(N_one_period + NSlopeD1 + 6);
    XData[5] = double(N_one_period + NSlope2 + NSlopeD1 + NSlopeD2);
    //printf("XData[0] = %f,  XData[1] = %f  XData[2] = %f  XData[3] = %f\n", XData[0], XData[1], XData[2], XData[3]);
    //getch();
	FData[0] = CurrhOld[N_one_period - NSlope1];
    FData[1] = 0;
    FData[2] = 0.02; //CurrhOld[N_one_period + 5];
    FData[3] = 0.055; //CurrhOld[N_one_period + 10];
    FData[4] = 0.1;
    FData[5] = CurrhOld[N_one_period + NSlope2];
    //printf("FData[0] = %f,  FData[1] = %f  FData[2] = %f  FData[3] = %f\n", FData[0], FData[1], FData[2], FData[3]);
    //getch();
    //printf("N = %d\n", N);
    //getch();
    for(i = 0; i <= NumberDelta; i++)
    {
      XVec[i] = N_one_period - NSlope1 + i;
      //printf("FData[2] = %f \n", FData[int(XVec[]);
      //getch();
    }
    
    //printf("Check\n");
    //getch();
    NSize = NumberDelta + 1;
    DCSIEZ(&NData, XData, FData, &NSize, XVec, NValue);
    FILE *FCur;
    FCur = fopen("NVALUE", "w+");
    printf("NumberDelta = %d", NumberDelta);
    getch();
	
	for(i = 0; i <= NumberDelta; i++)
    {
      
      printf("NValue[%d] = %f \n", i, NValue[i]);
      //fprintf(FCur, "NValue[%d]    %f XVec[%d]    %f\n", i, NValue[i], i, XVec[i]);
      //getch();
    }
    
    
    for(i = 0; i <= NSlopeD2 + NSlope2; i++)
    {
      //hOld[i] = NValue[i + NSlope1 + NSlopeD1];
       hOld[i] = NValue[i + NSlope1 + NSlopeD1];
      //printf("%d   %f\n", i, CurrhOld[i]);
      //getch();
    }
    
    //for(i = 0; i < N_one_period - NSlope1; i++)
    for(i = NSlopeD2 + NSlope2 + 1; i < N_one_period - NSlope1; i++)
    {
      hOld[i] = CurrhOld[i];
      fprintf(FCur, "%d    %f \n", i, CurrhOld[i]);
      //printf("%d   %f\n", i, CurrhOld[i]);
      //getch();
    }
    
	//hOld[N_one_period - NSlope1] = (CurrhOld[N_one_period - NSlope1 - 1] + NValue[1]) / 2;
    for(i = N_one_period - NSlope1 + 1; i <= N_one_period - NSlope1 + NumberDelta; i++)
	{
	  NCurrS = i - (N_one_period - NSlope1);
	  hOld[i] = NValue[NCurrS];
	  //fprintf(hO,"%d   %f\n", i, hOld[i]);	
	  //printf("%d   %f\n", i, CurrhOld[i]);
	  //getch();
	  
    }  
    //hOld[N_one_period - NSlope1 + NumberDelta] = (NValue[NumberDelta - 1] + CurrhOld[N_one_period + NSlope2 + 1]) / 2;
    //for(i = N_one_period - NSlope1 + NumberDelta + 1; i < N; i++)
    for(i = N_one_period - NSlope1 + NumberDelta + 1; i < N - NSlopeD1 - NSlope1; i++)
	{
	  hOld[i] = CurrhOld[i - NSlopeD1 - NSlopeD2];
	  //printf("i = %d   %f\n", i - NSlopeD1 - NSlopeD2, CurrhOld[i - NSlopeD1 - NSlopeD2]);
	  //getch();
	  //fprintf(hO,"%d   %f\n", i, hOld[i]);	
    }
    
    for(i = N - NSlopeD1 - NSlope1; i < N; i++)
	{
	  NCurrS = i - (N - NSlopeD1 - NSlope1);
	  hOld[i] = NValue[NCurrS];
	  //printf("i = %d   %f\n", i - NSlopeD1 - NSlopeD2, CurrhOld[i - NSlopeD1 - NSlopeD2]);
	  //getch();
	  //fprintf(hO,"%d   %f\n", i, hOld[i]);	
    }
    
    for(i = 0; i < N; i++)
      fprintf(hO,"%d   %f\n", i, hOld[i]);
    
    fclose(hO);
    fclose(FCur);
    //double Xvec[]
    
}


//initialisation of the hOld and hNew: read from file or generating profile
void Init ( double ang1, double ang2 )
{
  int i,j;
  double rand_max;
  double val;
  FILE *hO;
  FILE *htemp;
  FILE * profF;
  double temp;
  int num[100];
  double read[500];
  int k;
  int psd_calc_flaq;
  
  /* sawtooth parameters */
  //number of sawtooths
  int saw_num;
  //angles of sawtooths
  double saw_angle1,    // lest sawtooth angle
         saw_angle2;    // right sawtooth angle
  // length of the sawtooth basis
  double saw_length; 
  //number of points for the left sawtooth slope
  int saw_slope;
  //number of points per one sawtooth
  int saw_points; 
  //multipling factor = saw_points  / saw_length(nm)  
  double saw_factor;
  
  //profile choise flaq
  int profile_flaq;
  int read_flaq;
  int mult_flaq;
  int read_from_file = 0;
  int regular_profile_flaq; // flag for setting profile initial point to (0,0)
  
  std::fstream file;
  std::string str;
  char inputfile[20];
  int linesize = 512;
  char twolines[linesize];
  
  dx = 1.0; //Length / N; 
  dx2 = dx * dx;
  dx3 = dx * dx * dx;
  dx4 = dx * dx * dx * dx;
  
  Volume = 0.3;
  a_param = pow(Volume,0.33);
  Length =N * dx; 
  dx = 1.0; 
   
  for(i = 0; i < N; i++)
  {
    hOld[i] = 0.0;
    hNew[i] = 0.0;   
    hOld_linear[i] = 0.0;
    hNew_linear[i] = 0.0;
    length_array[i] = 0;
  } 
    
  profile_flaq = 2;
  read_flaq = 1;
  psd_calc_flaq = 1;
  mult_flaq = 0;
  regular_profile_flaq = 1;
  
  //reading profile from file 
  if(profile_flaq == 2)
  {
    
    if(read_flaq == 2)
    {
      cout << "Input file name" << endl;
      cin >> inputfile;
      cout << inputfile;
      file.open(inputfile,std::ios::in);
      getline(file,str);
      cout << str << endl;
      getline(file,str);
      cout << str << endl;
      getline(file,str);
      cout << str << endl;
      //getch();
    }
    
    if(read_flaq == 1)  
      //file.open("profile",std::ios::in);
      file.open("profile.dat",std::ios::in);
      
    i = 0;
  
    while(read_from_file != 1) //while(!file.eof())
    {   
      printf("%d  read...\n",i);
      //getch();
      file>> temp >> val;
      //Zr-Al system
      //hOld[i] = val * 200.0;
      
      //Mo-Si system
      if(mult_flaq == 1)
        hOld[i] = val * profile_mult;
      if(mult_flaq == 0)
        hOld[i] = val;
      if(linear_flaq == 1)
        hOld_linear[i] = val * profile_mult;
         
      i++;
      if(i == N)
        read_from_file = 1;
      
      
    }
    file.close(); 
    
    /* reading onw period of profile */
    /*
    for(i = 1000; i < N; i++)
      hOld[i] = hOld[i - 1000];
    */
    
    htemp = fopen("profile_copy","w+");
    //setting initial profile point to (0,0) if necessary
    if(regular_profile_flaq == 1)
    {
	  for(i = 0; i < N; i++)
	  {
	  	hNew[i] = hOld[i];
	  }
	  FullPeriodProfileData();
      for(i = 0; i < N; i++)
        fprintf(htemp,"%d  %1.10f\n",i,hNew[i]);
	}
	
    if(regular_profile_flaq != 1)
    {
      for(i = 0; i < N; i++)
      {
        fprintf(htemp,"%d  %1.10f\n",i,hOld[i]);
        if(psd_calc_flaq == 1)
          hNew[i] = hOld[i];
      }
      fclose(htemp);
    } 
  
   }
   
   
  //generating sawtooth profile
  if(profile_flaq == 3)
  {
      //hO = fopen("hOld","w+");
      // setting sawtooth ansamble parameters
      saw_num = 2;
      saw_length = double(N_one_period_W_B4C);
      //printf("saw_slope = %d\n", saw_slope);
      //getch();
      //printf("Check1\n");
      //getch();
      //saw_angle1 = 6.0;
      //saw_angle2 = 6.0; 
      saw_angle1 = ang1;
      saw_angle2 = ang2; 
      saw_points = N_one_period_W_B4C;
      saw_factor = saw_points / saw_length;
      saw_slope = (int)(saw_factor * saw_length * tan(3.14 * saw_angle2 / 180.0) / (tan(3.14 * saw_angle1 / 180.0) + tan(3.14 * saw_angle2 / 180.0)));
      //H0 = saw_slope * tan(3.14 * saw_angle1 / 180.0);
	  //printf("saw_slope = %d\n", saw_slope);
      //getch();
	  k = 1;
      hO = fopen("profileTriangle","w+");   
      //calculation of sawtooth discretisation for dx = 1 nm
      while(k <= saw_num)
      {
        for(i = N * (k-1) / saw_num; i <= N * (k-1) / saw_num + saw_slope; i++)
          hOld[i] = tan(3.14 * saw_angle1 / 180.0) * (i - N * (k - 1) / saw_num);
        for(i = N * (k-1) / saw_num + saw_slope + 1; i < N * k / saw_num;i++)
          hOld[i] = (N * k / saw_num - i) * tan(3.14 * saw_angle2 / 180.0);  
        k++;  
      }   
        
		//printf("Check1\n");
        //getch();
		//SubsRms();
        //SawToothInterpolationL(5,7,15,10);
        for(i = 0; i < N; i++)
          fprintf(hO,"%d   %f\n", i, hOld[i]);
        fclose(hO);
		
   }   
  
}
