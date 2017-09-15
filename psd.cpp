#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <conio.h>
#include <istream>
#include <fstream> 
#include <ios>
#include <string>

#include "equation.h"
#include "init.h"
#include "outdata.h"
#include "psd.h"
#include "imsl.h"

//double PSDfunction[200][2]; //Zr-Al
double PSDfunction[10000][2]; //Mo-Si

double tau0;

//top experimental profile
// double top_profile[100];  //Zr-Al
double top_profile[10000];  //Mo-Si

//profile of the one period of the sawtooth substrate
// double one_period[100]; // Zr-Al
double one_period[10000]; // Mo-Si

//profile of the one period of the sawtooth substrate for ".ggp" output 
//Zr-Al
double one_period_ggp[10000];

//final profile obtained with the theoretical model
double model_profile[10000];

//number of pints per one period of sawtooth profile
int N_one_period; 
//Zr-Al
int N_one_period_Zr_Al = 100;
//Mo-Si
int N_one_period_Mo_Si;
int N_one_period_W_B4C;

//multiplicator for .ggp profile 
double profile_mult;
double profile_mult_Zr_Al;
double profile_mult_Mo_Si;
double profile_mult_W_B4C;

double total_psd_thicknes;

//number of points per one period of sawtooth profile for ".ggp" output
int N_ggp;
//Zr-Al
int N_ggp_Zr_Al; 
//Mo-Si
int N_ggp_Mo_Si; 
int N_ggp_W_B4C;

//residials calculation data
//Zr-Al
//double residial[200];       // residial: residial[i] = hNew[i] - top_profile[i]
//Mo-Si
double residial[10000];       // residial: residial[i] = hNew[i] - top_profile[i]
double residial_sum;        //sum of squared residuals = Summ from 1 to N-1 of residial[i]^2 

//array of frequencies in PSD function
int length_array[10000];

void PSDCalculation ( void )
{
  int i,j,k,m;
  /*
  //Zr-Al and Mo-Si varuables
  int N_sum_partial_Zr_Al = 200;
  int N_sum_global_Zr_Al = 200;
  int N_sum_global_Mo_Si = 272;
  int N_sum_partial_Mo_Si = 272;
  int N_sum_partial;
  int N_sum_global;
  */
  
  //partial and global sums
  double sum_global[10000]; 
  double sum_partial[10000];
  double add_array_0[10000];
  double add_array_1[10000];
  double m_d, i_d;
  double N_d;
  
  FILE *PSDout;
  
  PSDout = fopen("psd_final", "w+");
  
  //setting parameters of the PSD function
  tau0 = 1.0;
  Length = N * tau0;
  N_d = (double)N;
  i_d = 0.0;
  
  
  //seting global and partial sums to zero
  for(i = 0; i < N; i++)
  {  
    sum_global[i] = 0;
    sum_partial[i] = 0;
    add_array_0[i] = 0; 
    add_array_1[i] = 0;
  }
  
  //calculating frequencies of the PSD function
  for(i = 0; i < N; i++)
  {
    m = -1 * N / 2 + i;
    PSDfunction[i][0] = m / Length;
    //i_d = i_d + 1.0;
    //m_d = -1.0 * N_d / 2.0 + i;
    //add_array_0[i] = m / Length;
    printf("%d calc part 1\n",i);
  }
  
  //calculating global sums
  for(k = 0; k < N; k++)
  {
    //seting partial sums to zero
    printf("%d calc part 2\n",k);
    /*
    for(i = 0; i < N; i++)
      sum_partial[i] = 0;
    */  
    m = -1 * N / 2 + k;  
      
    //calculating partial sums
    for(i = 0; i < N; i++)
      for(j = 0; j < N; j++)
      {
        //sum_global[k] = sum_global[k] + hNew[i] * hNew[j] * (cos(2 * 3.14 * m * i / N) * cos(2 * 3.14 * m * j / N) + sin(2 * 3.14 * m * i / N) * sin(2 * 3.14 * m * j / N));
        sum_global[k] = sum_global[k] + (hNew[i]-10.0) * (hNew[j]-10.0) * (cos(2 * 3.14 * m * i / N) * cos(2 * 3.14 * m * j / N) + sin(2 * 3.14 * m * i / N) * sin(2 * 3.14 * m * j / N));
        //sum_partial[i] = sum_partial[i] + hNew[i] * hNew[j] * (cos(2 * 3.14 * m * i / N) * cos(2 * 3.14 * m * j / N) + sin(2 * 3.14 * m * i / N) * sin(2 * 3.14 * m * j / N));
        //sum_partial[i] = sum_partial[i] + (hNew[i]-200.0) * (hNew[j]-200.0) * (cos(2 * 3.14 * m * i / N) * cos(2 * 3.14 * m * j / N) + sin(2 * 3.14 * m * i / N) * sin(2 * 3.14 * m * j / N));
        //printf("%f   %f\n",hNew[i], hNew[j]);
      }
    
    /*
    //calculating global sums
    for(i = 0; i < N; i++)  
      sum_global[k] = sum_global[k] + sum_partial[i];  
    */
  }  
  
  //calculating PSD function
  //printf("press any key\n");
  //getch();
  for(i = 0; i < N; i++)
  {
    if(i == 0 || i == N-1)
    {
      PSDfunction[i][1] = 1.0 * (tau0 / N) * sum_global[i]; //(tau0 / N) * (tau0 / N)
      //add_array_1[i] = 0.5 * (tau0 / N) * (tau0 / N) * sum_global[i];
    }
    else  
    {
      PSDfunction[i][1] = (tau0 / N) * sum_global[i]; 
      //add_array_1[i] = (tau0 / N) * (tau0 / N) * sum_global[i];
    }
    printf("%f   %f\n",PSDfunction[i][0], PSDfunction[i][1]);  
    
    fprintf(PSDout,"%f   %f\n",PSDfunction[i][0], PSDfunction[i][1]);  
    //fprintf(PSDout,"%f   %f\n",add_array_0[i], add_array_1[i]);  
      
  }
  
  printf("psd out");
  //getch();
 
  fclose(PSDout);

}


void PSDCalculationFFT ( void )
{
  int i,j,k;
  double input1[10000], input2[10000];
  double output1[10000], output2[10000];
  long int NNN = 1000;
  double freq;
  double factor, factor2;
  
  FILE *PSDout;
  
  
  
  for(i = 0; i < N; i++)
  {
    input1[i] = hNew[i] - 10.0;//deposited_thicknes
    //input2[i] = hNew[i]; 
  }  
  
  tau0 = 1.0;
  factor = tau0 / (double)N;
  factor2 = factor * factor;
  //printf("factor=%1.10f",factor);
  //getch();
  
  //fast furier transform
  DFFTRF(&NNN,input1,output1);
  
  PSDfunction[1][0] = 1.0 / N; 
  PSDfunction[1][1] = factor * output1[0] * output1[0]; 
  for(i = 2; i <= NNN / 2; i++)
  {
    freq = (double)i / N;
    //printf("%f\n",freq);
    //getch();
    PSDfunction[i][0] = freq;
    PSDfunction[i][1] = factor * (output1[2*i - 2] * output1[2*i - 2] + output1[2*i - 1] * output1[2*i - 1]);
  }  
  
  PSDout = fopen("psd_final", "w+");
  for(i = 1; i <= N / 2; i++)
    fprintf(PSDout,"%f   %f\n",PSDfunction[i][0], PSDfunction[i][1]);
  fclose(PSDout);  
}


void OnePeriodProfileData( void )
{
  int i, j, k;  
  
  double min_profile;         //minimal value of the profile height
  
  min_profile = hNew[0]; //hNew[0];
  j = 0;
  //printf("check5 output %f\n", min_profile);
  //getch();
  for(i = 0; i < N; i++)
  {
    //printf("%d  %f  %f\n",j, min_profile - hNew[i], min_profile);
    //getch(); 
    //printf("check hNew %f\n", hNew[i]); 
    if(min_profile > hNew[i])
    {
      //printf("minimal value is changed\n");
      //getch();
      min_profile = hNew[i];
      j = i;     
      //printf("check j output %d\n", j); 
      //getch();
      //getch();
    }  
  }
  
  //printf("check6 output\n"); 
  for(i = 0; i < N_ggp; i++)
  {
    k = j+i;
    //printf("check8 output %d\n", j); 
    //getch();
    if(k >= N)
      k = k - N;
    one_period_ggp[i] = hNew[k] - min_profile;
  }
  
  //Zr-Al
  //one_period_ggp[100]=0.0;
  //Mo-Si
  one_period_ggp[N_ggp]=0.0;
  for(i = 0; i <= N_ggp; i++)
    one_period_ggp[i] = one_period_ggp[i] / N_ggp;
    
    
}

void MacroProfileData( void )
{
	int i, j, k;  
	
	for(i = 0; i < N_ggp; i++)
	  one_period_ggp[i] = hNew[i] - hNew[0];
	
	one_period_ggp[N_ggp]=0.0;
    for(i = 0; i <= N_ggp; i++)
      one_period_ggp[i] = one_period_ggp[i] / N_ggp;
	
}

void ThicknessCalc( void )
{
   int i, j, k;

   ;
}


void FullPeriodProfileData( void )
{
  int i,j,k;  
  
  double min_profile;         //minimal value of the profile height
  FILE *FCh;
  
  min_profile = hNew[0]; //hNew[0];
  j = 0;
  for(i = 0; i < N; i++)
  {
    //printf("%d  %f  %f\n",j, min_profile - hNew[i], min_profile);
    //getch(); 
    
    if(min_profile > hNew[i])
    {
      //printf("minimal value is changed\n");
      //getch();
      min_profile = hNew[i];
      j = i;     
      
      //getch();
    }  
  }
  
  //FCh = fopen("full_period", "w+");
  for(i = 0; i < N; i++)
  {
  	k = j + i;
  	if(k >= N)
  	  k = k - N;
	hOld1[i] = hOld[k];
	hNew[i] = hOld1[i] - min_profile;
	//fprintf(FCh, "%d  %f\n", i, hOld1[i]);
  }
  
  //fclose(FCh);
}

void ProfileDeviation ( void )
{
  int i, j, k;
  
  double temp;
  double val;
  double min_profile;  //minimal value of the profile height
  
  
  
  int final_profile_flaq;


  std::fstream file1, file2;
  FILE *f1, *f2, *f3;

  file1.open("top",std::ios::in);
  printf("hello\n");
  
  final_profile_flaq = 2;
  
  i = 0;
  while(!file1.eof())
  {       
    file1 >> temp >> val;
    top_profile[i] = val * profile_mult / 2.0;
    i++;
  }
  
  file1.close(); 
  
  //reading final profile from the file
  if(final_profile_flaq = 2)
  {
    file2.open("out00000001.dat",std::ios::in);
    i = 0;
    while(!file2.eof())
    {   
      file2 >> temp >> val;
      model_profile[i] = val;
      i++;
    }
    for(i = 0; i < N; i++)
      hNew[i] = model_profile[i];
  
    file2.close(); 
  }  

  //testing
  f3 = fopen("final_model_profile","w+");
  for(i = 0; i < N; i++)
    fprintf(f3,"%d  %f\n",i,hNew[i]);
  fclose(f3);

  min_profile = hNew[0]; //hNew[0];
  for(i = 0; i < N; i++)
  {
    //printf("%d  %f  %f\n",j, min_profile - hNew[i], min_profile);
    //getch(); 
    ;
    if(min_profile > hNew[i] + 0.000000001)
    {
      printf("minimal value is changed\n");
      //getch();
      min_profile = hNew[i];
      j = i;     
      
      //getch();
    }  
  }
  
 
  
  for(i = 0; i < N_one_period; i++)
  {
    k = j+i;
    if(k > N)
      k = k - N;
    one_period[i] = hNew[k] - min_profile;
  }
  
  residial_sum = 0;  
  for(i = 0; i < N_one_period; i++)
  {
    residial[i] = one_period[i] - top_profile[i];
    residial_sum = residial_sum + residial[i] * residial[i];   
  }

  f1 = fopen("profile_data","w+");
  fprintf(f1,"%f\n",residial_sum);
  fclose(f1);
  
  f2 = fopen("model_profile","w+");
  for(i = 0; i < N_one_period; i++)
    fprintf(f2,"%d  %f\n",i,one_period[i]);
  fclose(f2);
  
}

//check calculation
void AutoCorrelationFunction( void )
{
  int i,j,k;
  double average;
  double AutoCorF[10000];
  FILE *f;
  
  average = 0.0;
  for(i = 0; i < N; i++)
  {
    average = average + hNew[i];
    AutoCorF[i] = 0.0;
  }
  
  average = average / N;
  for(i = 0; i < N; i++)
    for(j = 0; j < N; j++)
    {
      if((j+i) < N-1)
        AutoCorF[i] = AutoCorF[i] + hNew[j] * hNew[j+i+1];
      if((j+i) >= N-1)
      {
        k = j+i+1 - N;
        AutoCorF[i] = AutoCorF[i] + hNew[j] * hNew[k];
      }  
    }
  
   for(i = 0; i < N; i++)
     AutoCorF[i] = (AutoCorF[i] - average * average) / N;
  
  f = fopen("autocorr","w+");
  for(i = 0; i < N; i++)
    fprintf(f,"%d  %f\n",i,AutoCorF[i]);
    
  fclose(f);
}

void ProfileZshift( void )
{
	int i, j, k, l, ip2, ip3, ip1, im1, im2, im3;
	double p1;
	double S;
	double average_tg, average_angle;
	
	p1 = -1000.0;  
	S = 0; 	   
	i = 0;
	
	//looking for the first positive derivative
	while(p1 < 0)
	{
	  ip1 = i+1;
      ip2 = i+2;
      ip3 = i+3;
      im1 = i-1;
      im2 = i-2;
      im3 = i-3;
     
      if(i == 0)
      {          
        im3 = N-3;
        im2 = N-2;
        im1 = N-1;       
      }
      if( i == 1)
      {
        im3 = N-2;
        im2 = N - 1;       
      }
      if(i == 2)
        im3 = N-1;       
      if(i == N-3)
        ip3 = 0;
      if(i == N-2)
      {
        ip2 = 0;
        ip3 = 1;
      }    
      if(i == N-1)
      {
        ip1 = 0;
        ip2 = 1;
        ip3 = 2;
      } 
     
      p1 = (-1.0 * hOld[ip2] + 8.0 * hOld[ip1] - 8.0 * hOld[im1] + hOld[im2]) / (12.0 * dx);	  
      i++;
	}
	
	j = i - 1;
	i = j;
	
	//calculating mean angle, while derivative is positive, while on the working angle
	while(p1 > 0)
	{
	  ip1 = i+1;
      ip2 = i+2;
      ip3 = i+3;
      im1 = i-1;
      im2 = i-2;
      im3 = i-3;
     
      if(i == 0)
      {          
        im3 = N-3;
        im2 = N-2;
        im1 = N-1;       
      }
      if( i == 1)
      {
        im3 = N-2;
        im2 = N - 1;       
      }
      if(i == 2)
        im3 = N-1;       
      if(i == N-3)
        ip3 = 0;
      if(i == N-2)
      {
        ip2 = 0;
        ip3 = 1;
      }    
      if(i == N-1)
      {
        ip1 = 0;
        ip2 = 1;
        ip3 = 2;
      } 
     
      p1 = (-1.0 * hOld[ip2] + 8.0 * hOld[ip1] - 8.0 * hOld[im1] + hOld[im2]) / (12.0 * dx);
	  S = S + p1;	  
      i++;
	}
	
	k = i - 1;
	l = k - j;
	average_tg = S / (double)l;
	average_angle = atan(average_tg) * 180.0 / pi;
}
