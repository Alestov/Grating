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
#include "psd.h"


double pi = 3.1415926;
/* Zr-Al data */
/*
double hNew[200];
double hOld[200];
double dep[200];
//number of the spatial points
int N = 200;
*/

/* Mo-Si data */
double hNew[10000];
double hOld[10000];
double dep[10000];
double hNew_linear[10000];
double hOld_linear[10000];
double hOld1[10000];
//double ShiftedhOld[10000];
//number of the spatial points
int N;
//number of points in Zr-Al and Mo-Si material systems
int N_Zr_Al;
int N_Mo_Si;
int N_W_B4C;


/* parameters of the equation */
double nu1, nu2, nu3, nu4;
double lambda;
//end

//discretisation of time and spatial 
double dt = 0.00001;
double dx, dx2, dx3, dx4;
   
//atomic volume
double Volume;
//lattice parameter
double a_param;

//spatial length of the system
double Length;

/* equation discretization variables */
double p1, p2, p3, p4, p0;
double p2n2;  
double t1, t2, t3, t4, t0, tn;
/* equation discretization variables for linear equation */
double p1_linear, p2_linear, p3_linear, p4_linear, p0_linear, tn_linear;
//vkluchenie ili ne vkluchenie linejnogo rascheta
int linear_flaq;
//vybor tipa uravnenija: EW, KPZ, MKPZ, strict equation (EW,KPZ,MKPZ)
int nonlinearity_flaq;
  
//time inhomogenity of the deposition flux
int inhomogenity_time;
//length inhomogenity of the deposition flux
int inhomogenity_length;
//deposition fluctuation flux flaq
int inhomogenity_flaq;
int noise_length;
double max_noise;

//value of deposition inhomogenity
double inh_g;

//stationar deposition flux
double g;
//deposition flux angle
double alpha;
//flaq of the flux incident, 
int shadow_flaq;   //shadow_flaq = 1, deposition from the left to the right
                   //shadow_flaq = 2, deposition from the right to the left  

//profile angle
double profile_angle;

//derivatives data
double deriv_nonlinear[1000][6];
double deriv_linear[1000][6];
//linear equation difference scheme data
double t_linear[1000][6];
//nonlinear equation difference scheme data
double t_nonlinear[1000][6];

double pre_dep1, pre_dep2, pre_dep3;


double c_k = 1.0;

double H0;

FILE *fW;

//total thikness of the deposited layer with tilted deposition
double total_dep;
//average thikness of the deposited layer with tilted deposition
double average_dep;

// calculation of equation solution
void Solution( int time_stp  )
{
   int i,j, ip2, ip3, ip1, im1, im2, im3, k, l;
   FILE *f;
   double chisl1, chisl2, znam1, krivizna1, krivizna2, znam2;
   double depos_rate;
   
   double Deriv1[10000], Deriv2[10000], CurveR[10000], Chisl1[10000];
     
   //setting the spatial step discretisation
   //N = (int)(Length/a_param);
   //Length = 200.0;
   //Length = 272.0
   //f = fopen("calc","w+");   
   
   // nonlinear terms + linear exponents 
   for(i = 0; i < N; i++)
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
     
       //deposition assuming angle of the sawtooth
       
       // derivatives
       
	   p1 = (-1.0 * hOld[ip2] + 8.0 * hOld[ip1] - 8.0 * hOld[im1] + hOld[im2]) / (12.0 * dx);
       p2 = (-1.0 * hOld[ip2] + 16.0 * hOld[ip1] - 30.0 * hOld[i] + 16.0 * hOld[im1] - hOld[im2]) / (12.0 * dx2);       
       p3 = (hOld[im3] - 8.0 * hOld[im2] + 13.0 * hOld[im1] - 13.0 * hOld[ip1] + 8.0 * hOld[ip2] - hOld[ip3])/(8.0 * dx3);
       p4 = (-1.0 * hOld[im3] + 12.0 * hOld[im2] - 39.0 * hOld[im1] + 56.0 * hOld[i] - 39.0 * hOld[ip1] + 12.0 * hOld[ip2] - hOld[ip3]) / (6 * dx4);
       //p0 = sqrt( 1.0 + p1 * p1 );
       
	   
       //nonlinear term (second derivative from gradient square)
       //p2n2 = 2.0 * p1 * p3 + 2.0 * p2 * p2;
       /*
       p1 = (-1.0 * hOld1[ip2] + 8.0 * hOld1[ip1] - 8.0 * hOld1[im1] + hOld1[im2]) / (12.0 * dx);
       p2 = (-1.0 * hOld1[ip2] + 16.0 * hOld1[ip1] - 30.0 * hOld1[i] + 16.0 * hOld1[im1] - hOld1[im2]) / (12.0 * dx2);       
       p3 = (hOld1[im3] - 8.0 * hOld1[im2] + 13.0 * hOld1[im1] - 13.0 * hOld1[ip1] + 8.0 * hOld1[ip2] - hOld1[ip3])/(8.0 * dx3);
       p4 = (-1.0 * hOld1[im3] + 12.0 * hOld1[im2] - 39.0 * hOld1[im1] + 56.0 * hOld1[i] - 39.0 * hOld1[ip1] + 12.0 * hOld1[ip2] - hOld1[ip3]) / (6 * dx4);
       p0 = sqrt( 1.0 + p1 * p1 );
       */
     
     //Edward Wilkinson equation 
     if(nonlinearity_flaq == 4)
     {
       t0 = dep[i] * dt; 
       t2 = nu2 * p2 * dt;
       t4 = 0.0;
       tn = 0.0; 
       hNew[i] = hOld[i] + t0 + t2;
     }
    
     //modified Edward Wilkinson equation 
     if(nonlinearity_flaq == 5)
     {
       t0 = dep[i] * dt; 
       t2 = nu2 * p2 * dt;
       t4 = nu4 * p4 * dt;
       tn = 0.0; 
       hNew[i] = hOld[i] + t0 + t2 - t4;
     }
     
     //Edward Wilkinson equation with shadowing
     if(nonlinearity_flaq == 6)
     {
       //deposition from the left to right
       if(shadow_flaq == 1)
       {
         if(p1 >= 0.0 + 0.0000001)
         {
           profile_angle = atan(p1); 
           dep[i] = dep[i] * cos(pi / 2.0 - alpha - profile_angle);
         }    
       
         if(p1 < 0.0 + 0.0000001)
         {
           profile_angle = atan(fabs(p1)); 
           dep[i] = dep[i] * sin(alpha - profile_angle);
         }    
       }
       //deposition from the right to the left
       if(shadow_flaq == 2)
       {
         //see nonlinearity_flaq=7
         if(p1 < 0.0 + 0.0000001)
         {
           //g = 0.5 / cos(pi / 2.0 - alpha - profile_angle);
           profile_angle = atan(fabs(p1)); 
           dep[i] = dep[i] * cos(pi / 2.0 - alpha - profile_angle);
           
         }    
       
         if(p1 >= 0.0 + 0.0000001)
         {
           profile_angle = atan(fabs(p1)); 
           dep[i] = dep[i] * sin(alpha - profile_angle);
         }    
       }      
       t0 = dep[i] * dt; 
       t2 = nu2 * p2 * dt;
       t4 = 0.0;
       tn = 0.0; 
       hNew[i] = hOld[i] + t0 + t2;
     }
     
     
      //stogoe uravnenie dlja evolucii profilja (strict MKPZ)
     if(nonlinearity_flaq == 7)
     {
        //for the case of the slopes, deposition array should not be changed
	    depos_rate = dep[i];
        if(p1 < 0.0 + 0.0000001)
         {
           profile_angle = atan(fabs(p1)); 
           depos_rate = dep[i] * cos(pi / 2.0 - alpha - profile_angle);
           //printf("%d  %f  %f  \n", i,  g, dep[i]);
           //getch();
         }    
       
         if(p1 >= 0.0 + 0.0000001)
         {
           profile_angle = atan(p1); 
           depos_rate = dep[i] * sin(alpha - profile_angle);
         }    
         
       //pre_dep  
       //average_dep = average_dep + dep[i] * dt;
       t0 = depos_rate * dt; 
       t2 = 1.0 * nu2 * p2 * dt / (1.0 + p1 * p1);
       znam1 = (1.0 + p1 * p1) * (1.0 + p1 * p1) * (1.0 + p1 * p1);
       znam2 = (1.0 + p1 * p1) * (1.0 + p1 * p1) * (1.0 + p1 * p1) * (1.0 + p1 * p1);
       krivizna1 = (1 + p1 * p1);
       chisl1 = -1.0 * p4 * krivizna1 + 4.0 * p1 * p2 * p3 + 3 * p2 * p2 * p2;
       chisl2 = 6.0 * (-1.0 * p3 * krivizna1 + 3.0 * p1 * p2 * p2) * p1 * p2;
       t4 = nu4 * dt * (chisl1 / znam1 - chisl2 / znam2);
       Deriv1[i] = p1;
       Deriv2[i] = p2;
	   CurveR[i] = chisl1 / znam1;
	   Chisl1[i] = chisl1;
       //tn = lambda * p1 * p1;
	   //tn = depos_rate * sqrt(1.0 + p1 * p1) * dt;        
       //hNew[i] = hOld[i] + t0 + t2 + t4 + tn;
       //hNew[i] = hOld[i] + t0 + t2;
       //hNew[i] = hOld[i] + tn + t2;
       //t2 = nu2 * p2 * dt;
       //t4 = nu4 * p4 * dt;
       hNew[i] = hOld[i] + t0 + t4;
       
      
       //printf("%d  %f  %f  \n", i,  g, dep[i]);
       //getch();
       //hNew[i] = hOld[i] +  t0;
       
       //printf("%d  %f  %f  \n", i,  g, dep[i]);
       //getch();
       //calc add
       //hNew[i] = hOld1[i] + t0 / H0 + t2 + t4 + tn * H0;
       //hNew[i] = hOld1[i] + t0 / H0 + nu2 * p2 * dt + nu4 * p4 * dt + lambda * H0 * p1 * p1 * dt;
       //printf("t0 / H0 = %f\n nu2 * p2 * dt = %f \nnu4 * p4 * dt = %f \nlambda * H0 * p1 * p1 *dt = %f\n", t0 / H0, nu2 * p2 * dt, nu4 * p4 * dt, lambda * H0 * p1 * p1 * dt);
       //printf("hNew[%d] = %f hOld[i] = %f  t0 = %f  %f\n", i, hOld[i], hNew[i], t0, hNew[i] - hOld[i]);
	   //getch();
       
     }
     
     /* In double modified KPZ equation calculations stop due to -1.#IND
        double modified KPZ equations: KPZ - grad^2_h_; KPZ+grad^4_h_ 
        double KPZ is not used
     */
     
     //printf("hNew=%f\n",p1);
     /*
     fprintf(f,"i = %d\n",i);
     fprintf(f,"p0 = %1.15f\n",p0);
     fprintf(f,"p1 = %1.10f\n",p1);
     fprintf(f,"p2 = %1.10f\n",p2);
     fprintf(f,"p3 = %1.10f\n",p3);
     fprintf(f,"p4 = %1.10f\n",p4);
     fprintf(f,"p2n2 = %1.10f\n",p2n2);
     fprintf(f,"p1*p3= %1.10f\n",p1*p3);
     fprintf(f,"p2*p2 = %1.10f\n",p2*p2);
     fprintf(f,"t0 = %1.15f\n",t0);
     fprintf(f,"t2 = %1.10f\n",t2);
     fprintf(f,"t4 = %1.10f\n",t4);
     fprintf(f,"tn = %1.10f\n",tn);
     fprintf(f,"\n");  
     */   
   }
   
   if(time_stp == 100)
   {
      FILE *FchF, *FchF1;
      int ij;
      FchF = fopen("derivative","w+");
      FchF1 = fopen("profile_1","w+");
      for(ij = 0; ij < N; ij++)
      {
        fprintf(FchF,"%d  %f  %f  %f  %f\n",ij, Deriv1[ij], Deriv2[ij], CurveR[ij], Chisl1[ij]);
        fprintf(FchF1,"%d  %f\n", ij, hNew[ij]);
      }
      fclose(FchF);  
      fclose(FchF1);
	}
   
   
   //obnovlenie massiva profilja podlozhki   
   for(i = 0; i < N; i++)
   {
     //printf("Check for renewal %f\n", hOld[i]);
     //printf("hNew[0] = %1.10f  hOld[0] = %1.10f  t0 = %1.10f\n", hNew[0], hOld[0], t0);
	 hOld[i] = hNew[i];   
	 //printf("Check for renewal  %f   %f\n", hOld[i], hNew[i]);
	 //getch();
	 //hOld1[i] = hNew[i];    
     //hOld_linear[i] = hNew_linear[i];
   }  
   //fclose(f);
   
}

void WhiteNoiseDeposition( int time_step )
{
  int i,j,k;
  //number of points per current deposition inhomogenity
  int lengthPoints;
  //current point in the x-direction
  int counterPoints;
  int rand_dep_length;
  //distance between discretisation points
  double distX = 1.0;
  double Curr_Length;
  int rand_freq;
  double rand_max = (double)RAND_MAX; 
  double rand_noise_amplitude = (double)RAND_MAX; 
  double half_amplitude = (double)RAND_MAX / (2.0 * rand_noise_amplitude);
    
  double inh_g2;
  
  Length = N * distX;
  Curr_Length = Length;
  counterPoints = 0;
  //printf("%f\n",Curr_Length);
  //getch();
  
  
    while( Curr_Length > 0.0 + 0.00000001)
    {
      rand_dep_length = rand()%noise_length + 1;
      //printf("%d\n",rand_dep_length);
      //getch();
      length_array[rand_dep_length] = length_array[rand_dep_length] + 1;
      Curr_Length = Curr_Length - rand_dep_length;
      lengthPoints = rand_dep_length;
      //printf("%d\n",lengthPoints);
      //getch();
      
      //inh_g = max_noise * 2.0 * (rand()/rand_noise_amplitude - half_amplitude);
      inh_g = max_noise * (rand()/rand_noise_amplitude - half_amplitude);
      for(i = counterPoints; i < counterPoints + lengthPoints; i++)
      {
        dep[i] = g + inh_g * g;
        //printf("%d  %f  %f  %f\n", i,  g, inh_g, dep[i]);
        
      }
      counterPoints = counterPoints + lengthPoints;  
    }
  
  //printf("noise\n");
  //getch();
  /*
  for(i = 0; i < N; i++)
  {
  	printf("%d  %f\n", i, dep[i]);
  	getch(); 
  }
  */
}



