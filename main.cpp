/*
Program 
1. Solution of nonlinear kinetic equation (modified KPZ)
dh/dt= nu2 * d2h/dx2  - nu4 * d4h / dx4 + g(r,t) + lambda*d2((d2h)^2)

2. Numerical details
Central finite difference method, O(dx^4)
Periodic boundary conditions

3. Program version CDS.dx4.MKPZ.BLD.PSD

CDS - central differencies scheme 
dx4 with dx^4 accuracy
MKPZ - modified KPZ
BLD - bilayer deposition 
PSD - PSD function calculation
*/

#include <stdio.h>
#include <conio.h>
#include <time.h>
#include <conio.h>
#include <stdlib.h>
#include <math.h>

#include "imsl.h"

#include "equation.h"
#include "outdata.h"
#include "init.h"
#include "psd.h"

int fl;

int global_steps;

int main ( void )
{
  double check = 0.0;
  /* output parameters */
  double stepOut;
  double stepOut_first_layer;        // number of steps before first layer output
  double stepOut_second_layer;
  double stepOut_first_interlayer;
  double stepOut_second_interlayer;
  int stepOutflaq;
  double bilyaer_output;
  double double_layer_counter;
  double blaze_angle;
  double second_angle;
  double buffer_height;
  int buffer_stepout;
         
  double dj = 1.0;
  double c_temp1, c_temp2;
  int tott = 100;
  
  /* parameters of the deposited bilayer */
  int bilayer_flaq;                   // bilayer_flaq = 1 deposition of bilayers else without 
  double deposited_thicknes;          // thickness of the deposited layers
  double first_layer_thickness;       // thickness of the first layer in the bilayer
  double second_layer_thickness;      // thickness of the second layer in the bilayer
  double bilayer_thickness;           // thickness of the bilayer
  double bilayer_ratio;               // ratio between bilayer thickness and  first  layer thickness
  int bilayer_number;                 // number of bilayers
  int bilayer_counter;                // number of current depositiong bilayer 
    
  double inp;
  int flaq;
  int flaq_out = 1;
  int add_cur = 1;
  FILE *ff, *out, *adF, *ff1, *ff2;
   
  int i, k;
  int step = 1;
  int final_step;
  int slope_flaq;
  
  double t;
  double Total;
  double rand_max = 10.0 * (double)RAND_MAX; 
  double half_amplitude = (double)RAND_MAX / (2.0 * rand_max);
  double p1Left, p1Right;
  
  
  //initialisation of random number generator
  srand(time(NULL));  
  
  //out = fopen("output","w+");
  
  /* setting program parameters */
  // setting equation parameters
  //nu1 = nu3 = 0.0;
  //nu2 = 0.1;
  //nu4 =10.0;
  //lambda = 0.5; 
  //setting deposition flux parameters
  linear_flaq = 2;
  nonlinearity_flaq = 7;
  materials_flaq = 3;
  g = 0.1;
  shadow_flaq = 2;
  alpha =  75.0 * pi / 180.0;
  
  //Zr-Al parameters
  if(materials_flaq == 1)
  {
    N_Zr_Al = 200;
    profile_mult_Zr_Al = 200.0;
    N_one_period_Zr_Al = 100;
    N_ggp_Zr_Al = 100;
    N = N_Zr_Al;
    N_ggp = N_ggp_Zr_Al;
    N_one_period = N_one_period_Zr_Al;
    profile_mult = profile_mult_Zr_Al;
    /* Zr-Al multilyares */
    bilayer_flaq = 1;
    bilayer_ratio = 0.4;
    //for the mirrors
    //bilayer_thickness = 10.0;
    //for the sawtooth profile
    bilayer_thickness = 10.45;
    //thikness for the test
    //bilayer_thickness = 1.0;
    first_layer_thickness = bilayer_thickness * bilayer_ratio;
    second_layer_thickness = bilayer_thickness - first_layer_thickness;
    bilayer_number = 20;
    bilayer_counter = 1;
    double_layer_counter = 1.0;
    bilyaer_output = 1.0;
    /*
    blaze_angle = 6.0 * pi / 180;
    second_angle  = 64.0 * pi / 180;
    */
  }
  
  //Mo-Si parameters
  if(materials_flaq == 2)  
  {
    N_Mo_Si = 1360;
    N_one_period_Mo_Si = 136;
    profile_mult_Mo_Si = 1360.0; //272.0; //1360.0;
    N_ggp_Mo_Si = 136;
    N = N_Mo_Si; 
    N_ggp = N_ggp_Mo_Si;
    N_one_period = N_one_period_Mo_Si;
    profile_mult = profile_mult_Mo_Si;
    /* Mo-Si multilayers */
    bilayer_flaq = 1;
    bilayer_ratio = 0.45;
    bilayer_thickness = 6.9741;
    first_layer_thickness = bilayer_thickness * bilayer_ratio;
    second_layer_thickness = bilayer_thickness - first_layer_thickness;
    bilayer_number = 30;
    bilayer_counter = 1; 
    double_layer_counter = 1.0;
    bilyaer_output = 1.0;
  } 
  
  if(materials_flaq == 3)  
  {
    N_W_B4C = 400;
    N_one_period_W_B4C = 200;
    sawtooth_period = N_one_period_W_B4C;
    profile_mult_W_B4C = 400.0; //272.0; //1360.0;
	N = N_W_B4C; 
    N_ggp_W_B4C = N_one_period_W_B4C;
    N_ggp = N_ggp_W_B4C;
    N_one_period = N_one_period_W_B4C;
    profile_mult = profile_mult_W_B4C;
    bilayer_flaq = 1;
    bilayer_ratio = 0.5; //0.44; //0.5; //0.44;Mo/Be
    bilayer_thickness = 3.06; //5.88; //5.88; Mo/Be
    first_layer_thickness = bilayer_thickness * bilayer_ratio;
    second_layer_thickness = bilayer_thickness - first_layer_thickness;
    bilayer_number = 40;
    bilayer_counter = 1; 
    double_layer_counter = 1.0;
    bilyaer_output = 1.0;
    blaze_angle = 6;
    second_angle = 45;
  }  
  
  // setting time parameters
  t = 0.0;
  step = 1;
  Total = 100.0; 
  
  //setting stepOut values
  stepOut_first_layer = (double)(int)(first_layer_thickness / (g * dt));  
  stepOut_second_layer = (double)(int)(bilayer_thickness / (g * dt));
  stepOut = 100.0;  // setting number of steps before output, in t(sec) = stepOut * dt
  //printf("%f",stepOut_first_layer);
  //getch();
  
  
  //setting deposition flux inhomogenity 
  slope_flaq = 1;
  inhomogenity_flaq = 0;
  inhomogenity_time = 1000000;    // in t(sec.) = inhomogenity_time * dt
  //inhomogenity_length = 1;      // in l(nm.) = inhomogenity_length * dx
  /* end of setting program parameters */
  
  // auxiliary flaqs
  flaq = 1;        // calculation flaq
  flaq_out = 1;    // output flaq
  stepOutflaq = 2; // output conditions flaq
  
  //setting number of steps
  //testing program "setting of the final step"
  //final_step = 10000001; 
  final_step = bilayer_number * bilayer_thickness / (g * dt) + 5;
  total_psd_thicknes = bilayer_number * bilayer_thickness + 10 * g * dt;
  deposited_thicknes = 0.0;
  //deposition of buffer;
    
  
  //printf("%d",final_step);
  //getch();
  
  if(flaq == 2)
  {
    Init(blaze_angle, second_angle);
    //WhiteNoiseDeposition();
    //ProfileDeviation();
    //AutoCorrelationFunction();
    PSDCalculationFFT();
  }
  
  /* calculation of the profile evolution */
  if(flaq == 1)
  {
    /* initialization of the system */
    Init(blaze_angle, second_angle);
    adF = fopen("Check.dat", "w+");
    
    //Init();
    /*
	for(i = 0; i < N; i++)
    {
      hOld1[i] = hOld[i] / H0;
	}
	*/
    //buffer deposition
    /*
    buffer_flaq = 0;
    buffer_height = 10.0;
    buffer_stepout = (int)(buffer_height / (g * dt));
    if(buffer_flaq == 1)
    {
      while(step < buffer_stepout)
      {
         nu2 = 0.2; // * log(double(bilayer_counter)) / log(2.71);
         nu4 = 1.5;
         Solution(nu1,nu2,nu3,nu4,lambda);
         step = step + 1;
         //printf("%f\n",step * g * dt);
      }
      DegreeR((double)step,8.0);
      OutStep();
    }  
    */
    
    
    step = 1;
    while(deposited_thicknes < total_psd_thicknes) //while (t < Total) //while(step < 5)
    {
      //printf("hello");
      //printf("%d",step-1);
      //getch();
      
      //deposited_thicknes = g * dt * step * cos(blaze_angle);
      //printf("%f\n",deposited_thicknes);
      //getch();
      
      /* deposition of bilaters */
      if(bilayer_flaq == 1)
      {        
        
        //equation parameters for the first layer
        if(deposited_thicknes <= bilayer_thickness * (bilayer_counter-1) + first_layer_thickness + 2.0 * g * dt && deposited_thicknes > bilayer_thickness * (bilayer_counter-1) - 2.0 * g * dt)
        {
          /*
          if(bilayer_counter <= 7)
          {
             nu2 = 0.05;
             nu4 = 0.5; 
          }
          if(bilayer_counter > 7)
          {
             nu2 = 0.2;
             nu4 = 2.5; 
          }
          */
          nu2 = 0.15; // + 0.001 *bilayer_counter; //+ 0.0055 * bilayer_counter;1
          nu4 = 35.0; //2 + 0.5 * bilayer_counter;//0.5;
          lambda = 0.25; 
          
          /*
		  if(slope_flaq == 1)
		  {
            p1Left = 
            p1Right = 
		    
			if(p1Left > 0 && p1Right > 0)
			{
				
			}
			
			if(p1Left < 0 && p1Right < 0)
			if((p1Left < 0 && p1Right > 0) || (p1Left > 0 && p1Right < 0 ))
		    {
		  	
		    }
		  
		    if(p1Left < 0 && p1Right < 0) 
		    {
		  	
		    }
	      }*/
          //getch();
        }
        
        
        // equation parameters for the second layer
        if(deposited_thicknes > bilayer_thickness * (bilayer_counter-1) + first_layer_thickness - 2.0 * g * dt && deposited_thicknes <= bilayer_thickness * bilayer_counter + 2.0 * g * dt)
        {
          /*
          if(bilayer_counter <= 7)
          {
             nu2 = 0.05;
             nu4 = 1.0; 
          }
          if(bilayer_counter > 7)
          {
             nu2 = 0.2;
             nu4 = 3.0; 
          }
          */
          nu2 = 0.15; // + 0.001 *bilayer_counter; 
          //nu2 = 0.1; // + 0.0055 * bilayer_counter;
          nu4 = 25.0;// + 0.5 * bilayer_counter;
          //nu4 = 1.5;          
          lambda = 0.25; 
          //printf("%d  %f\n", bilayer_counter, nu4);
          //getch();
        }        
        
        //checking condition of next bilayer deposition
        if(deposited_thicknes > bilayer_thickness * bilayer_counter - 2.0 * g * dt && deposited_thicknes <= bilayer_thickness * bilayer_counter + 2.0 * g * dt)
        {
          bilayer_counter = bilayer_counter + 1;
          double_layer_counter = double_layer_counter + 1.0;
        }  
      }//end of setting parameters for the bilayer depostion 
      
     //deposition of layers (no bilyaer deposition)
     
      if(inhomogenity_flaq == 1)
      {
        noise_length = 1;
        max_noise = 0.15;
        
		/*if(step == 2)
        {
		  //printf("")
		  for(i = 0; i < N; i++)
		    fprintf(ff1, "%d  %1.10f\n", i, hNew[i]);   
	    }
        if(step == 10000)
        {
		  for(i = 0; i < N; i++)
		    fprintf(ff2, "%d  %1.10f\n", i, hNew[i]);   
	    }
		*/	    
		if((step-1) % inhomogenity_time == 0)
        {
          //printf("step = %d\n", step);
          //getch();
		  WhiteNoiseDeposition(step);
		  printf("step = %d  dep[0] = %f\n",step,  dep[0]);
          //printf("Noise %f  %f  %d\n", dep[5], dep[115], add_cur);
          //getch();
          add_cur++;
          if(add_cur == -1)
          {
            for(i = 0; i < N; i++)
              fprintf(adF,"%d  %f\n",i, hNew[i]);
          }
           
        }
      }
      
      /* deposition without noise from the source */
      if(inhomogenity_flaq == 0)
      {
         for(i = 0; i < N; i++)
           dep[i] = g;                   
      }
      
      deposited_thicknes = g * dt * step * cos(blaze_angle * pi / 180);
      /* evolution of the surface profile */  
      //average_dep = 0.0;
      //for(i = 0; i < N; i++)
      //{
       // hOld1[i] = hOld[i] / H0;
      //}
	  Solution(step); 
      //average_dep = average_dep / (double)N; 
      //total_dep = total_dep + average_dep;
      
      /* setting equation parameters for the different layers of the bilayer */
      //growth with normal deposition flux
      //deposited_thicknes = g * dt * step;
      
      //growth with tilted deposition flux
      //deposited_thicknes = total_dep;
      
      /*
      printf("average_dep=%1.15f\n",average_dep);
      getch();
      */
      //printf("p1=%f\n",p1);
              
      
      //printf("%f\n",c_k);
      //while(!getch())
           
      /* checking condition of output */
      // output through equal periods of time
      if(stepOutflaq == 1)
      {  
        //printf("%f   %f\n",c_k,stepOut);
        //printf("%f  calc...\n",c_k);
        if( c_k < (dj * stepOut + 1E-10) && c_k > ( dj * stepOut - 1E-10) && flaq_out == 1) //if(step == dj)
        {
          DegreeR(c_k,8.0);
          OutStep();
          dj = dj + 1.0;   
          //getch();
        }
        
      }  
      
      //printf("output=%f\n",dj * stepOut_first_layer + (dj - 1.0) * stepOut_second_layer);
      //getch();
      // output of the profiles of first and second bilayer     
      if(stepOutflaq == 2)
      {  
        //printf("%f\n",double_layer_counter);
        c_temp1 = stepOut_first_layer + (dj - 1.0) * stepOut_second_layer;
        c_temp2 = dj * stepOut_second_layer;
        //printf("bilayer deposition %f %f %f\n",dj, c_k, c_temp1);    
        /*
        if(step == (int)(stepOut_first_layer+1))
        {
          printf("%f  %f\n", dj,c_temp1);
          DegreeR(c_k,8.0);
          OutStep();
          //getch();
        }
        */
        
        //output of the first layer
        //if( c_k < c_temp1 + 0.0000000001 && c_k > c_temp1 - 0.0000000001)
        if(deposited_thicknes < bilayer_thickness * (bilayer_counter-1) + first_layer_thickness + 2 * g * dt && deposited_thicknes > bilayer_thickness * (bilayer_counter-1) + first_layer_thickness - 2 * g * dt)
        {
          //printf("first_layer %f\n",c_temp1);
          //getch();
          bilyaer_output = 2.0 * (double)bilayer_number - double_layer_counter * 2.0 + 2.0;      
          DegreeR(bilyaer_output,8.0);
          //printf("check1 output\n");
          OutStep();
          //printf("check2 output\n");
          //SawParams(blaze_angle, second_angle);
          //printf("output first layer\n");
          //getch();
          //dj = dj;
        }
       
                
        //output of the second layer
        //if( c_k < dj * stepOut_second_layer + 0.0000000001 && c_k > dj * stepOut_second_layer - 0.0000000001) 
        if(deposited_thicknes > bilayer_thickness * bilayer_counter - 2 * g * dt && deposited_thicknes <= bilayer_thickness * bilayer_counter + 2 * g * dt)
        {
          //printf("second layer\n");
          //getch();
          bilyaer_output = 2.0 * (double)bilayer_number - double_layer_counter * 2.0 + 1.0;
          DegreeR(bilyaer_output,8.0);
          OutStep();
          dj = dj + 1.0;   
          //SawParams(blaze_angle, second_angle);
        }
      }  
     
     //increasing of the time
      t = t + dt;
      // increase of the output counter
      c_k = c_k + 1.0;
      //printf("time=%3.10f\n",t); 
      
      printf("time=%3.10f  nu2=%f  nu4=%f lambda=%f h(time)=%1.10f nm \n",t, nu2, nu4, lambda, deposited_thicknes);
	  //deposited_thicknes=%1.10f",t,deposited_thicknes,);
      //printf("time=%3.10f  %f  %f\n",t, nu2, nu4);
      //printf("step=%d  %f\n",step, c_k);
      //printf("step=%d\n",step);
      //getch();
      step++;
      /*
	  if (step == 100000)
      {
      	for(i = 0; i < N; i++)
      	  fprintf(adF,"%d  %f  %f\n", i, hNew[i], hOld[i]);
	  }
	  */
      
    }
  
     /* calculation of surface profile integral data */
     //final output of the profile
     FinalOutput();
     //calculation of the profile PSD function
     //PSDCalculationFFT();
     
     /*
     fW = fopen("white noise","w+");
     for(i = 0; i < N; i++)
       fprintf(fW,"%d  %f\n",i,length_array[i]);
     fclose(fW);  
     */
     //calculation of model profile deviation from the experimnetal top surface profile
     //ProfileDeviation();
     
  }
  fclose(adF);
  
  
}
