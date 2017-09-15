#ifndef _H_PSD_H_
#define _H_PSD_H_

//extern double PSDfunction[200][2]; //Zr-Al 
extern double PSDfunction[10000][2]; //Mo-Si

extern double tau0;

//top experimental profile 
//extern double top_profile[100];  //Zr-Al
extern double top_profile[10000];  //Mo-Si

//profile of the one period of the sawtooth substrate
//extern double one_period[100];  //Zr-Al
extern double one_period[10000];  //Mo-Si

//profile of the one period of the sawtooth substrate for ".ggp" output 
//extern double one_period_ggp[101]; //Zr-Al
extern double one_period_ggp[10000]; //Mo-Si

//final profile obtained with the theoretical model
extern double model_profile[10000];

//number of pints per one period of sawtooth profile for ".ggp" output
extern int N_ggp;
extern int N_ggp_Zr_Al;
extern int N_ggp_Mo_Si;
extern int N_ggp_W_B4C;


//residials calculation data
//Zr-Al
//extern double residial[200];       // residial: residial[i] = hNew[i] - top_profile[i]
//Mo-Si
extern double residial[10000];       // residial: residial[i] = hNew[i] - top_profile[i]
extern double residial_sum;        //sum of squared residuals = Summ from 1 to N-1 of residial[i]^2 

//number of pints per one period of sawtooth profile
extern int N_one_period; 
//Zr-Al
extern int N_one_period_Zr_Al;
//Mo-Si
extern int N_one_period_Mo_Si;
extern int N_one_period_W_B4C;
//total deposited thickness for the PSD function
extern double total_psd_thicknes;

//multiplicator for .ggp profile 
extern double profile_mult;
extern double profile_mult_Zr_Al;
extern double profile_mult_Mo_Si;
extern double profile_mult_W_B4C;

//array of frequencies in PSD function
extern int length_array[10000];

void PSDCalculation ( void );

void PSDCalculationFFT ( void );

void OnePeriodProfileData( void );

void MacroProfileData( void );

void ThicknessCalc( void );

void FullPeriodProfileData( void );

void ProfileDeviation ( void );

void AutoCorrelationFunction( void );

void ProfileZshift( void );

#endif
