#ifndef _H_EQUATION_H_
#define _H_EQUATION_H_

//data arrays
//values of the height h(x,t) on timestep j+1
extern double hNew[10000];
//values of the height h(x,t) on timestep j
extern double hOld[10000];
//value of deposition flux in point x
extern double dep[10000];
//values of the hieght for linear equation
extern double hNew_linear[10000];
extern double hOld_linear[10000];
extern double hOld1[10000];
//extern double ShiftedhOld[10000];

//atomic volume
extern double Volume;
//lattice parameter
extern double a_param;
//number of the spatial points
extern int N;
//number of points in Zr-Al and Mo-Si material systems
extern int N_Zr_Al;
extern int N_Mo_Si;
extern int N_W_B4C;
//spatial length of the system
extern double Length;

/* parameters of the equation */
extern double nu1, nu2, nu3, nu4;
extern double lambda;
//end

/* discretisation of time and spatial */
extern double dt;
extern double dx, dx2, dx3, dx4;

/* equation discretization variables */
extern double p1, p2, p3, p4, p0;
extern double p2n2;  
extern double t1, t2, t3, t4, t0, tn;
/* equation discretization variables for linear equation */
extern double p1_linear, p2_linear, p3_linear, p4_linear, p0_linear, tn_linear;
//vkluchenie ili ne vkluchenie linejnogo rascheta
extern int linear_flaq;
//vybor tipa uravnenija: EW, KPZ, MKPZ, strict equation (EW,KPZ,MKPZ)
extern int nonlinearity_flaq;

//time inhomogenity of the deposition flux
extern int inhomogenity_time;
//length inhomogenity of the deposition flux
extern int inhomogenity_length;
//deposition fluctuation flux flaq
extern int inhomogenity_flaq;
extern int noise_length;
extern double max_noise;


//stationar deposition flux
extern double g;

//value of deposition inhomogenity
extern double inh_g;

extern double c_k;

//deposition flux angle
extern double alpha;
//profile angle
extern double profile_angle;

//derivatives data
extern double deriv_nonlinear[1000][6];
extern double deriv_linear[1000][6];
//linear equation difference scheme data
extern double t_linear[1000][6];
//nonlinear equation difference scheme data
extern double t_nonlinear[1000][6];

//flaq of the flux incident, 
extern int shadow_flaq;   //shadow_flaq = 1, deposition from the left to the right
                          //shadow_flaq = 2, deposition from the right to the left  

extern double pi;

//variables for weighted deposition
extern double pre_dep1, pre_dep2, pre_dep3;

extern double H0;

extern FILE *fW;

//total thikness of the deposited layer with tilted deposition
extern double total_dep;
//average thikness of the deposited layer with tilted deposition
extern double average_dep;

void Solution( int time_stp );

void WhiteNoiseDeposition( int time_step );

#endif
