/* ===============================================================

   nano-archimedes.c -- This file belongs to the nano-archimedes.

   nano-archimedes is a quantum mechanical simulator
   which implements the single-body Wigner Monte Carlo method.
   This code simulates a Gaussian wave-packet going towards
   a potential barrier.

   nano-archimedes is a GNU package.

   ===> You can extend easily this code to 2D and 3D space by following
        the details presented in this paper:

        J.M.~Sellier, M.~Nedjalkov, I.~Dimov, S.~Selberherr,
        A benchmark study of the Wigner Monte-Carlo method, Monte Carlo Methods and Applications,
        De Gruyter, DOI: 10.1515/mcma-2013-0018, (2014).

   ===> You can utilize this code for density functional theory (DFT) calculations
        by following the details reported in the this paper:

        J.M.~Sellier, I.~Dimov,
        A Wigner Monte Carlo Approach to Density Functional Theory,
       Journal of Computational Physics 270, pp. 265-277, (2014).

   ===> You can generalize this code to the time-dependent quantum many-body
        problem (ab-initio) by following the instructions reported in this
        paper:

        J.M.~Sellier, I.~Dimov,
        The Many-body Wigner Monte Carlo Method for Time-dependent Ab-initio Quantum Simulations,
        Journal of Computational Physics 273, pp. 589-597 (2014).

   =============================================================================
   If you use nano-archimedes in your research please consider citing one of the
   papers above.
   =============================================================================

   Copyright (C) 2004-2015 Jean Michel D. Sellier
   <jeanmichel.sellier@gmail.com>

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 3, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.

   ===============================================================
*/

// ===============================================================
// File Name     : nano-archimedes.c
// Author        : Jean Michel Sellier
// Revision      : 27 Jan. 2015, Beijing, China
// ===============================================================

// ===============================================================
// TO COMPILE TYPE IN A TERMINAL:
//
// > gcc nano-archimedes.c -Wall -lm -Ofast -o nano-archimedes
//
// ===============================================================

// ===============================================================
// TO RUN THE CODE IN A TERMINAL:
//
// > ./nano-archimedes
//
// ===============================================================

#include<stdio.h>
#include<unistd.h>
#include<stdlib.h>
#include<math.h>
#include<memory.h>
#include<time.h>
#ifdef	HAVE_STRING_H
#include<string.h>
#else
#include<strings.h>
#endif

// Preprocessor Definitions
#define NO 0
#define YES 1
#define NXM 1024
#define NKXM 512
#define NPMAX 4000000 // maximum number of super-particles

// definition of constants
const double Q=1.60217733e-19;    // Electron charge in absolute value (Coulomb)
const double HBAR=1.05457266e-34; // Reduced Planck constant (Joule*sec)
const double M=9.1093897e-31;     // free electron mass (Kg)
const double PI=3.141592654;      // Pi number
const double MSTAR=0.067;         // GaAs effective mass

// All integers here...
int NX;
int NKX;
int ITMAX,FINAL;
int INUM;
int *K;
int *W;
int DIST[NXM+1][2*NKXM+1];
int ISEED;
int *UPDATED;
int ANNIHILATION_FREQUENCY;

// All doubles here...
double FW1[NXM+1][2*NKXM+1];
double DENSX[NXM+1];
double DENSK[2*NKXM+1];
double PHI[NXM+1];
double LC;
double DX;
double DKX;
double TIME=0.;
double LX;
double DT;
double BKTQ,QH;
double *P;
double VW[NXM+1][2*NKXM+1];
double GAMMA[NXM+1];
double *PTIME;
double SIGMA_WAVE_PACKET;
double X0_WAVE_PACKET;
double K0_WAVE_PACKET;
double BARRIER_POTENTIAL;
double BARRIER_POSITION;
double BARRIER_WIDTH;

// All structures here...
time_t binarytime;
struct tm *nowtm;
char s[100];

// All files here...
FILE *fp;

#include "random.h"
#include "annihilation.h"
#include "distribution.h"
#include "density.h"
#include "save.h"
#include "config.h"
#include "kernel.h"
#include "gamma.h"
#include "wmc.h"

int main(void)
{
 int i;

 // The following parameters completely define
 // the simulation problem. For their meaning, see the comments below.

 INUM=20000;  // maximum number of particles in a phase-space cell for the initial distribution
 LX=200.e-9;  // total lenght of spatial domain
 LC=50.e-9;   // coherence lenght
 NX=200;      // number of cells in x-direction
 DT=0.01e-15; // time step
 ITMAX=4000;  // total number of time steps
 ANNIHILATION_FREQUENCY=100; // annihilation occurs every 100 time steps

 SIGMA_WAVE_PACKET=10.e-9;    // wave packet dispersion
 X0_WAVE_PACKET=LX/2-31.5e-9; // wave packet initial position

 BARRIER_POTENTIAL=-0.3;      // value of the potential barrier
 BARRIER_POSITION=0.5*LX;     // barrier center position
 BARRIER_WIDTH=6.e-9;         // barrier width

 // define random numbers generator seed
 ISEED=38467;

 // spatial cell lenght
 DX=LX/NX;

 // automatic calculation of NKX
 NKX=(int)(0.5*LC/DX);

 // pseudo-wave vector lenght
 DKX=PI/LC;
 K0_WAVE_PACKET=6.*DKX;       // wave packet initial wave vector

 // print the initial number of particles
 printf("\nMAXIMUM NUMBER OF PARTICLES ALLOWED = %d\n\n",NPMAX);

 // memory allocation
 K=malloc((NPMAX+1)*sizeof(*K));
 if(K==NULL){
  printf("Not enough memory to allocate\nint K[NPMAX+1]\n");
  exit(0);
 }
 W=malloc((NPMAX+1)*sizeof(*W));
 if(W==NULL){
  printf("Not enough memory to allocate\nint W[NPMAX+1]\n");
  exit(0);
 }
 UPDATED=malloc((NPMAX+1)*sizeof(*UPDATED));
 if(UPDATED==NULL){
  printf("Not enough memory to allocate\nint UPDATED[NPMAX+1]\n");
  exit(0);
 }
 P=malloc((NPMAX+1)*sizeof(*P));
 if(P==NULL){
  printf("Not enough memory to allocate\ndouble P[NPMAX+1]\n");
  exit(0);
 }
 PTIME=malloc((NPMAX+1)*sizeof(*PTIME));
 if(PTIME==NULL){
  printf("Not enough memory to allocate\ndouble PTIME[NPMAX+1]\n");
  exit(0);
 }

 // constant null potential
 for(i=1;i<=NX;i++) PHI[i]=0.;
 // defines the potential barrier
 for(i=1;i<=NX;i++){
  double pos;
  pos=(i-0.5)*DX;
  if(pos>=(BARRIER_POSITION-0.5*BARRIER_WIDTH) && pos<=(BARRIER_POSITION+0.5*BARRIER_WIDTH)) PHI[i]+=BARRIER_POTENTIAL;
 }

 // set gamma function to zero
 for(i=0;i<=NX+1;i++) GAMMA[i]=0.;

 // get initial time
 binarytime=time(NULL);
 nowtm=localtime(&binarytime);
 sprintf(s,"simulation started   : %s",asctime(nowtm));

 // initializations
 devconf();  // device configuration - distributes particles in the device
 density();
 save(0);    // save the initial configuration

 printf("\n");
 // updates the solution
 for(i=1;i<=ITMAX;i++){
  TIME+=DT;
  printf("%d of %d -- TIME=%g\n\n",i,ITMAX,TIME);
  if(i==1){
   printf("Calculating Wigner potential\n");
   kernel();
   printf("Calculating Gamma function\n");
   calculate_gamma();
  }
  printf("Evolving Wigner\n");
  WMC();
  printf("Calculating distribution function\n");
  distribution();
  printf("Calculating density in x- and k-space\n");
  density();
  if(((i%ANNIHILATION_FREQUENCY)==0) && (i!=ITMAX)){
   printf("Annihilation of particles\n");
   annihilation();
  }
  save(i); // save output every time step
 }
 printf("\n");

 printf("output files saved\n\n");

 // get final time and exit.
 binarytime=time(NULL);
 nowtm=localtime(&binarytime);

 printf("%s",s);
 printf("simulation ended     : %s\n",asctime(nowtm));

 // free the allocated memory
 free(K);
 free(W);
 free(UPDATED);
 free(P);
 free(PTIME);

 return(EXIT_SUCCESS); // Successfull exit
}

/* nano-archimedes.c ends here */

// ***********************************************************************
// ************************************************************************
