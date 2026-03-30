/* ===============================================================

   nano-archimedes.c -- This file belongs to the nano-archimedes.

   nano-archimedes is a quantum mechanical simulator
   which implements the signed particle formulation of quantum mechanics.

   nano-archimedes is a GNU package.

   The theory behind this code can be found here:

   1) J.M. Sellier,
      A signed particle formulation of non-relativistic quantum mechanics,
      Journal of Computational Phyisics (2015).

   2) J.M. Sellier, M. Nedjalkov, I. Dimov,
      An introduction to applied quantum mechanics in the Wigner Monte Carlo formalism,
      Physics Reports 577, pp. 1-34, (2015).

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
// Revision      : 10 Nov. 2015, Sofia, Bulgaria
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
#define ON 1
#define OFF 0
#define YES 1
#define NO 0
#define NPMAX 950000000 // maximum number of super-particles

// definition of constants
const double Q=1.60217733e-19;    // Electron charge in absolute value (Coulomb)
const double HBAR=1.05457266e-34; // Reduced Planck constant (Joule*sec)
const double M=9.1093897e-31;     // free electron mass (Kg)
const double PI=3.141592654;      // Pi number
const double MSTAR=0.067;          // Nedjalkov

// All integers here...
int NX;
int NY;
int NKX;
int NKY;
int ITMAX,FINAL;
int INUM;
int *KX;
int *KY;
int *W;
int ****DIST;
int ISEED;
int *UPDATED;

// All doubles here...
double **PHI;
double **RHO;
double **PSI;
double LC;
double DX;
double DY;
double DKX;
double DKY;
double TIME=0.;
double LX;
double LY;
double DT;
double BKTQ,QH;
double *PX;
double *PY;
double ****VW;
double **GAMMA;
double *PTIME;
double SIGMA_WAVE_PACKET;
double X0_WAVE_PACKET;
double Y0_WAVE_PACKET;
double KX0_WAVE_PACKET;
double KY0_WAVE_PACKET; 
double VBARRIER;

// All structures here...
time_t binarytime;
struct tm *nowtm;
char s[100];

// All files here...
FILE *fp;

#include "random.h"
#include "annihilation.h"
#include "distribution.h"
#include "save.h"
#include "config.h"
#include "wigner.h"
#include "gamma.h"
#include "emc.h"
#include "density.h"

int main(void)
{
 int i,j,k;

 INUM=512;    // define initial number of particles
 LX=150.e-9;     // x-direction total lenght of device
 LY=150.e-9;     // y-direction total lenght of device
 LC=100.e-9;      // coherence lenght
 NX=150;          // number of cells in x-direction
 NY=100;          // number of cells in y-direction
 NKX=32;         // number of cells in kx-direction
 NKY=32;         // number of cells in ky-direction
 DT=1.e-15;      // time step
 ITMAX=300;      // total number of time steps

 // against the wall in oblique direction
 SIGMA_WAVE_PACKET=10.e-9; // wave packet width
 X0_WAVE_PACKET=50.e-9;    // wave packet initial position in x-direction
 Y0_WAVE_PACKET=0.4*LY;    // wave packet initial position in y-direction
 KX0_WAVE_PACKET=2.5e8;    // wave packet initial Kx wave vector
 KY0_WAVE_PACKET=2.5e8;    // wave packet initial Ky wave vector

 VBARRIER=-0.05;

 // define random numbers generator seed
 ISEED=38467;

 // spatial cell lenght
 DX=LX/NX;
 DY=LY/NY;

 // pseudo-wave vector lenght
 DKX=PI/LC;
 DKY=PI/LC;

 // print the initial number of particles
 printf("\nMAXIMUM NUMBER OF PARTICLES ALLOWED = %d\n\n",NPMAX);

 // memory allocation
#include "alloc.h"

 // definition of a barrier on the upper half part of the domain
 printf("creating the potential profile..\n");
 // horizontal wall
 for(i=1;i<=NX;i++) for(j=1;j<=NY;j++){
  double x,y;
  x=(i-0.5)*DX;
  y=(j-0.5)*DY;
  PHI[i][j]=0.;
  if(x>=65.e-9) PHI[i][j]=VBARRIER;
 }
 printf("potential profile created.\n\n");

 // set gamma function to zero
 for(i=0;i<=NX;i++) for(j=0;j<=NY;j++) GAMMA[i][j]=0.;

 // get initial time
 binarytime=time(NULL);
 nowtm=localtime(&binarytime);
 sprintf(s,"simulation started   : %s",asctime(nowtm));

 // initializations
 devconf();               // device configuration - distributes particles in the device
 density();
 save(0);

 printf("\n");
 // updates the solution
 for(i=1;i<=ITMAX;i++){
  TIME+=DT;
  printf("%d of %d -- TIME=%g\n\n",i,ITMAX,TIME);
  if(i==1){
   printf("Calculating Wigner potential\n");
   wigner();
   printf("Calculating Gamma function\n");
   calculate_gamma();
  }
  printf("Evolving particles\n");
  EMC();
  printf("Calculating distribution function\n");
  distribution();
  printf("Calculating particle density\n");
  density();
  printf("Annihilation of particles\n");
  if(i!=ITMAX) annihilation();
  // eventually Poisson should go here
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
#include "free.h"

 return(EXIT_SUCCESS); // Successfull exit
}

/* main.c ends here */

// ***********************************************************************
// ************************************************************************
