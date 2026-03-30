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
#include<complex.h>
#include<omp.h>
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
#define MN3 3
#define NXM 256
#define NKXM 512
#define NPMAX 32000000 // maximum number of super-particles

// definition of constants
const double Q=1.60217733e-19;    // Electron charge in absolute value (Coulomb)
const double HBAR=1.05457266e-34; // Reduced Planck constant (Joule*sec)
const double ME=9.1093897e-31;     // free electron mass (Kg)
const double MP=1.6726217e-27;     // free proton mass (Kg)
const double PI=3.141592654;      // Pi number
const double EPS0=8.854187817e-12; // Permittivity of free space (F/m)

// All integers here...
int NX;
int NKX;
int ITMAX,FINAL;
int INUM;
int *K[2];
int *W;
int ****DIST;
int ISEED;
int *UPDATED;
int FTOT[NXM+1][NKXM+1];

// All doubles here...
double PHI[NXM+1];
double PHI0[NXM+1];
double RHO0[NXM+1];
double E[NXM+1];
double RHO[2][NXM+1];
double E[NXM+1];
double LC;
double DX;
double DKX;
double TIME=0.;
double LX;
double DT;
double BKTQ,QH;
double *P[2];
double ****VW;
double **GAMMA;
double *PTIME;
double SIGMA_WAVE_PACKET[2];
double X0_WAVE_PACKET[2];
double K0_WAVE_PACKET[2];

// All structures here...
time_t binarytime;
struct tm *nowtm;
char s[100];

// All files here...
FILE *fp;

#include "random.h"
#include "annihilation.h"
#include "density.h"
#include "hartree.h"
#include "distribution.h"
#include "save.h"
#include "config.h"
#include "wigner.h"
#include "gamma.h"
#include "emc.h"

int main(void)
{
 register int i,j,k;
 double a0;

 // Bohr radius
 a0=4.*PI*EPS0*HBAR*HBAR/(ME*Q*Q);


 INUM=1024;       // define initial number of particles
 LX=8.*a0;        // total lenght of device
 LC=LX/3.;        // coherence lenght
 NX=120;          // number of cells in x-direction
 DT=0.01e-18;     // time step
 ITMAX=1200;      // total number of time steps

 // Electron spread over the domain at rest
 SIGMA_WAVE_PACKET[0]=a0/3.;    // wave packet width
 X0_WAVE_PACKET[0]=0.5*LX; // wave packet initial position

 // Proton in the middle at rest
 SIGMA_WAVE_PACKET[1]=a0/12.; // wave packet width
 X0_WAVE_PACKET[1]=0.5*LX;   // wave packet initial position

 // define random numbers generator seed
 ISEED=38467;

 // spatial cell lenght
 DX=LX/NX;

 // automatic calculation of NKX
 NKX=(int)(0.5*LC/DX);

 // pseudo-wave vector lenght
 DKX=PI/LC;

 K0_WAVE_PACKET[0]=0.*DKX;   // wave packet initial wave vector
 K0_WAVE_PACKET[1]=0.*DKX;   // wave packet initial wave vector

 printf("\n\
===========================\n\
 K = [0.:%g:%g] 1/m\n\
===========================\n",DKX,NKX*DKX);

 // print the initial number of particles
 printf("\nMAXIMUM NUMBER OF PARTICLES ALLOWED = %d\n\n",NPMAX);

 // memory allocation
#include "alloc.h"
printf("\n");

 // constant electrostatic potential
 for(i=1;i<=NX;i++) PHI[i]=PHI0[i]=0.;

 // set gamma function to zero
 for(i=0;i<=NX;i++)
  for(j=0;j<=NX;j++)
   GAMMA[i][j]=0.;

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
  fflush(stdout);
  density();
  if((i%50)==0){
   hartree();
   wigner();
   calculate_gamma();
  }
  EMC();
  distribution();
  if((i%10)==0){
   printf("Annihilation of particles\n");
   if(i!=ITMAX) annihilation();
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
#include "free.h"

 return(EXIT_SUCCESS); // Successfull exit
}

/* main.c ends here */

// ***********************************************************************
// ************************************************************************
