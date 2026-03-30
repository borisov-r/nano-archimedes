/* ===============================================================

   config.h -- This file belongs to the nano-archimedes.

   nano-archimedes is a quantum mechanical simulator
   which implements the single-body Wigner Monte Carlo method.
   This code simulates a Gaussian wave-packet going towards
   a potential barrier.

   Copyright (C) 2004-2014 Jean Michel D. Sellier
   <jeanmichel.sellier@gmail.com>
   <jeanmichel.sellier@parallel.bas.bg>

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

// === initial configuration of the wave packet ===
// === Gaussian wave packet ===

void devconf(void){
 register int i,j;
 double epp,norm;
 double d_max;

 d_max=0.;

 // definition of the initial conditions
 for(i=1;i<=NX;i++)
  for(j=-NKX+1;j<NKX;j++){
   FW1[i][j+NKX-1]=exp(-pow(((i-0.5)*DX-X0_WAVE_PACKET)/SIGMA_WAVE_PACKET,2.))
                  *exp(-pow(((j*DKX)-K0_WAVE_PACKET)*SIGMA_WAVE_PACKET,2.));
  }

 // normalization of the initial conditions
 norm=0.;
 for(i=1;i<=NX;i++) for(j=-NKX+1;j<NKX;j++) norm+=FW1[i][j+NKX-1];
 norm*=DX*DKX;
 for(i=1;i<=NX;i++) for(j=-NKX+1;j<NKX;j++) FW1[i][j+NKX-1]/=norm;

 // calculates the EPP variable for the cloud in cell algorithm
 for(i=1;i<=NX;i++) for(j=0;j<2*NKX-1;j++) if(d_max<fabs(FW1[i][j])) d_max=fabs(FW1[i][j]);
 epp=d_max/INUM;

 // calculate initial distribution function
 printf("config() - calculating initial distribution\n");
 INUM=0;
 for(i=1;i<=NX;i++){
  for(j=0;j<2*NKX-1;j++){
   register int n;
   int local_number_of_particles;
   local_number_of_particles=(int)(fabs(FW1[i][j])/epp+0.5);
   // creates the new local particles in the (i,k)-th phase-space cell
   // the particles are uniformely distributed in space
   for(n=1;n<=local_number_of_particles;n++){
    int m;
    m=INUM+n-1;
    if(rnd()>0.5) P[m]=(i-0.5+0.5*rnd())*DX;
    else P[m]=(i-0.5-0.5*rnd())*DX;
    K[m]=j-NKX+1;
    if(FW1[i][j]>0) W[m]=+1;
    else W[m]=-1;
   }
   INUM+=local_number_of_particles;
  }
 }

 distribution();

 printf("Initial Number of Electron Super-particles = %d\n", INUM);
}

// ==============================================================
