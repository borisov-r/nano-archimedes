/* ===============================================================

   config.h -- This file belongs to the nano-archimedes.

   nano-archimedes is a quantum mechanical simulator
   which implements the signed particle formulation of quantum mechanics.

   Copyright (C) 2004-2015 Jean Michel D. Sellier
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


// initial configuration of electrons in device

double _Complex phi(int ind,double x){
 return cexp(-_Complex_I*0.5*pow((x-X0_WAVE_PACKET[ind])/SIGMA_WAVE_PACKET[ind],2.))*cexp(_Complex_I*K0_WAVE_PACKET[ind]*x);
}

void devconf(void){
 register int i,j,m,n;

// double entangled_sigma_x=0.5*SIGMA_WAVE_PACKET[0];
// double entangled_sigma_k=DKX;

 printf("\ndefining the initial distribution function\n");
 // definition of a wave packet distribution function
 for(i=1;i<=NX;i++)
  for(j=1;j<=NX;j++){
   for(m=-NKX+1;m<NKX;m++)
    for(n=-NKX+1;n<NKX;n++){
     // independent particles
     DIST[i][j][m+NKX-1][n+NKX-1]=(int)(INUM*exp(-pow(((i-0.5)*DX-X0_WAVE_PACKET[0])/SIGMA_WAVE_PACKET[0],2.))
                                            *exp(-pow(((j-0.5)*DX-X0_WAVE_PACKET[1])/SIGMA_WAVE_PACKET[1],2.))
                                            *exp(-pow(((m*DKX)-K0_WAVE_PACKET[0])*SIGMA_WAVE_PACKET[0],2.))
                                            *exp(-pow(((n*DKX)-K0_WAVE_PACKET[1])*SIGMA_WAVE_PACKET[1],2.)));
    }
  }
 annihilation();

 printf("Initial Number of Electron Super-particles = %d\n", INUM);

}

// ==============================================================
