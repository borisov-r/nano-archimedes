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

void devconf(void){

 // definition of a wave packet distribution function
 printf("\n\ndevconf() - calculation of the initial Wigner function..\n");
 register int i,j,l,m;
 for(i=1;i<=NX;i++)
  for(j=1;j<=NY;j++)
   for(l=-NKX+1;l<NKX;l++)
    for(m=-NKY+1;m<NKY;m++){
     DIST[i][j][l+NKX-1][m+NKY-1]=(int)(INUM*exp(-pow(((i-0.5)*DX-X0_WAVE_PACKET)/SIGMA_WAVE_PACKET,2.)
                                                 -pow(((j-0.5)*DY-Y0_WAVE_PACKET)/SIGMA_WAVE_PACKET,2.))
                                 *exp(-pow(((l*DKX)-KX0_WAVE_PACKET)*SIGMA_WAVE_PACKET,2.)
                                      -pow(((m*DKY)-KY0_WAVE_PACKET)*SIGMA_WAVE_PACKET,2.)));
   }
 printf("devconf() - calculation of the initial position of particles in phase space..\n");
 annihilation();

 printf("devconf() - Initial Number of Electron Super-particles = %d\n\n", INUM);

}

// ==============================================================
