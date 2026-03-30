/* ===============================================================

   distribution.h -- This file belongs to the nano-archimedes.

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

// calculates the distribution function according to the quantum signs

void distribution(void){
 register int n;
 int i,j,k;

 printf("\nCalculation of distribution function\n");
 printf("Number of particles = %d\n",INUM);

 for(i=0;i<=NX;i++) for(k=0;k<2*NKX-1;k++) DIST[i][k]=0;

 // cloud in cell algorithm
 for(n=0;n<INUM;n++){
  i=(int)(P[n]/DX)+1;
  k=K[n];
  if((0<i) && (i<=NX) && (-NKX<k) && (k<NKX)) DIST[i][k+NKX-1]+=W[n];
 }

 // stores the normalized quasi-distribution function
 double norm;
 norm=0.;
 for(i=1;i<=NX;i++)
  for(j=-NKX+1;j<NKX;j++)
   FW1[i][j+NKX-1]=(double)(DIST[i][j+NKX-1]);
 for(i=1;i<=NX;i++) for(j=-NKX+1;j<NKX;j++) norm+=FW1[i][j+NKX-1];
 norm*=DX*DKX;
  for(i=1;i<=NX;i++)
  for(j=-NKX+1;j<NKX;j++)
   FW1[i][j+NKX-1]/=norm;

 printf("end of distribution function calculation\n");
}
