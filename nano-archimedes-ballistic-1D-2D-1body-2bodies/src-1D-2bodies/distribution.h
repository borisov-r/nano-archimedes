/* ===============================================================

   distribution.h -- This file belongs to the nano-archimedes.

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

// calculates the distribution function according to the quantum weights

void distribution(void){
 register int n;
 int i[2],k[2];

 printf("\nCalculation of distribution function\n");
 printf("Number of particles = %d\n",INUM);
 fflush(stdout);

 for(i[0]=0;i[0]<=NX;i[0]++)
  for(i[1]=0;i[1]<=NX;i[1]++)
   for(k[0]=0;k[0]<2*NKX;k[0]++)
    for(k[1]=0;k[1]<2*NKX;k[1]++)
     DIST[i[0]][i[1]][k[0]][k[1]]=0;

 for(n=0;n<INUM;n++){
  i[0]=(int)(P[0][n]/DX)+1;
  i[1]=(int)(P[1][n]/DX)+1;
  k[0]=K[0][n];
  k[1]=K[1][n];
  if((0<i[0]) && (i[0]<=NX) && (-NKX<k[0]) && (k[0]<NKX) &&
     (0<i[1]) && (i[1]<=NX) && (-NKX<k[1]) && (k[1]<NKX))
   DIST[i[0]][i[1]][k[0]+NKX-1][k[1]+NKX-1]+=W[n];
 }

 printf("end of distribution function calculation\n");
 fflush(stdout);
}
