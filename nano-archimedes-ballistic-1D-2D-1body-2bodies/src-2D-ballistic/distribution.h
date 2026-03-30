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
 int i,j,kx,ky;

 printf("\nCalculation of distribution function\n");
 printf("Number of particles = %d\n",INUM);

 // set the distribution function to zero
 for(i=0;i<=NX;i++){
  for(j=0;j<=NY;j++){
   for(kx=0;kx<2*NKX+1;kx++){
    memset(DIST[i][j][kx],0,sizeof(DIST[i][j][kx][0])*(2*NKY+1));
   }
  }
 }

/*
 for(i=0;i<=NX;i++) for(j=0;j<=NY;j++)
  for(kx=0;kx<2*NKX-1;kx++) for(ky=0;ky<2*NKY-1;ky++)
   DIST[i][j][kx][ky]=0;*/

 for(n=0;n<INUM;n++){
  i=(int)(PX[n]/DX)+1;
  j=(int)(PY[n]/DY)+1;
  kx=KX[n];
  ky=KY[n];
  if((0<i) && (i<=NX) && (-NKX<kx) && (kx<NKX) &&
     (0<j) && (j<=NY) && (-NKY<ky) && (ky<NKY)) DIST[i][j][kx+NKX-1][ky+NKY-1]+=W[n];
 }

 printf("end of distribution function calculation\n");
}
