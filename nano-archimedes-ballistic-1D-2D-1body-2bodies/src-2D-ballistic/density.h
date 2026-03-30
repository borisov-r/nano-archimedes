/* ===============================================================

   density.h -- This file belongs to the nano-archimedes.

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

// calculate the density of particles in function of the position

void density(void){
 register int n;
 register int i,j;
 int kx,ky;

 // reset the density array
 for(i=0;i<=NX;i++) for(j=0;j<=NY;j++) RHO[i][j]=0.;
 for(i=0;i<2*NKX-1;i++) for(j=0;j<2*NKY-1;j++) PSI[i][j]=0.;

 // calculation of the real space density
 printf("calculation of the real space density.\n");
 for(n=0;n<INUM;n++){
  i=(int)(PX[n]/DX)+1;
  j=(int)(PY[n]/DY)+1;
  kx=KX[n];
  ky=KY[n];
  if((i>0 && i<=NX) && 
     (j>0 && j<=NY) &&
     (-NKX<kx && kx<NKX) &&
     (-NKY<ky && ky<NKY)) RHO[i][j]+=W[n];
 }

 // calculation of the phase space density
 printf("calculation of the phase space density.\n");
 for(n=0;n<INUM;n++){
  i=(int)(PX[n]/DX)+1;
  j=(int)(PY[n]/DY)+1;
  kx=KX[n];
  ky=KY[n];
  if((i>0 && i<=NX) &&
     (j>0 && j<=NY) &&
     (-NKX<kx && kx<NKX) &&
     (-NKY<ky && ky<NKY)) PSI[NKX-1+kx][NKY-1+ky]+=W[n];
 }
}
