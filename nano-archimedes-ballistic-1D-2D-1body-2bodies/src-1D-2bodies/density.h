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
 register int m,n;
 register int i,j;
 double sum;

 // reset the density array
 for(i=0;i<=NX;i++){
  RHO[0][i]=0.;
  RHO[1][i]=0.;
 }

 // calculation of the normalized density
 printf("\ncalculation of the particle densities\n");
 for(i=0;i<=NX;i++){
  for(j=0;j<=NX;j++)
  for(m=-NKX+1;m<NKX;m++)
   for(n=-NKX+1;n<NKX;n++)
    RHO[0][i]+=DIST[i][j][m+NKX-1][n+NKX-1];
 }

 for(j=0;j<=NX;j++){
  for(i=0;i<=NX;i++)
  for(m=-NKX+1;m<NKX;m++)
   for(n=-NKX+1;n<NKX;n++)
    RHO[1][j]+=DIST[i][j][m+NKX-1][n+NKX-1];
 }

 // normalization of the density
 int l;
 for(l=0;l<2;l++){
  sum=0.;
  for(i=1;i<=NX;i++) sum+=RHO[l][i];
  sum*=DX;
  for(i=1;i<=NX;i++) RHO[l][i]/=sum;
 }
}
