/* ===============================================================

   gamma.h -- This file belongs to the nano-archimedes.

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

// calculates gamma(x), i.e the number of couples created per unit of time
// see my manuscript on many-body Wigner MC method

void calculate_gamma(void){
 register int i,j,m,n;

 printf("\ncalculating the gamma function\n");

 #pragma omp parallel private(i,j,m,n)
 {
  #pragma omp for
  for(i=1;i<=NX;i++){
   printf("%d/%d - ",i,NX);
   fflush(stdout);
   for(j=1;j<=NX;j++){
   GAMMA[i][j]=0.;
   for(m=0;m<2*NKX;m++)
    for(n=0;n<2*NKX;n++)
     // use VW+
     if(VW[i][j][m][n]>0.) GAMMA[i][j]+=VW[i][j][m][n];
   }
  }
 }
 printf("\n");
 fflush(stdout);
}
