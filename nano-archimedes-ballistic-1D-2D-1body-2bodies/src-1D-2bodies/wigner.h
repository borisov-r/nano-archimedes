/* ===============================================================

   wigner.h -- This file belongs to the nano-archimedes.

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

// calculates the Wigner kernel/potential
// see my manuscript on many-body Wigner MC method.

// This is all we need to calculate the function gamma.

void wigner(void){
 register int i,j,m,n;
 register int l1,l2;

 printf("\ncalculating Wigner kernel\n");

 #pragma omp parallel private(i,j,l1,l2,m,n)
 {
  #pragma omp for
  for(i=1;i<=NX;i++){
   printf("%d/%d - ",i,NX);
   fflush(stdout);
   for(j=1;j<=NX;j++){
    for(m=0;m<2*NKX;m++)
     for(n=0;n<2*NKX;n++){
      VW[i][j][m][n]=0.;
      for(l1=-(int)(0.5*LC/DX)-1;l1<=(int)(0.5*LC/DX)+1;l1++)
       for(l2=-(int)(0.5*LC/DX)-1;l2<=(int)(0.5*LC/DX)+1;l2++){
        if((1<=(i+l1)) && ((i+l1)<=NX)
        && (1<=(i-l1)) && ((i-l1)<=NX)
        && (1<=(j+l2)) && ((j+l2)<=NX)
        && (1<=(j-l2)) && ((j-l2)<=NX))
         VW[i][j][m][n]+=sin(2.*((m-NKX+1)*DKX*(l1-0.5)*DX+(n-NKX+1)*DKX*(l2-0.5)*DX))
                        *(PHI[i+l1]+PHI[j+l2]-PHI[i-l1]-PHI[j-l2]);
      }
      VW[i][j][m][n]*=-2.*(-Q)*DX*DX/(HBAR*LC*LC);
     }
   }
  }
 }
 printf("\n");
 fflush(stdout);
}
