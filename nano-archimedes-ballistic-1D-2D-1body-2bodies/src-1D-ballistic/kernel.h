/* ===============================================================

   kernel.h -- This file belongs to the nano-archimedes.

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

// calculates the Wigner kernel

// we calculate only the Wigner kernel with only positive pseudo-wave vectors
// since this potential is antisymmetric. This is all we need to calculate the
// function gamma.

void kernel(void){
 register int i,j;
 for(i=1;i<=NX;i++){
  for(j=0;j<NKX;j++){
   VW[i][j]=0.;
   register int l;
   for(l=1;l<=(int)(0.5*LC/DX)+1;l++){
    if((1<=(i+l)) && ((i+l)<=NX)
     && (1<=(i-l)) && ((i-l)<=NX)) VW[i][j]+=sin(2.*j*DKX*(l-0.5)*DX)*(PHI[i+l]-PHI[i-l]);
   }
   VW[i][j]*=-2.*(-Q)*DX/(HBAR*LC);
  }
 }
}
