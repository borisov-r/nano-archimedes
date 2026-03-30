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

// calculates gamma(x,y), i.e the number of couples created per unit of time

void calculate_gamma(void){
 register int i,j,kx,ky;
 for(i=1;i<=NX;i++)
  for(j=1;j<=NY;j++){
   GAMMA[i][j]=0.;
   // Note that in 2D the Wigner potential is NOT anti-symmetric w.r.t. the k-space.
   // (see Sellier manuscript)
   for(kx=0;kx<2*NKX;kx++)
    for(ky=0;ky<2*NKY;ky++){
     // use VW+
     if(VW[i][j][kx][ky]>0.) GAMMA[i][j]+=VW[i][j][kx][ky];
    }
 }
}
