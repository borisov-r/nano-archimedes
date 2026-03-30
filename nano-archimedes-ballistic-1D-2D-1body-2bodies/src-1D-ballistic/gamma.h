/* ===============================================================

   gamma.h -- This file belongs to the nano-archimedes.

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

// calculates gamma(x), i.e the number of couples created per unit of time

void calculate_gamma(void){
 register int i,j;
 for(i=1;i<=NX;i++){
  GAMMA[i]=0.;
  // the implementation below holds taking into account the fact that
  // the Wigner potential is anti-symmetric w.r.t. the k-space.
  for(j=1;j<NKX;j++) GAMMA[i]+=fabs(VW[i][j]); // remember that VW[i][0]=0.;
 }
}
