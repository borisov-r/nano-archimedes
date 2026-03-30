/* ===============================================================

   density.h -- This file belongs to the nano-archimedes.

   nano-archimedes is a quantum mechanical simulator
   which implements the single-body Wigner Monte Carlo method.
   This code simulates a Gaussian wave-packet going towards
   a potential barrier.

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


// calculates probability density in x- and k-space

void density(void){
 register int i,j;

 // in x-space
 for(i=1;i<=NX;i++){
  double sum;
  sum=0.;
  for(j=-NKX+1;j<NKX;j++) sum+=FW1[i][j+NKX-1];
  DENSX[i]=sum*DKX;
 }
 DENSX[1]=DENSX[NX]=0.;
 
 // in k-space
 for(j=-NKX+1;j<NKX;j++){
  double sum;
  sum=0.;
  for(i=1;i<=NX;i++) sum+=FW1[i][j+NKX-1];
  DENSK[j+NKX-1]=sum*DX;
 }

}
