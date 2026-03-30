/* ===============================================================

   hartree.h -- This file belongs to the nano-archimedes.

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

// calculates the Hartree potentials due to the n-particles

void hartree(void){
 register int i;
 double sum;

 printf("\ncalculating the Hartree potential\n");

 for(i=1;i<=NX;i++){
  int n;
  double rho;
  rho=sum=0.;
  // calculation of the electron(s) potential
  {
   register int ii;
   double r;
   for(ii=1;ii<=NX;ii++){
    if(ii!=i){
     r=sqrt(pow((i-ii)*DX,2.));
     for(n=0;n<2;n++) sum+=RHO[n][ii]/r;
    }
   }
  }
  sum*=DX;
  for(n=0;n<2;n++) rho+=RHO[n][i]; // calculation of the actual total density
  PHI[i]=PHI0[i]-0.25*Q/(PI*EPS0)*sum; // potential barrier + hartree potential
 }

}

