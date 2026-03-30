/* ===============================================================

   annihilation.h -- This file belongs to the nano-archimedes.

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

// annihilates all unnecessary particles according to the
// previously calculated distribution function DIST[][]

void annihilation(void){
 register int i,j,m,n,l;

 printf("\n# of particles before annihilation = %d\n",INUM);
 fflush(stdout);

 // calculates the new array of particles
 INUM=0;
 for(i=1;i<=NX;i++)
  for(j=1;j<=NX;j++){
   for(m=0;m<2*NKX-1;m++)
    for(n=0;n<2*NKX-1;n++){
     int local_number_of_particles;
     local_number_of_particles=fabs(DIST[i][j][m][n]);
     // creates the new local particles in the (i,j,m,n)-th phase-space cell
     // the particles are uniformely distributed in space
     for(l=1;l<=local_number_of_particles;l++){
      int k;
      k=INUM+l;
      if(rnd()>0.5) P[0][k]=(i-0.5+0.5*rnd())*DX;
      else P[0][k]=(i-0.5-0.5*rnd())*DX;
      if(rnd()>0.5) P[1][k]=(j-0.5+0.5*rnd())*DX;
      else P[1][k]=(j-0.5-0.5*rnd())*DX;
      K[0][k]=m-NKX+1;
      K[1][k]=n-NKX+1;
      if(DIST[i][j][m][n]>0) W[k]=+1;
      else W[k]=-1;
     }
     INUM+=local_number_of_particles;
    }
  }
 printf("# of particles after annihilation  = %d\n\n",INUM);
 fflush(stdout);
}
