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
// previously calculated distribution function DIST[][][][]

void annihilation(void){
 register int i,j,kx,ky,n;
 printf("\n# of particles before annihilation = %d\n",INUM);
 // calculates the new array of particles
 INUM=0;
 for(i=1;i<=NX;i++)
  for(j=1;j<=NY;j++){
  for(kx=0;kx<2*NKX-1;kx++)
   for(ky=0;ky<2*NKY-1;ky++){
   int local_number_of_particles;
   local_number_of_particles=fabs(DIST[i][j][kx][ky]);
   // creates the new local particles in the (i,j,kx,ky)-th phase-space cell
   // the particles are uniformely distributed in space
   for(n=1;n<=local_number_of_particles;n++){
    int m;
    m=INUM+n;
    if(rnd()>0.5) PX[m]=(i-0.5+0.5*rnd())*DX;
    else PX[m]=(i-0.5-0.5*rnd())*DX;
    if(rnd()>0.5) PY[m]=(j-0.5+0.5*rnd())*DY;
    else PY[m]=(j-0.5-0.5*rnd())*DY;
    KX[m]=kx-NKX+1;
    KY[m]=ky-NKY+1;
    if(DIST[i][j][kx][ky]>0) W[m]=+1;
    else W[m]=-1;
   }
   INUM+=local_number_of_particles;
  }
 }
 printf("# of particles after annihilation  = %d\n\n",INUM);
}
