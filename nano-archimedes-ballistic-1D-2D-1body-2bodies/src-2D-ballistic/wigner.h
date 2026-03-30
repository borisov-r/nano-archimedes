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

// These calculations are useful to calculate the function gamma.

void wigner(void){
 register int i,j,kx,ky;
 for(i=1;i<=NX;i++){
  printf("%d/%d - ",i,NX);
  fflush(stdout);
  for(j=1;j<=NY;j++){
   for(kx=0;kx<2*NKX;kx++)
    for(ky=0;ky<2*NKY;ky++){
     VW[i][j][kx][ky]=0.;
     register int l,m;
     for(l=1;l<=(int)(0.5*LC/DX)+1;l++)
      for(m=-(int)(LC/DY)-1;m<=(int)(LC/DY)+1;m++){
       if((1<=(i+l)) && ((i+l)<=NX) &&
          (1<=(i-l)) && ((i-l)<=NX) &&
          (1<=(j+m)) && ((j+m)<=NY) &&
          (1<=(j-m)) && ((j-m)<=NY)) VW[i][j][kx][ky]+=sin(2.*(kx-NKX+1)*DKX*(l-0.5)*DX+2.*(ky-NKY+1)*DKY*(m-0.5)*DY)
                                                      *(PHI[i+l][j+m]-PHI[i-l][j-m]);
      }
      VW[i][j][kx][ky]*=-2.*(-Q)*DX*DY/(HBAR*LC*LC);
    }
  }
 }
 printf("\n");
}
