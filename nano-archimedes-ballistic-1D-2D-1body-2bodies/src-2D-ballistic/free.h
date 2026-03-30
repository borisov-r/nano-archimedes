/* ===============================================================

   free.h -- This file belongs to the nano-archimedes.

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

// free the memory from big arrays

free(KX);
free(KY);
free(W);
free(UPDATED);
free(PX);
free(PY);
free(PTIME);
for(i=0;i<=NX;i++){
 free(RHO[i]);
 free(PHI[i]);
 free(GAMMA[i]);
}
free(RHO);
free(PHI);
free(GAMMA);
for(i=0;i<2*NKX+1;i++){
 free(PSI[i]);
}
free(PSI);

for(i=0;i<=NX;i++){
 for(j=0;j<=NY;j++){
  for(k=0;k<2*NKX+1;k++){
   free(VW[i][j][k]);
   free(DIST[i][j][k]);
  }
  free(VW[i][j]);
  free(DIST[i][j]);
 }
 free(VW[i]);
 free(DIST[i]);
}
free(VW);
free(DIST);

