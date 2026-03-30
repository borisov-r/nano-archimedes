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

// free the allocated memory
printf("\nfreeing memory...\n");

free(K[0]);
free(K[1]);
free(W);
free(UPDATED);
free(P[0]);
free(P[1]);
free(PTIME);

for(i=0;i<=NX;i++){
 for(j=0;j<=NX;j++){
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

for(i=0;i<=NX;i++) free(GAMMA[i]);
free(GAMMA);

