/* ===============================================================

   alloc.h -- This file belongs to the nano-archimedes.

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

printf("allocating memory\n");

printf("allocating VW\n");
VW=(double ****)malloc((NX+1)*sizeof(double***));
if(VW==NULL){
 printf("Not enough memory to allocate\ndouble ***VW[NX+1]\n");
 exit(0);
}
for(i=0;i<=NX;i++){
 VW[i]=(double ***)malloc((NX+1)*sizeof(double**));
 if(VW[i]==NULL){
  printf("Not enough memory to allocate\ndouble **VW[%d][NX+1]\n",i);
  exit(0);
 }
 for(j=0;j<=NX;j++){
  VW[i][j]=(double **)malloc((2*NKX+1)*sizeof(double*));
  if(VW[i][j]==NULL){
   printf("Not enough memory to allocate\ndouble *VW[%d][%d][2*NKX+1]\n",i,j);
   exit(0);
  }
  for(k=0;k<2*NKX+1;k++){
   VW[i][j][k]=(double *)malloc((2*NKX+1)*sizeof(double));
   if(VW[i][j][k]==NULL){
    printf("Not enough memory to allocate\ndouble VW[%d][%d][%d][2*NKX+1]\n",i,j,k);
    exit(0);
   }
  }
 }
}

printf("allocating DIST\n");
DIST=(int ****)malloc((NX+1)*sizeof(int***));
if(DIST==NULL){
 printf("Not enough memory to allocate\nint ***DIST[NX+1]\n");
 exit(0);
}
for(i=0;i<=NX;i++){
 DIST[i]=(int ***)malloc((NX+1)*sizeof(int**));
 if(DIST[i]==NULL){
  printf("Not enough memory to allocate\nint **DIST[%d][NX+1]\n",i);
  exit(0);
 }
 for(j=0;j<=NX;j++){
  DIST[i][j]=(int **)malloc((2*NKX+1)*sizeof(int*));
  if(DIST[i][j]==NULL){
   printf("Not enough memory to allocate\nint *DIST[%d][%d][2*NKX+1]\n",i,j);
   exit(0);
  }
  for(k=0;k<2*NKX+1;k++){
   DIST[i][j][k]=(int *)malloc((2*NKX+1)*sizeof(int));
   if(DIST[i][j][k]==NULL){
    printf("Not enough memory to allocate\nint DIST[%d][%d][%d][2*NKX+1]\n",i,j,k);
    exit(0);
   }
  }
 }
}

printf("allocating GAMMA\n");
GAMMA=(double **)malloc((NX+1)*sizeof(double*));
if(GAMMA==NULL){
 printf("Not enough memory to allocate\ndouble *GAMMA[NX+1]\n");
 exit(0);
}
for(i=0;i<=NX;i++){
 GAMMA[i]=(double *)malloc((NX+1)*sizeof(double));
 if(GAMMA[i]==NULL){
  printf("Not enough memory to allocate\ndouble GAMMA[%d][NX+1]\n",i);
  exit(0);
 }
}

printf("allocating K\n");
for(i=0;i<2;i++){
 K[i]=malloc((NPMAX+1)*sizeof(*K[i]));
 if(K[i]==NULL){
  printf("Not enough memory to allocate\nint K[%d][NPMAX+1]\n",i);
  exit(0);
 }
}

printf("allocating W\n");
W=malloc((NPMAX+1)*sizeof(*W));
if(W==NULL){
 printf("Not enough memory to allocate\nint W[NPMAX+1]\n");
 exit(0);
}

printf("allocating UPDATED\n");
UPDATED=malloc((NPMAX+1)*sizeof(*UPDATED));
if(UPDATED==NULL){
 printf("Not enough memory to allocate\nint UPDATED[NPMAX+1]\n");
 exit(0);
}

printf("allocating P\n");
for(i=0;i<2;i++){
 P[i]=malloc((NPMAX+1)*sizeof(*P[i]));
 if(P[i]==NULL){
  printf("Not enough memory to allocate\ndouble P[%d][NPMAX+1]\n",i);
  exit(0);
 }
}

printf("allocating PTIME\n");
PTIME=malloc((NPMAX+1)*sizeof(*PTIME));
if(PTIME==NULL){
 printf("Not enough memory to allocate\ndouble PTIME[NPMAX+1]\n");
 exit(0);
}

