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

// memory allocation for big arrays

printf("allocating memory for VW[][][][]\n");
VW=(double ****)malloc((NX+1)*sizeof(double***));
if(VW==NULL){
 printf("Not enough memory to allocate\ndouble ***VW[NX+1]\n");
 exit(0);
}
for(i=0;i<=NX;i++){
 VW[i]=(double ***)malloc((NY+1)*sizeof(double**));
 if(VW[i]==NULL){
  printf("Not enough memory to allocate\ndouble **VW[%d][NY+1]\n",i);
  exit(0);
 }
 for(j=0;j<=NY;j++){
  VW[i][j]=(double **)malloc((2*NKX+1)*sizeof(double*));
  if(VW[i][j]==NULL){
   printf("Not enough memory to allocate\ndouble *VW[%d][%d][2*NKX+1]\n",i,j);
   exit(0);
  }
  for(k=0;k<2*NKX+1;k++){
   VW[i][j][k]=(double *)malloc((2*NKY+1)*sizeof(double));
   if(VW[i][j][k]==NULL){
    printf("Not enough memory to allocate\ndouble VW[%d][%d][%d][2*NKX+1]\n",i,j,k);
    exit(0);
   }
  }
 }
}

printf("allocating memory for DIST[][][][]\n");
DIST=(int ****)malloc((NX+1)*sizeof(int***));
if(DIST==NULL){
 printf("Not enough memory to allocate\nint ***DIST[NX+1]\n");
 exit(0);
}
for(i=0;i<=NX;i++){
 DIST[i]=(int ***)malloc((NY+1)*sizeof(int**));
 if(DIST[i]==NULL){
  printf("Not enough memory to allocate\nint **DIST[%d][NY+1]\n",i);
  exit(0);
 }
 for(j=0;j<=NY;j++){
  DIST[i][j]=(int **)malloc((2*NKX+1)*sizeof(int*));
  if(DIST[i][j]==NULL){
   printf("Not enough memory to allocate\nint *DIST[%d][%d][2*NKX+1]\n",i,j);
   exit(0);
  }
  for(k=0;k<2*NKX+1;k++){
   DIST[i][j][k]=(int *)malloc((2*NKY+1)*sizeof(int));
   if(DIST[i][j][k]==NULL){
    printf("Not enough memory to allocate\nint DIST[%d][%d][%d][2*NKX+1]\n",i,j,k);
    exit(0);
   }
  }
 }
}

printf("allocating memory for PHI[][]\n");
PHI=(double **)malloc((NX+1)*sizeof(double*));
if(PHI==NULL){
 printf("Not enough memory to allocate\ndouble *PHI[NX+1]\n");
 exit(0);
}
for(i=0;i<=NX;i++){
 PHI[i]=(double *)malloc((NY+1)*sizeof(double));
 if(PHI[i]==NULL){
  printf("Not enough memory to allocate\ndouble PHI[%d][NY+1]\n",i);
  exit(0);
 }
}

printf("allocating memory for RHO[][]\n");
RHO=(double **)malloc((NX+1)*sizeof(double*));
if(RHO==NULL){
 printf("Not enough memory to allocate\ndouble *RHO[NX+1]\n");
 exit(0);
}
for(i=0;i<=NX;i++){
 RHO[i]=(double *)malloc((NY+1)*sizeof(double));
 if(RHO[i]==NULL){
  printf("Not enough memory to allocate\ndouble RHO[%d][NY+1]\n",i);
  exit(0);
 }
}

printf("allocating memory for GAMMA[][]\n");
GAMMA=(double **)malloc((NX+1)*sizeof(double*));
if(GAMMA==NULL){
 printf("Not enough memory to allocate\ndouble *GAMMA[NX+1]\n");
 exit(0);
}
for(i=0;i<=NX;i++){
 GAMMA[i]=(double *)malloc((NY+1)*sizeof(double));
 if(GAMMA[i]==NULL){
  printf("Not enough memory to allocate\ndouble GAMMA[%d][NY+1]\n",i);
  exit(0);
 }
}

printf("allocating memory for PSI[][]\n");
PSI=(double **)malloc((2*NKX+1)*sizeof(double*));
if(PSI==NULL){
 printf("Not enough memory to allocate\ndouble *PSI[2*NKX+1]\n");
 exit(0);
}
for(i=0;i<2*NKX+1;i++){
 PSI[i]=(double *)malloc((2*NKY+1)*sizeof(double));
 if(PSI[i]==NULL){
  printf("Not enough memory to allocate\ndouble PSI[%d][2*NKY+1]\n",i);
  exit(0);
 }
}

printf("allocating memory for KX[]\n");
KX=malloc((NPMAX+1)*sizeof(*KX));
if(KX==NULL){
 printf("Not enough memory to allocate\nint KX[NPMAX+1]\n");
 exit(0);
}

printf("allocating memory for KY[]\n");
KY=malloc((NPMAX+1)*sizeof(*KY));
if(KY==NULL){
 printf("Not enough memory to allocate\nint KY[NPMAX+1]\n");
 exit(0);
}

printf("allocating memory for W[]\n");
W=malloc((NPMAX+1)*sizeof(*W));
if(W==NULL){
 printf("Not enough memory to allocate\nint W[NPMAX+1]\n");
 exit(0);
}

printf("allocating memory for UPDATED[]\n");
UPDATED=malloc((NPMAX+1)*sizeof(*UPDATED));
if(UPDATED==NULL){
 printf("Not enough memory to allocate\nint UPDATED[NPMAX+1]\n");
 exit(0);
}

printf("allocating memory for PX[]\n");
PX=malloc((NPMAX+1)*sizeof(*PX));
if(PX==NULL){
 printf("Not enough memory to allocate\ndouble PX[NPMAX+1]\n");
 exit(0);
}

printf("allocating memory for PY[]\n");
PY=malloc((NPMAX+1)*sizeof(*PY));
if(PY==NULL){
 printf("Not enough memory to allocate\ndouble PY[NPMAX+1]\n");
 exit(0);
}

printf("allocating memory for PTIME[]\n");
PTIME=malloc((NPMAX+1)*sizeof(*PTIME));
if(PTIME==NULL){
 printf("Not enough memory to allocate\ndouble PTIME[NPMAX+1]\n");
 exit(0);
}

printf("All memory has been allocated successfully.\n\n");
