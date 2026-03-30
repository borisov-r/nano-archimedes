/* ===============================================================

   save.h -- This file belongs to the nano-archimedes.

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

// save output files

void save(int ind)
{
 int i,j;
 FILE *fp;

 printf("\nsaving output files = %d..\n",ind);

 // save Wigner potential on 5 points
 fp=fopen("wigner_70x75.dat","w");
 for(i=0;i<2*NKX;i++){
  for(j=0;j<2*NKY;j++)
   fprintf(fp,"%g ",VW[(int)(70.e-9/DX)][(int)(NY/2)][i][j]);
  fprintf(fp,"\n");
 }
 fclose(fp);

 fp=fopen("wigner_70x37.dat","w");
 for(i=0;i<2*NKX;i++){
  for(j=0;j<2*NKY;j++)
   fprintf(fp,"%g ",VW[(int)(70.e-9/DX)][(int)(0.25*NY)][i][j]);
  fprintf(fp,"\n");
 }
 fclose(fp);

 fp=fopen("wigner_70x2.dat","w");
 for(i=0;i<2*NKX;i++){
  for(j=0;j<2*NKY;j++)
   fprintf(fp,"%g ",VW[(int)(70.e-9/DX)][2][i][j]);
  fprintf(fp,"\n");
 }
 fclose(fp);

 fp=fopen("wigner_70x148.dat","w");
 for(i=0;i<2*NKX;i++){
  for(j=0;j<2*NKY;j++)
   fprintf(fp,"%g ",VW[(int)(70.e-9/DX)][NY-2][i][j]);
  fprintf(fp,"\n");
 }
 fclose(fp);

 fp=fopen("wigner_70x108.dat","w");
 for(i=0;i<2*NKX;i++){
  for(j=0;j<2*NKY;j++)
   fprintf(fp,"%g ",VW[(int)(70.e-9/DX)][(int)(0.25*NY+NY/2)][i][j]);
  fprintf(fp,"\n");
 }
 fclose(fp);


 // saves potential
 fp=fopen("potential.dat","w");
 for(i=1;i<=NX;i++){
  for(j=1;j<=NY;j++) fprintf(fp,"%g ",PHI[i][j]);
  fprintf(fp,"\n");
 }
 fclose(fp);

 // saves gamma function
 fp=fopen("gamma.dat","w");
 for(i=1;i<=NX;i++){
  for(j=1;j<=NY;j++) fprintf(fp,"%g ",GAMMA[i][j]);
  fprintf(fp,"\n");
 }
 fclose(fp);

 char s[64];
/*
 // save distribution in some chosen (kx,ky) points
 sprintf(s,"distribution_%d.dat",ind);
 fp=fopen(s,"w");
 for(i=1;i<=NX;i++){
  for(j=-NKX+1;j<NKX;j++) fprintf(fp,"%d ",DIST[i][j+NKX-1]);
  fprintf(fp,"\n");
 }
 fclose(fp);*/

 // saves the electron density
 sprintf(s,"dens_%g_%d.dat",LC*1.e9,ind);
 fp=fopen(s,"w");
 for(i=1;i<=NX;i++){
  for(j=1;j<=NY;j++) fprintf(fp,"%g ",RHO[i][j]);
  fprintf(fp,"\n");
 }
 fclose(fp);

 // saves the phase-space density
 sprintf(s,"f_x_%g_%d.dat",LC*1.e9,ind);
 fp=fopen(s,"w");
 for(i=-NKX+1;i<NKX;i++){
  for(j=-NKY+1;j<NKY;j++) fprintf(fp,"%g ",PSI[i+NKX-1][j+NKY-1]);
  fprintf(fp,"\n");
 }
 fclose(fp);

 printf("output files saved successfully.\n\n");
}
