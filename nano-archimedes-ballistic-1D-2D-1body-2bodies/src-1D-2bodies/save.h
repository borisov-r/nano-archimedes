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
 int m,n;
 int sum;
 FILE *fp;

 char s[64];

 printf("\nsaving output files - %d\n",ind);
 fflush(stdout);

 if(ind==0){
  // saves potential
  fp=fopen("potential.dat","w");
  for(i=1;i<=NX;i++){
   fprintf(fp,"%g %g\n",(i-0.5)*DX,PHI[i]);
  }
  fclose(fp);

  // saves auxiliary files
  fp=fopen("x.dat","w");
  for(i=1;i<=NX;i++){
   fprintf(fp,"%g\n",(i-0.5)*DX*1.e9); // in nanometer
  }
  fclose(fp);

  fp=fopen("k.dat","w");
  for(i=-NKX+1;i<NKX;i++){
   fprintf(fp,"%g\n",i*DKX*1.e-9); // in 1/nm
  }
  fclose(fp);
 }

 // saves gamma function
 sprintf(s,"gamma_%d.dat",ind);
 fp=fopen(s,"w");
 for(i=1;i<=NX;i++){
  for(j=1;j<=NX;j++){
   fprintf(fp,"%g ",GAMMA[i][j]);
  }
  fprintf(fp,"\n");
 }
 fclose(fp);

 // saves the distribution in 2-space
 sprintf(s,"distribution_x_%d.dat",ind);
 fp=fopen(s,"w");
 for(i=1;i<=NX;i++){
  for(j=1;j<=NX;j++){
   sum=0;
   for(m=0;m<2*NKX;m++)
    for(n=0;n<2*NKX;n++){
     sum+=DIST[i][j][m][n];
    }
   fprintf(fp,"%d ",sum);
  }
  fprintf(fp,"\n");
 }
 fclose(fp);

 // saves the distribution in 2-wave space
 sprintf(s,"distribution_k_%d.dat",ind);
 fp=fopen(s,"w");
 for(i=1;i<2*NKX;i++){
  for(j=1;j<2*NKX;j++){
   sum=0;
   for(m=1;m<=NX;m++)
    for(n=1;n<=NX;n++){
     sum+=DIST[m][n][i][j];
    }
   fprintf(fp,"%d ",sum);
  }
  fprintf(fp,"\n");
 }
 fclose(fp);

 // saves the two distributions for each particle
 for(i=1;i<=NX;i++)
  for(m=1;m<2*NKX;m++)
   FTOT[i][m]=0;
 // particle #0
 sprintf(s,"distribution_xk_0_%d.dat",ind);
 fp=fopen(s,"w");
 for(i=1;i<=NX;i++){
  for(m=1;m<2*NKX;m++){
   sum=0;
   for(j=1;j<=NX;j++)
    for(n=1;n<2*NKX;n++){
     sum+=DIST[i][j][m][n];
    }
   fprintf(fp,"%d ",sum);
   FTOT[i][m]=sum;
  }
  fprintf(fp,"\n");
 }
 fclose(fp);
 // particle #1
 sprintf(s,"distribution_xk_1_%d.dat",ind);
 fp=fopen(s,"w");
 for(j=1;j<=NX;j++){
  for(n=1;n<2*NKX;n++){
   sum=0;
   for(i=1;i<=NX;i++)
    for(m=1;m<2*NKX;m++){
     sum+=DIST[i][j][m][n];
    }
   fprintf(fp,"%d ",sum);
   FTOT[j][n]+=sum;
  }
  fprintf(fp,"\n");
 }
 fclose(fp);

 // saves densities
 sprintf(s,"dens_%d.dat",ind);
 fp=fopen(s,"w");
 for(i=1;i<=NX;i++){
  int sum=0;
  for(m=1;m<2*NKX;m++) sum+=FTOT[i][m];
  fprintf(fp,"%g %g %g %d\n",(i-0.5)*DX,RHO[0][i],RHO[1][i],sum);
 }
 fclose(fp);

 // saves the induced dipole for the electron
 if(ind==0){
  for(i=1;i<=NX;i++) RHO0[i]=RHO[0][i];
  fp=fopen("induced_dipole.dat","w");
  fprintf(fp,"%g %g\n",ind*DT*1.e15,0.); // in attoseconds
  fclose(fp);
 } else {
  fp=fopen("induced_dipole.dat","a");
  double sum;
  sum=0.;
  for(i=1;i<=NX;i++) sum+=(i-0.5)*(RHO[0][i]-RHO0[i]);
  sum*=DX*DX;
  fprintf(fp,"%g %g\n",ind*DT*1.e15,sum);
  fclose(fp);
 }
}
