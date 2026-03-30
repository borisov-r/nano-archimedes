/* ===============================================================

   save.h -- This file belongs to the nano-archimedes.

   nano-archimedes is a quantum mechanical simulator
   which implements the single-body Wigner Monte Carlo method.
   This code simulates a Gaussian wave-packet going towards
   a potential barrier.

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

 if(ind==0 || ind==1){
  // saves potential
  fp=fopen("potential.dat","w");
  for(i=1;i<=NX;i++){
   fprintf(fp,"%g %g\n",(i-0.5)*DX,PHI[i]);
  }
  fclose(fp);

  // saves gamma function
  fp=fopen("gamma.dat","w");
  for(i=1;i<=NX;i++){
   fprintf(fp,"%g %g\n",(i-0.5)*DX,GAMMA[i]);
  }
  fclose(fp);

  // saves the coordinate axis values
  fp=fopen("x.dat","w");
  for(i=1;i<=NX;i++){
   fprintf(fp,"%g\n",(i-0.5)*DX);
  }
  fclose(fp);
  fp=fopen("k.dat","w");
  for(i=-NKX;i<NKX;i++){
   fprintf(fp,"%g\n",(i+0.5)*DKX);
  }
  fclose(fp);
 }

 char s[64];

 // saves normalized the Wigner quasi-distribution
 // ==============================================
 sprintf(s,"wigner_quasi_distribution_%d.dat",ind);
 fp=fopen(s,"w");
 for(i=1;i<=NX;i++){
  for(j=-NKX;j<NKX;j++){
   if(j==-NKX) fprintf(fp,"%g ",FW1[i][-NKX+1+NKX-1]);
   else fprintf(fp,"%g ",FW1[i][j+NKX-1]);
  }
  fprintf(fp,"\n");
 }
 fclose(fp);

 // saves the electron probability density in x-space
 sprintf(s,"wigner_probability_density_%d.dat",ind);
 fp=fopen(s,"w");
 for(i=1;i<=NX;i++) fprintf(fp,"%g %g\n",(i-0.5)*DX,DENSX[i]);
 fclose(fp);

 // saves the electron probability density in k-space
 sprintf(s,"wigner_k_space_probability_density_%d.dat",ind);
 fp=fopen(s,"w");
 for(i=-NKX;i<NKX;i++){
  if(i==-NKX) fprintf(fp,"%g %g\n",(i+0.5)*DKX,DENSK[0]);
  else fprintf(fp,"%g %g\n",(i+0.5)*DKX,DENSK[i+NKX-1]);
 }
 fclose(fp);
}
