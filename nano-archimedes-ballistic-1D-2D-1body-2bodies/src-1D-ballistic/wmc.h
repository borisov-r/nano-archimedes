/* ===============================================================

   wmc.h -- This file belongs to the nano-archimedes.

   nano-archimedes is a quantum mechanical simulator
   which implements the single-body Wigner Monte Carlo method.
   This code simulates a Gaussian wave-packet going towards
   a potential barrier.

   Copyright (C) 2004-2015 Jean Michel D. Sellier
   <jeanmichel.sellier@gmail.com>

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

// ensemble Wigner Monte Carlo method

void WMC(void){
 // evolution of the particles
 // and creation of (+,-) couples
 register int n;
 int i,j;
 int num;
 int all_particles_updated;
 int number_of_created_particles;
 int number_of_outside_particles;
 int created; // the purpose of this flag is only to make the routine faster
 double sum;
 double r;
 double p;
 double x0;
 double k0;
 double hmt;
 double time;
 double rdt;

 // initial settings
 number_of_outside_particles=0;
 all_particles_updated=NO;
 for(n=0;n<INUM;n++) UPDATED[n]=NO;
 for(n=0;n<INUM;n++) PTIME[n]=DT;

 for(;all_particles_updated==NO;){
  number_of_created_particles=0;
  // evolution and couples creation
  for(n=0;n<INUM;n++){
   if(UPDATED[n]==NO){
    hmt=HBAR/(MSTAR*M)*PTIME[n];
    // drift n-th particle
    x0=P[n];
    k0=K[n]*DKX;
    i=(int)(x0/DX)+1;
    // evolve position and wave vector of the n-th particle
    if((i>0 && i<=NX) && (-NKX<K[n] && K[n]<NKX)){
     P[n]=x0+hmt*k0;
     // reflective boundary conditions
     // corresponds to null wave packet on the boundaries
/*     if(P[n]<0.){
      P[n]=fabs(P[n]);
      k0*=-1.;
      K[n]*=-1;
     }
     if(P[n]>LX){
      P[n]=LX-fabs(P[n]-LX);
      k0*=-1.;
      K[n]*=-1;
     }*/
     // calculates the probability that the wave-vector
     // actually evolves according to the continuous dkx
     // check if a couple of (+,-) have to be created
     if(GAMMA[i]!=0.){
      for(time=0.;time<PTIME[n];){
       rdt=-log(rnd())/GAMMA[i];
       time+=rdt;
       if(time<PTIME[n]){
        created=NO;
        r=rnd();
        sum=0.;
        // random selection of the wave-vector
        for(j=0;((created==NO) && (j<NKX));j++){
         p=fabs(VW[i][j])/GAMMA[i];
         if((sum<=r) && (r<(sum+p))){
          number_of_created_particles+=2;
          num=INUM+number_of_created_particles;
          // select a random time interval when the creation happens
          // assign position
          P[num-2]=P[num-1]=x0+HBAR/(MSTAR*M)*time*k0;
          // assign wave-vector
          if(VW[i][j]>=0.){
           K[num-2]=K[n]+j;
           K[num-1]=K[n]-j;
          } else {
           K[num-2]=K[n]-j;
           K[num-1]=K[n]+j;
          }
          // assign quantum weight
          if(W[n]==1){
           W[num-2]=+1;
           W[num-1]=-1;
          } else {
           W[num-2]=-1;
           W[num-1]=+1;
          }
          // assign flag to evolve the particles at the next loop
          UPDATED[num-2]=UPDATED[num-1]=NO;
          // assign time
          PTIME[num-2]=PTIME[num-1]=PTIME[n]-time;
          // eventually ignore the just-created couple since
          // at least one of them is outside the device
          if(K[num-2]<=-NKX || K[num-2]>=NKX ||
             K[num-1]<=-NKX || K[num-1]>=NKX){
           num=INUM-2;
           number_of_created_particles-=2;
          }
          created=YES;
         }
         sum+=p;
        }
       }
      }
     }
    } else {
     number_of_outside_particles++;
    }
    UPDATED[n]=YES;
   }
  } // end of for(n=0;...
  INUM+=number_of_created_particles;
  printf("INUM = %d -- particles created = %d\n",INUM,number_of_created_particles);
  if(INUM>NPMAX){
   printf("Number of particles has exploded - please increase NPMAX and recompile\n");
   exit(0);
  }
  // checks if all particles have been updated
  int flag;
  flag=YES;
  for(n=0;n<INUM;n++) if(UPDATED[n]==NO) flag=NO;
  all_particles_updated=flag;
 }
 printf("-- number of particles outside = %d --\n",number_of_outside_particles);
}

