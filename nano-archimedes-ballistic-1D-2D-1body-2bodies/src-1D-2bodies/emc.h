/* ===============================================================

   emc.h -- This file belongs to the nano-archimedes.

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

// ensemble Monte Carlo method
// no need to parallelize this routine in 1D space,
// the total number of particles is relatively small.

void EMC(void){
 // evolution of the particles
 // and creation of (+,-) couples
 register int n;
 int i[2],j[2];
 int all_particles_updated;
 int number_of_created_particles;
 int number_of_outside_particles;
 double r;
 double x0[2];
 double k0[2];
 double hmt[2];

 printf("\nEvolving particles\n");
 fflush(stdout);

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
    hmt[0]=HBAR/(ME)*PTIME[n];
    hmt[1]=HBAR/(MP)*PTIME[n];
    // drift n-th particle
    x0[0]=P[0][n];
    x0[1]=P[1][n];
    k0[0]=K[0][n]*DKX;
    k0[1]=K[1][n]*DKX;
    i[0]=(int)(x0[0]/DX)+1;
    i[1]=(int)(x0[1]/DX)+1;
    // evolve position and wave vector of the n-th particle
    if((i[0]>0 && i[0]<=NX) && (-NKX<K[0][n] && K[0][n]<NKX) &&
       (i[1]>0 && i[1]<=NX) && (-NKX<K[1][n] && K[1][n]<NKX)){
     P[0][n]=x0[0]+hmt[0]*k0[0];
     P[1][n]=x0[1]+hmt[1]*k0[1];
     // reflective boundary conditions
     // corresponds to null wave packet on the boundaries
     if(P[0][n]<0.){
      K[0][n]*=-1;
      k0[0]*=-1;
      P[0][n]=fabs(P[0][n]);
      x0[0]=P[0][n];
     }
     if(P[0][n]>LX){
      K[0][n]*=-1;
      k0[0]*=-1;
      P[0][n]=LX-fabs(P[0][n]-LX);
      x0[0]=P[0][n];
     }
     if(P[1][n]<0.){
      K[1][n]*=-1;
      k0[1]*=-1;
      P[1][n]=fabs(P[1][n]);
      x0[1]=P[1][n];
     }
     if(P[1][n]>LX){
      K[1][n]*=-1;
      k0[1]*=-1;
      P[1][n]=LX-fabs(P[1][n]-LX);
      x0[1]=P[1][n];
     }
     // calculates the probability that the wave-vector
     // actually evolves according to the continuous dkx
     // check if a couple of (+,-) have to be created
     if(GAMMA[i[0]][i[1]]!=0.){
      double time;
      double rdt;
      for(time=0.;time<PTIME[n];){
       rdt=-log(rnd())/GAMMA[i[0]][i[1]];
       time+=rdt;
       if(time<PTIME[n]){
        int created; // the purpose of this flag is only to make the routine faster
        double sum;
        created=NO;
        r=rnd();
        sum=0.;
        // random selection of the wave-vector
        for(j[0]=0;((created==NO) && (j[0]<2*NKX));j[0]++)
        for(j[1]=0;((created==NO) && (j[1]<2*NKX));j[1]++){
         double p;
         p=fabs(VW[i[0]][i[1]][j[0]][j[1]])/GAMMA[i[0]][i[1]];
         if((sum<=r) && (r<(sum+p))){
          number_of_created_particles+=2;
          int num=INUM+number_of_created_particles;
          // select a random time interval when the creation happens
          // assign position
          P[0][num-2]=P[0][num-1]=x0[0]+HBAR/(ME)*time*k0[0];
          P[1][num-2]=P[1][num-1]=x0[1]+HBAR/(MP)*time*k0[1];
          // assign wave-vector
          if(VW[i[0]][i[1]][j[0]][j[1]]>=0.){
           K[0][num-2]=K[0][n]+j[0]-NKX+1;
           K[0][num-1]=K[0][n]-j[0]+NKX-1;
           K[1][num-2]=K[1][n]+j[1]-NKX+1;
           K[1][num-1]=K[1][n]-j[1]+NKX-1;
          } else {
           K[0][num-2]=K[0][n]-j[0]+NKX-1;
           K[0][num-1]=K[0][n]+j[0]-NKX+1;
           K[1][num-2]=K[1][n]-j[1]+NKX-1;
           K[1][num-1]=K[1][n]+j[1]-NKX+1;
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
          if(K[0][num-2]<=-NKX || K[0][num-2]>=NKX ||
             K[0][num-1]<=-NKX || K[0][num-1]>=NKX ||
             K[1][num-2]<=-NKX || K[1][num-2]>=NKX ||
             K[1][num-1]<=-NKX || K[1][num-1]>=NKX){
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
  }
  INUM+=number_of_created_particles;
  printf("INUM = %d -- particles created = %d\n",INUM,number_of_created_particles);
  fflush(stdout);
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
 fflush(stdout);
}

