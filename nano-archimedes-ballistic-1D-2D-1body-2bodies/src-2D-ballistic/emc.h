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

void EMC(void){
 // evolution of the particles
 // and creation of (+,-) couples
 register int n;
 int i,j;
 int kx,ky;
 int all_particles_updated;
 int number_of_created_particles;
 int number_of_outside_particles;
 double r;
 double x0;
 double y0;
 double kx0;
 double ky0;
 double hmt;

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
    x0=PX[n];
    y0=PY[n];
    kx0=KX[n]*DKX;
    ky0=KY[n]*DKY;
    i=(int)(x0/DX)+1;
    j=(int)(y0/DY)+1;
    // evolve position and wave vector of the n-th particle
    if((i>0 && i<=NX) && (-NKX<KX[n] && KX[n]<NKX) &&
       (j>0 && j<=NY) && (-NKY<KY[n] && KY[n]<NKY)){
     PX[n]=x0+hmt*kx0;
     PY[n]=y0+hmt*ky0;
     // calculates the probability that the wave-vector
     // actually evolves according to the continuous dkx
     // check if a couple of (+,-) have to be created
     if(GAMMA[i][j]!=0.){
      double time;
      double rdt;
      for(time=0.;time<PTIME[n];){
       rdt=-log(rnd())/GAMMA[i][j];
       time+=rdt;
       if(time<PTIME[n]){
        int created; // the purpose of this flag is only to make the routine faster
        double sum;
        created=NO;
        r=rnd();
        sum=0.;
        // random selection of the wave-vector
        for(kx=-NKX+1;((created==NO) && (kx<NKX));kx++)
         for(ky=-NKY+1;((created==NO) && (ky<NKY));ky++){
          double p;
//printf("here 1\n");
          if(VW[i][j][kx+NKX-1][ky+NKY-1]>0.) p=VW[i][j][kx+NKX-1][ky+NKY-1]/GAMMA[i][j];
          else p=0.;
//printf("here 2\n");
//printf("n = %d --- sum = %g  ---  p = %g\n",n,sum,p);
//if(n!=0) exit(0);
          if((sum<=r) && (r<(sum+p))){
           number_of_created_particles+=2;
           int num=INUM+number_of_created_particles;
           // select a random time interval when the creation happens
           // assign position
           PX[num-2]=PX[num-1]=x0+HBAR/(MSTAR*M)*time*kx0;
           PY[num-2]=PY[num-1]=y0+HBAR/(MSTAR*M)*time*ky0;
           // assign wave-vector
           if(VW[i][j][kx+NKX-1][ky+NKY-1]>=0.){
            KX[num-2]=KX[n]+kx;
            KX[num-1]=KX[n]-kx;
            KY[num-2]=KY[n]+ky;
            KY[num-1]=KY[n]-ky;
           } else {
            KX[num-2]=KX[n]-kx;
            KX[num-1]=KX[n]+kx;
            KY[num-2]=KY[n]-ky;
            KY[num-1]=KY[n]+ky;
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
           if(KX[num-2]<=-NKX || KX[num-2]>=NKX ||
              KX[num-1]<=-NKX || KX[num-1]>=NKX ||
              KY[num-2]<=-NKY || KY[num-2]>=NKY ||
              KY[num-1]<=-NKY || KY[num-1]>=NKY){
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

