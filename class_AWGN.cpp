#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "class_AWGN.h"

void AWGN::Gauss_initial()
{//initilization must be done before Gauss() and AWGN_noise() are used!.
		long seed_2; 
//	long*  seed_pointer;
//	double noise;
	
	int i;

	MBIG=0x0fffffff;
	MSEED=0x03c04ffd;
	MZ=0;
	
	seed_pointer=&seed_1;

	srand( (unsigned)time( NULL ) );
	seed_2=rand()*100/32767;
	for(i=0;i<seed_2;i++)
		seed_1=rand();
return ;
}

double AWGN::Crand(long *seed)
{//=> used in Gauss() function;

/*
	long int MBIG=0x0fffffff;
	long int MSEED=0x03c04ffd;
	int MZ=0;
	double FAC;
*/

   static int inext,inextp,iff=0;
   static long ma[56],mj,mk;
   int i,j,k;
   if(*seed<0 || iff==0)
      {
      iff=1;
      mj=MSEED-(*seed<0? -*seed:*seed);
      mj%=MBIG;
      ma[55]=mj;
      for(mk=1,i=1;i<=54;i++)
       {
         j=(21*i)%55;
         ma[j]=mk;
         mk=mj-mk;
         if(mk<MZ) mk+=MBIG;
         mj=ma[j];
         }
      for(k=1;k<=4;k++)
       for(i=1;i<=55;i++)
          {
            ma[i]-=ma[1+(i+30)%55];
            if(ma[i]<MZ) ma[i]+=MBIG;
            }
      inext=0;
      inextp=31;
      *seed=1;
      }
   if(++inext==56) inext=1;
   if(++inextp==56) inextp=1;
   mj=ma[inext]-ma[inextp];
   if(mj<MZ) mj+=MBIG;
   ma[inext]=mj;
	
   FAC=1.0/MBIG;
   return mj*FAC;
}

double AWGN::Gauss(long *seed)
{//return Gaussian noise with sigma=1 and mean=0;

//    double Crand(long *seed);
   static int        iflag=0;
   static double     gset;
   double   fac,  Rs,  v1,  v2;
   if(iflag==0)
      {
      do
       {
         v1=2.0*Crand(seed)-1.0;
         v2=2.0*Crand(seed)-1.0;
         Rs=v1*v1+v2*v2;
         }while(Rs==0.0 || Rs>=1.0);
      fac=sqrt(-2.0*log((float)Rs)/Rs);
      gset=v1*fac;
      iflag=1;

      return v2*fac;
      }
   else
    {
      iflag=0;
		
      return gset;
      }
}

double AWGN::AWGN_noise(double sigma, double mean)
{//
	double t;
	t=mean+sigma*Gauss(seed_pointer);
	return t;
}
