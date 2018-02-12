// RULKOV Global PARAMETRO DE ORDEM

#include<stdio.h>
#include<math.h>
#include<stdlib.h>

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)
#define PI acos(-1.0)
#define NP 500
#define N 100
#define nmax 90000
#define ntrans 50000


FILE *input;

float ran1(long *idum);

int main()
{
  long idum=-1;

  int n,i,cont[N+1],k[N+1],nk[NP][N+1],kmax[N+1],kinicial,kfinal,l,temp,
      nkk[NP][N+1],kkmax[N+1],um,dois,tres,quatro;
  double sigma,beta,xvelho[N+1],xnovo[N+1],ynovo[N+1],yvelho[N+1],
      epsilon,alfa[N+1],soma,yini[N+1],xini[N+1],phi[N+1],R,real,compl;
  

  input=fopen("/home/ewandson/C/rulkov/rulkovpadrao/rulkovacoplamento/parordem001.dat","wt");

  sigma=0.001;
  beta=0.001; 

  um = rand() % N; 
  dois = rand() % N;
  tres = rand() % N;
  quatro = rand() % N;   

  for(i=1;i<=N;i++)
    alfa[i]=5.0;

    for(i=1;i<=N;i++)
      while(alfa[i]>4.4)
        alfa[i]=4.1+ran1(& idum);
       
  alfa[um]= 4.1;
  alfa[dois] = 4.2;
  alfa[tres] = 4.3;
  alfa[quatro] = 4.4;

  for(i=1;i<=N;i++)
      alfa[i]=4.1; 

// CALCULO DOS COMECOS DOS BURSTS nk[k][i]


    for(i=1;i<=N;i++)
	xini[i]=pow(-1,1+(int) (10.0*rand()/(RAND_MAX+1.0)))*ran1(& idum);  
    for(i=1;i<=N;i++)
	yini[i]=pow(-1,1+(int) (10.0*rand()/(RAND_MAX+1.0)))*ran1(& idum);


  for(epsilon=0.01;epsilon<=0.01;epsilon=epsilon+0.001)
   {

    for(i=1;i<=N;i++)  
	xnovo[i]=xini[i];  
    for(i=1;i<=N;i++)  
	ynovo[i]=yini[i];

    for(i=1;i<=N;i++)
       {
        cont[i]=0.0;
        k[i]=0.0;
       }

  for(n=1;n<=nmax;n++) 
    {    

        for(i=1;i<=N;i++)
	{
	    xvelho[i]=xnovo[i];
	    yvelho[i]=ynovo[i];
	}


	/////////////////////RULKOV
      soma=0.0;
      for(i=1;i<=N;i++)
	soma=soma+xvelho[i];
  
      for(i=1;i<=N;i++)
	 {	  
	  xnovo[i]=alfa[i]/(1+xvelho[i]*xvelho[i])+yvelho[i]+
                   epsilon/N*soma;
                  
	  ynovo[i]=yvelho[i]-sigma*xvelho[i]-beta;         
	 }  
       
      for(i=1;i<=N;i++)
         { 
          if(yvelho[i]/ynovo[i]>1.0)         			  
	    cont[i]=cont[i]+1.0;  
     
          if(yvelho[i]/ynovo[i]<1.0 && cont[i]<20.0)	 	   
	      cont[i]=0.0;  
     
	  if(yvelho[i]/ynovo[i]<1.0 && cont[i]>20.0)
	    {      
	     if(n>ntrans)
	       {
	        k[i]=k[i]+1;		
	        nk[k[i]][i]=n; 
	       }	
	     cont[i]=0.0;  	  
	    }			    	
	 
	  kmax[i]=k[i];  
	 }
    }


 // CALCULO DO kinicial


    for(i=1;i<=N;i++)
       nkk[1][i]=nk[1][i];
    
    for(l=1;l<=N;l++)
       for(i=1;i<=N-1;i++)
          if(nkk[1][i+1]<=nkk[1][i])
	    {
             temp=nkk[1][i];
             nkk[1][i]=nkk[1][i+1];
             nkk[1][i+1]=temp;
	    }

    kinicial=nkk[1][N];


 // CALCULO DO kfinal


    for(i=1;i<=N;i++)
      {
       kkmax[i]=kmax[i];
       nkk[kkmax[i]][i]=nk[kmax[i]][i];
      }
    for(l=1;l<=N;l++)
      for(i=1;i<=N-1;i++)
	{
          if(nkk[kkmax[i+1]][i+1]<=nkk[kkmax[i]][i])
	    {   
             temp=nkk[kkmax[i]][i];
             nkk[kkmax[i]][i]=nkk[kkmax[i+1]][i+1];
             nkk[kkmax[i+1]][i+1]=temp;
	    }
	}

    kfinal=nkk[kkmax[1]][1];


// CALCULO DAS FREQUENCIAS


    for(i=1;i<=N;i++)  
      k[i]=1.0;
      
  for(n=kinicial+1;n<=kfinal;n++) 
     {    
      i=1;
      while(i<=N)
         {
	  if(n>nk[k[i]][i] && n<nk[kmax[i]][i])
            {
	     if(n<nk[k[i]+1][i])	  
		phi[i]=2*PI*k[i]+2*PI*(n-nk[k[i]][i])/(nk[k[i]+1][i]-
                       nk[k[i]][i]);

	     if(n>=nk[k[i]+1][i])
	        k[i]=k[i]+1;			  
	    }
	  i=i+1;
	 }


 // CALCULO DO PARAMETRO DE ORDEM


      if(n>kinicial+10 && n<kfinal-10)
	{
         real=0.0;
         compl=0.0;
         for(i=1;i<=N;i++)
	   {
	    real=real+cos(phi[i]);
	    compl=compl+sin(phi[i]);
	   }
	 real=real/N;
	 compl=compl/N;

	 R=sqrt(real*real+compl*compl);

	 fprintf(input,"%d %f\n",n,R);
	}
     }
   }

  fclose(input);
  return 0;    
}

float ran1(long *idum)
{
 int j;
 long k;
 static long iy=0;
 static long iv[NTAB];
 float temp;
 
 if(*idum<=0 || !iy)
   {
     if(-(*idum)<1) *idum=1;
     else *idum = -(*idum);
    for(j=NTAB+7;j>=0;j--)
      {
       k=(*idum)/IQ;
       *idum=IA*(*idum-k*IQ)-IR*k;
       if(*idum<0) *idum +=IM;
       if(j<NTAB) iv[j]=*idum;
      }
      iy=iv[0];
   }
   k=(*idum)/IQ;
   *idum=IA*(*idum-k*IQ)-IR*k;
   if(*idum<0) *idum += IM;
   j=iy/NDIV;
   iy=iv[j];
   iv[j]=*idum;
   if((temp=AM*iy)>RNMX) return RNMX;
   else return temp;
}
