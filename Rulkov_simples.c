// RULKOV - REDES SEM ESCALA  - PARAMETRO DE ORDEM

#include<stdio.h>
#include<math.h>
#include<stdlib.h>


//////////// PARAMETROS DO RAN1 ///////////////////

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

////////////////////////////////////////////////////


#define PI acos(-1.0)
#define nmax 80000      // MAXIMO DE ITERACOES
#define transiente 50000    // TRANSIENTE

#define NR_END 1
#define FREE_ARG char*


FILE *input;

int **matrix(long nrl, long nrh, long ncl, long nch);
void free_matrix(int **m, long nrl, long nrh, long ncl, long nch);
void nrerror(char error_text[]);
int inteiro(double x);

float ran1(long *idum);

int main()
{
  long idum=-1;

  int n,i,j;
  double sigma,beta,xvelho,xnovo,ynovo,
      yvelho,epsilon,alfa;
 
    
  input=fopen("teste.dat","wt");


    

  /////////////// RULKOV ///////////////////////////

  /////////////////  PARAMETROS ///////////////////

  sigma=0.001;
  beta=0.001; 
  alfa=4.1;


   /////condição inicial  
    xnovo=pow(-1,1+(int) (10.0*rand()/(RAND_MAX+1.0)))*ran1(& idum);  
     
    ynovo=pow(-1,1+(int) (10.0*rand()/(RAND_MAX+1.0)))*ran1(& idum);
 
 
  /////mapa de rulkov
   for(n=1;n<=nmax;n++) 
    {    
      xvelho=xnovo;
      yvelho=ynovo;
    	

     
       xnovo=alfa/(1+xvelho*xvelho)+yvelho;
	 
       ynovo=yvelho-sigma*xvelho-beta;    


     if(n>transiente)
       fprintf(input,"%d %9f  %9f\n",n,xnovo,ynovo);      
	      
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

int inteiro(double x)
{
  return(x);
}

int **matrix(long nrl, long nrh, long ncl, long nch)
{
   long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
   int **m;

   m=(int **) malloc((size_t)((nrow+NR_END)*sizeof(int*)));
   if (!m) nrerror("allocation failure 1 in matrix()");
   m += NR_END;
   m -= nrl;

   m[nrl]=(int *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(int)));
   if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
   m[nrl] += NR_END;
   m[nrl] -= ncl;

   for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

   return m;
}
  
void free_matrix(int **m, long nrl, long nrh, long ncl, long nch)
{
   free((FREE_ARG) (m[nrl]+ncl-NR_END));
   free((FREE_ARG) (m+nrl-NR_END));
}

void nrerror(char error_text[])
{
  fprintf(stderr,"Numerical Recipes run-time error...\n");
  fprintf(stderr,"%s\n",error_text);
  fprintf(stderr,"...now exiting to system...\n");
  exit(1);
}
  
