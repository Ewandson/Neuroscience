#include<math.h>
#include<stdio.h>
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

#define NR_END 1
#define FREE_ARG char*
#define N 3           // numero de equacoes
#define NHR 5       // quantidade de HR
#define numerocnl 5   // numero de conexoes nao locais
                      // NHR*NHR-3*NHR+2 numero maximo de numerocnl

float ran1(long *idum);

 FILE *o;

 void derivs(double y[],double df[],double **cnl);
 double *dvector(long nl,long nh);
 void free_dvector(double *v, long nl, long nh);
 void nrerror(char error_text[]);
 double **dmatrix(long nrl, long nrh, long ncl, long nch);
 void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch);


 int main(void)
{

 int n,j,i,loop,sum;
 double x[N*NHR+2],t,h,a[N*NHR+2],b[N*NHR+2],c[N*NHR+2],*df,*y,**cnl;
 long idum=(-2543251);

 y=dvector(1,N*NHR+1);
 df=dvector(1,N*NHR+1);
 cnl=dmatrix(1,NHR+1,1,NHR+1);

 // condicoes iniciais

 for(i=1;i<=N*NHR;i=i+3)
   {
    x[i]=-ran1(& idum);
    x[i+1]=0.0;
    x[i+2]=0.0;
   }

 t=0.0;
 h=0.001; // passo de integracao
 n=1500000;

    o=fopen ("/home/ewandson/doutorado/ML/V_t_SW.dat","wt");

    if (o==NULL)
    {
     puts ("erro no arquivo");
     getchar();
     exit(1);
    }


 ///////////// matriz de conexoes nao locais //////

 for(i=1;i<=NHR;i++)
   for(j=1;j<=NHR;j++)
     cnl[i][j]=0.0;

 sum=0.0;
 while(sum<numerocnl)
   {
    i=(int)NHR*ran1(& idum);
    j=(int)NHR*ran1(& idum);

    if(i==0)
      i=NHR;
    if(j==0)
      j=NHR;

    cnl[i][j]=1.0;

    // condicoes de contorno
    if(i==j)
      cnl[i][j]=0.0;
    cnl[i][i+1]=0.0;
    cnl[i][i-1]=0.0;

    sum=0.0;
    for(i=1;i<=NHR;i++)
      for(j=1;j<=NHR;j++)
	sum=sum+cnl[i][j];
   }

 ///////////// matriz de conexoes locais //////

for(i=1;i<=NHR;i++)
  {   cnl[i][i+1]=1.0;
      cnl[i][i-1]=1.0;
  }
  cnl[1][NHR]=1.0;
  cnl[NHR][1]=1.0;

 /////////////////////////////////////////////

 for(loop=1;loop<=n;loop++)
    {   
     t=t+h;

     for(i=1;i<=N*NHR;i++)
        y[i]=x[i];

 // ------------ Runge-Kutta --------------------         
     derivs(y,df,cnl);  
     for(i=1;i<=N*NHR;i++)
           {
            a[i]=h*df[i];
            y[i]=x[i]+a[i]/2.0;
           }
     derivs(y,df,cnl);
     for(i=1;i<=N*NHR;i++)
           {
            b[i]=h*df[i];
            y[i]=x[i]+b[i]/2.0;
	   }
     derivs(y,df,cnl);
     for(i=1;i<=N*NHR;i++)
           {
            c[i]=h*df[i];
            y[i]=x[i]+c[i]; 
	   }
     derivs(y,df,cnl);
     for(i=1;i<=N*NHR;i++)
          x[i]=x[i]+(a[i]+h*df[i])/6.0+(b[i]+c[i])/3.0;  
	 
  // -----------------------------------------------


     fprintf(o,"%6f",t);
     for(i=1;i<=N*NHR/1;i=i+3)
     fprintf(o,"%6f",x[i]);
     fprintf(o,"\n");

    }

   free_dvector(y,1,N*NHR+1);
   free_dvector(df,1,N*NHR+1); 
   free_dmatrix(cnl,1,NHR+1,1,NHR+1);
   return 0;

   fclose(o);
}

////////// equacoes Hindmarsh-Rose
void derivs(double y[],double df[],double **cnl)
{
 int i,j,contador;
 double gk,gl,gca,phi,Vca,Vk,Vl,V1,V2,V3,V4,m,soma[NHR+2],SOMA[N*NHR+2],epsilon,xx[NHR+2];

gk=2.0;
gl=0.5; 
gca=1.2;
phi=(1.0/3.0);
Vca=1.0;
Vk=-0.7;
Vl=-0.5;
V1=-0.01;
V2=0.15;
V3=0.1;
V4=0.05;
m=0.005;
epsilon=0.01; 


 for(i=1;i<=NHR;i++)
   soma[i]=0.0;    
 for(i=1;i<=N*NHR;i++)
   SOMA[i]=0.0;
 
 contador=0.0;
 for(i=1;i<=N*NHR;i=i+3)
   {
    contador=contador+1;
    xx[contador]=y[i];
   }

 for(i=1;i<=NHR;i++)
   for(j=1;j<=NHR;j++)
     soma[i]=soma[i]+cnl[i][j]*(xx[j]-xx[i]);

 contador=0.0;
 for(i=1;i<=N*NHR;i=i+3)
   {
    contador=contador+1;
    SOMA[i]=soma[contador];
   }
 
 for(i=1;i<=N*NHR;i=i+3)
   {


     df[i]=-y[i+2]+-gl*(y[i]-Vl)-gca*0.5*(1+tanh((y[i]-V1)/V2))*(y[i]-Vca)-gk*y[i+1]*(y[i]-Vk)+epsilon*SOMA[i];
     df[i+1]=(0.5*(1+tanh((y[i]-V3)/V4))-y[i+1])*phi*cosh((y[i]-V3)/(2.0*V4));
     df[i+2]=m*(0.2+y[i]);
   }
}
/////////////////////////////////////////

double *dvector(long nl,long nh)
{
   double *v;
   
   v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
   if (!v) nrerror("allocation failure in dvector()");
   return v-nl+NR_END;
}


void free_dvector(double *v, long nl, long nh)
{
   free((FREE_ARG) (v+nl-NR_END));
}

void nrerror(char error_text[])
{
  fprintf(stderr,"Numerical Recipes run-time error...\n");
  fprintf(stderr,"%s\n",error_text);
  fprintf(stderr,"...now exiting to system...\n");
  exit(1);
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

double **dmatrix(long nrl, long nrh, long ncl, long nch)
{
   long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
   double **m;

   m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double*)));
   if (!m) nrerror("allocation failure 1 in matrix()");
   m += NR_END;
   m -= nrl;

   m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
   if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
   m[nrl] += NR_END;
   m[nrl] -= ncl;

   for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

   return m;
}

void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch)
{
   free((FREE_ARG) (m[nrl]+ncl-NR_END));
   free((FREE_ARG) (m+nrl-NR_END));
}
