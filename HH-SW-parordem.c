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
#define PI acos(-1.0)
#define NP 100                //numero de picos para comparar
#define transiente 40000


#define NR_END 1
#define FREE_ARG char*
#define N 6           // numero de equacoes
#define NHH 10   // quantidade de HH
#define numerocnl 10  // numero de conexoes nao locais
                      // NHH*NHH-3*NHH+2 numero maximo de numerocnl

 float ran1(long *idum);

 FILE *o;

void derivs(double y[],double df[],double **cnl,double epsilon);
 double *dvector(long nl,long nh);
 void free_dvector(double *v, long nl, long nh);
 void nrerror(char error_text[]);
 double **dmatrix(long nrl, long nrh, long ncl, long nch);
 void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch);


 int main(void)
{

  int j,i,sum,k[N*NHH],cont[N*NHH+5],tn,t,nk[NP*NHH][N*NHH+2],kmax[N*NHH+5],
      kinicial,kfinal,l,temp,nkk[NP*NHH][N*NHH+5],kkmax[N*NHH+5],contR;

  double x[N*NHH+2],h,a[N*NHH+2],b[N*NHH+2],c[N*NHH+2],*df,*y,**cnl,epsilon,
   xvelho[N*NHH+2],xnovo[N*NHH+2],tempo,phi[N*NHH+5],R,real,compl,Rmedio;

 long idum=(-2543251);

 y=dvector(1,N*NHH+1);
 df=dvector(1,N*NHH+1);
 cnl=dmatrix(1,NHH+1,1,NHH+1);

 // condicoes iniciais

 for(i=1;i<=N*NHH;i=i+6)
   {
    x[i]=-10.0*ran1(& idum);
    x[i+1]=ran1(& idum);
    x[i+2]=ran1(& idum);
    x[i+3]=ran1(& idum);
    x[i+4]=ran1(& idum);
    x[i+5]=ran1(& idum);
   }

    o=fopen ("/home/ewandson/C/doutorado/ML/RR/parordemHH.dat","wt");

    if (o==NULL)
    {
     puts ("erro no arquivo");
     getchar();
     exit(1);
    }

 ///////////// matriz de conexoes nao locais //////

 for(i=1;i<=NHH;i++)
   for(j=1;j<=NHH;j++)
     cnl[i][j]=0.0;

 sum=0.0;
 while(sum<numerocnl)
   {
    i=(int)NHH*ran1(& idum);
    j=(int)NHH*ran1(& idum);

    if(i==0)
      i=NHH;
    if(j==0)
      j=NHH;

    cnl[i][j]=1.0;

    // condicoes de contorno
    if(i==j)
      cnl[i][j]=0.0;
    cnl[i][i+1]=0.0;
    cnl[i][i-1]=0.0;

    sum=0.0;
    for(i=1;i<=NHH;i++)
      for(j=1;j<=NHH;j++)
	sum=sum+cnl[i][j];
   }

 ////////// conexoes locais ///////////

 for(i=1;i<=NHH-1;i++)
     cnl[i][i+1]=1.0;
 for(i=2;i<=NHH;i++)
     cnl[i][i-1]=1.0;

 /////////////////////////////////////

for(epsilon=0.0;epsilon<=0.0;epsilon=epsilon+0.0001)
 {
   Rmedio=0.0;
   contR=0.0;
   tempo=0.0;
   h=0.1; // passo de integracao
   tn=500000;      // numero de iteracoes

 for(i=1;i<=N*NHH;i++)
   cont[i]=0.0;

  for(t=1;t<=tn;t++)
    {   
     tempo=tempo+h;

     for(i=1;i<=NHH;i++)
       {
	 y[i]=x[i];
	 xnovo[i]=1.0/x[5];
       }

     for(i=1;i<=NHH;i++)
       xvelho[i]=xnovo[i];



 // ------------ Runge-Kutta --------------------         
     derivs(y,df,cnl,epsilon);  
     for(i=1;i<=N*NHH;i++)
           {
            a[i]=h*df[i];
            y[i]=x[i]+a[i]/2.0;
           }
     derivs(y,df,cnl,epsilon);
     for(i=1;i<=N*NHH;i++)
           {
            b[i]=h*df[i];
            y[i]=x[i]+b[i]/2.0;
	   }
     derivs(y,df,cnl,epsilon);
     for(i=1;i<=N*NHH;i++)
           {
            c[i]=h*df[i];
            y[i]=x[i]+c[i]; 
	   }
     derivs(y,df,cnl,epsilon);
     for(i=1;i<=N*NHH;i++)
          x[i]=x[i]+(a[i]+h*df[i])/6.0+(b[i]+c[i])/3.0;


     
     ///////////////////picos////////////

     for(i=5;i<=N*NHH;i=i+6)
       {

     xvelho[i]=1.0/x[i]; 

     if(xvelho[i]/xnovo[i]>1.0) 
	cont[i]=cont[i]+1.0; 	
     
     if(xvelho[i]/xnovo[i]<1.0 && cont[i]<2000.0)	 	   
	      cont[i]=0.0; 
	  
	  if(xvelho[i]/xnovo[i]<1.0 && cont[i]>2000.0)
	  {      
	      if(tempo>transiente)
	     {
		  k[i]=k[i]+1;		
		  nk[k[i]][i]=tempo;		 

		  //    fprintf(o,"%d %f\n",nk[k],xvelho[5]);		 
	       }	
	      cont[i]=0.0;  	  
	  }     		
	  kmax[i]=k[i];
	   
       }
    }

  ////////// CALCULO DO kinicial //////////////////


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

/*

    for(i=1;i<=N*NHH;i++)  
      k[i]=1.0;
      
  for(tempo=kinicial+1;tempo<=kfinal;tempo=tempo+h) 
     {    
      i=1;
      while(i<=N)
         {
	  if(tempo>nk[k[i]][i] && tempo<nk[kmax[i]][i])
            {
	     if(tempo<nk[k[i]+1][i])	  
		phi[i]=2*PI*k[i]+2*PI*(tempo-nk[k[i]][i])/(nk[k[i]+1][i]-
                       nk[k[i]][i]);

	     if(tempo>=nk[k[i]+1][i])
	        k[i]=k[i]+1;			  
	    }
	  i=i+1;
	 }

      

 // CALCULO DO PARAMETRO DE ORDEM

      
      if(tempo>kinicial+10 && tempo<kfinal-10)
	{
         real=0.0;
         compl=0.0;
         for(i=1;i<=N*NHH;i++)
	   {
	    real=real+cos(phi[i]);
	    compl=compl+sin(phi[i]);
	   }
	 real=real/N*NHH;
	 compl=compl/N*NHH;

	 R=sqrt(real*real+compl*compl);
	 Rmedio=Rmedio+R;
	 contR=contR+1.0;
	}
     }
     fprintf(o,"%f %f\n",epsilon,Rmedio/contR);    */
   }

 


   free_dvector(y,1,N*NHH+1);
   free_dvector(df,1,N*NHH+1); 
   free_dmatrix(cnl,1,NHH+1,1,NHH+1);
   return 0;

   fclose(o);
}

////////// equacoes Hodgkin-Huxley
void derivs(double y[],double df[],double **cnl, double epsilon)
{
 int i,j,contador;
 double T,T0,phi,tauna,tauk,tausd,tausa,anainf,akinf,asdinf,eta,gamma,rho,
     gna,gk,gsd,gsa,gl,Vna,Vk,Vsd,Vsa,Vl,Cm,Ina,Ik,Isd,Isa,Il,soma[NHH+2],
   SOMA[N*NHH+2],xx[NHH+2],r[NHH+2],taur,taud,Vsyn;

 T=8.2;
 T0=20;
 phi=pow(3.0,(T-T0)/10.0);
 tauna=0.05;
 tauk=2.0;
 tausd=10.0;
 tausa=20.0;
 taur=0.50;
 taud=8.0;
 eta=0.012;
 gamma=0.17;
 rho=pow(1.3,(T-T0)/10.0);
 gna=1.5;
 gk=2.0;
 gsd=0.25;
 gsa=0.4;
 gl=0.1;
 Vna=50.0;
 Vk=-90.0;  
 Vsd=50.0;
 Vsa=-90.0;
 Vl=-60.0; 
 Vsyn=20.0;
 Cm=1.0;
 //epsilon=0.001; // coupling parameter

 for(i=1;i<=NHH;i++)
   soma[i]=0.0;    
 for(i=1;i<=N*NHH;i++)
   SOMA[i]=0.0;
 
 contador=0.0;
 for(i=1;i<=N*NHH;i=i+6)
   {
    contador=contador+1;
    xx[contador]=y[i];
   }

 contador=0.0;
 for(i=1;i<=N*NHH;i=i+6)
   {
    contador=contador+1;
    r[contador]=y[i+5];
   }

 for(i=1;i<=NHH;i++)
   for(j=1;j<=NHH;j++)
     soma[i]=soma[i]+cnl[i][j]*r[j]*(Vsyn-xx[i]);

 contador=0.0;
 for(i=1;i<=N*NHH;i=i+6)
   {
    contador=contador+1;
    SOMA[i]=soma[contador];
   }
 
 //// y[1] -> V
 //// y[2] -> ana
 //// y[3] -> ak
 //// y[4] -> asd
 //// y[5] -> asa

 for(i=1;i<=N*NHH;i=i+6)
   {
    anainf=1.0/(1.0+exp(-0.25*(y[i]+25.0))); 
    akinf=1.0/(1.0+exp(-0.25*(y[i]+25.0)));  
    asdinf=1.0/(1.0+exp(-0.09*(y[i]+40.0)));

    Ina=rho*gna*y[i+1]*(y[i]-Vna);
    Ik=rho*gk*y[i+2]*(y[i]-Vk);
    Isd=rho*gsd*y[i+3]*(y[i]-Vsd);
    Isa=rho*gsa*y[i+4]*(y[i]-Vsa);
    Il=gl*(y[i]-Vl);

    df[i]=(1.0/Cm)*(-Ina-Ik-Isd-Isa-Il)+epsilon*SOMA[i];
    df[i+1]=(phi/tauna)*(anainf-y[i+1]);
    df[i+2]=(phi/tauk)*(akinf-y[i+2]);
    df[i+3]=(phi/tausd)*(asdinf-y[i+3]);
    df[i+4]=(phi/tausa)*(-eta*Isd-gamma*y[i+4]);
    df[i+5]=(1.0/taur-1.0/taud)*1.0/(1.0+exp(-y[i]-20.0))*(1.0-y[i+5])
            -1.0/taud*y[i+5];
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
