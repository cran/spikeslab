////**********************************************************************
////**********************************************************************
////
////  SPIKE AND SLAB 1.0.1
////
////  Copyright 2010, Cleveland Clinic Foundation
////
////  This program is free software; you can redistribute it and/or
////  modify it under the terms of the GNU General Public License
////  as published by the Free Software Foundation; either version 2
////  of the License, or (at your option) any later version.
////
////  This program is distributed in the hope that it will be useful,
////  but WITHOUT ANY WARRANTY; without even the implied warranty of
////  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
////  GNU General Public License for more details.
////
////  You should have received a copy of the GNU General Public
////  License along with this program; if not, write to the Free
////  Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
////  Boston, MA  02110-1301, USA.
////
////  ----------------------------------------------------------------
////  Project Partially Funded By:
////    --------------------------------------------------------------
////    National Science Foundation, Grants DMS-0705037, DMS-0405675 and DMS-0405072
////
////    Hemant Ishwaran, Ph.D.
////    Dept of Quantitative Health Sciences/Wb4
////    Cleveland Clinic Foundation
////    9500 Euclid Avenue
////    Cleveland, OH 44195
////
////    email:  hemant.ishwaran@gmail.com
////    phone:  216-444-9932
////    URL:    www.bio.ri.ccf.org/Resume/Pages/Ishwaran/ishwaran.html
////
////
////	J. Sunil Rao, Ph.D.
////    Deparment of Biostatistics
////    University of Miami
////
////    email: rao.jsunil@gmail.com
////
////    --------------------------------------------------------------
////    Case Western Reserve University/Cleveland Clinic  
////    CTSA Grant:  XX1 RR000000, National Center for
////    Research Resources (NCRR), NIH
////
////  ----------------------------------------------------------------
////  Written by:
////    --------------------------------------------------------------
////    Hemant Ishwaran, Ph.D.
////    Dept of Quantitative Health Sciences/Wb4
////    Cleveland Clinic Foundation
////    9500 Euclid Avenue
////    Cleveland, OH 44195
////
////    email:  hemant.ishwaran@gmail.com
////    phone:  216-444-9932
////    URL:    www.bio.ri.ccf.org/Resume/Pages/Ishwaran/ishwaran.html
////
////  ----------------------------------------------------------------
////  Maintained by:
////    Udaya B. Kogalur, Ph.D.
////    Dept of Quantitative Health Sciences/Wb4
////    Cleveland Clinic Foundation
////    
////    Kogalur Shear Corporation
////    5425 Nestleway Drive, Suite L1
////    Clemmons, NC 27012
////
////    email:  kogalurshear@gmail.com
////    phone:  919-824-9825
////    URL:    www.kogalur-shear.com
////    --------------------------------------------------------------
////
////**********************************************************************
////**********************************************************************

#include <stdio.h>
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#define FREE_ARG char*
#define MAXLINE 10000
#define MCOL    30
#define NR_END 1
#define PI       3.141592654
#define SYSSMALL 1.0e-10
#define SYSLARGE 1.0e+99
static  double dsqrarg;
#define DSQR(a) ((dsqrarg=(a)) == 0.0 ? 0.0 : dsqrarg*dsqrarg)
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
double         **convert_dmatrix(double *a,long nrl,long nrh,long ncl,long nch);
float          **convert_matrix(float *a,long nrl,long nrh,long ncl,long nch);
unsigned char   *cvector(long nl,long nh);
double        ***d3tensor(long nrl,long nrh,long ncl,long nch,long ndl,long ndh);
double         **dmatrix(long nrl,long nrh,long ncl,long nch);
double         **dsubmatrix(double **a,long oldrl,long oldrh,long oldcl,long oldch,
	                   long newrl,long newcl);
double          *dvector(long nl,long nh);
float         ***f3tensor(long nrl,long nrh,long ncl,long nch,long ndl,long ndh);
int            **imatrix(long nrl,long nrh,long ncl,long nch);
int             *ivector(long nl,long nh);
unsigned long   *lvector(long nl,long nh);
float          **matrix(long nrl,long nrh,long ncl,long nch);
void             nrerror(char error_text[]);
float          **submatrix(float **a,long oldrl,long oldrh,long oldcl,long oldch,
	                  long newrl,long newcl);
float           *vector(long nl,long nh);
void free_convert_dmatrix(double **b,long nrl,long nrh,long ncl,long nch);
void free_convert_matrix(float **b,long nrl,long nrh,long ncl,long nch);
void free_cvector(unsigned char *v,long nl,long nh);
void free_d3tensor(double ***t,long nrl,long nrh,long ncl,long nch,
	           long ndl,long ndh);
void free_dmatrix(double **m,long nrl,long nrh,long ncl,long nch);
void free_dsubmatrix(double **b,long nrl,long nrh,long ncl,long nch);
void free_dvector(double *v,long nl,long nh);
void free_f3tensor(float ***t,long nrl,long nrh,long ncl,long nch,
	           long ndl,long ndh);
void free_imatrix(int **m,long nrl,long nrh,long ncl,long nch);
void free_ivector(int *v,long nl,long nh);
void free_lvector(unsigned long *v,long nl,long nh);
void free_matrix(float **m,long nrl,long nrh,long ncl,long nch);
void free_submatrix(float **b,long nrl,long nrh,long ncl,long nch);
void free_vector(float *v,long nl,long nh);
void fWmatrix(FILE *fp,double **a,int r,int c)
{
    int i,j;
    if (r==0 || c==0) return;
    for (i=1;i<=r;i++){
      for (j=1;j<=c;j++){     
         fprintf(fp,"%15.8f ",a[i][j]);
      }
      fprintf(fp,"\n");
    }
}
void fWvector(FILE *fp,double *a,int d)
{
    int i;
    if (d==0) return;
    for (i=1;i<=d;i++) fprintf(fp,"%15.8f ",a[i]);
    fprintf(fp,"\n");
}
double gammln(double xx)
{
	double x,y,tmp,ser;
	static double cof[6]={76.18009172947146,-86.50532032941677,
		24.01409824083091,-1.231739572450155,
		0.1208650973866179e-2,-0.5395239384953e-5};
	int j;
	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<=5;j++) ser += cof[j]/++y;
	return -tmp+log(2.5066282746310005*ser/x);
}
#define ITMAX 500
#define EPS 3.0e-7
void gser(double *gamser, double a, double x, double *gln)
{
	double gammln(double xx);
	void nrerror(char error_text[]);
	int n;
	double sum,del,ap;
	*gln=gammln(a);
	if (x <= 0.0) {
		if (x < 0.0) nrerror("x less than 0 in routine gser");
		*gamser=0.0;
		return;
	} else {
		ap=a;
		del=sum=1.0/a;
		for (n=1;n<=ITMAX;n++) {
			++ap;
			del *= x/ap;
			sum += del;
			if (fabs(del) < fabs(sum)*EPS) {
				*gamser=sum*exp(-x+a*log(x)-(*gln));
				return;
			}
		}
		return;
	}
}
#undef ITMAX
#undef EPS
#define ITMAX 500
#define EPS 3.0e-7
#define FPMIN 1.0e-30
void gcf(double *gammcf, double a, double x, double *gln)
{
	double gammln(double xx);
	void nrerror(char error_text[]);
	int i;
	double an,b,c,d,del,h;
	*gln=gammln(a);
	b=x+1.0-a;
	c=1.0/FPMIN;
	d=1.0/b;
	h=d;
	for (i=1;i<=ITMAX;i++) {
		an = -i*(i-a);
		b += 2.0;
		d=an*d+b;
		if (fabs(d) < FPMIN) d=FPMIN;
		c=b+an/c;
		if (fabs(c) < FPMIN) c=FPMIN;
		d=1.0/d;
		del=d*c;
		h *= del;
		if (fabs(del-1.0) < EPS) break;
	}
	*gammcf=exp(-x+a*log(x)-(*gln))*h;
}
#undef ITMAX
#undef EPS
#undef FPMIN
double gammp(double a, double x)
{
	void gcf(double *gammcf, double a, double x, double *gln);
	void gser(double *gamser, double a, double x, double *gln);
	void nrerror(char error_text[]);
	double gamser,gammcf,gln;
	if (x < 0.0 || a <= 0.0) nrerror("Invalid arguments in routine gammp");
	if (x < (a+1.0)) {
		gser(&gamser,a,x,&gln);
		return gamser;
	} else {
		gcf(&gammcf,a,x,&gln);
		return 1.0-gammcf;
	}
}
#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp;
#define M 7
void sort(int n,double *arr,int *indx, int *rank,
          double *uniq,int *nuniq,int *nclus)
{
	int i,indxt,ir=n,temp,j,k,l=1;
	int jstack=0,*istack;
	double a;
	istack=ivector(1,n);
	for (j=1;j<=n;j++) indx[j]=j;
	for (;;) {
	  if (ir-l < M) {
	    for (j=l+1;j<=ir;j++) {
	      indxt=indx[j];
	      a=arr[indxt];
	      for (i=j-1;i>=l;i--) {
		if (arr[indx[i]] <= a) break;
		indx[i+1]=indx[i];
	      }
	      indx[i+1]=indxt;
	    }
	    if (jstack == 0) break;
	    ir=istack[jstack--];
	    l=istack[jstack--];
	  } else {
	    k=(l+ir) >> 1;
	    SWAP(indx[k],indx[l+1]);
	    if (arr[indx[l]] > arr[indx[ir]]) {
	      SWAP(indx[l],indx[ir])
		}
	    if (arr[indx[l+1]] > arr[indx[ir]]) {
	      SWAP(indx[l+1],indx[ir])
		}
	    if (arr[indx[l]] > arr[indx[l+1]]) {
	      SWAP(indx[l],indx[l+1])
		}
	    i=l+1;
	    j=ir;
	    indxt=indx[l+1];
	    a=arr[indxt];
	    for (;;) {
	      do i++; while (arr[indx[i]] < a);
	      do j--; while (arr[indx[j]] > a);
	      if (j < i) break;
	      SWAP(indx[i],indx[j])
		}
	    indx[l+1]=indx[j];
	    indx[j]=indxt;
	    jstack += 2;
	    if (ir-i+1 >= j-l) {
	      istack[jstack]=ir;
	      istack[jstack-1]=i;
	      ir=j-1;
	    } else {
	      istack[jstack]=j-1;
	      istack[jstack-1]=l;
	      l=i;
	    }
	  }
	}
	for (j=1;j<=n;j++){
	  uniq[j]=0;
	  nuniq[j]=0;
	}
	(*nclus)=1;
	rank[indx[*nclus]]=(*nclus);
	uniq[*nclus]=arr[indx[*nclus]];
	for (j=2;j<=n;j++){
	  if(arr[indx[j]] > arr[indx[j-1]]){
	    (*nclus)++;
	    uniq[*nclus]=arr[indx[j]];
	  }
          rank[indx[j]]=(*nclus);
        }
	for (j=1;j<=n;j++) nuniq[rank[j]]++;
	free_ivector(istack,1,n);
}
#undef M
#undef SWAP
void chold(double **a,double **b,int dim)
{
  double **temp,*diag,sum;
  int i,j,k;
  temp=dmatrix(1,dim,1,dim);
  diag=dvector(1,dim);
  for (i=1;i<=dim;i++){
    for (j=i;j<=dim;j++)
      temp[i][j]=temp[j][i]=a[i][j];
  }
  for (i=1;i<=dim;i++){
    for (j=i;j<=dim;j++){
      for (sum=temp[i][j],k=i-1;k>=1;k--) sum -= temp[i][k]*temp[j][k];
      if (i == j) {
	if (sum <= 0.0) nrerror("choldc failed");
	diag[i]=sqrt(sum);
      } 
      else temp[j][i]=sum/diag[i];
    }
  }
  for (i=1;i<=dim;i++){
    b[i][i]=diag[i];
    for (j=1;j<=i-1;j++){
      b[i][j]=temp[i][j];
      b[j][i]=0.;
    }
  }
  free_dmatrix(temp,1,dim,1,dim);
  free_dvector(diag,1,dim);
}
void createMatrix(double **newX,double *X,int nrow,int ncol)
{
  int i,j;
  X=X-1;
  for (j=1;j<=ncol;j++){  
    for (i=1;i<=nrow;i++){  
      newX[i][j]=X[(j-1)*nrow+i];
    }
  }
}
void createIntMatrix(int **newX,int *X,int nrow,int ncol)
{
  int i,j;
  X=X-1;
  for (j=1;j<=ncol;j++){  
    for (i=1;i<=nrow;i++){  
      newX[i][j]=X[(j-1)*nrow+i];
    }
  }
}
void createVector(double *X,double **newX,int nrow,int ncol)
{
  int i,j;
  X=X-1;
  for (j=1;j<=ncol;j++){  
    for (i=1;i<=nrow;i++){  
      X[(j-1)*nrow+i]=newX[i][j];
    }
  }
}
void lubksb(double **a, int n, int *indx, double b[])
{
	int i,ii=0,ip,j;
	double sum;
	for (i=1;i<=n;i++) {
		ip=indx[i];
		sum=b[ip];
		b[ip]=b[i];
		if (ii)
			for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
		else if (sum) ii=i;
		b[i]=sum;
	}
	for (i=n;i>=1;i--) {
		sum=b[i];
		for (j=i+1;j<=n;j++) sum -= a[i][j]*b[j];
		b[i]=sum/a[i][i];
	}
}
#define TINY 1.0e-20;
void ludcmp(double **a, int n, int *indx, double *d)
{
	int i,imax,j,k;
	double big,dum,sum,temp;
	double *vv;
        imax = 0;  
	vv=dvector(1,n);
	*d=1.0;
	for (i=1;i<=n;i++) {
		big=0.0;
		for (j=1;j<=n;j++)
			if ((temp=fabs(a[i][j])) > big) big=temp;
		if (big == 0.0) nrerror("Singular matrix in routine ludcmp");
		vv[i]=1.0/big;
	}
	for (j=1;j<=n;j++) {
		for (i=1;i<j;i++) {
			sum=a[i][j];
			for (k=1;k<i;k++) sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
		}
		big=0.0;
		for (i=j;i<=n;i++) {
			sum=a[i][j];
			for (k=1;k<j;k++)
				sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
			if ( (dum=vv[i]*fabs(sum)) >= big) {
				big=dum;
				imax=i;
			}
		}
		if (j != imax) {
			for (k=1;k<=n;k++) {
				dum=a[imax][k];
				a[imax][k]=a[j][k];
				a[j][k]=dum;
			}
			*d = -(*d);
			vv[imax]=vv[j];
		}
		indx[j]=imax;
		if (a[j][j] == 0.0) a[j][j]=TINY;
		if (j != n) {
			dum=1.0/(a[j][j]);
			for (i=j+1;i<=n;i++) a[i][j] *= dum;
		}
	}
	free_dvector(vv,1,n);
}
#undef TINY
int matdot(double **a,double *b,double *c,int r1,int c1,int r2)
{
	int i,k,sdim=c1;
	for (i=1;i<=r1;i++){
	  c[i]=0.;
	  for (k=1;k<=sdim;k++){
	    c[i] += a[i][k] * b[k];
	  }
	}
	return (c1 != r2) ? -1 : 0;
}
void matinv(double **a,double **aInv,int N)
{
    int i,j,*indx;
    double **atemp,*col,d;
    indx=ivector(1,N);
    atemp=dmatrix(1,N,1,N);
    col=dvector(1,N);
    for (j=1;j<=N;j++){
      for (i=1;i<=N;i++) atemp[i][j]=a[i][j];
    }
    ludcmp(atemp,N,indx,&d);
    for (j=1;j<=N;j++){
      for (i=1;i<=N;i++) col[i]=0.;
      col[j]=1.0;
      lubksb(atemp,N,indx,col);
      for (i=1;i<=N;i++) aInv[i][j]=col[i];
    }
    free_ivector(indx,1,N);
    free_dmatrix(atemp,1,N,1,N);
    free_dvector(col,1,N);
}
void matinvDet(double **a,double **aInv,double *det,int N)
{
    int i,j,*indx;
    double **atemp,*col,d;
    indx=ivector(1,N);
    atemp=dmatrix(1,N,1,N);
    col=dvector(1,N);
    for (j=1;j<=N;j++){
      for (i=1;i<=N;i++) atemp[i][j]=a[i][j];
    }
    ludcmp(atemp,N,indx,&d);
    (*det)=d;
    for (j=1;j<=N;j++){      
      (*det) *= atemp[j][j];
      for (i=1;i<=N;i++) col[i]=0.;
      col[j]=1.0;
      lubksb(atemp,N,indx,col);
      for (i=1;i<=N;i++) aInv[i][j]=col[i];
    }
    (*det)=1/(*det);
    free_ivector(indx,1,N);
    free_dmatrix(atemp,1,N,1,N);
    free_dvector(col,1,N);
}
double vecdot(double *a,double *b,int d)
{
	int i;
	double c=0.;
	for (i=1;i<=d;i++)
	    c += a[i] * b[i];
	return c;
}
void vecadd(double *a,double *b,double *c,int d)
{
	int i;
	for (i=1;i<=d;i++)
	    c[i] = a[i] + b[i];
}
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)
double dran1(int *idum)
{
	int j;
	int k;
	static int iy=0;
	static int iv[NTAB];
	double temp;
	if (*idum <= 0 || !iy) {
		if (-(*idum) < 1) *idum=1;
		else *idum = -(*idum);
		for (j=NTAB+7;j>=0;j--) {
			k=(*idum)/IQ;
			*idum=IA*(*idum-k*IQ)-IR*k;
			if (*idum < 0) *idum += IM;
			if (j < NTAB) iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ;
	*idum=IA*(*idum-k*IQ)-IR*k;
	if (*idum < 0) *idum += IM;
	j=iy/NDIV;
	iy=iv[j];
	iv[j] = *idum;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}
#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX
double rdisc(double *supp,double *prob,int np,int *idum)
{
      double u=dran1(idum),l=0.,r=np;
      double cum=0.,*cumprob=dvector(1,np);
      int i,j,step;
      for (j=1;j<=np;j++){
	cum += prob[j];       
	cumprob[j]=cum;
      }
      for (j=1;j<=np;j++)  cumprob[j]=cumprob[j]/cum;
      for(;;){
	i=floor((l+r)/2);
	if (u>cumprob[i]){
	  l=i;
	  step=1;
        }
	else{
	  r=i;
          step=0;
        }
        if (l >= r-1){
	  free_dvector(cumprob,1,np);
          return supp[i+step];
	}
      }
}
double rGamma(double a,int *idum)
{
	double e=exp(1.),b=e/(a+e);
	double temp,u0,u1,u2,c1=(a-1.),c2=((a-1./(6.*a))/c1);
        double c3=(2./c1),c4=(2.*a/c1),c5=(1./DSQR(a)),w;
	if (a == 1){
            while((temp=dran1(idum)) != 0.0)         
                 return -log(temp);
	}
        else if (a < 1){
	    for(;;){
	      u0=dran1(idum);u1=dran1(idum);        
	      if ( u0 > b){
		temp=-1.*log((1-u0)/(a*b));
		if (u1 <= pow(temp,a-1.)) return temp;
              }
	      else{
		temp=pow(u0/b,1./a);                
		if (u1 <= exp(-1.*temp)) return temp;
              }
	    }
        }
        else{
	    for(;;){
	      do{
		u1=dran1(idum);u2=dran1(idum);
		if (a > 2.5) u1=u2+c5*(1.-1.86*u1);
	      } while (u1<0. || u1>1.);
	    w=c2*u2/u1;
            if ( (c3*u1+w+1/w <= c4) || 
                   (c3*log(u1)-log(w)+w < 1) ) return c1*w;
	    }
	}
        return 0.0;  
}
#define JMAX 150
#define XMAX 10000
#define SPI  1.772454
double rGammat(int *n3,double *x1,double *xacc,double *tgamma,int *seed)
{
	void nrerror(char error_text[]);
	int i,j;
	double dx,f,fmid,xmid,rtb,bb,u;
	x1=x1-1;tgamma=tgamma-1;
	for (i=1;i<=*n3;i++){
	  bb=gammp(0.5,x1[i])+exp(-1.*x1[i])/(sqrt(x1[i])*SPI);
	  u=dran1(seed);
	  f=-1.*u;
	  fmid=1.-u;
	  rtb = f < 0.0 ? (dx=XMAX-x1[i],x1[i]) : (dx=x1[i]-XMAX,XMAX);
	  for (j=1;j<=JMAX;j++){
	    xmid=rtb+(dx *= 0.5);
	    fmid=(gammp(0.5,xmid)+exp(-1.*xmid)/(sqrt(xmid)*SPI)-bb)/(1-bb)-u;
	    if (fmid <= 0.0) rtb=xmid;
	    if (fabs(dx) < (*xacc) || fmid == 0.0) break;
	  }
	  tgamma[i]=rtb;
	}
        return 0.0;  
}
#undef JMAX
#undef XMAX
#undef SPI
void rlgamma(double *a, double *plog,int *nsample,int *seed)
{
  double bound,x;
  int i;
  a=a-1;plog=plog-1;
  for (i=1;i<=*nsample;i++){
    bound=exp(-(1-a[i])*(pow(a[i],(a[i]/(1.-a[i])))));
    do{
    x=log(-log(1-dran1(seed)))/a[i];
    }while (dran1(seed)>=(exp(exp(x*a[i])-exp(x))*bound));
    plog[i]=x;
  }
}
double rnormal(int *idum)
{
	static int iset=0;
	static double gset;
	double fac,rsq,v1,v2;
	if  (iset == 0) {
		do {
			v1=2.0*dran1(idum)-1.0;
			v2=2.0*dran1(idum)-1.0;
			rsq=v1*v1+v2*v2;
		} while (rsq >= 1.0 || rsq == 0.0);
		fac=sqrt(-2.0*log(rsq)/rsq);
		gset=v1*fac;
		iset=1;
		return v2*fac;
	} else {
		iset=0;
		return gset;
	}
}
void rmnormal(double *mean,double **sigma,double *nvec,
	      int dim,int *idum)
{
  double **var,*zz;
  int i,j;
  var=dmatrix(1,dim,1,dim);
  zz=dvector(1,dim);
  for (i=1;i<=dim;i++){
    for (j=i;j<=dim;j++)
      var[i][j]=var[j][i]=sigma[i][j];
  }
  chold(var,var,dim);
  for (i=1;i<=dim;i++)
    zz[i]=rnormal(idum);
  matdot(var,zz,nvec,dim,dim,dim);
  vecadd(mean,nvec,nvec,dim);
  free_dmatrix(var,1,dim,1,dim);
  free_dvector(zz,1,dim);
}
void spikeSlabVar(double *beta,
               double *W,double *Waugment, 
               double *V,double *Vaugment, 
               int *prior,int *ncov,
               double *Wscl,double *Wshp,double *Vsmall,
               double *gg, double *ggaugment, int *seed)
{
  double scale,shape,bb,pp,pp1,pp2,Wk;
  int k;
  beta=beta-1;W=W-1;Waugment=Waugment-1;V=V-1;
  Vaugment=Vaugment-1;ggaugment=ggaugment-1;
  if ((*prior)==0) {
    for (k=1;k<=*ncov;k++) {
      bb=DSQR(beta[k])/2.;
      scale=(*Wscl)+bb/(V[k]*Vaugment[k]);
      shape=(*Wshp)+1/2.;
      W[k]=scale/rGamma(shape,seed);
      Wk=W[k]*Vaugment[k];
      pp1=(1-(*gg))*exp(-bb/((*Vsmall)*Wk))/sqrt((*Vsmall));
      pp2=(*gg)*exp(-bb/Wk);
      pp=pp1/(pp1+pp2);
      if (dran1(seed)<pp) {
	V[k]=(*Vsmall);
      }
      else {
	V[k]=1.;
      }
      Wk=W[k]*V[k];
      pp1=(1-(ggaugment[k]))*exp(-bb/((*Vsmall)*Wk))/sqrt((*Vsmall));
      pp2=(ggaugment[k])*exp(-bb/Wk);
      pp=pp1/(pp1+pp2);
      if (dran1(seed)<pp) {
	Vaugment[k]=(*Vsmall);
      }
      else {
	Vaugment[k]=1.;
      }
    }	
  }
  if ((*prior)==1) {
    for (k=1;k<=*ncov;k++) {
      bb=DSQR(beta[k])/2.;
      scale=(*Wscl)+bb/V[k];
      shape=(*Wshp)+1/2.;
      W[k]=scale/rGamma(shape,seed);
      pp1=(1-(*gg))*exp(-bb/((*Vsmall)*W[k]))/sqrt((*Vsmall));
      pp2=(*gg)*exp(-bb/W[k]);
      pp=pp1/(pp1+pp2);
      if (dran1(seed)<pp) {
	V[k]=(*Vsmall);
      }
      else {
	V[k]=1.;
      }
    }	
  }
  else if ((*prior)==2) {
    for (k=1;k<=*ncov;k++) {
      bb=DSQR(beta[k])/2.;
      if (V[k] == (*Vsmall)) {
        W[k]=(*Wscl)/rGamma((*Wshp),seed);
	pp1=(1-(*gg))*exp(-bb/Waugment[k])/sqrt(Waugment[k]);
	pp2=(*gg)*exp(-bb/W[k])/sqrt(W[k]);
        pp=pp1/(pp1+pp2);
	W[k]=Waugment[k];
      }
      else {
        scale=(*Wscl)+bb;
        shape=(*Wshp)+1/2.;
        W[k] = scale/rGamma(shape,seed);
	pp1=(1-(*gg))*exp(-bb/Waugment[k])/sqrt(Waugment[k]);
	pp2=(*gg)*exp(-bb/W[k])/sqrt(W[k]);
        pp=pp1/(pp1+pp2);
      }
      if (dran1(seed)<pp) {
	V[k]=(*Vsmall);
      }
      else {
	V[k]=1.;
      }
    }	
  }
  else {
    for (k=1;k<=*ncov;k++) {
      bb=DSQR(beta[k])/2.;
      pp1=(1-(*gg))*exp(-bb/((*Vsmall)*W[k]))/sqrt((*Vsmall));
      pp2=(*gg)*exp(-bb/W[k]);
      pp=pp1/(pp1+pp2);
      if (dran1(seed)<pp) {
	V[k]=(*Vsmall);
      }
      else {
	V[k]=1.;
      }
    }	
  }
}
double **convert_dmatrix(double *a,long nrl,long nrh,long ncl,long nch)
{
	long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1;
	double **m;
	m=(double **) malloc((size_t) ((nrow+NR_END)*sizeof(double*)));
	if (!m) nrerror("allocation failure in convert_dmatrix()");
	m += NR_END;
	m -= nrl;
	m[nrl]=a-ncl;
	for(i=1,j=nrl+1;i<nrow;i++,j++) m[j]=m[j-1]+ncol;
	return m;
}
float **convert_matrix(float *a,long nrl,long nrh,long ncl,long nch)
{
	long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1;
	float **m;
	m=(float **) malloc((size_t) ((nrow+NR_END)*sizeof(float*)));
	if (!m) nrerror("allocation failure in convert_matrix()");
	m += NR_END;
	m -= nrl;
	m[nrl]=a-ncl;
	for(i=1,j=nrl+1;i<nrow;i++,j++) m[j]=m[j-1]+ncol;
	return m;
}
unsigned char *cvector(long nl,long nh)
{
	unsigned char *v;
	v=(unsigned char *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(unsigned char)));
	if (!v) nrerror("allocation failure in cvector()");
	return v-nl+NR_END;
}
double ***d3tensor(long nrl,long nrh,long ncl,long nch,long ndl,long ndh)
{
	long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
	double ***t;
	t=(double ***) malloc((size_t)((nrow+NR_END)*sizeof(double**)));
	if (!t) nrerror("allocation failure 1 in d3tensor()");
	t += NR_END;
	t -= nrl;
	t[nrl]=(double **) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double*)));
	if (!t[nrl]) nrerror("allocation failure 2 in d3tensor()");
	t[nrl] += NR_END;
	t[nrl] -= ncl;
	t[nrl][ncl]=(double *) malloc((size_t)((nrow*ncol*ndep+NR_END)*sizeof(double)));
	if (!t[nrl][ncl]) nrerror("allocation failure 3 in d3tensor()");
	t[nrl][ncl] += NR_END;
	t[nrl][ncl] -= ndl;
	for(j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
	for(i=nrl+1;i<=nrh;i++) {
		t[i]=t[i-1]+ncol;
		t[i][ncl]=t[i-1][ncl]+ncol*ndep;
		for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
	}
	return t;
}
double **dmatrix(long nrl,long nrh,long ncl,long nch)
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	double **m;
	m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double*)));
	if (!m) nrerror("allocation failure 1 in dmatrix()");
	m += NR_END;
	m -= nrl;
	m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
	if (!m[nrl]) nrerror("allocation failure 2 in dmatrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;
	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
	return m;
}
double **dsubmatrix(double **a,long oldrl,long oldrh,long oldcl,long oldch,
	long newrl,long newcl)
{
	long i,j,nrow=oldrh-oldrl+1,ncol=oldcl-newcl;
	double **m;
	m=(double **) malloc((size_t) ((nrow+NR_END)*sizeof(double*)));
	if (!m) nrerror("allocation failure in dsubmatrix()");
	m += NR_END;
	m -= newrl;
	for(i=oldrl,j=newrl;i<=oldrh;i++,j++) m[j]=a[i]+ncol;
	return m;
}
double *dvector(long nl,long nh)
{
	double *v;
	v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
	if (!v) nrerror("allocation failure in dvector()");
	return v-nl+NR_END;
}
float ***f3tensor(long nrl,long nrh,long ncl,long nch,long ndl,long ndh)
{
	long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
	float ***t;
	t=(float ***) malloc((size_t)((nrow+NR_END)*sizeof(float**)));
	if (!t) nrerror("allocation failure 1 in f3tensor()");
	t += NR_END;
	t -= nrl;
	t[nrl]=(float **) malloc((size_t)((nrow*ncol+NR_END)*sizeof(float*)));
	if (!t[nrl]) nrerror("allocation failure 2 in f3tensor()");
	t[nrl] += NR_END;
	t[nrl] -= ncl;
	t[nrl][ncl]=(float *) malloc((size_t)((nrow*ncol*ndep+NR_END)*sizeof(float)));
	if (!t[nrl][ncl]) nrerror("allocation failure 3 in f3tensor()");
	t[nrl][ncl] += NR_END;
	t[nrl][ncl] -= ndl;
	for(j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
	for(i=nrl+1;i<=nrh;i++) {
		t[i]=t[i-1]+ncol;
		t[i][ncl]=t[i-1][ncl]+ncol*ndep;
		for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
	}
	return t;
}
int **imatrix(long nrl,long nrh,long ncl,long nch)
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	int **m;
	m=(int **) malloc((size_t)((nrow+NR_END)*sizeof(int*)));
	if (!m) nrerror("allocation failure 1 in imatrix()");
	m += NR_END;
	m -= nrl;
	m[nrl]=(int *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(int)));
	if (!m[nrl]) nrerror("allocation failure 2 in imatrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;
	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
	return m;
}
int *ivector(long nl,long nh)
{
	int *v;
	v=(int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int)));
	if (!v) nrerror("allocation failure in ivector()");
	return v-nl+NR_END;
}
unsigned long *lvector(long nl,long nh)
{
	unsigned long *v;
	v=(unsigned long *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(long)));
	if (!v) nrerror("allocation failure in lvector()");
	return v-nl+NR_END;
}
float **matrix(long nrl,long nrh,long ncl,long nch)
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	float **m;
	m=(float **) malloc((size_t)((nrow+NR_END)*sizeof(float*)));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;
	m[nrl]=(float *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(float)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;
	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
	return m;
}
void nrerror(char error_text[])
{
	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}
float **submatrix(float **a,long oldrl,long oldrh,long oldcl,long oldch,
	long newrl,long newcl)
{
	long i,j,nrow=oldrh-oldrl+1,ncol=oldcl-newcl;
	float **m;
	m=(float **) malloc((size_t) ((nrow+NR_END)*sizeof(float*)));
	if (!m) nrerror("allocation failure in submatrix()");
	m += NR_END;
	m -= newrl;
	for(i=oldrl,j=newrl;i<=oldrh;i++,j++) m[j]=a[i]+ncol;
	return m;
}
float *vector(long nl,long nh)
{
	float *v;
	v=(float *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(float)));
	if (!v) nrerror("allocation failure in vector()");
	return v-nl+NR_END;
}
void free_convert_dmatrix(double **b,long nrl,long nrh,long ncl,long nch)
{
	free((FREE_ARG) (b+nrl-NR_END));
}
void free_convert_matrix(float **b,long nrl,long nrh,long ncl,long nch)
{
	free((FREE_ARG) (b+nrl-NR_END));
}
void free_cvector(unsigned char *v,long nl,long nh)
{
	free((FREE_ARG) (v+nl-NR_END));
}
void free_d3tensor(double ***t,long nrl,long nrh,long ncl,long nch,
	long ndl,long ndh)
{
	free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END));
	free((FREE_ARG) (t[nrl]+ncl-NR_END));
	free((FREE_ARG) (t+nrl-NR_END));
}
void free_dmatrix(double **m,long nrl,long nrh,long ncl,long nch)
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}
void free_dsubmatrix(double **b,long nrl,long nrh,long ncl,long nch)
{
	free((FREE_ARG) (b+nrl-NR_END));
}
void free_dvector(double *v,long nl,long nh)
{
	free((FREE_ARG) (v+nl-NR_END));
}
void free_f3tensor(float ***t,long nrl,long nrh,long ncl,long nch,
	long ndl,long ndh)
{
	free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END));
	free((FREE_ARG) (t[nrl]+ncl-NR_END));
	free((FREE_ARG) (t+nrl-NR_END));
}
void free_imatrix(int **m,long nrl,long nrh,long ncl,long nch)
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}
void free_ivector(int *v,long nl,long nh)
{
	free((FREE_ARG) (v+nl-NR_END));
}
void free_lvector(unsigned long *v,long nl,long nh)
{
	free((FREE_ARG) (v+nl-NR_END));
}
void free_matrix(float **m,long nrl,long nrh,long ncl,long nch)
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}
void free_submatrix(float **b,long nrl,long nrh,long ncl,long nch)
{
	free((FREE_ARG) (b+nrl-NR_END));
}
void free_vector(float *v,long nl,long nh)
{
	free((FREE_ARG) (v+nl-NR_END));
}
