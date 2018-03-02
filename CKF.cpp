#include <vector>
#include <iostream>
#include <ctime>
#include <map>
#include <set>
#include <fstream>
#define _USE_MATH_DEFINES
#include <cmath>
#include <iomanip>
#include <string>
using namespace std;

#include "headerfile.h"

#ifndef M_PI
#define M_PI 3.1415926535897932384626
#endif


vector< vector<double> > MatrixInverse(const vector< vector<double> > & A, int k=0);
vector<double> Matrix_mul(const vector< vector<double> > & A, const vector<double> & X);
double vector_norm1(const vector<double> & v1, const vector<double> & v2);
double function_X(double x, double r1, double r2, double sigma);

extern string dataset_name;
extern string experiment_name;

class CKF{
public:
	CKF(){;}
	CKF(int _usern, int _itemn):usern(_usern),itemn(_itemn){
		K = 10;
		T0 = 86400;
		sigma=1.0;
		alpha_user=0;
		alpha_item=0;
		//select hyper-parameters in {0.5, 0.6, 0.7, ... , 2.5} to get best performance
		if(dataset_name=="ML100k" && experiment_name=="classical"){
			sigma=0.9;
		}
		if(dataset_name=="ML100k" && experiment_name=="time_order"){
			sigma=0.9;
		}
		if(dataset_name=="ML100k" && experiment_name=="time_order_online"){
			sigma=1.0;
		}
		if(dataset_name=="Epinions" && experiment_name=="classical"){
			sigma=2.1;
		}
		if(dataset_name=="Epinions" && experiment_name=="time_order"){
			sigma=1.1;
		}
		if(dataset_name=="Epinions" && experiment_name=="time_order_online"){
			sigma=1.4;
		}
		//
		pu=vector< vector<double> >(usern,vector<double>(K));
		pi=vector< vector<double> >(itemn,vector<double>(K));
		cov_user=vector< vector< vector<double> > >(usern,identityMatrix());
		cov_item=vector< vector< vector<double> > >(itemn,identityMatrix());
		t_user=vector<double>(usern,-1);
		t_item=vector<double>(itemn,-1);
		for(int m=0;m!=usern;++m)
			initialize(pu[m]);
		for(int m=0;m!=itemn;++m)
			initialize(pi[m]);
	}
	void initialize(vector<double> & p){
		for(int k=0;k!=K;++k)
			p[k]=random(-0.01,0.01);
	}
	vector< vector<double> > identityMatrix(){
		vector< vector<double> > cov(K,vector<double>(K,0));
		for(int k=0;k!=K;++k)
			cov[k][k]=1;
		return cov;
	}
	double calculate(int userid, int itemid){
		double r=0;
		for(int k=0;k!=K;++k)
			r+=pu[userid][k]*pi[itemid][k];
		return r;
	}
	void update(int userid, int itemid, double rating, double time){
		time /= T0;
		if(t_user[userid]!=-1){
			for(int k=0;k!=K;++k)
				cov_user[userid][k][k]+=(time-t_user[userid])*alpha_user;
		}
		if(t_item[itemid]!=-1){
			for(int k=0;k!=K;++k)
				cov_item[itemid][k][k]+=(time-t_item[itemid])*alpha_item;
		}
		vector< vector<double> > cov_user_1=MatrixInverse(cov_user[userid]);
		vector< vector<double> > cov_item_1=MatrixInverse(cov_item[itemid]);
		vector<double> cov_p_u=Matrix_mul(cov_user_1,pu[userid]),cov_p_i=Matrix_mul(cov_item_1,pi[itemid]);
		vector< vector<double> > & Cov_user=cov_user[userid], & Cov_item=cov_item[itemid];
		vector<double> & p_user=pu[userid], & p_item=pi[itemid];
		vector<double> P_user,P_item;
		double sigma2=sigma*sigma;
		double y=rating;
		for(int n=0;n!=100;++n){
			P_user=p_user;
			P_item=p_item;

			y=0;
			for(int k=0;k!=K;++k)
				y+=p_user[k]*p_item[k];
			y=function_X(y,rating-sigma/2,rating+sigma/2,sigma);

			Cov_user=cov_user_1;
			for(int k1=0;k1!=K;++k1)for(int k2=0;k2!=K;++k2)
				Cov_user[k1][k2]+=(p_item[k1]*p_item[k2]+Cov_item[k1][k2])/sigma2;
			Cov_user=MatrixInverse(Cov_user);
			p_user=cov_p_u;
			for(int k=0;k!=K;++k)
				p_user[k]+=y*p_item[k]/sigma2;
			p_user=Matrix_mul(Cov_user,p_user);

			Cov_item=cov_item_1;
			for(int k1=0;k1!=K;++k1)for(int k2=0;k2!=K;++k2)
				Cov_item[k1][k2]+=(p_user[k1]*p_user[k2]+Cov_user[k1][k2])/sigma2;
			Cov_item=MatrixInverse(Cov_item);
			p_item=cov_p_i;
			for(int k=0;k!=K;++k)
				p_item[k]+=y*p_user[k]/sigma2;
			p_item=Matrix_mul(Cov_item,p_item);

			if(vector_norm1(p_user,P_user)<0.00001 && vector_norm1(p_item,P_item)<0.00001)
				break;
		}
		t_user[userid]=time;
		t_item[itemid]=time;
	}
	int usern,itemn,K;
	double alpha_user,alpha_item,sigma,T0;
	vector< vector<double> > pu,pi;
	vector< vector< vector<double> > > cov_user,cov_item;
	vector<double> t_user,t_item;
};



void CKF_mainfunction(vector<Rating> & ratings, int usern, int itemn){

	CKF model1(usern,itemn);
	double rmse=0;
	double rmse1=0,num1=0,num_all=0;
	set<int> appeared_users, appeared_items;
	for(int m=0;m!=ratings.size();++m){
		if((ratings.size()/100)>0 && (m+1)%(ratings.size()/100)==0){
			cout<<".";
		}
		if((ratings.size()/10)>0 && (m+1)%(ratings.size()/10)==0){
			cout<<" "<<(m+1)/(ratings.size()/10)*10<<"%";
			cout<<" Current piece RMSE: "<<sqrt(rmse1/num1);
			cout<<endl;
			rmse1=num1=0;
		}
		double r;
		if(appeared_users.count(ratings[m].userid) && appeared_items.count(ratings[m].itemid)){
			r=model1.calculate(ratings[m].userid,ratings[m].itemid);
			rmse+=pow(ratings[m].rating-r,2);
			rmse1+=pow(ratings[m].rating-r,2);
			num1+=1;
			num_all+=1;
		}
		appeared_users.insert(ratings[m].userid);
		appeared_items.insert(ratings[m].itemid);
		model1.update(ratings[m].userid,ratings[m].itemid,ratings[m].rating,ratings[m].time);
	}
	cout << "Overall RMSE: " << endl;
	cout<<sqrt(rmse/num_all)<<endl;


}


void CKF_mainfunction_batch(vector<Rating> & rtraining, vector<Rating> & rtest, int usern, int itemn){

	CKF model1(usern,itemn);
	set<int> appeared_users, appeared_items;
	cout << "Training" << endl;
	for(int m=0;m!=rtraining.size();++m){
		if((rtraining.size()/100)>0 && (m+1)%(rtraining.size()/100)==0){
			cout<<".";
		}
		if((rtraining.size()/10)>0 && (m+1)%(rtraining.size()/10)==0){
			cout<<" "<<(m+1)/(rtraining.size()/10)*10<<"%"<<endl;
		}
		appeared_users.insert(rtraining[m].userid);
		appeared_items.insert(rtraining[m].itemid);
		model1.update(rtraining[m].userid,rtraining[m].itemid,rtraining[m].rating,rtraining[m].time);
	}
	cout << "Test" << endl;
	double rmse=0,num_all=0;
	for(int m=0;m!=rtest.size();++m){
		if(!appeared_users.count(rtest[m].userid) || !appeared_items.count(rtest[m].itemid))
			continue;
		double r=model1.calculate(rtest[m].userid,rtest[m].itemid);
		rmse+=pow(rtest[m].rating-r,2);
		num_all+=1;
	}
	
	cout<<num_all<<" ";
	cout << "RMSE: " ;
	cout<<sqrt(rmse/num_all)<<endl;

}



double function_phi(double x){
	return 1.0/sqrt(2*M_PI)*exp(-x*x/2);
}

double function_Phi(double x1, double x2){
	double y=0;
	int M=500;
	for(int m=0;m!=M;++m){
		y+=function_phi(x1+(x2-x1)*(m+0.5)/M);
	}
	return y*(x2-x1)/M;
}

double function_X(double x, double r1, double r2, double sigma){
	double e=(r1-x)/sigma;
	double f=(r2-x)/sigma;
	double X=x+sigma*(function_phi(e)-function_phi(f))/(function_Phi(e,f));
	return X;
}

vector<double> Matrix_mul(const vector< vector<double> > & A, const vector<double> & X){
	int N=A.size();
	vector<double> Y(N,0);
	for(int n=0;n!=N;++n)
		for(int k=0;k!=N;++k)
			Y[n]+=A[n][k]*X[k];
	return Y;
}

double vector_norm1(const vector<double> & v1, const vector<double> & v2){
	int K=v1.size();
	if(v2.size()!=K){
		return 0;
	}
	double d=0;
	for(int k=0;k!=K;++k)
		d+=fabs(v1[k]-v2[k]);
	return d/K;
}

void Matrix_mul(double * A, double * B, double * C, int N){
	for(int m=0;m!=N;++m){
		for(int n=0;n!=N;++n){
			C[m*N+n]=0;
			for(int k=0;k!=N;++k)
				C[m*N+n]+=A[m*N+k]*B[k*N+n];
		}
	}
}

int dluav(double a[],int m,int n,double u[],double v[],double eps,int ka);

vector< vector<double> > MatrixInverse(const vector< vector<double> > & A, int k){
	int N=A.size();
	if (A.empty() || A[0].size() != A.size()) {
		return vector< vector<double> >();
	}
	double * a=(double *)malloc(sizeof(double)*N*N);
	double * u=(double *)malloc(sizeof(double)*N*N);
	double * v=(double *)malloc(sizeof(double)*N*N);
	double eps=1e-8;
	for(int m=0;m!=N;++m)
		for(int n=0;n!=N;++n)
			a[m*N+n]=A[m][n];
	int r=dluav(a,N,N,u,v,eps,N+1);
	for(int m=0;m!=N;++m){
		for(int n=0;n!=m;++n){
			swap(u[m*N+n],u[n*N+m]);
			swap(v[m*N+n],v[n*N+m]);
		}
	}
	for(int m=0;m!=N;++m){
		if (fabs(a[m*N + m]) > 0)
			a[m*N + m] = 1.0 / a[m*N + m];
		else
			a[m*N + m] = 0;
	}
	for(int m=0;m!=N;++m)
		for(int n=0;n!=N;++n)
			v[m*N+n]=v[m*N+n]*a[n*N+n];
	Matrix_mul(v,u,a,N);
	vector< vector<double> > A1=A;
	for(int m=0;m!=N;++m)
		for(int n=0;n!=N;++n)
			A1[m][n]=a[m*N+n];
	free(a);free(u);free(v);
	return A1;
}

//The following code is from:
//Shiliang Xu and Erni Ma. Programs for Frequently-used Algorithms. Tsinghua University Press. 2013.

#define MAX_ITERA 60
#define MIN_DOUBLE (1e-30)

int dluav(double a[],int m,int n,double u[],double v[],double eps,int ka);
static void damul(double a[],double b[],int m,int n,int k,double c[]);
static void ppp(double a[],double e[],double s[],double v[],int m,int n);
static void sss(double fg[2],double cs[2]);

int dluav(double a[],int m,int n,double u[],double v[],double eps,int ka)
{
	int i,j,k,l,it,ll,kk,ix,iy,mm,nn,iz,ml,ks;
	double d,dd,t,sm,sml,eml,sk,ek,b,c,shh,fg[2],cs[2];
	double *s,*e,*w;
	s=(double*)malloc(ka*sizeof(double));
	e=(double*)malloc(ka*sizeof(double));
	w=(double*)malloc(ka*sizeof(double));
	for(i=1;i<=m;i++)
	{
		ix=(i-1)*m+i-1;
		u[ix]=0;
	}
	for(i=1;i<=n;i++)
	{
		iy=(i-1)*n+i-1;
		v[iy]=0;
	}
	it=MAX_ITERA;k=n;
	if(m-1<n)
		k=m-1;
	l=m;
	if(n-2<m) l=n-2;
	if(l<0) l=0;
	ll=k;
	if(l>k) ll=l;
	if(ll>=1)
	{
		for(kk=1;kk<=ll;kk++)
		{
			if(kk<=k)
			{
				d=0.0;
				for(i=kk;i<=m;i++)
				{
					ix=(i-1)*n+kk-1;d=d+a[ix]*a[ix];
				}
				s[kk-1]=sqrt(d);
				//if(s[kk-1]!=0.0)
				if(fabs(s[kk-1])>MIN_DOUBLE)
				{
					ix=(kk-1)*n+kk-1;
					//if(a[ix]!=0.0)
					if(fabs(a[ix])>MIN_DOUBLE)
					{
						s[kk-1]=fabs(s[kk-1]);
						if(a[ix]<0.0) s[kk-1]=-s[kk-1];
					}
					for(i=kk;i<=m;i++)
					{
						iy=(i-1)*n+kk-1;
						a[iy]=a[iy]/s[kk-1];
					}
					a[ix]=1.0+a[ix];
				}
				s[kk-1]=-s[kk-1];
			}
			if(n>=kk+1)
			{
				for(j=kk+1;j<=n;j++)
				{
					//if((kk<=k)&&(s[kk-1]!=0.0))
					if((kk<=k)&&(fabs(s[kk-1])>MIN_DOUBLE))
					{
						d=0.0;
						for(i=kk;i<=m;i++)
						{
							ix=(i-1)*n+kk-1;
							iy=(i-1)*n+j-1;
							d=d+a[ix]*a[iy];
						}
						d=-d/a[(kk-1)*n+kk-1];
						for(i=kk;i<=m;i++)
						{
							ix=(i-1)*n+j-1;
							iy=(i-1)*n+kk-1;
							a[ix]=a[ix]+d*a[iy];
						}
					}
					e[j-1]=a[(kk-1)*n+j-1];
				}
			}
			if(kk<=k)
			{
				for(i=kk;i<=m;i++)
				{
					ix=(i-1)*m+kk-1;iy=(i-1)*n+kk-1;
					u[ix]=a[iy];
				}
			}
			if(kk<=l)
			{
				d=0.0;
				for(i=kk+1;i<=n;i++)
					d=d+e[i-1]*e[i-1];
				e[kk-1]=sqrt(d);
				//if(e[kk-1]!=0.0)
				if(fabs(e[kk-1])>MIN_DOUBLE)
				{
					//if(e[kk]!=0.0)
					if(fabs(e[kk])>MIN_DOUBLE)
					{
						e[kk-1]=fabs(e[kk-1]);
						if(e[kk]<0.0)
							e[kk-1]=-e[kk-1];
					}
					for(i=kk+1;i<=n;i++)
						e[i-1]=e[i-1]/e[kk-1];
					e[kk]=1.0+e[kk];
				}
				e[kk-1]=-e[kk-1];
				//if((kk+1<=m)&&(e[kk-1]!=0.0))
				if((kk+1<=m)&&(fabs(e[kk-1])>MIN_DOUBLE))
				{
					for(i=kk+1;i<=m;i++) w[i-1]=0.0;
					for(j=kk+1;j<=n;j++)
						for(i=kk+1;i<=m;i++)
							w[i-1]=w[i-1]+e[j-1]*a[(i-1)*n+j-1];
					for(j=kk+1;j<=n;j++)
						for(i=kk+1;i<=m;i++)
						{
							ix=(i-1)*n+j-1;
							a[ix]=a[ix]-w[i-1]*e[j-1]/e[kk];
						}
				}
				for(i=kk+1;i<=n;i++)
					v[(i-1)*n+kk-1]=e[i-1];
			}
		}
	}
	mm=n;
	if(m+1<n) mm=m+1;
	if(k<n) s[k]=a[k*n+k];
	if(m<mm) s[mm-1]=0.0;
	if(l+1<mm) e[l]=a[l*n+mm-1];
	e[mm-1]=0.0;
	nn=m;
	if(m>n) nn=n;
	if(nn>=k+1)
	{
		for(j=k+1;j<=nn;j++)
		{
			for(i=1;i<=m;i++)
				u[(i-1)*m+j-1]=0.0;
			u[(j-1)*m+j-1]=1.0;
		}
	}
	if(k>=1)/////////////////////////////////
	{
		for(ll=1;ll<=k;ll++)
		{
			kk=k-ll+1;iz=(kk-1)*m+kk-1;
			//if(s[kk-1]!=0.0)
			if(fabs(s[kk-1])>MIN_DOUBLE)
			{
				if(nn>=kk+1)
					for(j=kk+1;j<=nn;j++)
					{
						d=0.0;
						for(i=kk;i<=m;i++)
						{
							ix=(i-1)*m+kk-1;
							iy=(i-1)*m+j-1;
							d=d+u[ix]*u[iy]/u[iz];
						}
						d=-d;
						for(i=kk;i<=m;i++)
						{
							ix=(i-1)*m+j-1;
							iy=(i-1)*m+kk-1;
							u[ix]=u[ix]+d*u[iy];
						}
					}
					for(i=kk;i<=m;i++)
					{
						ix=(i-1)*m+kk-1;
						u[ix]=-u[ix];
					}
					u[iz]=1.0+u[iz];
					if(kk-1>=1)//////////////////////////////////////
						for(i=1;i<=kk-1;i++)
							u[(i-1)*m+kk-1]=0.0;
			}
			else
			{
				for(i=1;i<=m;i++)
					u[(i-1)*m+kk-1]=0.0;
				u[(kk-1)*m+kk-1]=1.0;
			}
		}
	}
	for(ll=1;ll<=n;ll++)
	{
		kk=n-ll+1;iz=kk*n+kk-1;
		//if((kk<=l)&&(e[kk-1]!=0.0))/////////////////////////////
		if((kk<=l)&&(fabs(e[kk-1])>MIN_DOUBLE))
		{
			for(j=kk+1;j<=n;j++)
			{
				d=0.0;
				for(i=kk+1;i<=n;i++)
				{
				ix=(i-1)*n+kk-1;iy=(i-1)*n+j-1;
				d=d+v[ix]*v[iy]/v[iz];
				}
				d=-d;
				for(i=kk+1;i<=n;i++)
				{
					ix=(i-1)*n+j-1;iy=(i-1)*n+kk-1;
					v[ix]=v[ix]+d*v[iy];
				}
			}
		}
		for(i=1;i<=n;i++)
			v[(i-1)*n+kk-1]=0.0;
		v[iz-n]=1.0;
	}
	for(i=1;i<=m;i++)
		for(j=1;j<=n;j++)
			a[(i-1)*n+j-1]=0.0;
	ml=mm;
	it=MAX_ITERA;
	while(1==1)//////////////////////////////////
	{
		if(mm==0)
		{
			ppp(a,e,s,v,m,n);
			free(s);free(e);free(w);
			return l;
		}
		if(it==0)
		{
			ppp(a,e,s,v,m,n);
			free(s);free(e);free(w);
			return -1;
		}
		kk=mm-1;
		//while((kk!=0)&&(fabs(e[kk-1])!=0.0))
		while((kk!=0)&&(fabs(e[kk-1])>MIN_DOUBLE))
		{
			d=fabs(s[kk-1])+fabs(s[kk]);
			dd=fabs(e[kk-1]);
			if(dd>eps*d)
				kk=kk-1;
			else
				e[kk-1]=0.0;
		}
		if(kk==mm-1)
		{
			kk=kk+1;
			if(s[kk-1]<0.0)
			{
				s[kk-1]=-s[kk-1];
				for(i=1;i<=n;i++)
				{
					ix=(i-1)*n+kk-1;
					v[ix]=-v[ix];
				}
			}
			while((kk!=ml)&&(s[kk-1]<s[kk]))
			{
				d=s[kk-1];s[kk-1]=s[kk];s[kk]=d;
				if(kk<n)
					for(i=1;i<=n;i++)
					{
						ix=(i-1)*n+kk-1;iy=(i-1)*n+kk;
						d=v[ix];v[ix]=v[iy];v[iy]=d;
					}
					if(kk<m)
						for(i=1;i<=m;i++)
						{
							ix=(i-1)*m+kk-1;
							iy=(i-1)*m+kk;
							d=u[ix];u[ix]=u[iy];u[iy]=d;
						}
						kk=kk+1;
			}
			it=MAX_ITERA;
			mm=mm-1;
		}
		else
		{
			ks=mm;
			//while((ks>kk)&&(fabs(s[ks-1])!=0.0))
			while((ks>kk)&&(fabs(s[ks-1])>MIN_DOUBLE))
			{
				d=0.0;
				if(ks!=mm)
					d=d+fabs(e[ks-1]);
				if(ks!=kk+1) d=d+fabs(e[ks-2]);
				dd=fabs(s[ks-1]);
				if(dd>eps*d)
					ks=ks-1;
				else
					s[ks-1]=0.0;
			}
			if(ks==kk)
			{
				kk=kk+1;
				d=fabs(s[mm-1]);
				t=fabs(s[mm-2]);
				if(t>d)
					d=t;
				t=fabs(e[mm-2]);
				if(t>d)
					d=t;
				t=fabs(s[kk-1]);
				if(t>d)
					d=t;
				t=fabs(e[kk-1]);
				if(t>d)
					d=t;
				sm=s[mm-1]/d;sml=s[mm-2]/d;
				eml=e[mm-2]/d;
				sk=s[kk-1]/d;ek=e[kk-1]/d;
				b=((sml+sm)*(sml-sm)+eml*eml)/2.0;
				c=sm*eml;c=c*c;shh=0.0;
				//if((b!=0.0)||(c!=0.0))
				if((fabs(b)>MIN_DOUBLE)||(fabs(c)>MIN_DOUBLE))
				{
					shh=sqrt(b*b+c);
					if(b<0.0)
						shh=-shh;
					shh=c/(b+shh);
				}
				fg[0]=(sk+sm)*(sk-sm)-shh;
				fg[1]=sk*ek;
				for(i=kk;i<=mm-1;i++)
				{
					sss(fg,cs);
					if(i!=kk)
						e[i-2]=fg[0];
					fg[0]=cs[0]*s[i-1]+cs[1]*e[i-1];
					e[i-1]=cs[0]*e[i-1]-cs[1]*s[i-1];
					fg[1]=cs[1]*s[i];
					s[i]=cs[0]*s[i];
					//if((cs[0]!=1.0)||(cs[1]!=0.0))
					if((fabs(cs[0]-1.0)>MIN_DOUBLE)||(fabs(cs[1])>MIN_DOUBLE))
						for(j=1;j<=n;j++)
						{
							ix=(j-1)*n+i-1;
							iy=(j-1)*n+i;
							d=cs[0]*v[ix]+cs[1]*v[iy];
							v[iy]=-cs[1]*v[ix]+cs[0]*v[iy];
							v[ix]=d;
						}
					sss(fg,cs);
					s[i-1]=fg[0];
					fg[0]=cs[0]*e[i-1]+cs[1]*s[i];
					s[i]=-cs[1]*e[i-1]+cs[0]*s[i];
					fg[1]=cs[1]*e[i];
					e[i]=cs[0]*e[i];
					if(i<m)
						//if((cs[0]!=1.0)||(cs[1]!=0.0))
						if((fabs(cs[0]-1.0)>MIN_DOUBLE)||(fabs(cs[1])>MIN_DOUBLE))
							for(j=1;j<=m;j++)
							{
								ix=(j-1)*m+i-1;
								iy=(j-1)*m+i;
								d=cs[0]*u[ix]+cs[1]*u[iy];
								u[iy]=-cs[1]*u[ix]+cs[0]*u[iy];
								u[ix]=d;
							}
				}
				e[mm-2]=fg[0];
				it=it-1;
			}
			else
			{
				if(ks==mm)
				{
					kk=kk+1;
					fg[1]=e[mm-2];e[mm-2]=0.0;
					for(ll=kk;ll<=mm-1;ll++)
					{
						i=mm+kk-ll-1;
						fg[0]=s[i-1];
						sss(fg,cs);
						s[i-1]=fg[0];
						if(i!=kk)
						{
							fg[1]=-cs[1]*e[i-2];
							e[i-2]=cs[0]*e[i-2];
						}
						//if((cs[0]!=1.0)||(cs[1]!=0.0))
						if((fabs(cs[0]-1.0)>MIN_DOUBLE)||(fabs(cs[1])>MIN_DOUBLE))
							for(j=1;j<=n;j++)
							{
								ix=(j-1)*n+i-1;
								iy=(j-1)*n+mm-1;
								d=cs[0]*v[ix]+cs[1]*v[iy];
								v[iy]=-cs[1]*v[ix]+cs[0]*v[iy];
								v[ix]=d;
							}
					}
				}
				else
				{
					kk=ks+1;
					fg[1]=e[kk-2];
					e[kk-2]=0.0;
					for(i=kk;i<=mm;i++)
					{
						fg[0]=s[i-1];
						sss(fg,cs);
						s[i-1]=fg[0];
						fg[1]=-cs[1]*e[i-1];
						e[i-1]=cs[0]*e[i-1];
						//if((cs[0]!=1.0)||(cs[1]!=0.0))
						if((fabs(cs[0]-1.0)>MIN_DOUBLE)||(fabs(cs[1])>MIN_DOUBLE))
							for(j=1;j<=m;j++)
							{
								ix=(j-1)*m+i-1;
								iy=(j-1)*m+kk-2;
								d=cs[0]*u[ix]+cs[1]*u[iy];								
								u[iy]=-cs[1]*u[ix]+cs[0]*u[iy];
								u[ix]=d;
							}
					}
				}
			}
			}
			}
			free(s);free(e);free(w);
			return l;
}
static void ppp(double a[],double e[],double s[],double v[],int m,int n)
{
	int i,j,p,q;
	double d;
	if(m>=n)
		i=n;
	else
		i=m;
	for(j=1;j<=i-1;j++)
	{
		a[(j-1)*n+j-1]=s[j-1];
		a[(j-1)*n+j]=e[j-1];
	}
	a[(i-1)*n+i-1]=s[i-1];
	if(m<n)
		a[(i-1)*n+i]=e[i-1];
	for(i=1;i<=n-1;i++)
	for(j=i+1;j<=n;j++)
	{
		p=(i-1)*n+j-1;
		q=(j-1)*n+i-1;
		d=v[p];v[p]=v[q];v[q]=d;
	}
	return;
}
static void sss(double fg[2],double cs[2])
{
	double r,d;
	//if((fabs(fg[0])+fabs(fg[1]))==0.0)
	if((fabs(fg[0])+fabs(fg[1]))<MIN_DOUBLE)
	{
		cs[0]=1.0;cs[1]=0.0;d=0.0;
	}
	else
	{
		d=sqrt(fg[0]*fg[0]+fg[1]*fg[1]);
		if(fabs(fg[0])>fabs(fg[1]))
		{
			d=fabs(d);
			if(fg[0]<0.0)
				d=-d;
		}
		if(fabs(fg[1])>=fabs(fg[0]))
		{
			d=fabs(d);
			if(fg[1]<0.0)
				d=-d;
		}
		cs[0]=fg[0]/d;
		cs[1]=fg[1]/d;
	}
	r=1.0;
	if(fabs(fg[0])>fabs(fg[1]))
		r=cs[1];
	else
		//if(cs[0]!=0.0)
		if(fabs(cs[0])>MIN_DOUBLE)
			r=1.0/cs[0];
	fg[0]=d;
	fg[1]=r;
	return;
}

static void damul(double a[],double b[],int m,int n,int k,double c[])
{
	int i,j,l,u;
	for(i=0;i<=m-1;i++)
		for(j=0;j<=k-1;j++)
		{
			u=i*k+j;
			c[u]=0;
			for(l=0;l<=n-1;l++)
				c[u]=c[u]+a[i*n+l]*b[l*k+j];
		}
	return;
}


