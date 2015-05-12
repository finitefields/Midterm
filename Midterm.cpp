#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

int Generate_N(int p, int q, int r);
int Initial(double *x, double *y, int N);
int Print_Complex_Vector(double *x, double *y, int N);

int main(){
	int k, n, N, p, q, r;
	double *y_r, *y_i, *x_r, *x_i, w_r, w_i;
	clock_t t1, t2;
	printf("Please input p q r=");
	scanf("%d %d %d", &p, &q, &r);
	N = Generate_N(p, q, r);
	printf("N=2^%d 3^%d 5^%d = %d\n",p,q,r,N);		
	return 0;
}



int Generate_N(int p, int q, int r)
{
	int N = 1;
	for(int i=0;i<p;++i) N = N * 2;
	for(int i=0;i<q;++i) N = N * 3;
	for(int i=0;i<r;++i) N = N * 5;	
	return N;
}

int Initial(double *x, double *y, int N)
{
	for(int n=0;n<N;++n){
		x[n] = n; 
		y[n] = 0; 
	}
}
int Print_Complex_Vector(double *x, double *y, int N)
{
	for(int n=0;n<N;++n){
		if(y[n]>=0) printf("%d : %f +%f i\n", n, x[n], y[n]);
		else printf("%d : %f %f i\n", n, x[n], y[n]);
	}
}
