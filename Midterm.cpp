#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

int Generate_N(int p, int q, int r);
int Initial(double *x, double *y, int N);
int Print_Complex_Vector(double *x, double *y, int N);
int FFT_Multiple_of_two(double *x_r, double *x_i, double *y_r, double *y_i, int N);
int FFT_Multiple_of_three(double *x_r, double *x_i, double *y_r, double *y_i, int N);

int main(){
	int k, n, N, p, q, r;
	double *y_r, *y_i, *x_r, *x_i, w_r, w_i;
	clock_t t1, t2;
	printf("Please input p q r=");
	scanf("%d %d %d", &p, &q, &r);
	N = Generate_N(p, q, r);
	printf("N=2^%d 3^%d 5^%d = %d\n",p,q,r,N);	
	x_r = (double *) malloc(N*sizeof(double));
 	x_i = (double *) malloc(N*sizeof(double));
 	y_r = (double *) malloc(N*sizeof(double));
 	y_i = (double *) malloc(N*sizeof(double));
	Initial(x_r, x_i, N);	
	t1 = clock();
	FFT_Multiple_of_three(x_r, x_i, y_r, y_i, N);
	t2 = clock();
	printf("%f secs\n", 1.0*(t2-t1)/CLOCKS_PER_SEC);
	Print_Complex_Vector(y_r, y_i, N);
	
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

int FFT_Multiple_of_two(double *x_r, double *x_i, double *y_r, double *y_i, int N){
	if(N==1){
		y_r[0] = x_r[0];
		y_i[0] = x_i[0];
		return 0; 
	}
	
	double *u_r, *u_i, *v_r, *v_i, w_r, w_i;
		
	u_r = (double *) malloc(N*sizeof(double));
	u_i = (double *) malloc(N*sizeof(double));
	v_r = (double *) malloc(N*sizeof(double));
	v_i = (double *) malloc(N*sizeof(double));
	
	for(int n=0;n<N/2;n++){
		u_r[n] = x_r[2*n];
		u_i[n] = x_i[2*n];
		u_r[n+N/2] = x_r[2*n+1];
		u_i[n+N/2] = x_i[2*n+1];
	}
	
	FFT_Multiple_of_two(u_r, u_i, v_r, v_i, N/2);
	FFT_Multiple_of_two(u_r+N/2, u_i+N/2, v_r+N/2, v_i+N/2, N/2);
	
	double t_r,t_i;
	t_r = cos(-2*M_PI/N);
	t_i = sin(-2*M_PI/N);
	w_r = t_r;
	w_i = t_i;
		
	for(int k=0;k<N/2;++k){
		y_r[k] = v_r[k] + (w_r*v_r[k+N/2] - w_i*v_i[k+N/2]);
		y_i[k] = v_i[k] + (w_r*v_i[k+N/2] + w_i*v_r[k+N/2]);
		y_r[k+N/2] = v_r[k] - (w_r*v_r[k+N/2] - w_i*v_i[k+N/2]);
		y_i[k+N/2] = v_i[k] - (w_r*v_i[k+N/2] + w_i*v_r[k+N/2]);
		w_r = w_r*t_r - w_i*t_i;
		w_i = w_r*t_i + w_i*t_r; 	
	}
	
	free(u_r);
	free(u_i);
	free(v_r);
	free(v_i);
	
	return 0;
	
}

int FFT_Multiple_of_three(double *x_r, double *x_i, double *y_r, double *y_i, int N)
{
	if(N==1){
		y_r[0] = x_r[0];
		y_i[0] = x_i[0];
		return 0; 
	}
	
	double *u_r, *u_i , *v_r , *v_i , w_r , w_i , wr , wi ;
	
	u_r = (double *) malloc(N*sizeof(double));
	u_i = (double *) malloc(N*sizeof(double));
	v_r = (double *) malloc(N*sizeof(double));
	v_i = (double *) malloc(N*sizeof(double));
	
	for(int n=0;n<N/3;n++){
		u_r[n] = x_r[3*n];
		u_i[n] = x_i[3*n];
		u_r[n+N/3] = x_r[3*n+1];
		u_i[n+N/3] = x_i[3*n+1];
		u_r[n+2*N/3] = x_r[3*n+2];
		u_i[n+2*N/3] = x_i[3*n+2];
	}
	
	FFT_Multiple_of_three(u_r, u_i, v_r, v_i, N/3);
	FFT_Multiple_of_three(u_r+N/3, u_i+N/3, v_r+N/3, v_i+N/3, N/3);
	FFT_Multiple_of_three(u_r+2*N/3, u_i+2*N/3, v_r+2*N/3, v_i+2*N/3, N/3);
	
	double t_r , t_i , t1_r , t1_i ; 
	
	t1_r = cos(-2*M_PI/N);
	t1_i = sin(-2*M_PI/N);
	w_r = t1_r ;
	w_i = t1_i ;

	for(int k=0;k<N/3;++k){
		//w_r = cos(-k*2*M_PI/N);
		//w_i = sin(-k*2*M_PI/N);
		y_r[k] = v_r[k] + w_r*v_r[k+N/3] - w_i*v_i[k+N/3];
		y_i[k] = v_i[k] + w_r*v_i[k+N/3] + w_i*v_r[k+N/3];
		
		t_r = w_r;
		t_i = w_i;
		
		w_r = w_r * w_r - w_i * w_i ;
		w_i = 2 * w_r * w_i;
		
		//w_r = cos(-k*4*M_PI/N);
		//w_i = sin(-k*4*M_PI/N);
		
		y_r[k] += w_r*v_r[k+2*N/3] - w_i*v_i[k+2*N/3];
		y_i[k] += w_r*v_i[k+2*N/3] + w_i*v_r[k+2*N/3];
		
		w_r = cos(-(k+N/3)*2*M_PI/N);
		w_i = sin(-(k+N/3)*2*M_PI/N);
		y_r[k+N/3] = v_r[k] + w_r*v_r[k+N/3] - w_i*v_i[k+N/3];
		y_i[k+N/3] = v_i[k] + w_r*v_i[k+N/3] + w_i*v_r[k+N/3];
		w_r = cos(-(k+N/3)*4*M_PI/N);
		w_i = sin(-(k+N/3)*4*M_PI/N);
		y_r[k+N/3] += w_r*v_r[k+2*N/3] - w_i*v_i[k+2*N/3];
		y_i[k+N/3] += w_r*v_i[k+2*N/3] + w_i*v_r[k+2*N/3];
		
		w_r = cos(-(k+2*N/3)*2*M_PI/N);
		w_i = sin(-(k+2*N/3)*2*M_PI/N);
		y_r[k+2*N/3] = v_r[k] + w_r*v_r[k+N/3] - w_i*v_i[k+N/3];
		y_i[k+2*N/3] = v_i[k] + w_r*v_i[k+N/3] + w_i*v_r[k+N/3] ;
		w_r = cos(-(k+2*N/3)*4*M_PI/N);
		w_i = sin(-(k+2*N/3)*4*M_PI/N);
		y_r[k+2*N/3] += w_r*v_r[k+2*N/3] - w_i*v_i[k+2*N/3];
		y_i[k+2*N/3] += w_r*v_i[k+2*N/3] + w_i*v_r[k+2*N/3];
		
		w_r = t1_r * t_r + t1_i * t_i;
		w_i = t1_i * t_r + t1_r * t_r; 
		
	}
	
	free(u_r);
	free(u_i);
	free(v_r);
	free(v_i);
	
	return 0;
}
