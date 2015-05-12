#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

int Generate_N(int p, int q, int r);
int Initial(double *x, double *y, int N);
int Print_Complex_Vector(double *x, double *y, int N);
int FFT_Multiple_of_two(double *x_r, double *x_i, double *y_r, double *y_i, int N);
int FFT_Multiple_of_three(double *x_r, double *x_i, double *y_r, double *y_i, int N);
int FFT_Multiple_of_five(double *x_r, double *x_i, double *y_r, double *y_i, int N);
//int FFT_callfunction(double *x_r, double *x_i, double *y_r, double *y_i, int N , int p , int q, int r);
int FFT_2(double *x_r, double *x_i, double *y_r, double *y_i, int N);

int main(){
	int k, n, p, q, r, N ;
	double *y_r, *y_i, *x_r, *x_i, w_r, w_i;
	clock_t t1, t2,t3 ,t4;
	printf("Please input p q r=");
	scanf("%d %d %d", &p, &q, &r);
	N = Generate_N(p, q, r);
	printf("N=2^%d 3^%d 5^%d = %d\n",p,q,r,N);	
	x_r = (double *) malloc(N*sizeof(double));
 	x_i = (double *) malloc(N*sizeof(double));
 	y_r = (double *) malloc(N*sizeof(double));
 	y_i = (double *) malloc(N*sizeof(double));
	Initial(x_r, x_i, N);	
	//FFT_callfunction(x_r , x_i , y_r , y_i , N , p ,q , r);
	
	if (p == 0 && q == 0 && r == 0){
		return 0;
	}else if((N%5) == 0){
		t1 = clock();
		FFT_Multiple_of_five(x_r, x_i, y_r, y_i, N);
		t2 = clock();
		printf("%f secs\n", 1.0*(t2-t1)/CLOCKS_PER_SEC);	
		system("pause");
		//Print_Complex_Vector(y_r, y_i, N);	
	}else if((N%3) == 0){
		t1 = clock();
		FFT_Multiple_of_three(x_r, x_i, y_r, y_i, N);
		t2 = clock();
		printf("%f secs\n", 1.0*(t2-t1)/CLOCKS_PER_SEC);	
		system("pause");
		//Print_Complex_Vector(y_r, y_i, N);	
	}else{
		t1 = clock();
		FFT_Multiple_of_two(x_r, x_i, y_r, y_i, N);
		t2 = clock();
		printf("%f secs\n", 1.0*(t2-t1)/CLOCKS_PER_SEC);	
		system("pause");
		/*t3 = clock();
		FFT_2(x_r, x_i, y_r, y_i, N);
		t4 = clock();
		printf("%f secs\n", 1.0*(t4-t3)/CLOCKS_PER_SEC);
		*/
		//Print_Complex_Vector(y_r, y_i, N);	
	}
	
	Print_Complex_Vector(y_r, y_i, N);
	free(x_r);
	free(x_i);
	free(y_r);
	free(y_i);

	
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
	int k, n;
	double *u_r, *u_i, *v_r, *v_i, w_r, w_i;
		
	u_r = (double *) malloc(N*sizeof(double));
	u_i = (double *) malloc(N*sizeof(double));
	v_r = (double *) malloc(N*sizeof(double));
	v_i = (double *) malloc(N*sizeof(double));
	int N2=N/2;
	for(n=0;n<N2;n++){
		u_r[n] = x_r[2*n];
		u_i[n] = x_i[2*n];
		u_r[n+N/2] = x_r[2*n+1];
		u_i[n+N/2] = x_i[2*n+1];
	}
	FFT_Multiple_of_two(u_r, u_i, v_r, v_i, N/2);
	FFT_Multiple_of_two(u_r+N/2, u_i+N/2, v_r+N/2, v_i+N/2, N/2);

	double t_r , t_i , t1_r , t1_i ;
	t1_r = cos(-2*M_PI/N);
	t1_i = sin(-2*M_PI/N);
	w_r = 1;
	w_i = 0;	
	
	for(k=0;k<N2;++k){
		t_r = w_r;
		t_i = w_i;
		y_r[k] = v_r[k] + (w_r*v_r[k+N/2] - w_i*v_i[k+N/2]);
		y_i[k] = v_i[k] + (w_r*v_i[k+N/2] + w_i*v_r[k+N/2]);
		y_r[k+N/2] = v_r[k] - (w_r*v_r[k+N/2] - w_i*v_i[k+N/2]);
		y_i[k+N/2] = v_i[k] - (w_r*v_i[k+N/2] + w_i*v_r[k+N/2]);
		w_r = t1_r*t_r - t1_i*t_i;
		w_i = t1_r*t_i + t1_i*t_r; 	
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
	
	double *u_r, *u_i , *v_r , *v_i , w_r , w_i ;
	
	u_r = (double *) malloc(N*sizeof(double));
	u_i = (double *) malloc(N*sizeof(double));
	v_r = (double *) malloc(N*sizeof(double));
	v_i = (double *) malloc(N*sizeof(double));
	int N3=N/3;
	int n;
	for(n=0;n<N3;n++){
		u_r[n] = x_r[3*n];
		u_i[n] = x_i[3*n];
		u_r[n+N/3] = x_r[3*n+1];
		u_i[n+N/3] = x_i[3*n+1];
		u_r[n+2*N/3] = x_r[3*n+2];
		u_i[n+2*N/3] = x_i[3*n+2];
	}
	if((n%3) == 0){
		FFT_Multiple_of_three(u_r, u_i, v_r, v_i, N/3);
		FFT_Multiple_of_three(u_r+N/3, u_i+N/3, v_r+N/3, v_i+N/3, N/3);
		FFT_Multiple_of_three(u_r+2*N/3, u_i+2*N/3, v_r+2*N/3, v_i+2*N/3, N/3);
	}else{
		FFT_Multiple_of_two(u_r, u_i, v_r, v_i, N/3);
		FFT_Multiple_of_two(u_r+N/3, u_i+N/3, v_r+N/3, v_i+N/3, N/3);
		FFT_Multiple_of_two(u_r+2*N/3, u_i+2*N/3, v_r+2*N/3, v_i+2*N/3, N/3);
	}
	
	double t_r , t_i , t1_r , t1_i , p_r , p_i; 
	
	t1_r = cos(-2.0*M_PI/N);
	t1_i = sin(-2.0*M_PI/N);
	w_r = 1 ;
	w_i = 0 ;
    

	for(int k=0;k<N3;++k){
		
		//w_r = cos(-k*2*M_PI/N);
		//w_i = sin(-k*2*M_PI/N);
	
		y_r[k] = v_r[k] + w_r*v_r[k+N/3] - w_i*v_i[k+N/3];
		y_i[k] = v_i[k] + w_r*v_i[k+N/3] + w_i*v_r[k+N/3];
		
		t_r = w_r;
		t_i = w_i;
		
		p_r = t_r * t_r - t_i * t_i ;
		p_i = 2 * t_r * t_i;
		
		w_r = p_r;
		w_i = p_i;
	
		//w_r = cos(-k*4*M_PI/N);
		//w_i = sin(-k*4*M_PI/N);
		
		y_r[k] += w_r*v_r[k+2*N/3] - w_i*v_i[k+2*N/3];
		y_i[k] += w_r*v_i[k+2*N/3] + w_i*v_r[k+2*N/3];
		
		//w_r = cos(-(k+N/3)*2*M_PI/N) = cos( -k*2*M_PI/N + -2/3*M_PI ) = cos( -k*2*M_PI/N- 120 ) = cos(-k*2*M_PI/N)cos(-120) - sin(-k*2*M_PI/N)sin(-120)
		//w_i = sin(-(k+N/3)*2*M_PI/N) = sin( -k*2*M_PI/N + -2/3*M_PI ) = sin( -k*2*M_PI/N- 120 ) = sin(-k*2*M_PI/N)cos(-120) + cos(-k*2*M_PI/N)sin(-120)
		//w_r = cos(-(k+N/3)*2*M_PI/N);
		//w_i = sin(-(k+N/3)*2*M_PI/N);
		
		w_r = t_r*(-1.0/2.0)-t_i*(-1.732050807568877/2.0);
		w_i = t_i*(-1.0/2.0)+t_r*(-1.732050807568877/2.0);
		
		y_r[k+N/3] = v_r[k] + w_r*v_r[k+N/3] - w_i*v_i[k+N/3];
		y_i[k+N/3] = v_i[k] + w_r*v_i[k+N/3] + w_i*v_r[k+N/3];
		
		//w_r = cos(-(k+N/3)*4*M_PI/N) = cos( -k*4*M_PI/N + -4/3*M_PI ) = cos( -k*4*M_PI/N- 240 ) = cos(-k*4*M_PI/N)cos(-240) - sin(-k*4*M_PI/N)sin(-240)
		//w_i = sin(-(k+N/3)*4*M_PI/N) = sin( -k*4*M_PI/N + -4/3*M_PI ) = sin( -k*4*M_PI/N- 240 ) = sin(-k*4*M_PI/N)cos(-240) + cos(-k*4*M_PI/N)sin(-240)	
		//w_r = cos(-(k+N/3)*4*M_PI/N);
		//w_i = sin(-(k+N/3)*4*M_PI/N);
		
		w_r = p_r*(-1.0/2.0)-p_i*(1.732050807568877/2.0);
		w_i = p_i*(-1.0/2.0)+p_r*(1.732050807568877/2.0);
		
		y_r[k+N/3] += w_r*v_r[k+2*N/3] - w_i*v_i[k+2*N/3];
		y_i[k+N/3] += w_r*v_i[k+2*N/3] + w_i*v_r[k+2*N/3];
		
		//w_r = cos(-(k+2*N/3)*2*M_PI/N) = cos( -k*2*M_PI/N + -4/3*M_PI ) = cos( -k*2*M_PI/N- 240 ) = cos(-k*2*M_PI/N)cos(-240) - sin(-k*2*M_PI/N)sin(-240)
		//w_i = sin(-(k+2*N/3)*2*M_PI/N) = sin( -k*2*M_PI/N + -4/3*M_PI ) = sin( -k*2*M_PI/N- 240 ) = sin(-k*2*M_PI/N)cos(-240) + cos(-k*2*M_PI/N)sin(-240)
		//w_r = cos(-(k+2*N/3)*2*M_PI/N);
		//w_i = sin(-(k+2*N/3)*2*M_PI/N);
		
		w_r = t_r*(-1.0/2.0)-t_i*(1.732050807568877/2.0);
		w_i = t_i*(-1.0/2.0)+t_r*(1.732050807568877/2.0);
		
		y_r[k+2*N/3] = v_r[k] + w_r*v_r[k+N/3] - w_i*v_i[k+N/3];
		y_i[k+2*N/3] = v_i[k] + w_r*v_i[k+N/3] + w_i*v_r[k+N/3];
		
		//w_r = cos(-(k+2*N/3)*4*M_PI/N) = cos( -k*4*M_PI/N + -8/3*M_PI ) = cos( -k*4*M_PI/N- 480 ) = cos(-k*4*M_PI/N)cos(-120) - sin(-k*4*M_PI/N)sin(-120)
		//w_i = sin(-(k+2*N/3)*4*M_PI/N) = sin( -k*4*M_PI/N + -8/3*M_PI ) = sin( -k*4*M_PI/N- 480 ) = sin(-k*4*M_PI/N)cos(-120) + cos(-k*4*M_PI/N)sin(-120)	
		//w_r = cos(-(k+2*N/3)*4*M_PI/N);
		//w_i = sin(-(k+2*N/3)*4*M_PI/N);
		
		w_r = p_r*(-1.0/2.0)-p_i*(-1.732050807568877/2.0);
		w_i = p_i*(-1.0/2.0)+p_r*(-1.732050807568877/2.0);
		
		y_r[k+2*N/3] += w_r*v_r[k+2*N/3] - w_i*v_i[k+2*N/3];
		y_i[k+2*N/3] += w_r*v_i[k+2*N/3] + w_i*v_r[k+2*N/3];
			
		w_r = t1_r * t_r - t1_i * t_i;
		w_i = t1_i * t_r + t1_r * t_i; 
		
	}
	
	free(u_r);
	free(u_i);
	free(v_r);
	free(v_i);
	
	return 0;
}

int FFT_Multiple_of_five(double *x_r, double *x_i, double *y_r, double *y_i, int N)
{
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
	int N5=N/5;
	int n;
	for(n=0;n<N5;n++){
		u_r[n] = x_r[5*n];
		u_i[n] = x_i[5*n];
		u_r[n+N/5] = x_r[5*n+1];
		u_i[n+N/5] = x_i[5*n+1];
		u_r[n+2*N/5] = x_r[5*n+2];
		u_i[n+2*N/5] = x_i[5*n+2];
		u_r[n+3*N/5] = x_r[5*n+3];
		u_i[n+3*N/5] = x_i[5*n+3];
		u_r[n+4*N/5] = x_r[5*n+4];
		u_i[n+4*N/5] = x_i[5*n+4];
	}
	if((n%5) == 0){
		FFT_Multiple_of_five(u_r, u_i, v_r, v_i, N/5);
		FFT_Multiple_of_five(u_r+1*N/5, u_i+1*N/5, v_r+1*N/5, v_i+1*N/5, N/5);
		FFT_Multiple_of_five(u_r+2*N/5, u_i+2*N/5, v_r+2*N/5, v_i+2*N/5, N/5);
		FFT_Multiple_of_five(u_r+3*N/5, u_i+3*N/5, v_r+3*N/5, v_i+3*N/5, N/5);
		FFT_Multiple_of_five(u_r+4*N/5, u_i+4*N/5, v_r+4*N/5, v_i+4*N/5, N/5);
	}else if((n%3) == 0){
		FFT_Multiple_of_three(u_r, u_i, v_r, v_i, N/5);
		FFT_Multiple_of_three(u_r+1*N/5, u_i+1*N/5, v_r+1*N/5, v_i+1*N/5, N/5);
		FFT_Multiple_of_three(u_r+2*N/5, u_i+2*N/5, v_r+2*N/5, v_i+2*N/5, N/5);
		FFT_Multiple_of_three(u_r+3*N/5, u_i+3*N/5, v_r+3*N/5, v_i+3*N/5, N/5);
		FFT_Multiple_of_three(u_r+4*N/5, u_i+4*N/5, v_r+4*N/5, v_i+4*N/5, N/5);
	}else{
		FFT_Multiple_of_two(u_r, u_i, v_r, v_i, N/5);
		FFT_Multiple_of_two(u_r+1*N/5, u_i+1*N/5, v_r+1*N/5, v_i+1*N/5, N/5);
		FFT_Multiple_of_two(u_r+2*N/5, u_i+2*N/5, v_r+2*N/5, v_i+2*N/5, N/5);
		FFT_Multiple_of_two(u_r+3*N/5, u_i+3*N/5, v_r+3*N/5, v_i+3*N/5, N/5);
		FFT_Multiple_of_two(u_r+4*N/5, u_i+4*N/5, v_r+4*N/5, v_i+4*N/5, N/5);
	}
	
	double t_r , t_i , t1_r , t1_i , p_r , p_i , p6_r , p6_i , p8_r , p8_i ; 
	
	t1_r = cos(-2.0*M_PI/N);
	t1_i = sin(-2.0*M_PI/N);
	w_r = 1 ;
	w_i = 0 ;
	
	for(int k=0;k<N5;++k){
		
		//w_r = cos(-k*2*M_PI/N);
		//w_i = sin(-k*2*M_PI/N);
		
		y_r[k] = v_r[k] + w_r*v_r[k+N/5] - w_i*v_i[k+N/5];
		y_i[k] = v_i[k] + w_r*v_i[k+N/5] + w_i*v_r[k+N/5];
		
		//---------------------------------------------------------------
		t_r = w_r; // t_r = cos(-k*2*M_PI/N)
		t_i = w_i; // t_i = sin(-k*2*M_PI/N)
		
		p_r = t_r * t_r - t_i * t_i ;
		p_i = 2 * t_r * t_i;
		
		w_r = p_r; //p_r = cos(-k*4*M_PI/N);
		w_i = p_i; //p_i = sin(-k*4*M_PI/N);
		
		//w_r = cos(-k*4*M_PI/N);
		//w_i = sin(-k*4*M_PI/N);
		
		y_r[k] += w_r*v_r[k+2*N/5] - w_i*v_i[k+2*N/5]; 
		y_i[k] += w_r*v_i[k+2*N/5] + w_i*v_r[k+2*N/5];
		
		//---------------------------------------------------------------
		
		//w_r = cos(-k*6*M_PI/N) = cos(-k*2*M_PI/N)cos(-k*4*M_PI/N)-sin(-k*2*M_PI/N)sin(-k*4*M_PI/N)
		//w_i = sin(-k*6*M_PI/N) = sin(-k*2*M_PI/N)cos(-k*4*M_PI/N)+cos(-k*2*M_PI/N)sin(-k*4*M_PI/N)
		//w_r = cos(-k*6*M_PI/N);
		//w_i = sin(-k*6*M_PI/N);
		
		p6_r = t_r * p_r - t_i * p_i; // p6_r = cos(-k*6*M_PI/N);
		p6_i = t_i * p_r + t_r * p_i; // p6_i = sin(-k*6*M_PI/N);
		
		w_r = p6_r;
		w_i = p6_i;
		
		y_r[k] += w_r*v_r[k+3*N/5] - w_i*v_i[k+3*N/5];
		y_i[k] += w_r*v_i[k+3*N/5] + w_i*v_r[k+3*N/5];
		
		//---------------------------------------------------------------
		
		//w_r = cos(-k*8*M_PI/N) = cos(-k*4*M_PI/N)cos(-k*4*M_PI/N)-sin(-k*4*M_PI/N)sin(-k*4*M_PI/N)
		//w_i = sin(-k*8*M_PI/N) = sin(-k*4*M_PI/N)cos(-k*4*M_PI/N)+cos(-k*4*M_PI/N)sin(-k*4*M_PI/N)
		//w_r = cos(-k*8*M_PI/N);
		//w_i = sin(-k*8*M_PI/N);
		
		p8_r = p_r * p_r - p_i * p_i; // p8_r = cos(-k*8*M_PI/N);
		p8_i = p_i * p_r + p_r * p_i; // p8_i = sin(-k*8*M_PI/N);
		
		w_r = p8_r;
		w_i = p8_i;				
		
		y_r[k] += w_r*v_r[k+4*N/5] - w_i*v_i[k+4*N/5];
		y_i[k] += w_r*v_i[k+4*N/5] + w_i*v_r[k+4*N/5];
		
		//---------------------------------------------------------------
		
		//w_r = cos(-(k+N/5)*2*M_PI/N) = cos( -k*2*M_PI/N + -2/5*M_PI ) = cos( -k*2*M_PI/N- 72 ) = cos(-k*2*M_PI/N)cos(-72) - sin(-k*2*M_PI/N)sin(-72)
		//w_i = sin(-(k+N/5)*2*M_PI/N) = sin( -k*2*M_PI/N + -2/5*M_PI ) = sin( -k*2*M_PI/N- 72 ) = sin(-k*2*M_PI/N)cos(-72) + cos(-k*2*M_PI/N)sin(-72)	
		//w_r = cos(-(k+N/5)*2*M_PI/N);
		//w_i = sin(-(k+N/5)*2*M_PI/N);
		
		w_r = t_r*(0.309016994) - t_i*(-0.951056516);
		w_i = t_i*(0.309016994) + t_r*(-0.951056516);
		
		y_r[k+N/5] = v_r[k] + w_r*v_r[k+N/5] - w_i*v_i[k+N/5];
		y_i[k+N/5] = v_i[k] + w_r*v_i[k+N/5] + w_i*v_r[k+N/5];
		
		//---------------------------------------------------------------
		
		//w_r = cos(-(k+N/5)*4*M_PI/N) = cos( -k*4*M_PI/N + -4/5*M_PI ) = cos( -k*4*M_PI/N- 144 ) = cos(-k*4*M_PI/N)cos(-144) - sin(-k*4*M_PI/N)sin(-144)
		//w_i = sin(-(k+N/5)*4*M_PI/N) = sin( -k*4*M_PI/N + -4/5*M_PI ) = sin( -k*4*M_PI/N- 144 ) = sin(-k*4*M_PI/N)cos(-144) + cos(-k*4*M_PI/N)sin(-144)
		//w_r = cos(-(k+N/5)*4*M_PI/N);
		//w_i = sin(-(k+N/5)*4*M_PI/N);
		
		w_r = p_r*(-0.809016994) - p_i*(-0.587785252);
		w_i = p_i*(-0.809016994) + p_r*(-0.587785252);
		
		y_r[k+N/5] += w_r*v_r[k+2*N/5] - w_i*v_i[k+2*N/5]; 
		y_i[k+N/5] += w_r*v_i[k+2*N/5] + w_i*v_r[k+2*N/5];
	
		//---------------------------------------------------------------	

		//w_r = cos(-(k+N/5)*6*M_PI/N) = cos( -k*6*M_PI/N + -6/5*M_PI ) = cos( -k*6*M_PI/N- 144 ) = cos(-k*4*M_PI/N)cos(-216) - sin(-k*4*M_PI/N)sin(-216)
		//w_i = sin(-(k+N/5)*6*M_PI/N) = sin( -k*6*M_PI/N + -6/5*M_PI ) = sin( -k*6*M_PI/N- 144 ) = sin(-k*4*M_PI/N)cos(-216) + cos(-k*4*M_PI/N)sin(-216)		
		//w_r = cos(-(k+N/5)*6*M_PI/N);
		//w_i = sin(-(k+N/5)*6*M_PI/N);
		
		w_r = p6_r*(-0.809016994) - p6_i*(0.587785252);
		w_i = p6_i*(-0.809016994) + p6_r*(0.587785252);		
		
		y_r[k+N/5] += w_r*v_r[k+3*N/5] - w_i*v_i[k+3*N/5];
		y_i[k+N/5] += w_r*v_i[k+3*N/5] + w_i*v_r[k+3*N/5];

		//---------------------------------------------------------------
		
		//w_r = cos(-(k+N/5)*8*M_PI/N) = cos( -k*8*M_PI/N + -8/5*M_PI ) = cos( -k*8*M_PI/N- 288 ) = cos(-k*8*M_PI/N)cos(-288) - sin(-k*8*M_PI/N)sin(-288)
		//w_i = sin(-(k+N/5)*8*M_PI/N) = cos( -k*8*M_PI/N + -8/5*M_PI ) = cos( -k*8*M_PI/N- 288 ) = cos(-k*8*M_PI/N)cos(-288) - sin(-k*8*M_PI/N)sin(-288)	
		//w_r = cos(-(k+N/5)*8*M_PI/N);
		//w_i = sin(-(k+N/5)*8*M_PI/N);
		
		w_r = p8_r*(0.309016994) - p8_i*(0.951056516);
		w_i = p8_i*(0.309016994) + p8_r*(0.951056516);			
				
		y_r[k+N/5] += w_r*v_r[k+4*N/5] - w_i*v_i[k+4*N/5];
		y_i[k+N/5] += w_r*v_i[k+4*N/5] + w_i*v_r[k+4*N/5];

		//---------------------------------------------------------------
		
		//w_r = cos(-(k+2*N/5)*2*M_PI/N) = cos( -k*2*M_PI/N + -4/5*M_PI ) = cos( -k*2*M_PI/N- 144 ) = cos(-k*2*M_PI/N)cos(-144) - sin(-k*2*M_PI/N)sin(-144)
		//w_i = sin(-(k+2*N/5)*2*M_PI/N) = sin( -k*2*M_PI/N + -4/5*M_PI ) = sin( -k*2*M_PI/N- 144 ) = sin(-k*2*M_PI/N)cos(-144) + cos(-k*2*M_PI/N)sin(-144)				
		//w_r = cos(-(k+2*N/5)*2*M_PI/N);
		//w_i = sin(-(k+2*N/5)*2*M_PI/N);
		
		w_r = t_r*(-0.809016994) - t_i*(-0.587785252);
		w_i = t_i*(-0.809016994) + t_r*(-0.587785252);	
		
		y_r[k+2*N/5] = v_r[k] + w_r*v_r[k+N/5] - w_i*v_i[k+N/5];
		y_i[k+2*N/5] = v_i[k] + w_r*v_i[k+N/5] + w_i*v_r[k+N/5];

		//---------------------------------------------------------------

		//w_r = cos(-(k+2*N/5)*4*M_PI/N) = cos( -k*4*M_PI/N + -8/5*M_PI ) = cos( -k*4*M_PI/N- 288) = cos(-k*4*M_PI/N)cos(-288) - sin(-k*4*M_PI/N)sin(-288)
		//w_i = sin(-(k+2*N/5)*4*M_PI/N) = sin( -k*4*M_PI/N + -8/5*M_PI ) = sin( -k*4*M_PI/N- 288) = sin(-k*4*M_PI/N)cos(-288) + cos(-k*4*M_PI/N)sin(-288)	
		//w_r = cos(-(k+2*N/5)*4*M_PI/N);
		//w_i = sin(-(k+2*N/5)*4*M_PI/N);
		
		w_r = p_r*(0.309016994) - p_i*(0.951056516);
		w_i = p_i*(0.309016994) + p_r*(0.951056516);			
		
		y_r[k+2*N/5] += w_r*v_r[k+2*N/5] - w_i*v_i[k+2*N/5]; 
		y_i[k+2*N/5] += w_r*v_i[k+2*N/5] + w_i*v_r[k+2*N/5];

		//---------------------------------------------------------------

		//w_r = cos(-(k+2*N/5)*6*M_PI/N) = cos( -k*6*M_PI/N + 8/5*M_PI ) = cos( -k*2*M_PI/N+ 288 ) = cos(-k*2*M_PI/N)cos(288) - sin(-k*2*M_PI/N)sin(288)
		//w_i = sin(-(k+2*N/5)*6*M_PI/N) = sin( -k*6*M_PI/N + 8/5*M_PI ) = sin( -k*2*M_PI/N+ 288 ) = sin(-k*2*M_PI/N)cos(288) + cos(-k*2*M_PI/N)sin(288)
		//w_r = cos(-(k+2*N/5)*6*M_PI/N);
		//w_i = sin(-(k+2*N/5)*6*M_PI/N);
		
		w_r = p6_r*(0.309016994) - p6_i*(-0.951056516);
		w_i = p6_i*(0.309016994) + p6_r*(-0.951056516);
		
		y_r[k+2*N/5] += w_r*v_r[k+3*N/5] - w_i*v_i[k+3*N/5];
		y_i[k+2*N/5] += w_r*v_i[k+3*N/5] + w_i*v_r[k+3*N/5];

		//---------------------------------------------------------------	

		//w_r = cos(-(k+2*N/5)*8*M_PI/N) = cos( -k*8*M_PI/N + 4/5*M_PI ) = cos( -k*2*M_PI/N+ 144 ) = cos(-k*2*M_PI/N)cos(144) - sin(-k*2*M_PI/N)sin(144)
		//w_i = sin(-(k+2*N/5)*8*M_PI/N) = sin( -k*8*M_PI/N + 4/5*M_PI ) = sin( -k*2*M_PI/N+ 144 ) = sin(-k*2*M_PI/N)cos(144) + cos(-k*2*M_PI/N)sin(144)	
		//w_r = cos(-(k+2*N/5)*8*M_PI/N);
		//w_i = sin(-(k+2*N/5)*8*M_PI/N);
		
		w_r = p8_r*(-0.809016994) - p8_i*(0.587785252);
		w_i = p8_i*(-0.809016994) + p8_r*(0.587785252);	
		
		y_r[k+2*N/5] += w_r*v_r[k+4*N/5] - w_i*v_i[k+4*N/5];
		y_i[k+2*N/5] += w_r*v_i[k+4*N/5] + w_i*v_r[k+4*N/5];

		//---------------------------------------------------------------
		
		//w_r = cos(-(k+3*N/5)*2*M_PI/N);
		//w_i = sin(-(k+3*N/5)*2*M_PI/N);
		
		w_r = t_r*(-0.809016994) - t_i*(0.587785252);
		w_i = t_i*(-0.809016994) + t_r*(0.587785252);	
	
		y_r[k+3*N/5] = v_r[k] + w_r*v_r[k+N/5] - w_i*v_i[k+N/5];
		y_i[k+3*N/5] = v_i[k] + w_r*v_i[k+N/5] + w_i*v_r[k+N/5];

		//---------------------------------------------------------------

		//w_r = cos(-(k+3*N/5)*4*M_PI/N);
		//w_i = sin(-(k+3*N/5)*4*M_PI/N);
		
		w_r = p_r*(0.309016994) - p_i*(-0.951056516);
		w_i = p_i*(0.309016994) + p_r*(-0.951056516);
		
		y_r[k+3*N/5] += w_r*v_r[k+2*N/5] - w_i*v_i[k+2*N/5]; 
		y_i[k+3*N/5] += w_r*v_i[k+2*N/5] + w_i*v_r[k+2*N/5];
	
		//---------------------------------------------------------------	

		//w_r = cos(-(k+3*N/5)*6*M_PI/N);
		//w_i = sin(-(k+3*N/5)*6*M_PI/N);
		
		w_r = p6_r*(0.309016994) - p6_i*(0.951056516);
		w_i = p6_i*(0.309016994) + p6_r*(0.951056516);
		
		y_r[k+3*N/5] += w_r*v_r[k+3*N/5] - w_i*v_i[k+3*N/5];
		y_i[k+3*N/5] += w_r*v_i[k+3*N/5] + w_i*v_r[k+3*N/5];

		//---------------------------------------------------------------

		//w_r = cos(-(k+3*N/5)*8*M_PI/N);
		//w_i = sin(-(k+3*N/5)*8*M_PI/N);
			
		w_r = p8_r*(-0.809016994) - p8_i*(-0.587785252);
		w_i = p8_i*(-0.809016994) + p8_r*(-0.587785252);
		
		y_r[k+3*N/5] += w_r*v_r[k+4*N/5] - w_i*v_i[k+4*N/5];
		y_i[k+3*N/5] += w_r*v_i[k+4*N/5] + w_i*v_r[k+4*N/5];

		//---------------------------------------------------------------
		
		//w_r = cos(-(k+4*N/5)*2*M_PI/N);
		//w_i = sin(-(k+4*N/5)*2*M_PI/N);
			
		w_r = t_r*(0.309016994) - t_i*(0.951056516);
		w_i = t_i*(0.309016994) + t_r*(0.951056516);
		
		y_r[k+4*N/5] = v_r[k] + w_r*v_r[k+N/5] - w_i*v_i[k+N/5];
		y_i[k+4*N/5] = v_i[k] + w_r*v_i[k+N/5] + w_i*v_r[k+N/5];
		
		//---------------------------------------------------------------	

		//w_r = cos(-(k+4*N/5)*4*M_PI/N);
		//w_i = sin(-(k+4*N/5)*4*M_PI/N);
			
		w_r = p_r*(-0.809016994) - p_i*(0.587785252);
		w_i = p_i*(-0.809016994) + p_r*(0.587785252);
		
		y_r[k+4*N/5] += w_r*v_r[k+2*N/5] - w_i*v_i[k+2*N/5]; 
		y_i[k+4*N/5] += w_r*v_i[k+2*N/5] + w_i*v_r[k+2*N/5];

		//---------------------------------------------------------------

		//w_r = cos(-(k+4*N/5)*6*M_PI/N);
		//w_i = sin(-(k+4*N/5)*6*M_PI/N);

		w_r = p6_r*(-0.809016994) - p6_i*(-0.587785252);
		w_i = p6_i*(-0.809016994) + p6_r*(-0.587785252);		
		
		y_r[k+4*N/5] += w_r*v_r[k+3*N/5] - w_i*v_i[k+3*N/5];
		y_i[k+4*N/5] += w_r*v_i[k+3*N/5] + w_i*v_r[k+3*N/5];
	
		//---------------------------------------------------------------	
	
		//w_r = cos(-(k+4*N/5)*8*M_PI/N);
		//w_i = sin(-(k+4*N/5)*8*M_PI/N);
	
		w_r = p8_r*(0.309016994) - p8_i*(-0.951056516);
		w_i = p8_i*(0.309016994) + p8_r*(-0.951056516);
		
		y_r[k+4*N/5] += w_r*v_r[k+4*N/5] - w_i*v_i[k+4*N/5];
		y_i[k+4*N/5] += w_r*v_i[k+4*N/5] + w_i*v_r[k+4*N/5];
		
		//---------------------------------------------------------------	
		
		w_r = t1_r * t_r - t1_i * t_i;
		w_i = t1_i * t_r + t1_r * t_i; 
	}
	
	free(u_r);
	free(u_i);
	free(v_r);
	free(v_i);
	
	return 0;
}
/*
int FFT_callfunction(double *x_r, double *x_i, double *y_r, double *y_i, int N , int p , int q, int r){
	clock_t t1, t2;
	if(p == 0 && q == 0 && r == 0){
		return 0;
	}else if((N%5) == 0){
		t1 = clock();
		FFT_Multiple_of_five(x_r, x_i, y_r, y_i, N);
		t2 = clock();
		printf("%f secs\n", 1.0*(t2-t1)/CLOCKS_PER_SEC);
	}else if((N%3) == 0){
		t1 = clock();
		FFT_Multiple_of_three(x_r, x_i, y_r, y_i, N);
		t2 = clock();
		printf("%f secs\n", 1.0*(t2-t1)/CLOCKS_PER_SEC);
	}else{
		t1 = clock();
		FFT_Multiple_of_two(x_r, x_i, y_r, y_i, N);
		t2 = clock();
		printf("%f secs\n", 1.0*(t2-t1)/CLOCKS_PER_SEC);
	}
	free(x_r);
	free(x_i);
	free(y_r);
	free(y_i);
	return 0;	
}
*/
