/*************************************************************************/
/*                                                                       */
/*   SNU-RT Benchmark Suite for Worst Case Timing Analysis               */
/*   =====================================================               */
/*                              Collected and Modified by S.-S. Lim      */
/*                                           sslim@archi.snu.ac.kr       */
/*                                         Real-Time Research Group      */
/*                                        Seoul National University      */
/*                                                                       */
/*                                                                       */
/*        < Features > - restrictions for our experimental environment   */
/*                                                                       */
/*          1. Completely structured.                                    */
/*               - There are no unconditional jumps.                     */
/*               - There are no exit from loop bodies.                   */
/*                 (There are no 'break' or 'return' in loop bodies)     */
/*          2. No 'switch' statements.                                   */
/*          3. No 'do..while' statements.                                */
/*          4. Expressions are restricted.                               */
/*               - There are no multiple expressions joined by 'or',     */
/*                'and' operations.                                      */
/*          5. No library calls.                                         */
/*               - All the functions needed are implemented in the       */
/*                 source file.                                          */
/*                                                                       */
/*                                                                       */
/*************************************************************************/
/*                                                                       */
/*  FILE: fft1.c                                                         */
/*  SOURCE : Turbo C Programming for Engineering by Hyun Soon Ahn        */
/*                                                                       */
/*  DESCRIPTION :                                                        */
/*                                                                       */
/*     FFT using Cooly-Turkey algorithm.                                 */
/*     There are two inputs, ar[] and ai[]. ar[] is real number parts    */
/*     of input array and the ai[] is imaginary number parts of input.   */
/*     The function fft1 process FFT or inverse FFT according to the    .*/
/*     parameter flag. (FFT with flag=0, inverse FFT with flag=1).       */
/*                                                                       */
/*                                                                       */
/*  REMARK :                                                             */
/*                                                                       */
/*  EXECUTION TIME :                                                     */
/*                                                                       */
/*                                                                       */
/*************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <stdint.h>


using namespace std;

//#define DEBUG
#define PI 3.14159
#define M_PI 3.14159
#define NUM 8


float* fft1(int n, int flag);
int reverse(int n, int flag);


float ar[NUM];
float ai[NUM] = {0.,  };

static float fabs(float n)
{
  float f;

  if (n >= 0) f = n;
  else f = -n;
  return f;
}

static float sin(float rad)
{
  app = (rad-(rad*rad*rad/6)+(rad*rad*rad*rad*rad/120)-(rad*rad*rad*rad*rad*rad*rad/5040)+(rad*rad*rad*rad*rad*rad*rad*rad*rad/362880));
  return(app);
}


static float cos(float rad)
{
	  return (1-(rad*rad/2)+(rad*rad*rad*rad/24)-(rad*rad*rad*rad*rad*rad/720)+(rad*rad*rad*rad*rad*rad*rad*rad/40320));
}


int main()
{


	int  i, n = NUM, flag, chkerr;

	float *temp;

	temp = (float*)malloc(sizeof(float)*(2*NUM));

	for(int p = 0; p < 1000; p++){


	/* ar  */
	for(i = 0; i < n; i++)
	  //ar[i] = cos(2*M_PI*i/n);
		ar[i] = i;

	/* forward fft */
	flag = 0;
	temp = fft1(n, flag);


	for(i=0; i<NUM; i++){
		ai[i] = temp[NUM + i];
		ar[i] = temp[i];
	}


	reverse(n, flag);

	//printf("after reverse 1\n");
	/*for(i = 0; i< NUM; i++){
			printf("%f ", ar[i]);
	}
	printf("\n");
	for(i = 0; i< NUM; i++){
		printf("%f ", ai[i]);
	}
	printf("\n");
	*/


	/* inverse fft */
	flag = 1;
	temp = fft1(n, flag);

	for(i=0; i<NUM; i++){
		ai[i] = temp[NUM + i];
		ar[i] = temp[i];
	}

	reverse(n, flag);



	/*for(i = 0; i< NUM; i++){
			printf("%f ", ar[i]);
	}
	printf("\n");
	for(i = 0; i< NUM; i++){
		printf("%f ", ai[i]);
	}

	printf("\n");
*/
	//printf("check err = %d\n", chkerr);

	}

	return 0;

}



float* fft1(int n, int flag)
{
	int i;
	//printf("in fft1\n");
	float ar_hw[2*NUM];
	float ai_hw[NUM] = {0.,  };

	for(int run = 0; run < 100000; run++){

	for(i=0; i<NUM; i++){
				ar_hw[i] = ar[i];
				ai_hw[i] = ai[i];
			}

	int j, k, it, xp, xp2, j1, j2, iter;
	float sign, w, wr, wi, dr1, dr2, di1, di2, tr, ti, arg;

	//printf("after init\n");

	 iter = 3;
	 j = 1;
#ifdef DEBUG
	printf("iter=%d\n",iter);
	printf("n=%d\n",n);
#endif
	 for(i = 0; i < iter; i++)
	   j *= 2;

	 sign = ((flag == 1) ? 1.0 : -1.0);
	 xp2 = n;
	 for(it = 0; it < iter; it++)
	 {
			 xp = xp2;
			 xp2 /= 2;
			 w = PI / xp2;

			 for(k = 0; k < xp2; k++)
			 {
					 arg = k * w;
					 wr = cos(arg);
					 wi = sign * sin(arg);
					 i = k - xp;
					 for(j = xp; j <= n; j += xp)
					 {
							 j1 = j + i;
							 j2 = j1 + xp2;
							 dr1 = ar_hw[j1];
							 dr2 = ar_hw[j2];
							 di1 = ai_hw[j1];
							 di2 = ai_hw[j2];
							 tr = dr1 - dr2;
							 ti = di1 - di2;
							 ar_hw[j1] = dr1 + dr2;
							 ai_hw[j1] = di1 + di2;
							 ar_hw[j2] = tr * wr - ti * wi;
							 ai_hw[j2] = ti * wr + tr * wi;
					 }
			 }
	 }
	}
	// printf("after computation\n");

 	 for(i = 0; i < NUM; i++){
		 ar_hw[NUM+i] = ai_hw[i];
		 //printf("%f", ai_hw[i]);
	 }

return ar_hw;


}

int reverse(int n, int flag){

	 /*  Digit Reverse Counter  */

	int j1, j2, j, i, k, w;
	float tr, ti;

	j1 = n / 2;
		 j2 = n - 1;
		 j = 1;
	#ifdef DEBUG
		printf("j2=%d\n",j2);
	#endif
		 for(i = 1; i <= j2; i++)
		 {
				 if(i < j)
				 {
						tr = ar[j-1];
						ti = ai[j-1];
						ar[j-1] = ar[i-1];
						ai[j-1] = ai[i-1];
						ar[i-1] = tr;
						ai[i-1] = ti;
				 }
				 k = j1;
				 while(k < j)
				 {
						j -= k;
						k /= 2;
				 }
				 j += k;
		 }
		 if(flag == 0) return(0);
		 w = n;
		 for(i = 0; i < n; i++)
		 {
				 ar[i] /= w;
				 ai[i] /= w;
		 }

		 return 0;
}
