/*
 * compute_ab.c
 *
 * 	See below
 *
 * Copyright (C) 2009 Tom Richardson <tomr@qualcomm.com>
 *
 * Wrapped into matlab: Yury Polyanskiy <polyanskiy@gmail.com>
 */

static char * help_msg = "\n\
 *  function [a b] = compute_ab(k, R) \n\
 *   This function computes the LDPC scaling parameters (a,b) for a given\n\
 *   data size k and rate R. \n\
 *\n\
 *   The usage of (a,b) is: WER(SNR) = Q(a*(10*log10(SNR)-b));\n\
 *   Note: limitations: 250 <= k <= 10000, 25/74 <= R <= 25/32\n";


#include "mex.h"
#include "matrix.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

double Qdata[] = {
1.64741,4.37417,32,25,10,
1.7823,4.39283,32,25,20,
2.73602,4.23319,32,25,40,
4.62843,4.22278,32,25,100,
6.40456,4.11633,32,25,200,
10.5233,4.11227,32,25,400,
1.77178,3.83057,34,25,10,
2.35475,3.72985,34,25,20,
3.12746,3.54719,34,25,40,
4.865,3.48064,34,25,100,
7.83115,3.4784,34,25,200,
11.6305,3.4492,34,25,400,
1.78154,3.10672,36,25,10,
2.2443,3.01614,36,25,20,
3.20881,3.04292,36,25,40,
5.27143,3.00499,36,25,100,
7.72024,2.98124,36,25,200,
10.8001,2.95346,36,25,400,
1.91845,2.70442,38,25,10,
2.29064,2.60278,38,25,20,
3.30915,2.57697,38,25,40,
5.33661,2.55343,38,25,100,
7.38134,2.52523,38,25,200,
10.0935,2.49549,38,25,400,
1.81146,2.31754,40,25,10,
2.34017,2.23708,40,25,20,
3.09574,2.13417,40,25,40,
5.4329,2.0844,40,25,100,
7.88705,2.06081,40,25,200,
10.9385,2.00164,40,25,400,
1.77572,1.91878,42,25,10,
2.39381,1.84634,42,25,20,
3.22597,1.73751,42,25,40,
5.36644,1.71465,42,25,100,
6.78258,1.64867,42,25,200,
10.4735,1.66151,42,25,400,
1.89117,1.59822,44,25,10,
2.53221,1.53728,44,25,20,
3.16796,1.40632,44,25,40,
5.37753,1.4069,44,25,100,
7.75319,1.36131,44,25,200,
11.2761,1.36596,44,25,400,
1.80817,1.23207,46,25,10,
2.49785,1.19559,46,25,20,
3.51756,1.16054,46,25,40,
5.57619,1.10975,46,25,100,
7.7055,1.06553,46,25,200,
10.8983,1.05042,46,25,400,
1.97638,1.03702,48,25,10,
2.38761,0.879989,48,25,20,
3.46179,0.857227,48,25,40,
5.82853,0.825244,48,25,100,
8.09585,0.783093,48,25,200,
12.0788,0.770085,48,25,400,
1.86676,0.754931,50,25,10,
2.58492,0.689422,50,25,20,
3.48458,0.608158,50,25,40,
5.99078,0.534573,50,25,100,
8.0175,0.488617,50,25,200,
11.8174,0.469132,50,25,400,
1.81518,0.487696,52,25,10,
2.43107,0.385302,52,25,20,
3.28604,0.258776,52,25,40,
5.51571,0.224051,52,25,100,
8.90916,0.216781,52,25,200,
12.5695,0.1848,52,25,400,
1.82939,0.236548,54,25,10,
2.51242,0.145907,54,25,20,
3.44233,0.0442195,54,25,40,
5.75171,0.0109963,54,25,100,
8.52972,-0.0375244,54,25,200,
12.5976,-0.0521179,54,25,400,
1.79078,-0.0217274,56,25,10,
2.47329,-0.0930272,56,25,20,
3.49786,-0.183913,56,25,40,
5.91739,-0.231223,56,25,100,
8.39529,-0.262209,56,25,200,
12.5628,-0.28216,56,25,400,
1.88899,-0.156296,58,25,10,
2.46691,-0.304364,58,25,20,
3.69912,-0.3562,58,25,40,
6.40493,-0.4022,58,25,100,
8.37809,-0.476277,58,25,200,
11.9779,-0.507924,58,25,400,
1.72595,-0.419539,60,25,10,
2.5904,-0.472487,60,25,20,
3.63603,-0.5761,60,25,40,
5.93633,-0.647782,60,25,100,
8.68311,-0.681829,60,25,200,
13.1691,-0.704235,60,25,400,
2.00796,-0.469754,62,25,10,
2.05649,-0.756879,62,25,20,
3.4608,-0.761793,62,25,40,
6.09034,-0.850691,62,25,100,
8.8754,-0.867034,62,25,200,
12.9577,-0.903798,62,25,400,
1.86804,-0.747658,64,25,10,
2.70444,-0.831951,64,25,20,
3.29555,-1.00959,64,25,40,
5.85237,-1.03588,64,25,100,
8.4459,-1.08069,64,25,200,
12.8935,-1.10355,64,25,400,
1.89483,-0.943391,66,25,10,
2.54966,-1.04475,66,25,20,
3.67842,-1.11365,66,25,40,
5.89775,-1.2075,66,25,100,
9.54893,-1.22467,66,25,200,
14.1153,-1.25599,66,25,400,
1.96293,-1.04485,68,25,10,
2.60835,-1.20094,68,25,20,
3.77648,-1.26874,68,25,40,
6.55415,-1.35325,68,25,100,
9.17148,-1.39756,68,25,200,
13.8715,-1.4355,68,25,400,
1.90733,-1.22753,70,25,10,
2.35522,-1.38685,70,25,20,
3.68205,-1.44182,70,25,40,
6.29001,-1.52611,70,25,100,
9.16421,-1.57885,70,25,200,
14.2846,-1.59307,70,25,400,
1.9763,-1.31647,72,25,10,
2.51301,-1.52869,72,25,20,
3.67062,-1.58159,72,25,40,
6.4886,-1.64355,72,25,100,
9.39234,-1.71058,72,25,200,
13.2467,-1.73076,72,25,400,
1.85956,-1.53614,74,25,10,
2.66185,-1.66955,74,25,20,
3.44503,-1.78944,74,25,40,
6.16869,-1.80959,74,25,100,
9.66082,-1.85092,74,25,200,
11.9337,-1.8857,74,25,400,
};

/*
 * fit parameters to get estimate of Q function parameters
 * for given info block lenght and rate
 * must be in range  10*25 < k < 400*25,  25/74 < r < 25/32
 */
int getQcoeff(double infoblockLength,double rate,double *a,double *b){
	int nest, k1;
	double asample1,asample2,lsample1,lsample2,a1scaled, a2scaled, alpha, alb, blb, aub, bub;
	int nSampleLength = 6; /* number of lengths sampled  */

	if(rate>=25.0/32.0) return -1;
	if(rate<=25.0/74.0) return -1;
	
	if(infoblockLength<10*25.0) return -1;
	if(infoblockLength>400*25.0) return -1;

	nest = (int)floor(25.0/rate);
	nest-=nest%2;


	/* linear interpolate parameters */
	
	k1 = nSampleLength*((nest-32)/2); 
	if(infoblockLength>20*25) k1++; 
	if(infoblockLength>40*25) k1++; 
	if(infoblockLength>100*25) k1++; 
	if(infoblockLength>200*25) k1++;
	/* increment according to sampled lengths */
	
	/* get parameters by interpolation (with scaling for a) */
	asample1 = Qdata[5*k1];     lsample1 = Qdata[5*k1+4];
	asample2 = Qdata[5*(k1+1)]; lsample2 = Qdata[5*(k1+1)+4];
	a1scaled = asample1/sqrt(lsample1);
	a2scaled = asample2/sqrt(lsample2);
	/* interpolate linearly */
	alpha = ((infoblockLength/25.0) - lsample1)/(lsample2-lsample1);
	alb = (alpha*a2scaled + (1-alpha)*a1scaled) *sqrt(infoblockLength/25.0);
	blb = alpha*Qdata[5*(k1+1)+1] + (1-alpha)*Qdata[5*(k1)+1];

	/* now upper length bounds */
	k1+=nSampleLength;

	asample1 = Qdata[5*k1];     lsample1 = Qdata[5*k1+4];
	asample2 = Qdata[5*(k1+1)]; lsample2 = Qdata[5*(k1+1)+4];
	a1scaled = asample1/sqrt(lsample1);
	a2scaled = asample2/sqrt(lsample2);
	/* interpolate linearly */
	alpha = ((infoblockLength/25.0) - lsample1)/(lsample2-lsample1);
	aub = (alpha*a2scaled + (1-alpha)*a1scaled) *sqrt(infoblockLength/25.0);
	bub = alpha*Qdata[5*(k1+1)+1] + (1-alpha)*Qdata[5*(k1)+1];

	/* now interpolate over rate */
	alpha = (rate - 25.0/(nest+2)) /( (25.0/nest) - (25.0/(nest+2)) );
	*a = alb*alpha+(1-alpha)*aub;
	*b = blb*alpha+(1-alpha)*bub;
	return 0;
}

/* simple Q function approximation (3% error bound) */
double Q(double x){
	if(x<0) return 1.0-Q(-x);
	double xsq=x*x;
	return exp(-xsq/2.0)/(1.64*x+sqrt(0.76*xsq+4));
}

/* estimated FER given a,b,snr */
double Qest(double a,double b,double snr){
	return Q(a*(snr-b));
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double k;
	double R = 1/2;
	double a,b;
	
	if( (nrhs > 2) || !nrhs) goto usage;

	/*
	 * The first param is integer, but matlab works with doubles by default
	 */
	if(!mxIsDouble(prhs[0])) goto nondouble;
	
	if(nrhs == 2){
		if(!mxIsDouble(prhs[1])) goto nondouble;
		R = *mxGetPr(prhs[1]);
	}
	k = *mxGetPr(prhs[0]);

#if DEBUG>1
	mexPrintf(" +++ compute_ab: nlhs = %d, k = %g,  R= %.3g ", nlhs, k, R);
#endif

	if(getQcoeff(k,R,&a,&b)){
		/* 
		 * Note: my code relies on compute_ab() returning NAN for
		 * unsupported rates
		 */
		mexPrintf("Warning: getQcoeff() failed. k,R out of range?\n");
		a = b = NAN;
	}

	if(nlhs == 0){
		int i;
		mexPrintf(" +++ compute_ab: a = %.6g, b = %.6g\n", a, b);
		for(i=0;i<50;i++)
			mexPrintf("      SNR = %g dB:\t\tWER %g\n",-2.0+0.2*i,Qest(a,b,-2.0+0.2*i));
	}
	if(nlhs > 0){
		plhs[0] = mxCreateNumericMatrix(1,1, mxDOUBLE_CLASS, mxREAL);
		if(!plhs[0]) goto oom;
		*mxGetPr(plhs[0]) = a;
	}
	if(nlhs > 1){
		plhs[1] = mxCreateNumericMatrix(1,1, mxDOUBLE_CLASS, mxREAL);
		if(!plhs[1]) goto oom;
		*mxGetPr(plhs[1]) = b;
	}
	return;
oom:
	mexPrintf("ERR: OOM\n");
	return;

nondouble:
	mexPrintf("ERR: nondouble arguments\n");
	goto helpmsg;
usage:
	mexPrintf("ERR: usage\n");
helpmsg:
	mexPrintf("%s\n", help_msg);
	return;
}
