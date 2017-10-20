
#pragma OPENCL EXTENSION cl_amd_fp64 : enable

#define real double
#define real_in float

#warning CL_VERSION_1_2

#ifdef CL_VERSION_2_0   
#define HASCL
#endif

#ifdef CL_VERSION_1_2  
   #define HASCL
#endif

#ifdef CL_VERSION_1_1  
#define HASCL
#endif


  #ifndef  HASCL
    #define global
    #define local
    #define kernel

    #include <math.h>
    #include <stdio.h>
    #include <stdlib.h>

    #ifndef max
    #define max(a,b) ((a) > (b) ? (a) : (b))
    #endif

#ifndef min
#define min(a,b) ((a) < (b) ? (a) : (b))
#endif


  #endif


#define Integra( r, ff , x, x1 , x2 , dx )  for(x = (x1) ; x <= (x2) ;   ){ r = r + (dx) * ( ff )/6.0 ; x+=(dx)/2.0 ; r = r + 4*(dx) * ( ff )/6.0; x+=(dx)/2.0 ;r = r + (dx) * ( ff )/6.0; };



typedef struct ModalParameterIc
{
	real a, x, s;
} ModalParameterIc ;

typedef struct ModalParameter
{
	ModalParameterIc p[4];
	real result;
} ModalParameter;


#define H (7.0)

real_in gaussian_integrate_range(real_in x, real_in s, real_in x1, real_in x2)
{
	real_in w1 = 0.5*erf(0.707107 * (x1 - x) / s);
	real_in w2 = 0.5*erf(0.707107 * (x2 - x) / s);
	return (w2 - w1);

}

real_in gaussian(real_in  x, real_in x0, real_in sigma)
{
	if (fabs(x - x0) < sigma * H)
		return ( ((real_in)1.0) / (sigma* sqrt(2.0 * 3.1415928))) * exp(-0.5*pow((x - x0) / sigma, ((real_in)2.0)));
	return 0.0f; 
}

real_in cov_gauss_num(real_in x1, real_in q1, real_in x2, real_in q2)
{

	if (q2 < 0.1*q1) return gaussian(x2, x1, q1);
	else if (q1 < 0.1*q2) return gaussian(x1, x2, q2);


	real_in r = 1.0;
	real_in u = 0;

	real_in a1 = max(x1 - 5*q1, x2 - 5*q2);
	real_in a2 = min(x1 + 5*q1, x2 + 5*q2);
	if (a2<=a1) return 0.0;

	real_in dx = (a2 - a1) / 8.0;
	if (dx < 1e-6) return 0;
	//Integra(r, (gaussian(u,x1,q1)*gaussian(u, x2, q2)) , u, a1, a2, dx);
	Integra(r, ( exp(-0.5*pow((u - x1) / q1, ((real_in)2.0))) * exp(-0.5*pow((u - x2) / q2, ((real_in)2.0)))), u, a1, a2, dx);

	r = r *  (((real_in)1.0) / (sqrt(2.0 * 3.1415928* q1*q1)));
	r = r *  (((real_in)1.0) / (sqrt(2.0 * 3.1415928* q2*q2)));

	return r;
}


double conv_math(double x1,double s1,double x2,double s2, double a,double b)
{
	 
	double c1 = (s1*s2) / exp( pow(x1 - x2, 2.0) /	(2.0*(pow(s1, 2.0) + pow(s2, 2.0))));
	double c2 = erf((pow(s2, 2.0)*(a - x1) + pow(s1, 2.0)*(a - x2)) / (sqrt(2.0)*s1*s2*sqrt(pow(s1, 2.0) + pow(s2, 2.0))));
	double c3 = -erf((pow(s2, 2.0)*(b - x1) + pow(s1, 2.0)*(b - x2)) /	(sqrt(2.0)*s1*s2*sqrt(pow(s1, 2.0) + pow(s2, 2.0))));
	double c4 = 2.0 * sqrt(2.0 * M_PI)*sqrt(pow(s1, 2.0))*sqrt(pow(s2, 2.0))*	sqrt(pow(s1, 2.0) + pow(s2, 2.0));

	return -c1*(c2 + c3) / c4;
}
real_in cov_gauss(real_in x1, real_in q1, real_in x2, real_in q2)
{
	//real_in a = cov_gauss_num(x1, q1, x2, q2);
	real_in b = conv_math(x1, q1, x2, q2,-10,10);
	//printf("[%f %f %f %f ] diff %g %g \n", x1, q1, x2, q2, a, b);
	return b;
	if (q2 < 1e-12) 
	{
		// delta dirac ?
		if (x1 == x2) return gaussian(x2, x1, q1);
		return 1e-12;
	}

	//if (q2 < 0.05*q1) return gaussian(x2, x1, q1);
	//else if (q1 < 0.05*q2) return gaussian(x1, x2, q2);
	

	// calcula a covolucao de duas gaussianas
	real_in a1, a2;
	real_in ss;
	//float u;
	a1 = max(x1 - H*q1, x2 - H*q2);
	a2 = min(x1 + H*q1, x2 + H*q2);
	if (a2<a1) return 0.0;


	real_in f1 = (exp(-pow(x1 - x2,  ((real_in)2.0)) / (2.0*(pow(q1, ((real_in)2.0)) + pow(q2,  ((real_in)2.0))))));
	real_in f2 = sqrt((1.0f / (q1*q1)) + (1.0f / (q2*q2)));
	real_in f3 = 15.7496*q1*q2 * f2;

	double gsum = f1 / (f3);
	ss = gsum;
	 
	
	return ss;
}




real ModalLikehod_i(ModalParameterIc mip,   real m1,   real m2 )
{
    if ( mip.x + H * mip.s < m1) return 0.0;
	if ( mip.x - H * mip.s > m2) return 0.0;
	return  mip.a *gaussian_integrate_range(mip.x, mip.s, m1, m2);
}

real_in ModalLikehod_i_np(real_in x1, real_in q1, real_in x2, real_in q2)
{



	return gaussian(x2, x1, q1);
	 
	//real_in p1 = x2 - H * q2;
	//real_in p5 = x2 + H * q2;
 
	real_in dx = 0.05;

	if (fabs(x2 - x1) > 3*q1) return 0.0;
	if (fabs(x2 - x1) > q1) return 0.16;
	return 0.68;
	
	
	return  0.5*(1.0 - fabs(erf(0.707107 *  (x2 - x1) / q1))) ;




	real_in p1 = x2 - dx;
	real_in p5 = x2 + dx;	
	return  0.5*(erf(0.707107 * (p5 - x1) / q1) - erf(0.707107 * (p1 - x1) / q1))   ;


	//if (q2 < q1*0.1) return 0.5*fabs(erf(0.707107 * (p5 - x1) / q1) - erf(0.707107 * (p1 - x1) / q1));

	//if (p1 < x1 - H* q1) return 0.0;
	//if (p5 > x1 + H* q1) return 0.0;


	real_in p2 = x2 - q2;
	//real_in p3 = x2 ;
	real_in p4 = x2 + q2;


	real_in w1 = erf(0.707107 * (p1 - x1) / q1);
	real_in w2 = erf(0.707107 * (p2 - x1) / q1);
	//real_in w3 = 0.5*erf(0.707107 * (p3 - x1) / q1);
	real_in w4 = erf(0.707107 * (p4 - x1) / q1);
	real_in w5 = erf(0.707107 * (p5 - x1) / q1);

	//return 0.16*(w2 - w1) + 0.34*(w3 - w2) + 0.34*(w4 - w3) + 0.16*(w5 - w4);
	return 0.5* (0.16*(w2 - w1) + 0.68*(w4 - w2) + 0.16*(w5 - w4));



	//return gaussian_integrate_range(x1, q1, x2 - q2, x2 + q2);

	//return gaussian(x2, x1, q1);

	//if (q2 < 0.1*q1)
	//{
	//	return gaussian(x2, x1, q1);
	//}
	//
	//real_in a1 = max(x1 - 5 * q1, x2 - 5 * q2);
	//real_in a2 = min(x1 + 5 * q1, x2 + 5 * q2);
	//if (a2 <= a1) return 0.0;

	////gaussian_integrate_range(x1, q1, x2 - 2*q2, x2-q2)
	////gaussian_integrate_range(x1, q1, x2 - q2, x2 )
	////gaussian_integrate_range(x1, q1, x2     , x2 + q2)
	////gaussian_integrate_range(x1, q1, x2 + q2, x2 + 2*q2)
	//return 0;

}

real_in integra_np(ModalParameter mp, const int modal_order, real_in u1, real_in u2)
{
	double uu;
	double du = (u2 - u1) / 50;
	double logProd = 0.0;
	real_in noise = 1e-10; 

	real_in term = 0.0;


 


	for (real_in u = u1; u <= u2; u += du)
	{
		
		if (modal_order >= 0) term += mp.p[0].a* ModalLikehod_i_np(mp.p[0].x, mp.p[0].s, u, du);
		if (modal_order >= 1) term += mp.p[1].a* ModalLikehod_i_np(mp.p[1].x, mp.p[1].s, u, du);
		if (modal_order >= 2) term += mp.p[2].a* ModalLikehod_i_np(mp.p[2].x, mp.p[2].s, u, du);
		if (modal_order >= 3) term += mp.p[3].a* ModalLikehod_i_np(mp.p[3].x, mp.p[3].s, u, du);

		//if (modal_order >= 0) term += du* mp.p[0].a* gaussian(u, mp.p[0].x, mp.p[0].s);
		//if (modal_order >= 1) term += du* mp.p[1].a* gaussian(u, mp.p[1].x, mp.p[1].s);
		//if (modal_order >= 2) term += du* mp.p[2].a* gaussian(u, mp.p[2].x, mp.p[2].s);
		//if (modal_order >= 3) term += du* mp.p[3].a* gaussian(u, mp.p[3].x, mp.p[3].s);


		//logProd = logProd + log(max(noise , du*term ));
		//if (logProd < -99999) return -99999;

	}
	return term/50.0;

	return logProd;
}



real ModalLikehod_h(ModalParameter mp, global real *mass_h, global real *count_h, const int num_h, const int modal_order)
{
	real termo1 = 1.0 ;
	real logProd = 0.0;
	real noise = 1e-10;

	termo1 =  integra_np(mp, modal_order, -1.0, 3.0);
	for (int i = 0; i < num_h-1; ++i)
	{
		real m1 = mass_h[i];
		real m2 = mass_h[i+1];
		real c = count_h[i];
		real termo2 = 0.0;

		if (modal_order >= 0)  termo2 += c * ModalLikehod_i(mp.p[0], m1, m2);
		if (modal_order >= 1)  termo2 += c * ModalLikehod_i(mp.p[1], m1, m2);
		if (modal_order >= 2)  termo2 += c * ModalLikehod_i(mp.p[2], m1, m2);
		if (modal_order >= 3)  termo2 += c * ModalLikehod_i(mp.p[3], m1, m2);


		logProd = logProd + log(noise + termo2);

		//if (isnan(logProd) || !(isfinite(logProd) || !(isfinite(termo2)) || isnan(termo2)))
		//{
		//	printf("@2 Nan %g %g  \n", logProd, termo2);

		//}
	}
	return     (-termo1 + logProd);
}


real np(real  uu, ModalParameter mp, const int modal_order)
{
	double result = 0.0;
	//if ((a + gg)>1) return 1e-40;
	if (modal_order >= 0) result += mp.p[0].a * gaussian(uu, mp.p[0].x, mp.p[0].s);
	if (modal_order >= 1) result += mp.p[1].a * gaussian(uu, mp.p[1].x, mp.p[1].s);
	if (modal_order >= 2) result += mp.p[2].a * gaussian(uu, mp.p[2].x, mp.p[2].s);
	if (modal_order >= 3) result += mp.p[3].a * gaussian(uu, mp.p[3].x, mp.p[3].s);
	return 1e-12 + result;
}

real_in prop_relative(ModalParameter mp, const int modal_order, real_in u1, real_in u2)
{

	real_in result = 0.0;

	if (modal_order >= 0) result += mp.p[0].a * 4 * mp.p[0].s;
	if (modal_order >= 1) result += mp.p[1].a * 4 * mp.p[1].s;
	if (modal_order >= 2) result += mp.p[2].a * 4 * mp.p[2].s;
	if (modal_order >= 3) result += mp.p[3].a * 4 * mp.p[3].s;
	return result;


	//real_in result = 0.0;
	//if (modal_order >= 0) result += mp.p[0].a * 2 * mp.p[0].s;
	//if (modal_order >= 1) result += mp.p[1].a * 2 * mp.p[1].s;
	//if (modal_order >= 2) result += mp.p[2].a * 2 * mp.p[2].s;
	//if (modal_order >= 3) result += mp.p[3].a * 2 * mp.p[3].s;
	//return result;
}

  real_in ModalLikehod_np(bool is_valid,    ModalParameter mp, global real *mass_v, global real *sigma_v,   const int num_v, const int modal_order)
{ 
	 
	// double prod;
	real_in soma, termo1,   u;
	soma = 0.0; 
	termo1 = 0.0; 
	real_in logProd = 0.0;
	
	real_in prod = 1.0;
	real_in noise = 1e-99;
	real_in log_noise = -9999;

	logProd = 0.0;

	termo1 = 1.0; // prop_relative(mp, modal_order, 0.0, 2.0);
	//termo1 = integra_np(mp, modal_order, 0.0, 2.0);
	
	real_in log_term1 = 0.0;

	if (termo1 < 0.000001) return -99999;
	 

	//termo1 = 1.0;
	//if (isnan(termo1) || !(isfinite(termo1) || !(isfinite(termo1)) || isnan(termo1)))
	//	{
	//		printf("@2 Nan %g  \n", termo1);
	//	}
	 
	
	 
	for (int i = 0; i < num_v; i+=1)
	{
	//	int j = get_global_id(1);
	//	int j0 = get_local_id(1);

	//	int i = j;


		 

		real_in mi = mass_v[i];
		real_in si = sigma_v[i];

		//if (mi < 0.001) continue;
		//printf("@3  %d %d  \n", get_global_id(0),i);

		//if (i < num_v && is_valid)
		
		if (mi > 0.001)			
		{
			si = max(((real_in)0.001f), si);
			{
				//if (prod<1e-20) break;
				real_in termo2 = 0.0;
				//if (modal_order >= 0) termo2 += mp.p[0].a * cov_gauss(mp.p[0].x, mp.p[0].s, mi, si );
				//if (modal_order >= 1) termo2 += mp.p[1].a * cov_gauss(mp.p[1].x, mp.p[1].s, mi, si );
				//if (modal_order >= 2) termo2 += mp.p[2].a * cov_gauss(mp.p[2].x, mp.p[2].s, mi, si );
				//if (modal_order >= 3) termo2 += mp.p[3].a * cov_gauss(mp.p[3].x, mp.p[3].s, mi, si );

 
                 termo2 += mp.p[0].a * ModalLikehod_i_np( mp.p[0].x, mp.p[0].s, mi ,si);

#ifdef MODAL1
				 termo2 += mp.p[1].a * ModalLikehod_i_np( mp.p[1].x, mp.p[1].s, mi, si);
#endif

#ifdef MODAL2
				 termo2 += mp.p[2].a * ModalLikehod_i_np( mp.p[2].x, mp.p[2].s, mi, si);  
#endif

#ifdef MODAL3
				 termo2 += mp.p[3].a * ModalLikehod_i_np( mp.p[3].x, mp.p[3].s, mi, si);
#endif

				 //if (termo2 > termo1) printf("@LL %g %g \n", termo2, termo1);

				 
				 logProd = logProd + log(max( noise, termo2));
			 

				
			}
		}
	}


	//barrier(CLK_LOCAL_MEM_FENCE);
	//if (is_valid)
	//{
	//	if (j0 == 0)
	//	{
	//		for (int j = 0; j < get_local_size(1); ++j)
	//		{
	//			logProd += logProd_local[j];
	//		}
	//	}
	//}
	 
	//if (isnan(termo1) || !(isfinite(termo1) || !(isfinite(logProd)) || isnan(logProd)))
	//{
	//	printf("@2 Nan %g %g \n", termo1, logProd);
	//}
 
return     (   logProd );


}

 kernel void ModalLikehodList_h(global const ModalParameter *inputParams, global real *output,   const int nArgs  , global real *mass_h, global real *count_h, const int num_h, const int modal_order )
{
	
	int tid = get_global_id(0);
	//if ( tid ==0) printf((__constant char *)" id \n");
	__local  real_in logProd_local[64];
	if (tid < nArgs)
	{
		real res = ModalLikehod_h(  inputParams[tid] , mass_h, mass_h, num_h, modal_order);
		output[tid] = res;
	}
}


 kernel void ModalLikehodList(global const ModalParameter *inputParams, global real *output, const int nArgs, global real *mass_v, global real *sigma_v, const int num_v,  const int modal_order)
 {

	 int tid = get_global_id(0);
 
	 //if ( tid ==0) printf((__constant char *)" id \n");
	 
	 if (tid < nArgs) 
	 {
		 real res = ModalLikehod_np(tid < nArgs ,   inputParams[tid], mass_v, sigma_v, num_v, modal_order);
		 output[tid] = res;
	 }
	 else
	 {
		 output[tid] = -99999;
	 }
 }

 


 kernel void makeHistogram( global real *mass_v, global real *sigma_v, const int num_v, global real *mass_h, global real *count_h, const int num_h )
 {

	 int tid = get_global_id(0);	 
	 if (tid < num_v)
	 {
		 for (int i = 0; i < num_h -1 ; ++i)
		 {
			 real n = gaussian_integrate_range(mass_v[tid], sigma_v[tid], mass_h[i], mass_h[i + 1]);
			 count_h[i] += n;
		 }
	 }
 }

kernel void vectorAdd(global const real *inputA, global const real *inputB, global real *output, const real x)
{
	 
	output[get_global_id(0)] = inputA[get_global_id(0)] + inputB[get_global_id(0)] + x;
}