// WhiteDwarfV6.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <vector>
#include <algorithm>
#include <functional>
#include "WhiteDwarfV6.h"

#define real double

typedef struct ModalParameterIc
{
	real a, x, s;
	ModalParameterIc(): a(0), x(0), s(1.0)
	{}
}ModalParameterIc;


 

 

typedef struct ModalParameter
{
	//std::vector<ModalParameterIc> p_vector;
	ModalParameterIc p[4];
	real result;
	ModalParameter(): result(0)
	{}
}ModalParameter;


class ParameterValue
{
	
};

class ParameterRange:ParameterValue
{
public:
	float x1, x2, dx;
	ParameterRange(float _x1, float _x2, int  _n): x1(_x1), x2(_x2), dx(fabs(_x2 - _x1) / _n){ } 
	ParameterRange(double _x1, double _x2, int  _n) : x1(_x1), x2(_x2), dx( fabs(_x2-_x1)/_n ) {}
};




typedef struct ModalParameterRange
{
	ParameterRange a, x, s;
	ModalParameterRange(ParameterRange _a , ParameterRange _x, ParameterRange _s) : a(_a), x(_x), s(_s)
	{}

	void expand(const std::function<void(float , float, float )>& f  ) const
	{
		for( float _a = a.x2; _a>= a.x1; _a -= a.dx  )
			for (float _x = x.x1; _x <= x.x2; _x += x.dx)
				for (float _s = s.x1; _s <= s.x2; _s += s.dx)
				{
					f(_a, _x, _s);
				}
	}

	 ModalParameterRange  refine_around(ModalParameterIc p) const
	 {
		 double a_1 = 0.5*(a.x1 + p.a);
		 double a_2 = 0.5*(a.x2 + p.a);

		 double x_1 = 0.5*(x.x1 + p.x);
		 double x_2 = 0.5*(x.x2 + p.x);

		 double s_1 = 0.5*(s.x1 + p.s);
		 double s_2 = 0.5*(s.x2 + p.s);

		 int na = (a.x2 - a.x1) / a.dx;
		 int nx = (x.x2 - x.x1) / x.dx;
		 int ns = (s.x2 - s.x1) / s.dx;

		 na = std::min(na, 5);
		 nx = std::min(nx, 5);
		 ns = std::min(ns, 5);

		 if (na % 2 == 1) na++;
		 if (nx % 2 == 1) nx++;
		 if (ns % 2 == 1) ns++;
			 

		 printf("%f %f  [%f]  -> %f %f \n", a.x1, a.x2, p.a, a_1, a_2);
		 printf("%f %f  [%f]  -> %f %f \n", x.x1, x.x2, p.x, x_1, x_2);
		 printf("%f %f  [%f]  -> %f %f \n", s.x1, s.x2, p.s, s_1, s_2);
		 printf("\n");
		 return ModalParameterRange({ a_1,a_2, na }, { x_1,x_2, nx }, { s_1,s_2, ns });

	}
}ModalParameterRange;



typedef struct ModalParameterRegion
{
public:
	std::vector<ModalParameterRange> p;	
	ModalParameterRegion(ModalParameterRange m1) : p({ m1 }) 	{}	
	ModalParameterRegion(ModalParameterRange m1, ModalParameterRange m2) : p({ m1,m2 })  {}
	ModalParameterRegion(ModalParameterRange m1, ModalParameterRange m2, ModalParameterRange m3) : p({ m1,m2,m3 })  {}
	ModalParameterRegion(ModalParameterRange m1, ModalParameterRange m2, ModalParameterRange m3, ModalParameterRange m4) : p({ m1,m2,m3,m4 })  {}

}ModalParameterRegion;

void dump_region(ModalParameterRegion r)
{
	for (auto &g : r.p)
	{
		printf("%f %f [%f] ", g.a.x1, g.a.x2, g.a.dx);
		printf("%f %f [%f] ", g.x.x1, g.x.x2, g.x.dx);
		printf("%f %f [%f]\n", g.s.x1, g.s.x2, g.s.dx);
	}
}

ModalParameterRegion refine_around(ModalParameterRegion base, ModalParameter center)
{
	int modal_order = base.p.size();

	ModalParameterRange m1 = base.p[0].refine_around(center.p[0]);

	if (modal_order > 1)
	{
		ModalParameterRange m2 = base.p[1].refine_around(center.p[1]);

		if (modal_order > 2)
		{

			ModalParameterRange m3 = base.p[2].refine_around(center.p[2]);

			if (modal_order > 3)
			{
				ModalParameterRange m4 = base.p[3].refine_around(center.p[3]);
				return ModalParameterRegion(m1, m2, m3, m4);
			}
			return ModalParameterRegion(m1, m2, m3);
		}
		return ModalParameterRegion(m1, m2);
	}
	return ModalParameterRegion(m1);

}


//=======================================================

void dump( const ModalParameter& mp)
{
	printf("[%f %f %f] ", mp.p[0].a, mp.p[0].x, mp.p[0].s);
	if (mp.p[1].a>0) printf("[%f %f %f] ", mp.p[1].a, mp.p[1].x, mp.p[1].s);
	if (mp.p[2].a>0)printf("[%f %f %f] ", mp.p[2].a, mp.p[2].x, mp.p[2].s);
	if (mp.p[3].a>0)printf("[%f %f %f] ", mp.p[3].a, mp.p[3].x, mp.p[3].s);
	printf("= %g \n", mp.result );
}
void dump_f( FILE* f, const ModalParameter& mp)
{
	fprintf(f," %2.5f %2.5f %2.5f ",  mp.p[0].a, mp.p[0].x, mp.p[0].s);
	if (mp.p[1].a > 0) fprintf(f, "  %2.5f %2.5f %2.5f ", mp.p[1].a, mp.p[1].x, mp.p[1].s);
	if (mp.p[2].a > 0)fprintf(f, "  %2.5f %2.5f %2.5f ", mp.p[2].a, mp.p[2].x, mp.p[2].s);
	if (mp.p[3].a > 0)fprintf(f, "  %2.5f %2.5f %2.5f ", mp.p[3].a, mp.p[3].x, mp.p[3].s);
	fprintf(f,"  %g \n", mp.result);
}

void computeParams(std::vector<ModalParameter>& points, int modalOrder);


bool same_covergence(const ModalParameter& mp, int modal_order)
{
	for(int i = 0;i  <= modal_order-1;++i)
		for (int j = i+1 ; j <= modal_order; ++j)
		{
			if (mp.p[i].x == mp.p[j].x) 
				//if (mp.p[i].s == mp.p[j].s)
				{
					return  true;
				}
		}

	return false;
}

struct ModalParameterBuffer
{
	std::vector<ModalParameter>  v_parameters_to_process;
	std::vector<ModalParameter>  v_parameters_filled;

	ModalParameter best_parameter;
	double  max_likehood = -99999;
	
	void save_progress_log()
	{
		FILE* f = fopen("progress.dat", "a");
		for (auto &mpj : v_parameters_filled)
		{
			if (mpj.result > max_likehood - 5)
			{
				dump_f(f, mpj);
				if (mpj.result > max_likehood)
				{
					max_likehood = std::max(mpj.result, max_likehood);
					best_parameter = mpj;
				}
			}
		}
		fclose(f);
	}


	void save_likehood(const char* detail_filename)
	{
		FILE* f = fopen(detail_filename, "a");
		for (auto &mpj : v_parameters_filled)
		{
				dump_f(f, mpj);
		}
		fclose(f);
	}

	void process_buffer(int modal_order, const char* detail_filename ,bool force_out =false )
	{
		 
			computeParams(v_parameters_to_process, modal_order);
			v_parameters_filled.insert(v_parameters_filled.end(), v_parameters_to_process.begin(), v_parameters_to_process.end());
			v_parameters_to_process.clear();

			 
			if(v_parameters_filled.size()> 1024 * 200 || force_out)
			{
				 
				save_progress_log( );

				if (detail_filename!=nullptr) save_likehood(detail_filename);


				v_parameters_filled.clear();
			}

	}

	void add(ModalParameter mp, int modal_order , const char* detail_filename   )
	{
		//if (mp.p[0].x < 0.09) return;
		//if (mp.p[0].x > 0.20) return;

		//if (mp.p[0].a < 0.6) return;
		//if (mp.p[0].s > 0.10) return;

		//if (mp.p[1].s < 0.20) return;

		 

		if (modal_order >= 1) if (fabs(mp.p[1].x - mp.p[0].x) <0.001) return;
		if (modal_order >= 2) if (fabs(mp.p[2].x - mp.p[0].x) <0.001) return;
		if (modal_order >= 2) if (fabs(mp.p[2].x - mp.p[1].x) <0.001) return;

		if (modal_order >= 3) if (fabs(mp.p[3].x - mp.p[0].x) <0.001) return;
		if (modal_order >= 3) if (fabs(mp.p[3].x - mp.p[1].x) <0.001) return;
		if (modal_order >= 3) if (fabs(mp.p[3].x - mp.p[2].x) <0.001) return;


		

		if (modal_order >= 1) if ( mp.p[1].a > mp.p[0].a) return;

		if (modal_order >= 2) if (mp.p[2].a > mp.p[0].a) return;
		if (modal_order >= 2) if (mp.p[2].a > mp.p[1].a) return;

		if (modal_order >= 3) if (mp.p[3].a > mp.p[0].a) return;
		if (modal_order >= 3) if (mp.p[3].a > mp.p[1].a) return;
		if (modal_order >= 3) if (mp.p[3].a > mp.p[2].a) return;
		
		//normalize
		real aTotal = 0.0;
		for (int i = 0; i <= modal_order; ++i) aTotal += mp.p[i].a;

		
		if (aTotal < FLT_EPSILON) return;

		 
		for (int i = 0; i <= modal_order; ++i) mp.p[i].a = mp.p[i].a / aTotal;		
	 

		v_parameters_to_process.push_back(mp);
		//dump(mp); 

		//if (int64_t(sizeof(ModalParameter))  * int64_t(v_parameters_to_process.size()) > int64_t(100 * 1024 * 1024)) // 100Mb per batch
		if ( int64_t(v_parameters_to_process.size()) > 1024*10) // 100Mb per batch
		{
			process_buffer(modal_order, detail_filename);
		}

		if ( int64_t(sizeof(ModalParameter))  * int64_t( v_parameters_filled.size()) > int64_t(1 * 1024) * int64_t(1024 * 1024) ) // 2GB max
		{
			printf("too many params num points to compute \n");
			exit(0);
		}

	}

	std::vector<ModalParameter>  get_results(int modal_order, const char* detail_filename)
	{
		if (v_parameters_to_process.empty() == false)
		{
			process_buffer(modal_order, detail_filename, true);
		}
		return v_parameters_filled;
	}

	ModalParameter get_best_parameter(const int modal_order, const char* detail_filename)
	{
		if (v_parameters_to_process.empty() == false)
		{
			process_buffer(modal_order,detail_filename, true);
		}
		return best_parameter;
	}
};
 
struct Parameters
{
	real m_min; 
	real m_max;
	real dm;

	real s_min;
	real s_max;
	real ds;

	real da;
};

 

void fill_range(ModalParameterBuffer  &v_modalpoints, ModalParameter  mp, int modalOrder, int orderMax, real amax,  real da, const Parameters& p, const char* detail_filename)
{ 
	for (real a = amax; a >= 0 ; a -= p.da)
	{
		if (a < 0) return;
		if ((modalOrder > 0) && (mp.p[modalOrder - 1].a < a)) return; //nao prossege pois ha o valor anterior eh maior que o atual

		if (modalOrder == orderMax)
		{
			for (real s = p.s_min; s <= p.s_max; s += p.ds)
				for (real m = p.m_min; m <= p.m_max; m += p.dm)
				{
					mp.p[modalOrder].a = a;
					mp.p[modalOrder].x = m;
					mp.p[modalOrder].s = s*s;
					v_modalpoints.add(mp, orderMax, detail_filename);
				}
		}
		else
		{


			for (real s = p.s_min; s <= p.s_max; s += p.ds)
				for (real m = p.m_min; m <= p.m_max; m += p.dm)
				{
					mp.p[modalOrder].a = a;
					mp.p[modalOrder].x = m;
					mp.p[modalOrder].s = s*s;

					real amax_next = a;
					fill_range(v_modalpoints, mp, modalOrder + 1, orderMax, amax_next,da, p, detail_filename);
				}

		}
	}

}

double amp_total(const ModalParameter&  mp , int orderMax  )
{
	double acc = 0.0;
	for (auto i = 0; i <= orderMax; ++i) acc += mp.p[i].a;
	return acc;
}
 
void fill_range_var(ModalParameterBuffer  &v_modalpoints, ModalParameter  mp, int modalOrder, int orderMax, real aMax,   const ModalParameterRegion& prange, const char* detail_filename)
{
	prange.p[modalOrder].expand( [&](float a,float m,float s)
	    {
		   if (modalOrder == orderMax)
		   {
			   mp.p[modalOrder].a = std::max( a, 0.0f );
			   mp.p[modalOrder].x = m;
			   mp.p[modalOrder].s = s; 
			   v_modalpoints.add(mp, orderMax, detail_filename); 
		   }
		   else 
		   {   
			   real aNext = std::max(a, 0.0f);
		       mp.p[modalOrder].a = aNext;
			   mp.p[modalOrder].x = m;
			   mp.p[modalOrder].s = s;
			   fill_range_var(v_modalpoints, mp, modalOrder + 1, orderMax, aNext,    prange,detail_filename);
		   }

	    }
	);


}

 

 

std::vector<ModalParameter>  process_data(   int orderMax,  const Parameters& p , const char* detail_filename)
{


	ModalParameter mp_zero;
	ModalParameterBuffer  v_modalpoints;
	
	fill_range(v_modalpoints, mp_zero, 0, std::min(orderMax,3),  1.0,  p.da, p, detail_filename );
	int64_t psize = 0;
	psize += v_modalpoints.v_parameters_to_process.size();
	printf("has %llu data \n", psize);
	return v_modalpoints.get_results(std::min(orderMax, 3),detail_filename);

}

 ModalParameter   process_data(  const ModalParameterRegion& p, const char* detail_filename)
{

	int orderMax = p.p.size()-1;
	ModalParameter mp_zero;
	ModalParameterBuffer  v_modalpoints;

	fill_range_var(v_modalpoints, mp_zero, 0, std::min(orderMax, 3), 1.0,   p, detail_filename);
	int64_t psize = 0;
	psize += v_modalpoints.v_parameters_to_process.size();
	printf("has %llu data \n", psize);
	return v_modalpoints.get_best_parameter(std::min(orderMax, 3),detail_filename);

}

 ParameterRange expand_border(ParameterRange a, int nborder , float xmin, float xmax)
 {
	 ParameterRange e = ParameterRange(a.x1,a.x2, 2 );
	 e.x1 = a.x1 - nborder * a.dx;
	 e.x2 = a.x2 + nborder * a.dx;
	 e.x1 = std::max(a.x1, xmin);
	 e.x2 = std::min(a.x2, xmax);
	 return e;
}

 ParameterRange expand_around_ii(ParameterRange a , float xcenter, int nborder, float xmin, float xmax)
 {
	 ParameterRange e = ParameterRange(a.x1, a.x2, 2);
	 e.x1 = xcenter - 5 * nborder * a.dx;
	 e.x2 = xcenter + 5 * nborder * a.dx;
	 e.x1 = std::max(a.x1, xmin);
	 e.x2 = std::min(a.x2, xmax);
	 e.dx = 5 * a.dx;
	 return e;
 }


 ModalParameterRange expand_border(ModalParameterRange r, int nborder)
 {
	 return  ModalParameterRange(
		 expand_border(r.a, nborder, 0.0, 1.0),
		 expand_border(r.x, nborder, 0.0, 2.5),
		 expand_border(r.s, nborder, 0.001, 2.0));
 }
 ModalParameterRange expand_around_i(ModalParameterRange r, ModalParameterIc center, int nborder)
 {
	 return  ModalParameterRange(
		 expand_around_ii(r.a, center.a, nborder, 0.0, 1.0),
		 expand_around_ii(r.x, center.x, nborder, 0.0, 2.5),
		 expand_around_ii(r.s, center.s, nborder, 0.001, 2.0));
 }
 
 ModalParameterRegion expand_around(ModalParameterRegion region, ModalParameter center )
 {
	 int orderMax = region.p.size() - 1;

	 ModalParameterRegion e = ModalParameterRegion(region);	 //copia, mas depois sobrescreve
	 for (int j = 0; j <= 3 ;++j)
	 {
		 if (orderMax >= j)
		 {
			 e.p[j + 1] = expand_around_i(region.p[j + 1], center.p[j + 1], 4); // determina uma regiao de 5 dx em torno do centro
		 }
	 }
	 
	 return e;
 }

 ParameterRange expand(float x, float dx, float n, float xmin, float xmax)
 {
	 float x1 = x - n*dx;
	 float x2 = x + n*dx;
	 x1 = std::max(x1, xmin);
	 x2 = std::min(x2, xmax);
	 return  ParameterRange(x1, x2, n);
 }

 ModalParameterRegion expand_around(ModalParameter center, int order, float da, float dx, float ds)
 {
	 int n = 4;
	 if (order < 0)
	 {
		 throw "unable to do";
	 }
	 
	 
		 ModalParameterRange m0 = ModalParameterRange(
			 expand(center.p[0].a, da, n, 0.001f, 1.0f ),
			 expand(center.p[0].x , dx,n, 0.001f, 2.0f),
			 expand(center.p[0].s , ds,n, 0.001f, 1.0f));
		 if (order >= 1)
		 {
			 ModalParameterRange m1 = ModalParameterRange(
				 expand(center.p[1].a, da, n, 0.001f, 1.0f),
				 expand(center.p[1].x, dx, n, 0.001f, 2.0f),
				 expand(center.p[1].s, ds, n, 0.001f, 1.0f));
			 if (order >= 2)
			 {
				 ModalParameterRange m2 = ModalParameterRange(
					 expand(center.p[2].a, da, n, 0.001f, 1.0f),
					 expand(center.p[2].x, dx, n, 0.001f, 2.0f),
					 expand(center.p[2].s, ds, n, 0.001f, 1.0f));
				 if (order >= 3)
				 {
					 ModalParameterRange m3 = ModalParameterRange(
						 expand(center.p[3].a, da, n, 0.001f, 1.0f),
						 expand(center.p[3].x, dx, n, 0.001f, 2.0f),
						 expand(center.p[3].s, ds, n, 0.001f, 1.0f));
					 return ModalParameterRegion(m0, m1, m2, m3);
				 }
				 return ModalParameterRegion(m0, m1, m2);
			 }
			 return ModalParameterRegion(m0, m1);
		 }
		 return ModalParameterRegion(m0);
 
 }

int main()
{
	FILE* f = fopen("Lk.dat", "w");
	fclose(f);
	FILE* fs = fopen("progress.dat", "w");
	fclose(fs);
	{

//		auto m1 = ModalParameterRange({ 0.85,0.95, 6 }, { 0.13,0.18,6 }, { 0.005,0.03,8 });
//		auto m2 = ModalParameterRange({ 0.06,0.1,7 }, { 0.45,0.65,5 }, { 0.08,0.20,4 });
//      auto m3 = ModalParameterRange({ 0.02,0.08,5 }, { 0.1,1.7,12 }, { 0.01,0.8,8});

		auto m1 = ModalParameterRange({ 0.80,0.90, 5 }, { 0.13,0.18,4 }, { 0.005,0.03,4 });
		auto m2 = ModalParameterRange({ 0.0,0.1,4 },  { 0.45,0.65,4 },   { 0.01,0.8,4 });
		auto m3 = ModalParameterRange({ 0.0,0.1,5 }, { 0.0, 1.7, 5 },    { 0.01,0.5, 5});
		auto m4 = ModalParameterRange({ 0.0,0.05,6 }, { 0.1, 1.9, 6 },   { 0.05,0.2, 6 });
		

		ModalParameterRegion region( m1,m2,m3 ,m4 );
		ModalParameter best_parameter;

		for (int j = 5; j >0; --j)
		{
			best_parameter = process_data(region, nullptr);
			printf("============================\n");
			dump(best_parameter);
			if (j !=0 ) region = refine_around(region, best_parameter);
		}

		
		// regiao de maxima likehood ja identificada, faca um relatorio em torno dessa regiao
		
		ModalParameterRegion region_analise = expand_around(  best_parameter , region.p.size(), 0.05 , 0.02,0.02 );
		dump_region(region_analise);
		process_data(region_analise, "Lk.dat");

		return 0;

	}


	Parameters params;
	params.da = 0.02;
	

	params.m_min = 0.0;
	params.m_max = 1.0;
	params.dm = (params.m_max - params.m_min)/ 20.0;


	params.s_min = sqrt( 0.01);
	params.s_max = sqrt( 1.0 );
	params.ds = (params.s_max - params.s_min) / 12.0;
	

	auto v_data = process_data(2, params, "Lk.dat");
	 
	 
	printf("============================\n");
	 

    return 0;
}

