/*
 *  Utilities.cpp
 *
 *  Created by Gert Thijs on 27/03/06.
 *  Copyright 2006 Silicos NV. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 *     + Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     + Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY SILICOS NV AND CONTRIBUTORS ``AS IS'' AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL SILICOS NV AND CONTRIBUTORS BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. 
 *
 */

#include "Utilities.h"
#include <float.h>


double 
SiMath::randD(double a, double b)
{
	double d(a);
	d += (b-a) * ((double)rand()/RAND_MAX);
	return d;
}


int 
SiMath::randI(int a, int b)
{
	int d(a);
	d += (int)((b-a) * ((double)rand()/RAND_MAX));
	return d;
}


double 
SiMath::lnGamma(double x)
{
	// Variables
	double coef[7] = {
		0.0, 
		76.18009173, 
		-86.50532033, 
		24.01409822, 
		-1.231739516,
		0.120858003e-2, 
		-0.536382e-5};
	double stp(2.50662827465);
	double xx(x - 1.0);
	double ser(1.0);
	double tmp = xx + 5.5;
	tmp = (xx + 0.5) * log(tmp) - tmp;
	for (unsigned int j = 1; j <= 6; j++)
	{
		xx += 1.0;
		ser += coef[j] / xx;
	}
	return (tmp + log(stp * ser));
}


double
SiMath::incompGamma(double a, double x)
{
	if ( x < 0.0 || a <= 0.0  ) 
		SiMathError("a or x out of scope for incomplete gamma function");

	if ( x == 0.0 )
		return 0.0;

	const double dmin = DBL_MIN/DBL_EPSILON;
	
	if ( x < a+1.0 )
	{
		double sum(1.0/a), del(1.0/a), anext(a);
		for ( unsigned int n=0; n < 100; ++n )
		{
			++anext;
			del *= x/anext;
			sum += del;
			if ( fabs(del) < fabs(sum) * DBL_EPSILON )
			{
				return sum * exp( a * log(x) - lnGamma(a) - x);
			}
		}
	}
	else
	{
		double del(0.0), anext(a), c(1.0/dmin), b(x + 1.0-a);
		double d(1.0/b);
		double h(d);
		
		for ( unsigned int n=1; n < 100; ++n )
		{
			anext = -n * (n-a);
			b += 2.0;
			d = anext * d + b;
			if ( fabs(d) < dmin) d = dmin;
			c = b + anext/c;
			if ( fabs(c) < dmin ) c = dmin;
			del = c / d;
			h *= del;
			if ( fabs(del - 1.0) <= DBL_EPSILON )
			{	
				return 1.0 - h * exp(a * log(x) - lnGamma(a) - x);
			}
		}
	}
	throw(SiMathError("Unable to complete computation of incomplete gamma function"));
}


double
SiMath::beta(double x, double y)
{
	return exp(lnGamma(x) + lnGamma(y) - lnGamma(x+y));
}


double 
SiMath::incompBeta(double a, double b, double x)
{
	double bt(0.0);
	if ( x < 0.0 || x > 1.0 ) 
		SiMathError("x out of scope for incomplete beta function");
	if ( x != 0.0 && x != 1.0 )
		bt = exp(lnGamma(a+b) - lnGamma(a) - lnGamma(b) +
						 a * log(x) + b * log(1.0 -x));
	
	if ( x < (a+1.0)/(a+b+2.0) )
		return bt * cfBeta(a,b,x)/a;
	else
		return 1.0 - bt * cfBeta(b,a,1.0-x)/b;
}

double 
SiMath::cfBeta(double a, double b, double x)
{
	// Variables
	unsigned int m2;
	double aa;
	double c(1.0);
	double d;
	double del;
	double h;
	double qab(a + b);
	double qam(a - 1.0);
	double qap(a + 1.0);
	
	const double dmin = DBL_MIN/DBL_EPSILON;
	d = 1.0 - qab * x / qap;
	if (fabs(d) < dmin) d = dmin;
	d = 1.0 / d;
	h = d;
	
	for (unsigned int i = 1; i <= 100; ++i )
	{
		m2 = 2 * i;
		aa = i * (b - i) * x / ((qam + m2) * (a + m2));
		d = 1.0 + aa * d;
		if (fabs(d) < dmin ) d = dmin;
		c = 1.0 + aa / c;
		if (fabs(c) < dmin ) c = dmin;
		d = 1.0 / d;
		h *= d * c;
		
		aa = -(a + i) * (qab + i) * x / ((qap + m2) * (a + m2));
		d = 1.0 + aa * d;
		if ( fabs(d) < dmin ) d = dmin;
		c = 1.0 + aa / c;
		if ( fabs(c) < dmin ) c = dmin;
		d = 1.0 / d;
		del = d * c;
		h *= del;
		if ( fabs(del - 1.0) < DBL_EPSILON )
			return h;
	}
	throw(SiMathError("Unable to complete computation of incomplete beta function"));
}


double 
SiMath::FDist(double x, double df1, double df2)
{
	double z(df2 / (df2 + df1 * x));
	double a(df2 / 2.0);
	double b(df1 / 2.0);
	
	if (z < 0 || z > 1)
		SiMathError("Bad argument x in FDist.");

	return incompBeta(a,b,z);
}



void 
SiMath::fitSigmoid( const Vector & values, const Vector & targets, 
										 double & A, double & B)
{
	if ( values.size() != targets.size() )
		throw(SiMathError("Unable to fit sigmoid if values and targets differ in length"));
	
	unsigned int l = values.size();

	// compute prior
	double prior1=0, prior0 = 0;
	for (unsigned int i=0; i<l; ++i)
	{	
		if (targets[i] > 0) 
			prior1 += 1;
		else 
			prior0 += 1;
	}
	
	// Initial Point and Initial Fun Value
	A = 0.0;
	B = log((prior0+1.0)/(prior1+1.0));

	double hiTarget=(prior1+1.0)/(prior1+2.0);
	double loTarget=1/(prior0+2.0);            
	
	int maxIterations=100; 	 // Maximal number of iterations
	double min_step=1e-10;	 // Minimal step taken in line search
	double sigma=1e-3;	     // For numerically strict PD of Hessian
	double eps=1e-5;
	
	// local variables
	double fApB,p,q,h11,h22,h21,g1,g2,det,dA,dB,gd,stepsize;
	double newA,newB,newf,d1,d2;
	
	double fval = 0.0;
	double * t = new double[l];
	for (unsigned int i=0; i<l; ++i)
	{
		if (targets[i]>0) 
			t[i]=hiTarget;
		else 
			t[i]=loTarget;
		fApB = values[i]*A+B;
		if (fApB>=0)
			fval += t[i]*fApB + log(1+exp(-fApB));
		else
			fval += (t[i] - 1)*fApB +log(1+exp(fApB));
	}
	
	// do the iterations
	unsigned int iter=0;
	while ( iter<maxIterations )
	{
		// Update Gradient and Hessian (use H' = H + sigma I)
		h11=sigma; // numerically ensures strict PD
		h22=sigma;
		h21=0.0;g1=0.0;g2=0.0;
		for (unsigned int i=0; i<l; ++i)
		{
			fApB = values[i]*A+B;
			if (fApB >= 0)
			{
				p=exp(-fApB)/(1.0+exp(-fApB));
				q=1.0/(1.0+exp(-fApB));
			}
			else
			{
				p=1.0/(1.0+exp(fApB));
				q=exp(fApB)/(1.0+exp(fApB));
			}
			d2 = p*q;
			h11 += values[i] * values[i]*d2;
			h22 += d2;
			h21 += values[i]*d2;
			d1 = t[i]-p;
			g1 += values[i]*d1;
			g2 += d1;
		}
		
		// Stopping Criteria
		if (fabs(g1)<eps && fabs(g2)<eps)
			break;
		
		// Finding Newton direction: -inv(H') * g
		det = h11*h22-h21*h21;
		dA = (h21*g2 - h22*g1);
		dB = (h21*g1 - h11*g2);   
		gd = (g1*dA + g2*dB)/det;  // removed /det from dA and dB, now only one / needed 
		
		stepsize = 1; 		// Line Search
		while (stepsize >= min_step)
		{
			newA = A + stepsize * dA;
			newB = B + stepsize * dB;
			
			// New function value
			newf = 0.0;
			for (unsigned int i=0;i<l;++i)
			{
				fApB = values[i]*newA+newB;
				if (fApB >= 0)
					newf += t[i]*fApB + log(1+exp(-fApB));
				else
					newf += (t[i] - 1)*fApB +log(1+exp(fApB));
			}
			// Check sufficient decrease
			if (newf<fval+0.0001*stepsize*gd)
			{
				A=newA;B=newB;fval=newf;
				break;
			}
			else
			{	
				stepsize = stepsize / 2.0;
			}
		}
		
		if (stepsize < min_step)
		{
			std::cerr << "WARNING: Line search fails in two-class probability estimates." << std::endl;
			break;
		}
		
		++iter;
	}
	
	if (iter>=maxIterations)
		std::cerr << "Reached maximal iterations in two-class probability estimates" << std::endl;
	
	delete[] t;
}


double 
SiMath::sigmoid(double x, double A, double B)
{
	double fApB = x*A+B;
	if (fApB >= 0)
		return exp(-fApB)/(1.0+exp(-fApB));
	else
		return 1.0/(1+exp(fApB)) ;
}



double 
SiMath::powi(double x, unsigned int n)
{
	double tmp(x), v(1.0);
	
	for(unsigned int i=n; i>0; i/=2)
	{
		if( i%2 == 1) 
			v *= tmp;
		tmp = tmp * tmp;
	}
	return v;
}


double 
SiMath::triangle(const double & a, const double & b )
{
	double A(fabs(a)), B(fabs(b));
	if ( A > B )
	{
		return A * sqrt(1.0 + (B/A)*(B/A));
	}
	else if ( B == 0 )
	{
		return 0;
	}
		
	return B * sqrt(1.0 + (A/B)*(A/B));			
}


unsigned int 
SiMath::nchoosek(unsigned int k, unsigned int n){
	unsigned int c;
	
	if (k > n)
		return 0;
	
	if (k > n/2)
		k = n-k; // faster
	
	unsigned int accum = 1;
	for (unsigned int i = 1; i <= k; i++){
		accum = accum * (n-k+i) / i;		
	}
	
	return accum + 0.5; // avoid rounding error
}

