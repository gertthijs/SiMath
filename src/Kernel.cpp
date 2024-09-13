/*
 *  Kernel.cpp
 *
 *  Created by Gert Thijs on 04/04/06.
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

#include "Kernel.h"
#include "Utilities.h"

//------------------------------------------------------------------------------------------------------------------------------------
// Constructor
SiMath::Kernel::Kernel( SiMath::SVMProblem & prob, const SiMath::SVMParameters & param )
	: 
	_kernelType(param.kernelType), 
	_degree(param.degree),
	_gamma(param.gamma), 
	_coef0(param.coef0), 
	_dataPointer(&prob),
	_xSquare(NULL),
	_kernelFunction(NULL)
{
	switch(_kernelType)
	{
		case LINEAR:
			_kernelFunction = &Kernel::_linearKernel;
			break;
		case POLY:
			_kernelFunction = &Kernel::_polyKernel;
			break;
		case RBF:
			_kernelFunction = &Kernel::_rbfKernel;
			break;
		case SIGMOID:
			_kernelFunction = &Kernel::_sigmoidKernel;
			break;
	}
	
	if( _kernelType == RBF)
	{
		_xSquare = new double[prob.nbrData];
		for(int i=0;i<prob.nbrData;++i){
			_xSquare[i] = _dot(prob.nbrVar, prob.x[i], prob.x[i]);
			// _xSquare[prob.nbrData+i] = _xSquare[i];
		}
	}	
}

 
SiMath::Kernel::~Kernel()
{
	// reset data pointer, do not delete data here 
	_dataPointer = NULL; 
	_kernelFunction = NULL;
	if ( _xSquare != NULL )
		delete[] _xSquare;
}

// dot-product
double 
SiMath::Kernel::_dot(unsigned int l, const double * px, const double * py)
{
	double sum(0);
	for ( unsigned int i=0; i<l; ++i )
		sum += px[i] * py[i];
	return sum;
}


// single kernel evaluation
double 
SiMath::Kernel::function(unsigned int l, const double * px, const double * py, const SiMath::SVMParameters & param)
{
	double d(0.0);
	switch(param.kernelType)
	{
		case LINEAR:
			return _dot(l,px,py);
		case POLY:
			return SiMath::powi(param.gamma*_dot(l,px,py) + param.coef0, param.degree);
		case RBF:
		{
			for ( unsigned int i=0; i<l; ++i ){
				d += (px[i] - py[i]) * (px[i] - py[i]);				
			}
			return exp(-param.gamma * d);
		}
		case SIGMOID:
			return tanh(param.gamma * _dot(l,px,py) + param.coef0);
		case SPLINE:
			d = _dot(l,px,py);
			if (	d <= 1 ) 
				return (2/3)-d*d + 0.5*d*d*d;
			else if ( d <= 2 ) 
				return (2-d)*(2-d)*(2-d)/6;
			else 
				return 0;			
		default:
			return 0;	/* Unreachable */
	}
}


double 
SiMath::Kernel::_linearKernel(unsigned int i, unsigned  int j) const 
{ 
	return _dot(_dataPointer->nbrVar,_dataPointer->x[i],_dataPointer->x[j]);
}

double 
SiMath::Kernel::_polyKernel(unsigned int i, unsigned int j) const 
{ 
	return powi(_gamma*_dot(_dataPointer->nbrVar,_dataPointer->x[i],_dataPointer->x[j])+_coef0,_degree); 
};

double 
SiMath::Kernel::_rbfKernel(unsigned int i, unsigned int j) const 
{ 
	return exp(-_gamma*(_xSquare[i]+_xSquare[j]-2*_dot(_dataPointer->nbrVar,_dataPointer->x[i],_dataPointer->x[j]))); 
};

double 
SiMath::Kernel::_sigmoidKernel(unsigned int i, unsigned int j) const 
{ 
	return tanh(_gamma*_dot(_dataPointer->nbrVar,_dataPointer->x[i],_dataPointer->x[j])+_coef0); 
};

double 
SiMath::Kernel::_splineKernel(unsigned int i, unsigned int j) const
{
  double d = _dot(_dataPointer->nbrVar,_dataPointer->x[i],_dataPointer->x[j]);
  if (	d <= 1 ) 
		return (2/3)-d*d + 0.5*d*d*d;
  else if ( d <= 2 ) 
		return (2-d)*(2-d)*(2-d)/6;
  else 
		return 0;
}



//------------------------------------------------------------------------------------------------------------------------------------
// Implementation of kernel matrices for various problem formulations
//------------------------------------------------------------------------------------------------------------------------------------

// ------ Support Vector Classification -------
SiMath::svcQMatrix::svcQMatrix (SVMProblem & prob, const SVMParameters & param) : 
	Kernel(prob, param), 
	_index(prob.nbrData),
	_QD(NULL),
	_cache(NULL)
{
	// set up a new cache 
	_cache = new Cache(prob.nbrData,(int)(param.cacheSize*(1<<20)));
	
	// setup indices
	_QD = new double[prob.nbrData];
	for( int i=0; i<prob.nbrData; ++i)
	{	
		_index[i] = i;
		_QD[i]= (double)(this->*_kernelFunction)(i,i);
	}
}

double *
SiMath::svcQMatrix::getQ(unsigned int i, unsigned int len) const
{
	double * data = NULL; // pointer to the kernel cache
	int start(0);
	unsigned int realI = _index[i];
	if( (start = _cache->getData(i,&data,len)) < len )
	{
		double d = _dataPointer->y[realI];
		for(int j=start;j<len;++j)
		{	
			data[j] = (double)(d  * _dataPointer->y[_index[j]] * (this->*_kernelFunction)(realI,_index[j]));
		}
	}
	return data;
}

void 
SiMath::svcQMatrix::swapIndex(const unsigned int i, const unsigned int j)
{
	std::swap(_index[i],_index[j]);
	if( _xSquare != 0 )
		std::swap(_xSquare[i],_xSquare[j]);

	_cache->swapIndex(i,j);
	std::swap(_QD[i],_QD[j]);
}

void 
SiMath::svcQMatrix::printIndex(){
	std::cerr << "_index:";
	for ( unsigned int i=0; i<_index.size(); ++i){
		std::cerr << "(" << i << "," << _index[i] << ")";
	}
	std::cerr << std::endl;
}

SiMath::svcQMatrix::~svcQMatrix()
{
	delete _cache;
	delete[] _QD;
}


//---------------------------------------------------------------------------------------------------------------------------------------
//------ One Class Problem -------
SiMath::oneClassQMatrix::oneClassQMatrix(SiMath::SVMProblem & prob, const SiMath::SVMParameters & param) : 
	Kernel(prob, param),
	_index(prob.nbrData),
	_QD(NULL),
	_cache(NULL)
{
	_cache = new Cache(prob.nbrData,(int)(param.cacheSize*(1<<20)));
	_QD = new double[prob.nbrData];
	for( int i=0;i<prob.nbrData; ++i )
	{	
		_index[i] = i;
		_QD[i]= (double)(this->*_kernelFunction)(i,i);
	}
}

double *
SiMath::oneClassQMatrix::getQ(unsigned int i, unsigned int len) const
{
	double * data(0); // pointer to the kernel cache
	int start;
	unsigned int real_i = _index[i];
	if( (start = _cache->getData(i, &data, len)) < len)
	{
		for(int j=start;j<len;++j)
		{	
			data[j] = (double)(this->*_kernelFunction)(real_i,_index[j]);
		}
	}
	return data;
}

void 
SiMath::oneClassQMatrix::swapIndex(unsigned int i, unsigned int j)
{
	std::swap(_index[i],_index[j]);
	if( _xSquare != 0 )
		std::swap(_xSquare[i],_xSquare[j]);
	_cache->swapIndex(i,j);
	std::swap(_QD[i],_QD[j]);
}

void 
SiMath::oneClassQMatrix::printIndex(){
	std::cerr << "_index:";
	for ( unsigned int i=0; i<_index.size(); ++i){
		std::cerr << "(" << i << "," << _index[i] << ")";
	}
	std::cerr << std::endl;
}

SiMath::oneClassQMatrix::~oneClassQMatrix()
{
	delete _cache;
	delete[] _QD;
}



//---------------------------------------------------------------------------------------------------------------------------------------
//------ Support Vector Regression -------
SiMath::svrQMatrix::svrQMatrix(SiMath::SVMProblem & prob, const SiMath::SVMParameters & param)	: 
	Kernel(prob, param),
	_l(prob.nbrData),
	_index(2*prob.nbrData),
	_QD(NULL),
	_sign(NULL),
	_nextBuffer(0),
	_cache(NULL)
{
	_cache = new Cache(_l,(int)(param.cacheSize*(1<<20)));
	_QD = new double[2*_l];
	_sign = new signed char[2*_l];
	for(int k=0;k<_l;k++)
	{
		_sign[k] = 1;
		_sign[k+_l] = -1;
		_index[k] = k;
		_index[k+_l] = k;
		_QD[k]= (double)(this->*_kernelFunction)(k,k);
		_QD[k+_l] = _QD[k];
	}
	_buffer[0] = new double[2*_l];
	_buffer[1] = new double[2*_l];
}

void
SiMath::svrQMatrix::swapIndex(unsigned int i, unsigned int j)
{
	std::swap(_index[i],_index[j]);
	std::swap(_sign[i],_sign[j]);
//	if( _xSquare != 0 )
//		std::swap(_xSquare[i],_xSquare[j]);
	std::swap(_QD[i],_QD[j]);
}


double *
SiMath::svrQMatrix::getQ(unsigned int i, unsigned int len) const
{
	double * data = NULL; // pointer to acccess the kernel cache
	int realI = _index[i];
	if( _cache->getData(realI, &data, _l) < _l )
	{
		for(int j=0;j<_l;++j){
			data[j] = (double)(this->*_kernelFunction)(realI,j);
		}
	}
	
	// reorder and copy
	double * buf = _buffer[_nextBuffer];
	_nextBuffer = 1 - _nextBuffer;
	signed char si = _sign[i];
	for(int j=0;j<len;++j){
		buf[j] = si * _sign[j] * data[_index[j]];
	}
	return buf;
}

void 
SiMath::svrQMatrix::printIndex(){
	std::cerr << "_index:";
	for ( unsigned int i=0; i<_index.size(); ++i){
		std::cerr << "(" << i << "," << _index[i] << ")";
	}
	std::cerr << std::endl;
}


SiMath::svrQMatrix::~svrQMatrix()
{
	delete _cache;
	delete[] _sign;
	delete[] _buffer[0];
	delete[] _buffer[1];
	delete[] _QD;
}


