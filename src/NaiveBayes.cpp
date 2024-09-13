/*
 *  NaiveBayes.cpp
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

#include "NaiveBayes.h"
#include "MatrixOperations.h"

SiMath::NaiveBayes::NaiveBayes(const SiMath::NaiveBayesParameters & p) :
	_nbrClasses(p.nbrClasses),
	_nbrVar(p.dimension),
	_classes(p.nbrClasses)
	//_classTrained(2),
	//_classMean(0,0,0.0),
	//_classSDev(0,0,0.0),
	//_classPrior(2,0.5)
{
	if ( _nbrClasses < 2 )
		throw(NaiveBayesError("There should be at least 2 classes in a naive Bayesian classifier."));
		
	for ( unsigned int i=0; i<_nbrClasses; ++i)
	{
		_classes[i].mean.resize(_nbrVar);
		_classes[i].mean = 0.0;
		_classes[i].stdev.resize(_nbrVar);
		_classes[i].stdev = 1.0;
		_classes[i].prior = 1.0/_nbrClasses;
		_classes[i].trained = false;
	}
}

/****
SiMath::NaiveBayes::NaiveBayes(unsigned int n, unsigned int v) :
	_nbrClasses(n),
	_nbrVar(v),
	_classTrained(n),
	_classMean(n,v,0.0),
	_classSDev(n,v,0.0),
	_classPrior(n,1/n)
{
}
*******/

SiMath::NaiveBayes::~NaiveBayes()
{
}

void
SiMath::NaiveBayes::train(unsigned int c, const SiMath::Matrix & data) 
{
	if ( c >= _nbrClasses )
		throw(NaiveBayesError("Class label is outside the scope. Unable to train class mean and sdev"));

	if ( data.nbrColumns() != _nbrVar )
		throw(NaiveBayesError("Dimension of data do not match model dimensions."));
	
	_classes[c].mean = SiMath::meanOfAllColumns(data);
	_classes[c].stdev = SiMath::stDevOfAllColumns(data);
	_classes[c].trained = true;

	return;
}


void 
SiMath::NaiveBayes::setModel(unsigned int c, const Vector & m, const Vector & sd)
{
	if ( c < 0 || c >= _nbrClasses )
		throw(NaiveBayesError("Trying to set parameters of class outside the range of defined classes."));

	if ( m.size() != sd.size() || m.size() != _nbrVar )
		throw(NaiveBayesError("Incompatible dimensions of mean and standard deviation of class and class definition."));
	
	for ( unsigned int i=0; i<m.size(); ++i )
	{
		_classes[c].mean[i] = m[i];
		_classes[c].stdev[i] = sd[i];
	}
	_classes[c].trained = true;
	return;
}

void
SiMath::NaiveBayes::setClassPrior(unsigned int c, double d)
{
	if ( c < 0 || c >= _nbrClasses )
		throw(NaiveBayesError("Trying to set parameters of class outside the range of defined classes."));

	if ( d < 0 || d > 1 )
		throw(NaiveBayesError("Wrong value to define class prior."));
	
	_classes[c].prior = d;
	return;
}


SiMath::Matrix
SiMath::NaiveBayes::classify(const SiMath::Matrix & data)
{
	
	bool check = true;
	for ( unsigned int i=0; i<_nbrClasses; ++i )
	{
		check &= _classes[i].trained;
	}
	if ( !check )
		throw(NaiveBayesError("Not all classes have been trained. Unable to classify data."));
					
	if ( data.nbrColumns() != _nbrVar )
		throw(NaiveBayesError("Dimension of input data incompatible with defined class model."));

	SiMath::Matrix v(data.nbrRows(),_nbrClasses,0.0);
	double s1 =  log(sqrt(2 * PI));
	// double s2 = 2 * PI * PI;
	for ( unsigned int i=0; i<data.nbrRows(); ++i)
	{
		for ( unsigned int c=0; c<_nbrClasses; ++c )
		{
			v[i][c] = log(_classes[c].prior);
			for ( unsigned int j=0; j<_nbrVar; ++j )
			{
				double mm = _classes[c].mean[j];
				double ss = _classes[c].stdev[j];
				double d = data[i][j] - mm;
				//v[i][c] += -(d*d) / s2 - log(ss) - s1;
				//std::cerr << -(d*d) / (2*ss*ss) - log(ss) - s1 << " ";
				v[i][c] += -(d*d) / (2*ss*ss) - log(ss) - s1;
			}
		}
	}
	// std::cerr << std::endl;
	return v;
}


SiMath::Vector
SiMath::NaiveBayes::classify(const Vector & data)
{
	if ( data.size() != _nbrVar )
		throw(NaiveBayesError("Dimension of input data incompatible with defined class model."));
	
	Vector v(_nbrClasses,0.0);
	double s1 =  log(sqrt(2 * PI));
	// double s2 = 2 * PI * PI;

	for ( unsigned int c=0; c<_nbrClasses; ++c )
	{
		v[c] = log(_classes[c].prior);
		for ( unsigned int j=0; j<_nbrVar; ++j )
		{
			double mm = _classes[c].mean[j];
			double ss = _classes[c].stdev[j];
			double d = data[j] - mm;
			// std::cerr << -(d*d) / (2*ss*ss) - log(ss) - s1 << " ";
			v[c] += -(d*d) / (2*ss*ss) - log(ss) - s1;
		}
	}
	// std::cerr << v << std::endl;
	return v;
}


void 
SiMath::NaiveBayes::setNbrOfClasses(unsigned int n)
{
	if ( _nbrClasses != n )
		_classes.resize(n);
	
	_nbrClasses = n;
	for ( unsigned int c=0; c<_nbrClasses; ++c )
	{	
		// classes are no longer trained
		_classes[c].prior = 1.0/_nbrClasses;
		_classes[c].trained = false;
		
		// change in number of variables changes local matrices
		_classes[c].mean = 0.0;
		_classes[c].stdev = 1.0;			
	}
}


/******
void 
SiMath::NaiveBayes::setNbrOfVariables(unsigned int v)
{
	if ( _nbrVar != v )
		_nbrVar = v;
		
	for ( unsigned int c=0; c<_nbrClasses; ++c )
	{	
		// classes are no longer trained
		_classes[c].prior = 1.0/_nbrClasses;
		_classes[c].trained = false;
			
		// change in number of variables changes local matrices
		_classes[c].mean.reset(_nbrVar);
		_classes[c].mean = 0.0;
		_classes[c].stdev.reset(_nbrVar);	
		_classes[c].stdev = 1.0;			
	}
}
********/

bool 
SiMath::NaiveBayes::isClassTrained(unsigned int c)
{
	if ( c < 0 || c >= _nbrClasses )
		return false;
	return _classes[c].trained;
}

double 
SiMath::NaiveBayes::getClassPrior(unsigned int c)
{
	if ( c < 0 || c >= _nbrClasses )
		throw(NaiveBayesError("Trying to get class prior for undefined class number"));
	return _classes[c].prior;
}


SiMath::Vector 
SiMath::NaiveBayes::getClassMean(unsigned int c)
{
	if ( c < 0 || c >= _nbrClasses )
		throw(NaiveBayesError("Trying to get class mean for undefined class number"));
	return _classes[c].mean;
}

SiMath::Vector 
SiMath::NaiveBayes::getClassStDev(unsigned int c)
{
	if ( c < 0 || c >= _nbrClasses )
		throw(NaiveBayesError("Trying to get class st.dev. for undefined class number"));
	return _classes[c].stdev;
}
