/*
 *  SVM.cpp
 *
 *  Created by Gert Thijs.
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

#include "SVM.h"
#include "SVMSolver.h"
#include "Utilities.h"

// constructor
SiMath::SVM::SVM(SiMath::SVMParameters & param) :
	_param(param),
	_nbrSV(0),
	_nbrVar(0),
	_supVecs(0,0),
	_svCoef(0,0),
	_probA(0),
	_probB(0),
	_decFunctions(0),
	_nbrClasses(2),
	_nbrClassSV(0),
	_classLabels(),
	_dataLabels(0)
{ 
	_classLabels.clear();
	if ( _param.svmType == ONE_CLASS || _param.svmType == EPSILON_SVR || _param.svmType == NU_SVR)
	{
		// regression or one-class-svm (set classes to two)
		_nbrClasses = 2;
	}
}


SiMath::SVM::SVM(const SiMath::SVM & src ) : 
	_param(src._param),
	_nbrSV(src._nbrSV),
	_nbrVar(src._nbrVar)
{ 
	throw(SiMath::SVMError("Copy constructor not yet implemented"));
}

// SVM Model class
SiMath::SVM::~SVM()
{
}


bool
SiMath::SVM::checkProbabilityModel()
{
	return ( (_param.svmType == C_SVC || _param.svmType == NU_SVC) && _probA.size() !=0 && _probB.size() != 0 ) || 
	( (_param.svmType == EPSILON_SVR || _param.svmType == NU_SVR) && _probA.size() != 0 );
}


SiMath::Vector 
SiMath::SVM::predict(const SiMath::Matrix & data)
{
	if ( data.nbrColumns() != _nbrVar )
	{
		throw(SiMath::SVMError("Number of variables in data matrix does not match number of variables in model."));
	}
	unsigned int l(data.nbrRows());
	SiMath::Vector v(l);
	
	for ( unsigned int i=0; i<l; ++i )
		v[i] = _predict(data[i]);

	return v;
}


SiMath::Vector 
SiMath::SVM::predictProbability(const SiMath::Matrix & data)
{
	if ( data.nbrColumns() != _nbrVar )
	{
		throw(SiMath::SVMError("Number of variables in data matrix does not match number of variables in model."));
	}
	unsigned int l(data.nbrRows());
	SiMath::Vector v(l);
	
	for ( unsigned int i=0; i<l; ++i )
		v[i] = _predictProbability(data[i]);
	
	return v;
}


void 
SiMath::SVM::train(SiMath::SVMProblem & prob)
{
	// initialize data dependent data members
	_nbrVar = prob.nbrVar;
	
	// loop iterator variables
	unsigned int i(0), j(0), k(0);
	
	// check the type of svm 
	if( _param.svmType == ONE_CLASS || _param.svmType == EPSILON_SVR || _param.svmType == NU_SVR)
	{
		// regression or one-class-svm (set classes to two)
		_nbrClasses = 2;		
		_classLabels.clear();
		_nbrClassSV.clear();
		_dataLabels.clear();
		_probA.reset(0); 
		_probB.reset(0);

		if( _param.probability && ( _param.svmType == EPSILON_SVR || _param.svmType == NU_SVR) )
		{
			_probA.reset(1);
			_probA[0] = _svrProbability(prob,_param);
		}

		// get the decision function
		_decFunctions.resize(1);
		_trainOne(prob,_param,0,0,_decFunctions[0]);
				

		// count the number of support vectors
		_nbrSV = 0;
		for( i=0; i<prob.nbrData; ++i)
		{	
			if ( fabs(_decFunctions[0].alpha[i]) > 0 )
				++_nbrSV;
		}

		// only one set of support vector coefficients
		_supVecs.reset(_nbrSV,_nbrVar);
		_svCoef.reset(1,_nbrSV);
		
		unsigned int idx(0);
		for(i=0;i<prob.nbrData;++i)
		{	
			if(fabs(_decFunctions[0].alpha[i]) > 0)
			{
				for ( j=0; j<_nbrVar; ++j ) {
					_supVecs[idx][j] = prob.x[i][j]; // copy data row
				}
				_svCoef[0][idx] = _decFunctions[0].alpha[i]; // updata coefficient
				++idx;
			}		
		}
	}
	else
	{
		// classification
		_nbrSV = prob.nbrData;
		_dataLabels.resize(_nbrSV);

		// functions to hold starting positions of new class labels 
		std::vector<unsigned int> classStart;
		std::vector<unsigned int> classCount;
	
		// group training data of the same class
		_groupClasses(prob,classStart,classCount);
		
		// calculate weighted C
		std::vector<double> weightedC(_nbrClasses);
		for( i=0; i<_nbrClasses; ++i)
			weightedC[i] = _param.C;
		
		for ( std::map<int,double>::iterator mi=_param.weights.begin(); mi != _param.weights.end(); ++mi )
		{	
			if ( _classLabels.count(mi->first) == 0 ){
				std::cerr << "WARNING: class label " << mi->first << " specified in weight is not found in data set" << std::endl;
			}else{
				weightedC[_classLabels[mi->first]] *= mi->second;
			}
		}

		// train k*(k-1)/2 models
		std::vector<bool> nonZero(_nbrSV);
		for( i=0; i<_nbrSV; ++i)
			nonZero[i] = false;

		// update and reset decision functions
		unsigned int nComb = _nbrClasses * (_nbrClasses - 1)/2;
		_decFunctions.resize(nComb);
		if (_param.probability)
		{
			_probA.reset(nComb);
			_probB.reset(nComb);
		}
		else
		{
			_probA.reset(0);
			_probB.reset(0);
		}

		int p = 0;
		// loop over all combination of classes
		for( i=0; i<_nbrClasses; ++i)
		{
			for( j=i+1; j<_nbrClasses; ++j)
			{
				unsigned int si = classStart[i], sj = classStart[j];
				unsigned int ci = classCount[i], cj = classCount[j];
				// std::cerr << si << "," << ci << " <--> " << sj << "," << cj << std::endl;

				// create a new sub problem by comparing class i and j
				SiMath::SVMProblem subProb;
				subProb.nbrData = ci+cj;
				subProb.nbrVar = prob.nbrVar;
				subProb.x.reset(subProb.nbrData,subProb.nbrVar);
				subProb.y.reset(subProb.nbrData);
				for( k=0; k<ci; k++)
				{
					for ( unsigned int jj=0; jj<subProb.nbrVar; ++jj) {
						subProb.x[k][jj] = prob.x[_dataLabels[si+k]][jj];
					}
					subProb.y[k] = +1;
				}
				for( k=0; k<cj; k++)
				{
					for ( unsigned int jj=0; jj<subProb.nbrVar; ++jj) {
						subProb.x[ci+k][jj] = prob.x[_dataLabels[sj+k]][jj];
					}
					subProb.y[ci+k] = -1;
				}

				if(_param.probability)
					_binarySVCProbability(subProb,_param,weightedC[i],weightedC[j],_probA[p],_probB[p]);

				_trainOne(subProb,_param,weightedC[i],weightedC[j],_decFunctions[p]);
				
				for( k=0; k<ci; k++)
				{	
					if( !nonZero[si+k] && fabs(_decFunctions[p].alpha[k]) > 0)
						nonZero[si+k] = true;
				}
				for( k=0; k<cj; k++)
				{	
					if( !nonZero[sj+k] && fabs(_decFunctions[p].alpha[ci+k]) > 0)
						nonZero[sj+k] = true;
				}
				// keep track on the number of combinations already processed.
				// std::cerr << p << ":" << _decFunctions[p].alpha << std::endl;
				++p;
			}
		}
		// build output
		std::vector<int> nzCount(_nbrClasses);
		_nbrClassSV.resize(_nbrClasses);

		// reset number of support vectors
		_nbrSV = 0;
		for( i=0; i< _nbrClasses; ++i)
		{
			int nSV = 0;
			for( j=0; j<classCount[i]; ++j)
			{	
				if(nonZero[classStart[i]+j])
				{	
					++nSV;
					++_nbrSV;
				}
			}
			_nbrClassSV[i] = nSV;
			nzCount[i] = nSV;
			// std::cerr << i << ":" << nSV << " ";
		}
		
		std::cerr << "Total number of support vectors = " << _nbrSV << std::endl;

		_supVecs.reset(_nbrSV,_nbrVar);
		p = 0;
		for( i=0; i<prob.nbrData; ++i)
		{	
			if(nonZero[i]) 
			{	
				unsigned int c = _dataLabels[i];
				for ( j=0; j<_nbrVar; ++j ) {
					_supVecs[p][j] = prob.x[c][j]; // copy data row
				}
				++p;
			}
		}

		std::vector<int> nzStart(_nbrClasses);
		nzStart[0] = 0;
		for( i=1; i<_nbrClasses; ++i)
			nzStart[i] = nzStart[i-1]+nzCount[i-1];

		// reset coefficients
		_svCoef.reset(_nbrClasses-1,_nbrSV);
		
		p = 0;
		for( i=0; i<_nbrClasses; ++i)
		{
			for( j=i+1; j<_nbrClasses; ++j)
			{
				// classifier (i,j): coefficients with
				// i are in sv_coef[j-1][nz_start[i]...],
				// j are in sv_coef[i][nz_start[j]...]

				int si = classStart[i];
				int sj = classStart[j];
				int ci = classCount[i];
				int cj = classCount[j];
				
				int q = nzStart[i];
				for( k=0;k<ci;k++)
				{	
					if(nonZero[si+k])
					{	
						_svCoef[j-1][q] = _decFunctions[p].alpha[k];					
						++q;
					}
				}
				
				q = nzStart[j];
				for(k=0;k<cj;k++)
				{	
					if(nonZero[sj+k])
					{	
						_svCoef[i][q] = _decFunctions[p].alpha[ci+k];					
						++q;
					}
				}
				++p;
			}
		}
	}


	//std::cerr << "Support Vectors: " << _nbrSV << std::endl; 
	//	std::cerr << _supVecs << std::endl;
	//std::cerr << "Support Vector Coefficients:" << std::endl; 
	//	std::cerr << _svCoef << std::endl;
	
	return;
}


// Stratified cross validation
SiMath::Vector 
SiMath::SVM::crossValidation(SiMath::SVMProblem & prob, int nFold)
{
	int i;
		
	std::vector<int> foldStart(nFold);
	int l = prob.nbrData;
	SiMath::Vector target(l);
	
	
	// stratified cv may not give leave-one-out rate
	// Each class to l folds -> some folds may have zero elements
	if( ( _param.svmType == C_SVC || _param.svmType == NU_SVC) && nFold < l)
	{
		// group classes together (labels stored in _dataLabels
		std::vector<unsigned int> classStart;
		std::vector<unsigned int> classCount;
		_groupClasses(prob,classStart,classCount);
		
		// random shuffle and then data grouped by fold using the array perm
		std::vector<unsigned int> foldCount(nFold);
		int c;
		// int *index = Malloc(int,l);
		std::vector<int> index(_nbrSV);
		for(i=0; i<l; ++i)
			index[i]=_dataLabels[i];
		
		for ( c=0; c<_nbrClasses; ++c ) 
		{	
			for( i=0; i<classCount[c]; ++i)
			{
				int j = i+rand()%(classCount[c]-i);
				std::swap(index[classStart[c]+j],index[classStart[c]+i]);
			}
		}
		
		for(i=0;i<nFold;++i)
		{
			foldCount[i] = 0;
			for ( c=0; c<_nbrClasses; ++c )
				foldCount[i]+=(i+1)*classCount[c]/nFold-i*classCount[c]/nFold;
		}
				
		foldStart[0]=0;
		for (i=1; i<=nFold; ++i)
			foldStart[i] = foldStart[i-1]+foldCount[i-1];
		
		for ( c=0; c<_nbrClasses; ++c )
		{	
			for( i=0; i<nFold; ++i)
			{
				int begin = classStart[c]+i*classCount[c]/nFold;
				int end = classStart[c]+(i+1)*classCount[c]/nFold;
				for(int j=begin;j<end;++j)
				{
					_dataLabels[foldStart[i]] = index[j];
					foldStart[i]++;
				}
			}
		}
		foldStart[0]=0;
		for (i=1;i<=nFold;++i)
			foldStart[i] = foldStart[i-1]+foldCount[i-1];
	}
	else
	{
		for(i=0;i<l;++i) 
			_dataLabels[i] = i;
		
		for(i=0;i<l;++i)
		{
			int j = i+rand()%(l-i);
			std::swap(_dataLabels[i],_dataLabels[j]);
		}
		for(i=0;i<=nFold;++i)
			foldStart[i] = i*l/nFold;
	}
	
	
	for(i=0;i<nFold;++i)
	{
		int begin = foldStart[i];
		int end = foldStart[i+1];
		int j,k;
		SiMath::SVMProblem subProb;
		
		subProb.nbrData = l-(end-begin);
		subProb.nbrVar = prob.nbrVar;
		
		subProb.x.reset(subProb.nbrData,subProb.nbrVar);
		subProb.y.reset(subProb.nbrData);
		
		k=0;
		for(j=0;j<begin;++j)
		{
			for ( unsigned int jj=0; jj<_nbrVar; ++jj )
				subProb.x[k][jj] = prob.x[_dataLabels[j]][jj];
			subProb.y[k] = prob.y[_dataLabels[j]];
			++k;
		}
		for(j=end;j<l;++j)
		{
			for ( unsigned int jj=0; jj<_nbrVar; ++jj )
				subProb.x[k][jj] = prob.x[_dataLabels[j]][jj];
			subProb.y[k] = prob.y[_dataLabels[j]];
			++k;
		}
		// store model
		SVM submodel(_param);
		submodel.train(subProb);
		
		if(_param.probability && ( _param.svmType == C_SVC || _param.svmType == NU_SVC) )
		{
			for(j=begin;j<end;++j)
			{	
				unsigned int idx(_dataLabels[j]);
				target[idx] = _predictProbability(prob.x[idx]);
			}
		}
		else
		{	
			for(j=begin;j<end;++j)
			{	
				unsigned int idx(_dataLabels[j]);
				target[idx] = _predict(prob.x[idx]);
			}
		}
	}		
	
	return target;
}



void
SiMath::SVM::_trainOne(SiMath::SVMProblem & prob, const SiMath::SVMParameters & param, double Cp, double Cn, SiMath::DecisionFunction & f)
{
	SiMath::SolutionInfo si;
	f.alpha.clear();
	f.alpha.resize(prob.nbrData);
	f.rho = 0.0;
	switch(param.svmType)
	{
		case C_SVC:
			std::cerr << "C_SVC" << std::endl;
			si = solveCSVC(prob,param,f.alpha,Cp,Cn);
			break;
		case NU_SVC:
			std::cerr << "NU_SVC" << std::endl;
			si = solveNuSVC(prob,param,f.alpha);
			break;
		case ONE_CLASS:
			std::cerr << "ONE_CLASS" << std::endl;
			si = solveOneClass(prob,param,f.alpha);
			break;
		case EPSILON_SVR:
			std::cerr << "EPSILON_SVR" << std::endl;
			si = solveEpsilonSVR(prob,param,f.alpha);
			break;
		case NU_SVR:
			std::cerr << "NU_SVR" << std::endl;
			si = solveNuSVR(prob,param,f.alpha);
			break;
	}
	//std::cerr << "resulting alpha: " << f.alpha << std::endl;
	//std::cerr << "obj = " << si.obj << " rho = " << si.rho << std::endl;
	
	// output SVs
	
	int nSV = 0;
	int nBSV = 0;
	for(int i=0; i<prob.nbrData; ++i)
	{
		if(fabs(f.alpha[i]) > 0)
		{
			++nSV;
			if(prob.y[i] > 0)
			{
				if(fabs(f.alpha[i]) >= si.pUpperBound)
					++nBSV;
			}
			else
			{
				if(fabs(f.alpha[i]) >= si.nUpperBound)
					++nBSV;
			}
		}
	}
	
	std::cerr << "Training done: " << std::endl;
	std::cerr << " => number of Support Vectors = " << nSV << std::endl;
	std::cerr << " => number of Bounded Support Vectors = " << nBSV << std::endl;
	
	f.rho = si.rho;
	return;
}



// process the input data to find the class labels
void 
SiMath::SVM::_groupClasses(SiMath::SVMProblem & prob, std::vector<unsigned int> & classStart, std::vector<unsigned int> & classCount)
{
	classStart.clear();
	classCount.clear();
	std::vector<int> dataLabel(_nbrSV,0);
	unsigned int counter = 0;
	for ( unsigned int i=0; i<_nbrSV; ++i)
	{
		int l = (int)prob.y[i];
		if ( _classLabels.count(l) == 0 )  // check if class label already exists
		{	
			_classLabels[l] = counter;
			classCount.push_back(0);
			++counter;
		}
		unsigned int ll(_classLabels[l]);
		dataLabel[i] = ll;
		classCount[ll] += 1;
	}
	// we have the number of classes 
	_nbrClasses = counter;
	
	classStart.push_back(0);
	for ( unsigned int i=1; i<_nbrClasses; ++i )
		classStart.push_back(classStart[i-1] + classCount[i-1]);
	for ( unsigned int i=0; i<_nbrSV; ++i)
	{
		_dataLabels[classStart[dataLabel[i]]] = i;
		classStart[dataLabel[i]] += 1; 
	}
	classStart[0] = 0;
	for ( unsigned int i=1; i<_nbrClasses; ++i )
		classStart[i] = classStart[i-1] + classCount[i-1];
	
	return;
}



// Cross-validation decision values for probability estimates
void 
SiMath::SVM::_binarySVCProbability(SVMProblem & prob, const SVMParameters & param,
																	double Cp, double Cn, double & probA, double & probB)
{
	unsigned int i;
	int nbrFold = 5; // HARD CODED !!!!!!!!!!!!!
	int l = prob.nbrData;
	std::vector<unsigned int> permutation(l);
	
	SiMath::Vector predValues(l,0.0);
	
	// random shuffle the indices from 0 to l
	for(i=0; i<l; ++i) 
		permutation[i]=i;
	srand(rand());
	random_shuffle(permutation.begin(),permutation.end());
	
	for(i=0; i<nbrFold; ++i)
	{
		int begin = i*l/nbrFold;
		int end = (i+1)*l/nbrFold;
		int j,k;
		
		// create a new problem
		SVMProblem subProb;
		subProb.nbrData = l-(end-begin);
		subProb.nbrVar = prob.nbrVar;
		subProb.x.reset(subProb.nbrData,subProb.nbrVar);
		subProb.y.reset(subProb.nbrData);
		
		// copy pointers to data rows
		k=0;
		for(j=0;j<begin;++j)
		{
			for ( unsigned int jj=0; jj<prob.nbrVar; ++jj)
				subProb.x[k][jj] = prob.x[permutation[j]][jj];
			subProb.y[k] = prob.y[permutation[j]];
			++k;
		}
		for(j=end;j<l;++j)
		{
			for ( unsigned int jj=0; jj<prob.nbrVar; ++jj)
				subProb.x[k][jj] = prob.x[permutation[j]][jj];
			subProb.y[k] = prob.y[permutation[j]];
			++k;
		}
		
		int p_count=0,n_count=0;
		for(j=0;j<k;++j)
			if(subProb.y[j]>0)
				p_count++;
			else
				n_count++;
		
		if(p_count==0 && n_count==0)
			for(j=begin;j<end;++j)
				predValues[permutation[j]] = 0;
		else if(p_count > 0 && n_count == 0)
			for(j=begin;j<end;++j)
				predValues[permutation[j]] = 1;
		else if(p_count == 0 && n_count > 0)
			for(j=begin;j<end;++j)
				predValues[permutation[j]] = -1;
		else
		{
			SVMParameters subparam = param;
			subparam.probability=0;
			subparam.C=1.0;
			
			// update weights of classes
			subparam.weights.clear();
			subparam.weights[1] = Cp;
			subparam.weights[-1] = Cn;
			
			// create a new model
			SVM submodel(subparam);
			
			// train the model
			submodel.train(subProb);
			
			for(j=begin;j<end;++j)
			{
				predValues[permutation[j]] = submodel._predictValue(prob.x[permutation[j]]); 
				// ensure +1 -1 order; reason not using CV subroutine
				predValues[permutation[j]] *= submodel._classLabels[0];
			}
		}
	}		
	SiMath::fitSigmoid(predValues, prob.y, probA, probB);

	return;
}


double 
SiMath::SVM::_predict(const double * xRow)
{
	if( _param.svmType == ONE_CLASS || _param.svmType == EPSILON_SVR || _param.svmType == NU_SVR)
	{
		double * coef0 = _svCoef[0];
		double sum(0.0);
		double d(0.0);
		for( int i=0; i<_nbrSV; ++i)
		{	
			d = Kernel::function(_nbrVar,xRow,_supVecs[i],_param);
			sum += (coef0[i] * d);
		}
		sum -= _decFunctions[0].rho;
		
		if(_param.svmType == ONE_CLASS)
			return (sum>0) ? 1 : -1; // return class label
		else
			return sum;              // return value
	}
	else
	{
		int i;
		SiMath::Vector decValues = _predictValues(xRow);
		// std::cerr << "prediction: " << decValues << std::endl;
		
		std::vector<int> vote(_nbrClasses,0);
		int pos=0;
		for(i=0;i<_nbrClasses;++i)
		{	
			for(int j=i+1;j<_nbrClasses;++j)
			{
				if(decValues[pos] > 0)
					++vote[i];
				else
					++vote[j];
				++pos;
			}
		}	
		// std::cerr << "votes: ";
		//for ( i=0; i<vote.size(); ++i )
			// std::cerr << vote[i];
		//std::cerr << std::endl;
		
		int maxVote = 0;
		for( i=1;i<_nbrClasses;++i)
		{	
			if(vote[i] > vote[maxVote])
				maxVote = i;			
		}

		//for ( i=0; i<_nbrClasses; ++i)
		//
		for ( std::map<int,unsigned int>::iterator mi=_classLabels.begin(); mi!=_classLabels.end(); ++mi)
		{	
			if ( mi->second == maxVote )
				return mi->first;
		}
		// if we end up here, we have not found the class label
		throw(SVMError("Class label not found"));
	}
}


double 
SiMath::SVM::_predictValue(const double * x)
{
	double * coef = _svCoef[0];
	double sum = 0;
	for(int i=0;i<_nbrSV;++i)
		sum += coef[i] * Kernel::function(_nbrVar,x,_supVecs[i],_param);
	sum -= _decFunctions[0].rho;
	return sum;
}


SiMath::Vector 
SiMath::SVM::_predictValues(const double * x)
{
	int i;
	SiMath::Vector dv(_nbrClasses * (_nbrClasses-1)/2);

	std::vector<double> kvalue(_nbrSV);
	for(i=0; i<_nbrSV; ++i)
	{	
		kvalue[i] = Kernel::function(_nbrVar, x, _supVecs[i], _param);
		// std::cerr << " " << i << ":" << kvalue[i];
	}
	// std::cerr << std::endl;
	
	std::vector<unsigned int> start(_nbrClasses);
	start[0] = 0;
	for(i=1; i<_nbrClasses; ++i)
	{	
		start[i] = start[i-1]+_nbrClassSV[i-1];
		// std::cerr << "start " << i << " = " << start[i] << std::endl;
	}
	// std::cerr << _svCoef << std::endl;
	int p=0;
	for( i=0; i<_nbrClasses; ++i)
	{	
		for( int j=i+1; j<_nbrClasses; ++j)
		{
			double sum = 0;
			int si = start[i];
			int sj = start[j];
			int ci = _nbrClassSV[i];
			int cj = _nbrClassSV[j];
			int k;
			double * coef1 = _svCoef[j-1];
			double * coef2 = _svCoef[i];
			for(k=0;k<ci;k++)
				sum += coef1[si+k] * kvalue[si+k];
			// std::cerr << "sum = " << sum << " ";
			for(k=0;k<cj;k++)
				sum += coef2[sj+k] * kvalue[sj+k];
			//std::cerr << "sum = " << sum << " ";
			sum -= _decFunctions[p].rho;
			// std::cerr << "sum = " << sum << " ";
			dv[p] = sum;
			++p;
		}
	}
	// std::cerr << dv << std::endl;
	return dv;
}



double 
SiMath::SVM::_predictProbability(const double * xRow)
{
	if ( (_param.svmType == C_SVC || _param.svmType == NU_SVC) && _probA.size() && _probB.size() )
	{
		int i;
		SiMath::Vector dv = _predictValues(xRow);
		
		double min_prob=1e-7;

		unsigned int counter(0);
		SiMath::Matrix pairwiseProb(_nbrClasses,_nbrClasses,0.0);
		for(i=0;i<_nbrClasses;++i)
		{	
			for(int j=i+1;j<_nbrClasses;++j)
			{
				double d = max(min_prob,SiMath::sigmoid(dv[counter],_probA[counter],_probB[counter]));
				pairwiseProb[i][j] = min(d,1-min_prob);
				pairwiseProb[j][i] = 1-pairwiseProb[i][j];
				counter++;
			}
		}
	
		SiMath::Vector probEstimates = _multiclassProbability(pairwiseProb);
		
		int prob_max_idx = 0;
		for( i=1; i<_nbrClasses; ++i)
		{	
			if(probEstimates[i] > probEstimates[prob_max_idx])
				prob_max_idx = i;
		}
		
		// return class labels
		for ( std::map<int,unsigned int>::iterator mi=_classLabels.begin(); mi!=_classLabels.end(); ++mi)
		{	
			if ( mi->second == prob_max_idx )
				return mi->first;
		}
		return 0.0;
	}
	else 
	{
		return _predict(xRow);
	}
}



SiMath::Vector
SiMath::SVM::_multiclassProbability(SiMath::Matrix & r)
{
	int t,j;
	int iter = 0, max_iter=100;
	SiMath::Matrix Q(_nbrClasses,_nbrClasses,0.0);
	SiMath::Vector Qp(_nbrClasses,0.0);
	double pQp, eps=0.005/_nbrClasses;

	SiMath::Vector p(_nbrClasses,1.0/_nbrClasses);
	for (t=0;t<_nbrClasses;t++)
	{
		for (j=0;j<t;++j)
		{
			Q[t][t]+=r[j][t]*r[j][t];
			Q[t][j]=Q[j][t];
		}
		for (j=t+1;j<_nbrClasses;++j)
		{
			Q[t][t]+=r[j][t]*r[j][t];
			Q[t][j]=-r[j][t]*r[t][j];
		}
	}
	for (iter=0;iter<max_iter;iter++)
	{
		// stopping condition, recalculate QP,pQP for numerical accuracy
		pQp=0;
		for (t=0;t<_nbrClasses;t++)
		{
			Qp[t]=0;
			for (j=0;j<_nbrClasses;++j)
				Qp[t]+=Q[t][j]*p[j];
			pQp+=p[t]*Qp[t];
		}
		double max_error=0;
		for (t=0;t<_nbrClasses;t++)
		{
			double error=fabs(Qp[t]-pQp);
			if (error>max_error)
				max_error=error;
		}
		if (max_error<eps) break;
		
		for (t=0;t<_nbrClasses;t++)
		{
			double diff=(-Qp[t]+pQp)/Q[t][t];
			p[t]+=diff;
			pQp=(pQp+diff*(diff*Q[t][t]+2*Qp[t]))/(1+diff)/(1+diff);
			for (j=0;j<_nbrClasses;++j)
			{
				Qp[j]=(Qp[j]+diff*Q[t][j])/(1+diff);
				p[j]/=(1+diff);
			}
		}
	}
	if (iter>=max_iter)
		std::cerr << "SVM::_multiclassProbability(): maximum iterations exceeded." << std::endl;
	
	return p;
}


// Return parameter of a Laplace distribution 
double 
SiMath::SVM::_svrProbability(SiMath::SVMProblem & prob, const SiMath::SVMParameters & param)
{
	int i;
	int nFold = 5; // hard-coded
	int l = prob.nbrData;
	double mae = 0;
	
	// store old parameter setting
	bool oldProb = _param.probability;
	_param.probability = false;
	SiMath::Vector ymv = crossValidation(prob, nFold);

	// reset old parameters
	_param.probability = oldProb;
	
	for(i=0;i<l;++i)
	{
		ymv[i] = prob.y[i] - ymv[i];
		mae += fabs(ymv[i]);
	}		
	mae /= l;
	
	double std=sqrt(2*mae*mae);
	int count=0;
	mae=0;
	for(i=0;i<l;++i)
	{	
		if (fabs(ymv[i]) > 5*std) 
			count=count+1;
		else 
			mae+=fabs(ymv[i]);
	}
	mae /= (l-count);
	
	std::cerr << "Prob. model for test data: target value = predicted value + z,\nz: Laplace distribution e^(-|z|/sigma)/(2sigma),sigma= " << mae << std::endl;
	
	return mae;
}




//------------------------------------------------------------------------------------------------------------------------------------
// construct and solve various formulations
//------------------------------------------------------------------------------------------------------------------------------------

SiMath::SolutionInfo 
SiMath::SVM::solveCSVC(SiMath::SVMProblem & prob, const SiMath::SVMParameters & param, std::vector<double> & alpha,  double Cp, double Cn)
{
	SiMath::SolutionInfo si;
	
	// initial values
	unsigned int l = prob.nbrData;
	SiMath::Vector minusOnes(l,-1.0);
	std::vector<int> y(l);
	for( unsigned int i=0; i<l; ++i)
	{
		alpha[i] = 0;
		minusOnes[i] = -1;
		y[i] = (prob.y[i] > 0) ? +1 : -1;
	}

	// solve the equations
	Solver s;
	SiMath::svcQMatrix QM(prob, param);
	s.solve(l, QM, minusOnes, y, alpha, Cp, Cn, param.eps, si, param.shrinking);
		
	// some info
	if (Cp==Cn)
	{	
		double alphaSum(0);
		for( unsigned int i=0; i<l; ++i)
			alphaSum += alpha[i];
		std::cerr << "nu = " << alphaSum/(Cp * prob.nbrData) << std::endl;
	}

	// update alpha's
	for( unsigned int i=0; i<l; ++i )
		alpha[i] *= y[i];

	return si;
}


SiMath::SolutionInfo 
SiMath::SVM::solveNuSVC(SiMath::SVMProblem & prob, const SiMath::SVMParameters & param, std::vector<double> & alpha)
{
	SiMath::SolutionInfo si;
	int i;
	int l = prob.nbrData;
	double nu = param.nu;
	
	std::vector<int> y(l);
	
	for( i=0; i<l; ++i)
		y[i] = (prob.y[i] > 0) ? +1 : -1;
	
	double sum_pos = nu*l/2;
	double sum_neg = nu*l/2;
	
	for(i=0;i<l;++i)
	{	
		if(y[i] == +1)
		{
			alpha[i] = min(1.0,sum_pos);
			sum_pos -= alpha[i];
		}
		else
		{
			alpha[i] = min(1.0,sum_neg);
			sum_neg -= alpha[i];
		}
	}
	
	Vector zeros(l,0.0);
	
	nuSolver s;
	SiMath::svcQMatrix QM(prob, param);
	
	s.solve(l, QM, zeros, y, alpha, 1.0, 1.0, param.eps, si,  param.shrinking);
	double r = si.r;
	
	std::cerr << "C = " << 1/r << std::endl;
	
	for(i=0;i<l;++i)
		alpha[i] *= y[i]/r;
	
	si.rho /= r;
	si.obj /= (r*r);
	si.pUpperBound = 1/r;
	si.nUpperBound = 1/r;
	
	return si;
}


SiMath::SolutionInfo 
SiMath::SVM::solveOneClass(SiMath::SVMProblem & prob, const SiMath::SVMParameters & param, std::vector<double> & alpha)
{
	SiMath::SolutionInfo si;
	unsigned int l = prob.nbrData;
	SiMath::Vector zeros(l,0.0);
	std::vector<int> ones(l,1);
	
	unsigned int n = (unsigned int)(param.nu * l);	// # of alpha's at upper bound
	for(unsigned int i=0;i<n;++i)
	{	
		alpha[i] = 1;
	}

	if( n<l )
	{	
		alpha[n] = param.nu * l - n;
		for(unsigned int i=n+1;i<l;++i)
			alpha[i] = 0;		
	}
	
	Solver s;
	SiMath::oneClassQMatrix QM(prob,param);
	s.solve(l, QM, zeros, ones, alpha, 1.0, 1.0, param.eps, si, param.shrinking);
	
	return si;
}


SiMath::SolutionInfo 
SiMath::SVM::solveEpsilonSVR(SiMath::SVMProblem & prob, const SiMath::SVMParameters & param, std::vector<double> & alpha)
{
	SiMath::SolutionInfo si;
	unsigned int l = prob.nbrData;
	std::vector<double> alpha2(2*l,0.0);	
	std::vector<double> linearTerm(2*l,0.0);
	std::vector<int> y(2*l);
	
	// initialisation
	for(unsigned i=0;i<l;++i)
	{
		alpha2[i] = 0;
		linearTerm[i] = param.p - prob.y[i];
		y[i] = 1;
		
		alpha2[i+l] = 0;
		linearTerm[i+l] = param.p + prob.y[i];
		y[i+l] = -1;
	}
	
	Solver s;
	SiMath::svrQMatrix QM(prob,param);
	s.solve(2*l, QM, linearTerm, y, alpha2, param.C, param.C, param.eps, si, param.shrinking);
	
	double alphaSum = 0;
	for(unsigned int i=0;i<l;++i)
	{
		alpha[i] = alpha2[i] - alpha2[i+l];
		alphaSum += fabs(alpha[i]);
	}
	std::cerr << "nu = " << alphaSum/(param.C*l) << std::endl;
	
	return si;
}


SiMath::SolutionInfo 
SiMath::SVM::solveNuSVR(SiMath::SVMProblem & prob, const SiMath::SVMParameters & param, std::vector<double> & alpha)
{
	SiMath::SolutionInfo si;
	
	int l = prob.nbrData;
	double C = param.C;
	std::vector<double> alpha2(2*l,0.0);
	std::vector<double> linearTerm(2*l,0.0);
	std::vector<int> y(2*l);
	int i;
	
	double sum = C * param.nu * l / 2;
	for(i=0;i<l;++i)
	{
		alpha2[i] = alpha2[i+l] = min(sum,C);
		sum -= alpha2[i];
		
		linearTerm[i] = - prob.y[i];
		y[i] = 1;
		
		linearTerm[i+l] = prob.y[i];
		y[i+l] = -1;
	}
	
	nuSolver s;
	SiMath::svrQMatrix QM(prob,param);
	s.solve(2*l, QM, linearTerm, y, alpha2, C, C, param.eps, si, param.shrinking);
	
	std::cerr << "epsilon = " << -si.r << std::endl;
			 
	for(i=0;i<l;++i)
		alpha[i] = alpha2[i] - alpha2[i+l];
			 
	return si;
}


std::vector<int>
SiMath::SVM::getClassCount()
{
	return _nbrClassSV;
}


std::vector<int> 
SiMath::SVM::getClassLabels()
{
	std::vector<int> l(_nbrClasses);
	for ( std::map<int,unsigned int>::iterator mi=_classLabels.begin(); mi!=_classLabels.end(); ++mi)
	{	
		l[mi->second] = mi->first;
		std::cerr << mi->second << ":" << l[mi->second] << " "; 
	}
	std::cerr << std::endl;
	return l;
}

double
SiMath::SVM::getSVRProbability()
{
	if ( (_param.svmType == EPSILON_SVR || _param.svmType == NU_SVR) && _probA.size() )
	{	
		return _probA[0];
	}
	else
	{
		std::cerr << "WARNING: " << SVMTypeTable[_param.svmType] << " doesn't contain information for SVR probability inference." << std::endl;
		return 0.0;
	}
}


void
SiMath::SVM::setSupportVectors(const Matrix & sv, const Matrix & svc)
{
	if ( sv.nbrRows() != svc.nbrColumns() )
		throw(SVMError("Incompatible dimension of support vectors and coefficients."));

	_nbrSV = sv.nbrRows();
	_nbrVar = sv.nbrColumns();
	_nbrClasses = svc.nbrRows() + 1;
	
	if ( _nbrClasses == 2 )
	{
		_nbrClassSV.resize(1);
		_nbrClassSV[0] = _nbrSV;
	}

	// copy the matrices
	_supVecs = sv;
	_svCoef = svc;
	
}

void 
SiMath::SVM::setNbrOfClasses(unsigned int n)
{
	if ( n < 2 )
		throw(SVMError("Number of classes should be at least 2 (also for regression and one class problems)"));
	
	_nbrClasses = n;
	// clear data members related to number of classes
	_nbrClassSV.clear();
	_classLabels.clear();
}

void 
SiMath::SVM::setClassLabels(std::vector<int> & cl)
{
	if ( cl.size() != _nbrClasses )
		throw(SVMError("Number of class labels does not corresponds to number of classes in model."));
	
	_classLabels.clear();
	for ( unsigned int i=0; i<_nbrClasses; ++i)
		_classLabels[cl[i]] = i;
}

void 
SiMath::SVM::setClassCount(std::vector<int> & cn)
{
	if ( cn.size() != _nbrClasses )
		throw(SVMError("Number of class counts does not corresponds to number of classes in model."));
	
	_nbrClassSV.clear();
	for ( unsigned int i=0; i<_nbrClasses; ++i)
		_nbrClassSV.push_back(cn[i]);
}


void 
SiMath::SVM::setClassLabels(const std::vector<int> & cl, const std::vector<int> cn)
{
	if ( cl.size() != cn.size() || cl.size() != _nbrClasses || cn.size() != _nbrClasses )
		throw(SVMError("Number of class labels does not corresponds to number of classes in model."));
	
	_nbrClasses = cl.size();
	
	_classLabels.clear();  // clear map
	_nbrClassSV.clear();   // clear count vector
	for (unsigned int i=0; i<_nbrClasses; ++i)
	{	
		_classLabels[cl[i]] = i;      // add cl[i[],i to map
		_nbrClassSV.push_back(cn[i]);
	}
	
	return;
}


void 
SiMath::SVM::setPairWiseProbabilities(Vector & pA, Vector & pB)
{
	if( _param.svmType == ONE_CLASS || _param.svmType == EPSILON_SVR || _param.svmType == NU_SVR)
	{	
		if ( pA.size() != 1 && pB.size() != 0 )
		{
			throw(SVMError("Incompatible pairwise probabilities"));
		}
		else
		{
			_probA = pA;
			_probB.reset(0);
		}
	}
	else 
	{
		
		if ( pA.size() != pB.size() || pA.size() != _nbrClasses*(_nbrClasses-1)/2)
		{
			throw(SVMError("Incompatible pairwise probabilities"));
		}
		else
		{
			_probA = pA;
			_probB = pB;
		}
	}
	
	return;
}



SiMath::DecisionFunction 
SiMath::SVM::getDecisionFunction(unsigned int i) const
{
	if ( i<0 || i>_decFunctions.size() )
		throw(SVMError("Index out of bounds for decision functions."));
	
	return _decFunctions[i];
}


void 
SiMath::SVM::setDecisionFunction(unsigned int i, double rho, std::vector<double> & alpha)
{
	if ( i<0 || i>_decFunctions.size() )
		throw(SVMError("Index out of bounds for decision functions."));

	if ( i >= _decFunctions.size() )
		_decFunctions.resize(i+1);
	
	_decFunctions[i].rho = rho;
	_decFunctions[i].alpha = alpha;

}



