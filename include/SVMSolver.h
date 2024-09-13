/*
 *  SVMSolver.h
 *
 *  Created by Gert Thijs on 07/04/06.
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

#ifndef __SVMSOLVER_H__
#define __SVMSOLVER_H__

#include "Definitions.h"
#include "Vector.h"
#include "Kernel.h"
#include "SVM.h"

namespace SiMath
{
	
	/**
		\class Solver
		\brief Class to solve the dual constraint minimisation problem in SVM training
		
		Adapted from the code in libsvm version 2.8.1 (see also http://www.csie.ntu.edu.tw/~cjlin/libsvm/ )
	 */
	class Solver 
	{
		public:
			/** 
				\brief Default constructor
				*/
			Solver();
			
			
			virtual ~Solver();
		
			void solve(int l, Kernel & Q, const Vector & b, const std::vector<int> & y,
								 std::vector<double> & alpha, double Cp, double Cn, double eps,
								 SolutionInfo & si, int shrinking);
			
	
		protected:
			int _activeSize;
			std::vector<unsigned int> _activeSet;
			std::vector<int> _y;
			Vector _b;
			Kernel * _Q;
			const double * _QD;
			double _eps;
			double _Cp;
			double _Cn;
			Vector _G;       ///< storage of the gradient of the objective function
			Vector _Gbar;  ///< storage of the gradient of the objective function, with free variables as 0

			int _l;
			bool _unshrinked;	// XXX
		
			double _getC(int i) { return (_y[i] > 0) ? _Cp : _Cn; };
		
	
			// alpha status
		std::vector<double> _alpha;
			std::vector<AlphaStatus> _alphaStatus;	// LOWER_BOUND, UPPER_BOUND, FREE
	
			void _updateAlphaStatus(int i);
			bool _isUpperBound(unsigned int i) const { return _alphaStatus[i] == UPPER_BOUND; }
			bool _isLowerBound(unsigned int i) const { return _alphaStatus[i] == LOWER_BOUND; }
			bool _isFree(unsigned int i) const { return _alphaStatus[i] == FREE; }
	
	
			void _swapIndex(unsigned int i, unsigned int j);
			void _reconstructGradient();
			virtual int _selectWorkingSet(int &i, int &j);
			virtual int _maxViolatingPair(int &i, int &j);
			virtual double _calculateRho();
			virtual void _doShrinking();
			virtual bool _beShrunk(int i, double Gmax1, double Gmax2);
	};

//
// Solver for nu-svm classification and regression
//
// additional constraint: e^T \alpha = constant
//
	class nuSolver : public Solver
	{
		public:
			nuSolver() : Solver() {};
		
		private:
			int _selectWorkingSet(int &i, int &j);
			double _calculateRho(SolutionInfo & si);
			void _doShrinking();
		bool _beShrunk(int i, double Gmax1, double Gmax2, double Gmax3, double Gmax4);
	};

	
	
};


#endif __SVMSOLVER_H__
