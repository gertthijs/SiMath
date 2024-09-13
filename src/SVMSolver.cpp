/*
 *  SVMSolver.cpp
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

#include "SVMSolver.h"

//------------------------------------------------------------------------------------------------------------------------------------
// Solver Implementation
//------------------------------------------------------------------------------------------------------------------------------------
SiMath::Solver::Solver() :
	_activeSize(0),
	_activeSet(0),
	_y(0),
	_b(0),
	_Q(NULL),
	_QD(NULL),
	_eps(0.001),
	_Cp(1),
	_Cn(1),
	_G(),
	_Gbar(),
	_l(0),
	_unshrinked(false)
{ }

SiMath::Solver::~Solver()
{
	if ( _Q != NULL )
		_Q = NULL;
	
	//if ( _G != NULL )
		//delete[] _G;
	//if ( _Gbar != NULL )
		//delete[] _Gbar;
}

void 
SiMath::Solver::_swapIndex(unsigned int i, unsigned int j)
{
	if ( i == j ) return;
	
	_Q->swapIndex(i,j);
	std::swap(_y[i],_y[j]);
	std::swap(_G[i],_G[j]);	
	std::swap(_alphaStatus[i],_alphaStatus[j]);
	std::swap(_alpha[i],_alpha[j]);
	std::swap(_b[i],_b[j]);
	std::swap(_activeSet[i],_activeSet[j]);
	std::swap(_Gbar[i],_Gbar[j]);
}


void
SiMath::Solver::_updateAlphaStatus(int i)
{
	if ( _alpha[i] >= _getC(i))
	{
		_alphaStatus[i] = UPPER_BOUND;
	}
	else if(_alpha[i] <= 0)
	{
		_alphaStatus[i] = LOWER_BOUND;
	}
	else 
	{
		_alphaStatus[i] = FREE;
	}
	return;
}


void 
SiMath::Solver::_reconstructGradient()
{
	// reconstruct inactive elements of _G from _Gbar and free variables
	if(_activeSize == _l) 
		return;

	unsigned int i;
	for( i=_activeSize; i<_l; ++i)
	{	
		_G[i] = _Gbar[i] + _b[i];
	}

	int nbrFree(0);
	for( i=0; i<_activeSize;++i)
	{	
		if( _isFree(i) ) ++nbrFree;
	}
	
	if ( nbrFree*_l > 2 * _activeSize * (_l-_activeSize) ){
		for( i=_activeSize; i<_l; ++i ){
			const double * Q_i = _Q->getQ(i,_activeSize);
			for (int j=0; j<_activeSize; ++j ){
				if ( _isFree(j) ) _G[i] += _alpha[j] * Q_i[j];
			}
		} 
	}else{
		for( i=0; i<_activeSize;++i)
		{	
			if( _isFree(i) )
			{
				const double * Q_i = _Q->getQ(i,_l);
				for(int j=_activeSize; j<_l; ++j)
					_G[j] += _alpha[i] * Q_i[j];
			}
		}
	}
	
}




void 
SiMath::Solver::solve( int l, SiMath::Kernel & Q, const SiMath::Vector & b, const std::vector<int> & y,
											std::vector<double> & alpha, double Cp, double Cn, double eps, SiMath::SolutionInfo & si, int shrinking)
{
	_l = l;
	_Q = &Q;
	_QD = Q.getQD();
	// clone(_b, b, _l);
	//	clone(_y, y, _l);
	//	clone(_alpha,alpha, _l);
	_b = b;
	_y = y;
	_alpha = alpha;
	_Cp = Cp;
	_Cn = Cn;
	_eps = eps;
	_unshrinked = false;

	// initialize alpha_status
	{		
		_alphaStatus.resize(l);
		for(int i=0; i<_l; ++i)
			_updateAlphaStatus(i);
	}

	// initialize active set (for shrinking)
	{
		_activeSet.resize(_l);
		for( int i=0; i<_l; ++i )
			_activeSet[i] = i;
		_activeSize = _l;
	}

	// initialize gradient
	{
		_G.resize(_l);
		_Gbar.resize(_l);
		for( int i=0; i<_l; ++i)
		{
			_G[i] = b[i];
			_Gbar[i] = 0;
		}
		for( int i=0; i<_l; ++i)
		{	
			if(!_isLowerBound(i))
			{
				const double * Qi = Q.getQ(i,_l);
				double alpha_i = _alpha[i];
				for( int j=0; j<_l; ++j)
					_G[j] += alpha_i*Qi[j];

				if( _isUpperBound(i) )
					for( int j=0; j<_l; ++j )
						_Gbar[j] += _getC(i) * Qi[j];
			}
		}
	}

	// optimization step

	int iter = 0;
	int counter = min(_l,1000)+1;
	std::cerr << "Optimisation step: " << std::endl;
	while(1)
	{
		// show progress and do shrinking
		if(--counter == 0)
		{
			counter = min(_l,1000);
			if( shrinking ) _doShrinking();
			std::cerr << ".";
		}

		int i(0),j(0);		
		if( _selectWorkingSet(i,j) != 0 )
		{
			// reconstruct the whole gradient
			_reconstructGradient();
			// reset active set size and check
			_activeSize = l;
			std::cerr << "*";
			if( _selectWorkingSet(i,j)!=0)
				break;
			else
				counter = 1;	// do shrinking next iteration
		}
		++iter;
		
		// update alpha[i] and alpha[j], handle bounds carefully
		
		const double *Q_i = Q.getQ(i, _activeSize);
		const double *Q_j = Q.getQ(j, _activeSize);
		
		double C_i = _getC(i);
		double C_j = _getC(j);

		double iOldAlpha = _alpha[i];
		double jOldAlpha = _alpha[j];

		if( _y[i] != _y[j] )
		{
			double quad_coef = Q_i[i]+Q_j[j]+2*Q_i[j];
			if (quad_coef <= 0)
				quad_coef = TAU;
			double delta = (-_G[i]-_G[j])/quad_coef;
			double diff = _alpha[i] - _alpha[j];
			_alpha[i] += delta;
			_alpha[j] += delta;
			
			if(diff > 0)
			{
				if( _alpha[j] < 0 )
				{
					_alpha[j] = 0;
					_alpha[i] = diff;
				}
			}
			else
			{
				if( _alpha[i] < 0 )
				{
					_alpha[i] = 0;
					_alpha[j] = -diff;
				}
			}
			if(diff > C_i - C_j)
			{
				if(_alpha[i] > C_i)
				{
					_alpha[i] = C_i;
					_alpha[j] = C_i - diff;
				}
			}
			else
			{
				if(_alpha[j] > C_j)
				{
					_alpha[j] = C_j;
					_alpha[i] = C_j + diff;
				}
			}
			
		}
		else
		{
			double quad_coef = Q_i[i]+Q_j[j]-2*Q_i[j];
			if (quad_coef <= 0)
				quad_coef = TAU;
			double delta = (_G[i]-_G[j])/quad_coef;
			double sum = _alpha[i] + _alpha[j];
			_alpha[i] -= delta;
			_alpha[j] += delta;

			if(sum > C_i)
			{
				if(_alpha[i] > C_i)
				{
					_alpha[i] = C_i;
					_alpha[j] = sum - C_i;
				}
			}
			else
			{
				if(_alpha[j] < 0)
				{
					_alpha[j] = 0;
					_alpha[i] = sum;
				}
			}
			if(sum > C_j)
			{
				if(_alpha[j] > C_j)
				{
					_alpha[j] = C_j;
					_alpha[i] = sum - C_j;
				}
			}
			else
			{
				if(_alpha[i] < 0)
				{
					_alpha[i] = 0;
					_alpha[j] = sum;
				}
			}
		}
		
		// update G
		double iDeltaAlpha = _alpha[i] - iOldAlpha; // change of i
		double jDeltaAlpha = _alpha[j] - jOldAlpha; // change of j
		
		for( int k=0; k<_activeSize; k++ )
		{
			_G[k] += Q_i[k]*iDeltaAlpha + Q_j[k]*jDeltaAlpha;
		}
		
		// update _alphaStatus and _Gbar
		{
			bool ui = _isUpperBound(i);
			bool uj = _isUpperBound(j);
			
			_updateAlphaStatus(i);
			_updateAlphaStatus(j);
			
			int k;
			if( ui != _isUpperBound(i) )
			{
				Q_i = Q.getQ(i,l);
				if(ui)
				{	
					for(k=0;k<_l;++k)
						_Gbar[k] -= C_i * Q_i[k];					
				}
				else
				{	
					for(k=0;k<_l;++k)
						_Gbar[k] += C_i * Q_i[k];					
				}
			}

			if(uj != _isUpperBound(j))
			{
				Q_j = Q.getQ(j,_l);
				if(uj)
				{	
					for(k=0;k<_l;++k)
						_Gbar[k] -= C_j * Q_j[k];					
				}
				else
				{	
					for(k=0;k<_l;++k)
						_Gbar[k] += C_j * Q_j[k];					
				}
			}
		}
	}

	// std::cerr << "Solver::_alpha: " << _alpha;
	// calculate rho
	si.rho = _calculateRho();

	// calculate objective value
	{
		double v = 0;
		int i;
		for(i=0;i<l;++i)
			v += _alpha[i] * (_G[i] + _b[i]);

		si.obj = v/2;
	}

	// put back the solution
	{
		//std::cerr << "active set: " << std::endl;
		for(int i=0;i<l;++i)
		{
			//std::cerr << " " << i << ":" << _activeSet[i];
			alpha[_activeSet[i]] = _alpha[i];
		}
		//std::cerr << std::endl;
	}
	
	
	si.pUpperBound = Cp;
	si.nUpperBound = Cn;

	std::cerr << "positive upper bound: " << Cp << std::endl;
	std::cerr << "negative upper bound: " << Cn << std::endl;
	
	std::cerr << "Optimization finished after " << iter << " iterations." << std::endl;

	// cleanup local variable
	_G.clear();
	_Gbar.clear();
}

//----------------------------------------------------------------------------------------------------------
// return 1 if already optimal, return 0 otherwise
int 
SiMath::Solver::_selectWorkingSet(int & iOut, int & jOut)
{
	// return i,j such that
	// i: maximizes -y_i * grad(f)_i, i in I_up(\alpha)
	// j: minimizes the decrease of obj value
	//    (if quadratic coefficeint <= 0, replace it with tau)
	//    -y_j*grad(f)_j < -y_i*grad(f)_i, j in I_low(\alpha)
	
	double Gmax = -INF;
	double Gmax2 = -INF;
	int indexGMax = -1;
	int indexGMin = -1;
	double obj_diff_min = INF;

	for( int t=0; t<_activeSize; t++)
	{	
		if( _y[t]==+1)	
		{
			if( !_isUpperBound(t) && -_G[t] >= Gmax)
			{
				Gmax = -_G[t];
				indexGMax = t;
			}
		}
		else
		{
			if( !_isLowerBound(t) && _G[t] >= Gmax)
			{
				Gmax = _G[t];
				indexGMax = t;
			}
		}
	}
				
	int i = indexGMax;
	const double *Q_i = NULL;
	if(i != -1) // NULL Q_i not accessed: Gmax=-INF if i=-1
		Q_i = _Q->getQ(i,_activeSize);

	for( int j=0; j<_activeSize; ++j )
	{
		if( _y[j]==+1 )
		{
			if ( !_isLowerBound(j) )
			{
				double grad_diff(Gmax+_G[j]);
				if ( _G[j] >= Gmax2 )
					Gmax2 = _G[j];
				if (grad_diff > 0)
				{
					double obj_diff; 
					double quad_coef(Q_i[i]+_QD[j]-2*_y[i]*Q_i[j]);
					if (quad_coef > 0)
						obj_diff = -(grad_diff*grad_diff)/quad_coef;
					else
						obj_diff = -(grad_diff*grad_diff)/TAU;

					if (obj_diff <= obj_diff_min)
					{
						indexGMin=j;
						obj_diff_min = obj_diff;
					}
				}
			}
		}
		else
		{
			if ( !_isUpperBound(j) )
			{
				double grad_diff(Gmax - _G[j]);
				if ( -_G[j] >= Gmax2 )
					Gmax2 = -_G[j];
				if (grad_diff > 0 )
				{
					double obj_diff; 
					double quad_coef(Q_i[i]+ _QD[j]+2*_y[i]*Q_i[j]);
					if (quad_coef > 0)
						obj_diff = -(grad_diff*grad_diff)/quad_coef;
					else
						obj_diff = -(grad_diff*grad_diff)/TAU;

					if (obj_diff <= obj_diff_min)
					{
						indexGMin=j;
						obj_diff_min = obj_diff;
					}
				}
			}
		}
	}

	if(Gmax + Gmax2 < _eps)
 		return 1;

	iOut = indexGMax;
	jOut = indexGMin;
	return 0;
}

//-------------------------------------------------------------------------------------------------------
// return 1 if already optimal, return 0 otherwise
int 
SiMath::Solver::_maxViolatingPair(int & iOut, int & jOut)
{
	// return i,j: maximal violating pair

	double Gmax1 = -INF;		// max { -y_i * grad(f)_i | i in I_up(\alpha) }
	int Gmax1_idx = -1;

	double Gmax2 = -INF;		// max { y_i * grad(f)_i | i in I_low(\alpha) }
	int Gmax2_idx = -1;

	for( int i=0; i<_activeSize; ++i)
	{
		if( _y[i]==+1 )	// y = +1
		{
			if( !_isUpperBound(i) )	// d = +1
			{
				if(-_G[i] >= Gmax1)
				{
					Gmax1 = -_G[i];
					Gmax1_idx = i;
				}
			}
			if( !_isLowerBound(i) )	// d = -1
			{
				if(_G[i] >= Gmax2)
				{
					Gmax2 = _G[i];
					Gmax2_idx = i;
				}
			}
		}
		else		// y = -1
		{
			if( !_isUpperBound(i) )	// d = +1
			{
				if( -_G[i] >= Gmax2 )
				{
					Gmax2 = -_G[i];
					Gmax2_idx = i;
				}
			}
			if( ! _isLowerBound(i) )	// d = -1
			{
				if( _G[i] >= Gmax1 )
				{
					Gmax1 = _G[i];
					Gmax1_idx = i;
				}
			}
		}
	}

	if(Gmax1+Gmax2 < _eps)
 		return 1;

	iOut = Gmax1_idx;
	jOut = Gmax2_idx;
	return 0;
}


//--------------------------------------------------------------------------------------------
bool
SiMath::Solver::_beShrunk(int i, double Gmax1, double Gmax2)
{
	if ( _isUpperBound(i) ){
		if ( _y[i] == + 1 )
			return ( -_G[i] > Gmax1 );
		else
			return ( -_G[i] > Gmax2 );
	}else if ( _isLowerBound(i) ){
		if ( _y[i] == + 1 )
			return ( _G[i] > Gmax2 );
		else
			return ( _G[i] > Gmax1 );
	}else{
		return false;
	}	
}

void 
SiMath::Solver::_doShrinking()
{
	int i,j,k;

	double Gm1 = -INF; 
	double Gm2 = -INF;

	// find maximal violating point
	for( k=0; k<_activeSize; k++ )
	{
		if( _y[k]==+1 ){
			if( !_isUpperBound(k) )
			{
				if( -_G[k] >= Gm1)
					Gm1 =  -_G[k];
			}
			if( !_isLowerBound(k) )
			{
				if( _G[k] >= Gm2)
					Gm2 =  _G[k];
			}
		}else{
			if( !_isUpperBound(k) )
			{
				if( -_G[k] >= Gm2)
					Gm2 =  -_G[k];
			}
			if( !_isLowerBound(k) )
			{
				if( _G[k] >= Gm1)
					Gm1 =  _G[k];
			}
		}
	}

	// unshrink, check all variables again before final iterations
	if( _unshrinked == false && (Gm1 + Gm2) <= _eps*10) {
		_unshrinked = true;
		_reconstructGradient();
		_activeSize = _l;
		std::cerr << "*";
	}
		
	for( k=0; k<_activeSize; ++k )
	{
		_activeSize--;
		if ( ! _beShrunk(_activeSize,Gm1,Gm2) ){
			_swapIndex(k, _activeSize);
			break;
		}
	}

}



//----------------------------------------------------------------------
double 
SiMath::Solver::_calculateRho()
{
	int nbrFree = 0;
	double ub = INF, lb = -INF, sumFree = 0;
	for( unsigned int i=0; i<_activeSize; ++i )
	{
		double yG = _y[i]*_G[i];

		if( _isLowerBound(i) )
		{
			if( _y[i] > 0)
				ub = min(ub,yG);
			else
				lb = max(lb,yG);
		}
		else if( _isUpperBound(i) )
		{
			if( _y[i] < 0 )
				ub = min(ub,yG);
			else
				lb = max(lb,yG);
		}
		else
		{
			++nbrFree;
			sumFree += yG;
		}
	}

	// return rho
	if( nbrFree>0)
		return sumFree/nbrFree;
	else
		return (ub+lb)/2;
}




//----------------------------------------------------------------------
// SOLVER-NU 
//----------------------------------------------------------------------

// return 1 if already optimal, return 0 otherwise
int 
SiMath::nuSolver::_selectWorkingSet(int &iOut, int &jOut)
{
	// return i,j such that y_i = y_j and
	// i: maximizes -y_i * grad(f)_i, i in I_up(\alpha)
	// j: minimizes the decrease of obj value
	//    (if quadratic coefficeint <= 0, replace it with tau)
	//    -y_j*grad(f)_j < -y_i*grad(f)_i, j in I_low(\alpha)

	double Gmaxp = -INF;
	double Gmaxp2 = -INF;
	int indexGMaxP = -1;

	double Gmaxn = -INF;
	double Gmaxn2 = -INF;
	int indexGMaxN = -1;

	int indexGMin = -1;
	double obj_diff_min = INF;

	for( int t=0; t<_activeSize; ++t )
	{	
		if(_y[t]==+1)
		{
			if( !_isUpperBound(t) && -_G[t] >= Gmaxp )
			{
				Gmaxp = -_G[t];
				indexGMaxP = t;
			}
		}
		else
		{
			if( !_isLowerBound(t) && _G[t] >= Gmaxn)
			{
				Gmaxn = _G[t];
				indexGMaxN = t;
			}
		}
	}
	int ip = indexGMaxP;
	int in = indexGMaxN;
	
	double * Q_ip = NULL;
	double * Q_in = NULL;
	if(ip != -1) // NULL Q_ip not accessed: Gmaxp=-INF if ip=-1
		Q_ip = _Q->getQ(ip, _activeSize);
	if(in != -1)
		Q_in = _Q->getQ(in, _activeSize);

	for(int j=0;j<_activeSize; ++j)
	{
		if( _y[j]==+1 )
		{
			if (!_isLowerBound(j) )	
			{
				double grad_diff(Gmaxp+_G[j]);
				if ( _G[j] >= Gmaxp2 )
					Gmaxp2 = _G[j];
				if ( grad_diff > 0 )
				{
					double obj_diff; 
					double quad_coef(Q_ip[ip]+ _QD[j]-2*Q_ip[j]);
					if (quad_coef > 0)
						obj_diff = -(grad_diff*grad_diff)/quad_coef;
					else
						obj_diff = -(grad_diff*grad_diff)/TAU;

					if (obj_diff <= obj_diff_min)
					{
						indexGMin=j;
						obj_diff_min = obj_diff;
					}
				}
			}
		}
		else
		{
			if ( !_isUpperBound(j) )
			{
				double grad_diff(Gmaxn-_G[j]);
				if ( -_G[j] >= Gmaxn2 )
					Gmaxn2 = -_G[j];
				if (grad_diff > 0)
				{
					double obj_diff; 
					double quad_coef(Q_in[in] + _QD[j]-2*Q_in[j]);
					if (quad_coef > 0)
						obj_diff = -(grad_diff*grad_diff)/quad_coef;
					else
						obj_diff = -(grad_diff*grad_diff)/TAU;

					if (obj_diff <= obj_diff_min)
					{
						indexGMin=j;
						obj_diff_min = obj_diff;
					}
				}
			}
		}
	}
	
	if( max(Gmaxp+Gmaxp2,Gmaxn+Gmaxn2) < _eps )
 		return 1;

	if (_y[indexGMin] == +1)
		iOut = indexGMaxP;
	else
		iOut = indexGMaxN;
	jOut = indexGMin;

	return 0;
}


bool
SiMath::nuSolver::_beShrunk(int i, double Gmax1, double Gmax2, double Gmax3, double Gmax4)
{
	if ( _isUpperBound(i) ){
		if ( _y[i] == +1 )
			return ( -_G[i] > Gmax1 );
		else
			return ( -_G[i] > Gmax4 );
	}else if ( _isLowerBound(i) ){
		if ( _y[i] == +1 )
			return ( _G[i] > Gmax2 );
		else
			return ( _G[i] > Gmax3 );
	}else{
		return false;
	}	
}



void 
SiMath::nuSolver::_doShrinking()
{
	double Gmax1 = -INF;	// max { -y_i * grad(f)_i | y_i = +1, i in I_up(\alpha) }
	double Gmax2 = -INF;	// max { y_i * grad(f)_i | y_i = +1, i in I_low(\alpha) }
	double Gmax3 = -INF;	// max { -y_i * grad(f)_i | y_i = -1, i in I_up(\alpha) }
	double Gmax4 = -INF;	// max { y_i * grad(f)_i | y_i = -1, i in I_low(\alpha) }

	// find maximal violating pair first
	int k;
	for( k=0; k<_activeSize; k++)
	{
		if( !_isUpperBound(k) )
		{
			if( _y[k]==+1 ){
				if(-_G[k] > Gmax1) Gmax1 = -_G[k];
			}else{
				if(-_G[k] > Gmax4) Gmax4 = -_G[k];	
			}
		}
		if(!_isLowerBound(k))
		{
			if( _y[k]==+1 ){	
				if(_G[k] > Gmax2) Gmax2 = _G[k];
			}else{
				if(_G[k] > Gmax3){
					Gmax3 = _G[k];
				} 	
			}
		}
	}
	
	if ( _unshrinked == false && max(Gmax1+Gmax2,Gmax3+Gmax4) <= _eps*10){
		_unshrinked = true;
		_reconstructGradient();
		_activeSize = _l;
	}
	
	for( k=0; k<_activeSize; ++k )
	{
		if ( _beShrunk(k, Gmax1,Gmax2,Gmax3,Gmax4) ){
			--_activeSize;
			while ( _activeSize > k ){
				if ( !_beShrunk(_activeSize, Gmax1,Gmax2,Gmax3,Gmax4) ){
					_swapIndex(k,_activeSize);
					break;
				}
				_activeSize--;
			}
		}
	}

}


double 
SiMath::nuSolver::_calculateRho(SiMath::SolutionInfo & si)
{
	int nbrFree1(0), nbrFree2(0);
	double ub1 = INF, 
		ub2 = INF,  // upper
		lb1 = -INF, 
		lb2 = -INF; // lower
	double sumFree1 = 0, sumFree2 = 0;

	for( int i=0; i<_activeSize; ++i)
	{
		if( _y[i]==+1 )
		{
			if( _isLowerBound(i) )
				ub1 = min(ub1,_G[i]);
			else if( _isUpperBound(i) )
				lb1 = max(lb1,_G[i]);
			else
			{
				++nbrFree1;
				sumFree1 += _G[i];
			}
		}
		else
		{
			if( _isLowerBound(i) )
				ub2 = min(ub2,_G[i]);
			else if( _isUpperBound(i) )
				lb2 = max(lb2,_G[i]);
			else
			{
				++nbrFree2;
				sumFree2 += _G[i];
			}
		}
	}

	double r1,r2;
	if(nbrFree1 > 0)
		r1 = sumFree1/nbrFree1;
	else
		r1 = (ub1+lb1)/2;
	
	if(nbrFree2 > 0)
		r2 = sumFree2/nbrFree2;
	else
		r2 = (ub2+lb2)/2;
	
	si.r = (r1+r2)/2;
	return (r1-r2)/2;
}




