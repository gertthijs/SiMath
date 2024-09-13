/*
 *  SVD.cpp
 *
 *  Created by Gert Thijs on 9/08/06.
 *  Copyright 2006 Silicos. All rights reserved.
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

#include "SVD.h"
#include "Utilities.h"


SiMath::SVD::SVD (const Matrix & Aorig, bool bU, bool bV) :
	_m(Aorig.nbrRows()),
	_n(Aorig.nbrColumns()),
	_U(),
	_V(),
	_S(0),
	_computeV(bV),
	_computeU(bU)
{
	// dimensionality of the problem
	int nu = min(_m,_n);
	int nct = min(_m-1,_n);
	int nrt = max(0,min(_n-2,_m));

	// define the dimensions of the internal matrices and vetors
	_S.reset(min(_m+1,_n)); 
	
	if ( _computeU )
		_U.reset(_m,nu);

	if ( _computeV )
		_V.reset(_n,_n);

	// local working vectors
	Vector e(_n);
	Vector work(_m);
	
	// make a copy of A to do the computations on
	Matrix Acopy(Aorig);

	// loop indices
	int i=0, j=0, k=0;
	
	// Reduce A to bidiagonal form, storing the diagonal elements
	// in _S and the super-diagonal elements in e.
	
	for (k = 0; k < max(nct,nrt); k++)
	{
		if (k < nct) 
		{
			// Compute the transformation for the k-th column and place the k-th diagonal in _S[k].
			_S[k] = 0;
			for (i = k; i < _m; i++) {
				_S[k] = triangle(_S[k],Acopy[i][k]);
			}
			if (_S[k] != 0.0) {
				if (Acopy[k][k] < 0.0) {
					_S[k] = -_S[k];
				}
				for (i = k; i < _m; i++) {
					Acopy[i][k] /= _S[k];
				}
				Acopy[k][k] += 1.0;
			}
			_S[k] = -_S[k];
		}
		for (j = k+1; j < _n; j++) 
		{
			if ((k < nct) && (_S[k] != 0.0))
			{
				// Apply the transformation to Acopy
				double t = 0;
				for (i = k; i < _m; i++) {
					t += Acopy[i][k]*Acopy[i][j];
				}
				t = -t/Acopy[k][k];
				for (i = k; i < _m; i++) {
					Acopy[i][j] += t*Acopy[i][k];
				}
			}
			
			// Place the k-th row of A into e for the subsequent calculation of the row transformation.
			e[j] = Acopy[k][j];
		}
		
		// Place the transformation in _U for subsequent back multiplication.
		if ( _computeU & (k < nct) ) 
		{			
			for (i = k; i < _m; i++) 
			{
				_U[i][k] = Acopy[i][k];
			}
		}
		
		if ( k < nrt )
		{
			// Compute the k-th row transformation and place the k-th super-diagonal in e[k].
			// Compute 2-norm without under/overflow.
			e[k] = 0.0;
			for (i = k+1; i < _n; i++) {
				e[k] = triangle(e[k],e[i]);
			}
			if (e[k] != 0.0) 
			{
				if (e[k+1] < 0.0) { // switch sign
					e[k] = -e[k]; 
				}
				for (i = k+1; i < _n; i++) { // scale
					e[i] /= e[k];
				}
				e[k+1] += 1.0;
			}
			e[k] = -e[k]; 
			if ((k+1 < _m) & (e[k] != 0.0))
			{
				// Apply the transformation.
				
				for (i = k+1; i < _m; i++) {
					work[i] = 0.0;
				}
				for (j = k+1; j < _n; j++) {
					for (i = k+1; i < _m; i++) {
						work[i] += e[j]*Acopy[i][j];
					}
				}
				for (j = k+1; j < _n; j++) {
					double t = -e[j]/e[k+1];
					for (i = k+1; i < _m; i++) {
						Acopy[i][j] += t*work[i];
					}
				}
			}
			
			
			// Place the transformation in _V for subsequent back multiplication.			
			if ( _computeV ) 
			{
				for (i = k+1; i < _n; i++) 
				{
					_V[i][k] = e[i];
				}
			}
		}
	}
	
	// Set up the final bidiagonal matrix of order p.
	int p = min(_n,_m+1);
	if (nct < _n) {
		_S[nct] = Acopy[nct][nct];
	}
	if (_m < p) {
		_S[p-1] = 0.0;
	}
	if (nrt+1 < p) {
		e[nrt] = Acopy[nrt][p-1];
	}
	e[p-1] = 0.0;
	
	// If required, generate U.
	if ( _computeU ) 
	{
		for (j = nct; j < nu; j++) {
			for (i = 0; i < _m; i++) {
				_U[i][j] = 0.0;
			}
			_U[j][j] = 1.0;
		}
		for (k = nct-1; k >= 0; k--) 
		{
			if (_S[k] != 0.0) 
			{
				for (j = k+1; j < nu; j++) 
				{
					double t = 0;
					for (i = k; i < _m; i++) {
						t += _U[i][k]*_U[i][j];
					}
					t = -t/_U[k][k];
					for (i = k; i < _m; i++) {
						_U[i][j] += t*_U[i][k];
					}
				}
				for (i = k; i < _m; i++ ) {
					_U[i][k] = -_U[i][k];
				}
				_U[k][k] = 1.0 + _U[k][k];
				for (i = 0; i < k-1; i++) {
					_U[i][k] = 0.0;
				}
			} 
			else
			{
				for (i = 0; i < _m; i++) {
					_U[i][k] = 0.0;
				}
				_U[k][k] = 1.0;
			}
		}
	}
	
	// If required, generate _V.
	if ( _computeV ) 
	{
		for (k = _n-1; k >= 0; k--) 
		{
			if ((k < nrt) & (e[k] != 0.0)) 
			{
				for (j = k+1; j < nu; j++) {
					double t = 0;
					for (i = k+1; i < _n; i++) {
						t += _V[i][k]*_V[i][j];
					}
					t = -t/_V[k+1][k];
					for (i = k+1; i < _n; i++) {
						_V[i][j] += t*_V[i][k];
					}
				}
			}
			for (i = 0; i < _n; i++) {
				_V[i][k] = 0.0;
			}
			_V[k][k] = 1.0;
		}
	}
	
	// Main iteration loop for the singular values.
	int pp = p-1;
	int iter = 0;
	double eps = pow(2.0,-52.0);
	while (p > 0) {
		k=0;
		unsigned int mode=0;
		
		// Here is where a test for too many iterations would go.
		// This section of the program inspects for negligible elements in the s and e arrays.  
		// On completion the variables mode and k are set as follows.
		
		// mode = 1     if s(p) and e[k-1] are negligible and k<p
		// mode = 2     if s(k) is negligible and k<p
		// mode = 3     if e[k-1] is negligible, k<p, and s(k), ..., s(p) are not negligible (qr step).
		// mode = 4     if e(p-1) is negligible (convergence).
		for (k = p-2; k >= -1; k--)
		{
			if (k == -1) {
				break;
			}
			if (fabs(e[k]) <= eps*(fabs(_S[k]) + fabs(_S[k+1]))) {
				e[k] = 0.0;
				break;
			}
		}
		if ( k == p-2 ) 
		{
			mode = 4;
		} 
		else 
		{
			int ks(p-1); // start from ks == p-1
			for ( ; ks >= k; ks--) {
				if (ks == k) {
					break;
				}
				double t = ( (ks != p) ? fabs(e[ks]) : 0.0) + ( (ks != k+1) ? fabs(e[ks-1]) : 0.0);
				if (fabs(_S[ks]) <= eps*t)  
				{
					_S[ks] = 0.0;
					break;
				}
			}
			if (ks == k) {
				mode = 3;
			} else if (ks == p-1) {
				mode = 1;
			} else {
				mode = 2;
				k = ks;
			}
		}
		k++;
				
		// Perform the task indicated by the selected mode.		
		switch ( mode ) {
			
			case 1: 
			{ 			// Deflate negligible _S[p]
				double f = e[p-2];
				e[p-2] = 0.0;
				for (j = p-2; j >= k; j--) 
				{
					double t = SiMath::triangle(_S[j],f);
					double cs = _S[j]/t;
					double sn = f/t;
					_S[j] = t;
					if (j != k) {
						f = -sn*e[j-1];
						e[j-1] = cs*e[j-1];
					}
					
					// update V 
					if ( _computeV ) 
					{
						for (i = 0; i < _n; i++) 
						{
							t = cs*_V[i][j] + sn*_V[i][p-1];
							_V[i][p-1] = -sn*_V[i][j] + cs*_V[i][p-1];
							_V[i][j] = t;
						}
					}
				}
			}
			break; // end case 1
				
			case 2: 
			{ // Split at negligible _S[k]
				double f = e[k-1];
				e[k-1] = 0.0;
				for (j = k; j < p; j++) 
				{
					double t = triangle(_S[j],f);
					double cs = _S[j]/t;
					double sn = f/t;
					_S[j] = t;
					f = -sn*e[j];
					e[j] = cs*e[j];
					
					if ( _computeU ) 
					{
						for (i = 0; i < _m; i++) {
							t = cs*_U[i][j] + sn*_U[i][k-1];
							_U[i][k-1] = -sn*_U[i][j] + cs*_U[i][k-1];
							_U[i][j] = t;
						}
					}
				}
			}
			break; // end case 2 
				
			case 3: 
			{ // Perform one qr step.
				
				// Calculate the shift.
				double scale = max(max(max(max(fabs(_S[p-1]),fabs(_S[p-2])),fabs(e[p-2])),fabs(_S[k])),fabs(e[k]));
				double sp = _S[p-1]/scale;
				double spm1 = _S[p-2]/scale;
				double epm1 = e[p-2]/scale;
				double sk = _S[k]/scale;
				double ek = e[k]/scale;
				double b = ((spm1 + sp)*(spm1 - sp) + epm1*epm1)/2.0;
				double c = (sp*epm1)*(sp*epm1);
				double shift = 0.0;
				if ((b != 0.0) || (c != 0.0)) {
					shift = sqrt(b*b + c);
					if (b < 0.0) {
						shift = -shift;
					}
					shift = c/(b + shift);
				}
				double f = (sk + sp)*(sk - sp) + shift;
				double g = sk*ek;
				
				// Chase zeros.
				
				for (j = k; j < p-1; j++) 
				{
					double t = SiMath::triangle(f,g);
					double cs = f/t;
					double sn = g/t;
					if (j != k) {
						e[j-1] = t;
					}
					f = cs*_S[j] + sn*e[j];
					e[j] = cs*e[j] - sn*_S[j];
					g = sn*_S[j+1];
					_S[j+1] = cs*_S[j+1];
					
					if ( _computeV ) 
					{
						for (i = 0; i < _n; i++) 
						{
							t = cs*_V[i][j] + sn*_V[i][j+1];
							_V[i][j+1] = -sn*_V[i][j] + cs*_V[i][j+1];
							_V[i][j] = t;
						}
					}
					t = SiMath::triangle(f,g);
					cs = f/t;
					sn = g/t;
					_S[j] = t;
					f = cs*e[j] + sn*_S[j+1];
					_S[j+1] = -sn*e[j] + cs*_S[j+1];
					g = sn*e[j+1];
					e[j+1] = cs*e[j+1];
					
					if ( _computeU && (j < _m-1) ) 
					{
						for (i = 0; i < _m; i++) 
						{
							t = cs*_U[i][j] + sn*_U[i][j+1];
							_U[i][j+1] = -sn*_U[i][j] + cs*_U[i][j+1];
							_U[i][j] = t;
						}
					}
				}
				e[p-2] = f;
				iter++;
			}
				break; // end case 3
				
				// convergence step				
			case 4: 
			{ 
				
				// Make the singular values positive.
				if (_S[k] <= 0.0) 
				{
					_S[k] = (_S[k] < 0.0 ) ? -_S[k] : 0.0;

					if ( _computeV )
					{
						for (i = 0; i <= pp; i++) {
							_V[i][k] = -_V[i][k];
						}
					}
				}
				
				// Order the singular values.
				while (k < pp) 
				{
					if (_S[k] >= _S[k+1])
						break;

					// swap values and columns if necessary
					_S.swap(k,k+1);
					
					if ( _computeV && (k < _n-1)) 
						_V.swapColumns(k,k+1);
					
					if ( _computeU && (k < _m-1) )
						_U.swapColumns(k,k+1);

					k++;
				}
				iter = 0;
				p--;
			}
				break; // end case 4
		} 
	}
}


SiMath::Matrix 
SiMath::SVD::getSingularMatrix() 
{
	unsigned int n = _S.size();
	Matrix A(n,n,0.0);
	// set diagonal elements
	for (int i = 0; i < n; i++) 
	{
		A[i][i] = _S[i];
	}
	
	return A;
}

int 
SiMath::SVD::rank () 
{
	double eps = pow(2.0,-52.0);
	double tol = max(_m,_n)*_S[0]*eps;
	int r = 0;
	for (int i = 0; i < _S.size(); i++) {
		if (_S[i] > tol) {
			r++;
		}
	}
	return r;
}
