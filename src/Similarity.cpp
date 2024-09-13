/*
 *  Similarity.cpp

 *
 *  Created by Gert Thijs on 24/03/06.
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

#include "Similarity.h"

//------------------------------------------------------------------------------------------
// Vector based implementations
//------------------------------------------------------------------------------------------

// Euclidean distance
double
SiMath::similarityByEuclidean(const Vector & d1, const Vector & d2)
{
	// Variables
	unsigned int n(d1.size());
	double sumRS(0.0);
	
	// Euclidean similarity
	for (unsigned int i = 0; i < n; ++i)
	{
		sumRS += (d1[i] - d2[i]) * (d1[i] - d2[i]);
	}
	
	// Return
	return sqrt(sumRS);
}


// weighted Euclidean distance
double
SiMath::similarityByEuclidean(const Vector & d1, const Vector & d2, const Vector & weight)
{
	// Variables
	unsigned int n(d1.size());
	double sumRS(0.0);
	
	// Euclidean similarity
	for (unsigned int i = 0; i < n; ++i)
	{
		sumRS += (d1[i] - d2[i]) * (d1[i] - d2[i]) * weight[i];
	}
	
	// Return
	return sqrt(sumRS);
}


// normalised Euclidean distance
double
SiMath::similarityByNormalisedEuclidean(const Vector & d1, const Vector & d2)
{
	// Variables
	unsigned int n(d1.size());
	double sumRS(0.0);
	double sumR(0.0);
	double sumS(0.0);
	
	// Euclidean similarity
	for (unsigned int i = 0; i < n; ++i)
	{
		sumRS += (d1[i] - d2[i]) * (d1[i] - d2[i]);
		sumR += d1[i] * d1[i];
		sumS += d2[i] * d2[i];
	}
	
	if (sqrt(sumR) + sqrt(sumS) != 0.0)
	{
		return sqrt(sumRS) / (sqrt(sumR) + sqrt(sumS));
	}
	else
	{
		return 0.0;
	}
}


// weighted normalised euclidean distance
double
SiMath::similarityByNormalisedEuclidean(const Vector & d1, const Vector & d2, const Vector & weight)
{
	// Variables
	unsigned int n(d1.size());
	double sumRS(0.0);
	double sumR(0.0);
	double sumS(0.0);
	
	// Euclidean similarity
	for (unsigned int i = 0; i < n; ++i)
	{
		sumRS += (d1[i] - d2[i]) * (d1[i] - d2[i]) * weight[i];
		sumR += d1[i] * d1[i] * weight[i];
		sumS += d2[i] * d2[i] * weight[i];
	}
	
	if (sqrt(sumR) + sqrt(sumS) != 0.0)
	{
		return sqrt(sumRS) / (sqrt(sumR) + sqrt(sumS));
	}
	else
	{
		return 0.0;
	}
}


double
SiMath::similarityByCosine(const Vector & d1, const Vector & d2, const Vector & weight)
{
	// Variables
	unsigned int n(d1.size());
	double sumRS(0.0);
	double sumR(0.0);
	double sumS(0.0);
	
	// Cosine similarity
	for (unsigned int i = 0; i < n; ++i)
	{
		sumRS += d1[i] * d2[i] * weight[i];
		sumR  += d1[i] * d1[i] * weight[i];
		sumS  += d2[i] * d2[i] * weight[i];
	}
	
	if (sumR * sumS != 0.0)
	{
		return sumRS / sqrt(sumR * sumS);
	}
	else
	{
		return -1.0;
	}
}


double
SiMath::similarityByCosine(const Vector & d1, const Vector & d2)
{
	// Variables
	unsigned int n(d1.size());
	double sumRS(0.0);
	double sumR(0.0);
	double sumS(0.0);
	
	// Cosine similarity
	for (unsigned int i = 0; i < n; ++i)
	{
		sumRS += d1[i] * d2[i];
		sumR  += d1[i] * d1[i];
		sumS  += d2[i] * d2[i];
	}
	
	if (sumR * sumS != 0.0)
	{
		return sumRS / sqrt(sumR * sumS);
	}
	else
	{
		return -1.0;
	}
}


double
SiMath::similarityByTanimoto(const Vector & d1, const Vector & d2, const Vector & weight)
{
	// Variables
	unsigned int n(d1.size());
	double sumRS(0.0);
	double sumR(0.0);
	double sumS(0.0);
	
	// Tanimoto similarity
	for (unsigned int i = 0; i < n; ++i)
	{
		sumRS += d1[i] * d2[i] * weight[i];
		sumR  += d1[i] * d1[i] * weight[i];
		sumS  += d2[i] * d2[i] * weight[i];
	}
	
	if (sumR + sumS != sumRS)
	{
		return sumRS / (sumR + sumS - sumRS);
	}
	else
	{
		return 0.0;
	}
}


double
SiMath::similarityByTanimoto(const Vector & d1, const Vector & d2)
{
	// Variables
	unsigned int n(d1.size());
	double sumRS(0.0);
	double sumR(0.0);
	double sumS(0.0);
	
	// Tanimoto similarity
	for (unsigned int i = 0; i < n; ++i)
	{
		sumRS += d1[i] * d2[i];
		sumR  += d1[i] * d1[i];
		sumS  += d2[i] * d2[i];
	}
	
	if (sumR + sumS != sumRS)
	{
		return sumRS / (sumR + sumS - sumRS);
	}
	else
	{
		return 0.0;
	}
}


double
SiMath::similarityByDice(const Vector & d1, const Vector & d2, const Vector & weight)
{
	// Variables
	unsigned int n(d1.size());
	double sumRS(0.0);
	double sumR(0.0);
	double sumS(0.0);
	
	// Dice similarity
	for (unsigned int i = 0; i < n; ++i)
	{
		sumRS += d1[i] * d2[i] * weight[i];
		sumR  += d1[i] * d1[i] * weight[i];
		sumS  += d2[i] * d2[i] * weight[i];
	}
	
	if ( sqrt(sumR) != sqrt(sumS) )
	{
		return (2.0 * sumRS) / ( sqrt(sumR) + sqrt(sumS) );
	}
	else
	{
		return 0.0;
	}
}


double
SiMath::similarityByDice(const Vector & d1, const Vector & d2)
{
	// Variables
	unsigned int n(d1.size());
	double sumRS(0.0);
	double sumR(0.0);
	double sumS(0.0);
	
	// Dice similarity
	for (unsigned int i = 0; i < n; ++i)
	{
		sumRS += d1[i] * d2[i];
		sumR  += d1[i] * d1[i];
		sumS  += d2[i] * d2[i];
	}
	
	if (sqrt(sumR) != sqrt(sumS))
	{
		return (2.0 * sumRS) / (sqrt(sumR) + sqrt(sumS));
	}
	else
	{
		return 0.0;
	}
}


double
SiMath::similarityByPearson(const Vector & d1, const Vector & d2, const Vector & weight)
{
	// Variables
	unsigned int n(d1.size());
	double sumAB(0.0);
	double sumAA(0.0);
	double sumBB(0.0);
	double avgA(0.0);
	double avgB(0.0);
	
	// Pearson similarity
	// Averages
	for (unsigned int i = 0; i < n; ++i)
	{
		avgA += d1[i];
		avgB += d2[i];
	}
	avgA /= n;
	avgB /= n;
	
	for (unsigned int i = 0; i < n; ++i)
	{
		sumAB += (d1[i] - avgA) * (d2[i] - avgB) * weight[i];
		sumAA += (d1[i] - avgA) * (d1[i] - avgA) * weight[i];
		sumBB += (d2[i] - avgB) * (d2[i] - avgB) * weight[i];
	}
	
	if (sumAA * sumBB != 0.0)
	{
		return sumAB / sqrt(sumAA * sumBB);
	}
	else
	{
		return 0.0;
	}
}

double
SiMath::similarityByPearson(const Vector & d1, const Vector & d2)
{
   // Variables
	unsigned int n(d1.size());
	double sumSqrtX = 0.0;
	double sumSqrtY = 0.0;
	double sumXY = 0.0;
	double mX = d1[0];
	double mY = d2[0];
		
	for (unsigned int i=1; i<n; ++i ){
		double sweep = (i - 1.0) / i;
		double dX = d1[i] - mX;
		double dY = d2[i] - mY;
				
		sumSqrtX += dX * dX * sweep;
		sumSqrtY += dY * dY * sweep;
		sumXY += dX * dY * sweep;
		
		mX += dY / i;
		mY += dY / i;
	}
	
	// normalise
	sumSqrtX = sqrt( sumSqrtX / n );
	sumSqrtY = sqrt( sumSqrtY / n );
	
	if ( sumSqrtX == 0 || sumSqrtY == 0 )
		return 0;
	else
		return sumXY/(n * sumSqrtX * sumSqrtY );	
}

//-------------------------------------------------------------------------------------------------------------------
// array-based implentation
//-------------------------------------------------------------------------------------------------------------------
double
SiMath::similarityByEuclidean(unsigned int n, const double * d1, const double * d2, const double * weight)
{
	// Variables
	double sumRS(0.0);
	
	// Euclidean similarity
	for (unsigned int i = 0; i < n; ++i)
	{
		sumRS += (d1[i] - d2[i]) * (d1[i] - d2[i]) * weight[i];
	}
	
	// Return
	return sqrt(sumRS);
}

double
SiMath::similarityByEuclidean(const unsigned int n, const double * d1, const double * d2)
{
	// Variables
	double sumRS(0.0);
	
	// Euclidean similarity
	for (unsigned int i = 0; i < n; ++i)
	{
		sumRS += (d1[i] - d2[i]) * (d1[i] - d2[i]);
	}
	
	// Return
	return sqrt(sumRS);
}


double
SiMath::similarityByNormalisedEuclidean(unsigned int n, const double * d1, const double * d2)
{
	// Variables
	double sumRS(0.0);
	double sumR(0.0);
	double sumS(0.0);
	
	// Euclidean similarity
	for (unsigned int i = 0; i < n; ++i)
	{
		sumRS += (d1[i] - d2[i]) * (d1[i] - d2[i]);
		sumR += d1[i] * d1[i];
		sumS += d2[i] * d2[i];
	}
	
	if (sqrt(sumR) + sqrt(sumS) != 0.0)
	{
		return sqrt(sumRS) / (sqrt(sumR) + sqrt(sumS));
	}
	else
	{
		return 0.0;
	}
}


double
SiMath::similarityByNormalisedEuclidean(unsigned int n, const double * d1, const double * d2, const double * weight)
{
	// Variables
	double sumRS(0.0);
	double sumR(0.0);
	double sumS(0.0);
	
	// Euclidean similarity
	for (unsigned int i = 0; i < n; ++i)
	{
		sumRS += (d1[i] - d2[i]) * (d1[i] - d2[i]) * weight[i];
		sumR += d1[i] * d1[i] * weight[i];
		sumS += d2[i] * d2[i] * weight[i];
	}
	
	if (sqrt(sumR) + sqrt(sumS) != 0.0)
	{
		return sqrt(sumRS) / (sqrt(sumR) + sqrt(sumS));
	}
	else
	{
		return 0.0;
	}
}


double
SiMath::similarityByCosine(unsigned int n, const double * d1, const double * d2)
{
	// Variables
	double sumRS(0.0);
	double sumR(0.0);
	double sumS(0.0);
	
	// Cosine similarity
	for (unsigned int i = 0; i < n; ++i)
	{
		sumRS += d1[i] * d2[i];
		sumR  += d1[i] * d1[i];
		sumS  += d2[i] * d2[i];
	}
	
	if (sumR * sumS != 0.0)
	{
		return sumRS / sqrt(sumR * sumS);
	}
	else
	{
		return -1.0;
	}
}


double
SiMath::similarityByCosine(unsigned int n, const double * d1, const double * d2, const double * weight)
{
	// Variables
	double sumRS(0.0);
	double sumR(0.0);
	double sumS(0.0);
	
	// Cosine similarity
	for (unsigned int i = 0; i < n; ++i)
	{
		sumRS += d1[i] * d2[i] * weight[i];
		sumR  += d1[i] * d1[i] * weight[i];
		sumS  += d2[i] * d2[i] * weight[i];
	}
	
	if (sumR * sumS != 0.0)
	{
		return sumRS / sqrt(sumR * sumS);
	}
	else
	{
		return -1.0;
	}
}


double
SiMath::similarityByTanimoto(unsigned int n, const double * d1, const double * d2, const double * weight)
{
	// Variables
	double sumRS(0.0);
	double sumR(0.0);
	double sumS(0.0);
	
	// Tanimoto similarity
	for (unsigned int i = 0; i < n; ++i)
	{
		sumRS += d1[i] * d2[i] * weight[i];
		sumR  += d1[i] * d1[i] * weight[i];
		sumS  += d2[i] * d2[i] * weight[i];
	}
	
	if (sumR + sumS != sumRS)
	{
		return sumRS / (sumR + sumS - sumRS);
	}
	else
	{
		return 0.0;
	}
}


double
SiMath::similarityByTanimoto(unsigned int n, const double * d1, const double * d2)
{
	// Variables
	double sumRS(0.0);
	double sumR(0.0);
	double sumS(0.0);
	
	// Tanimoto similarity
	for (unsigned int i = 0; i < n; ++i)
	{
		sumRS += d1[i] * d2[i];
		sumR  += d1[i] * d1[i];
		sumS  += d2[i] * d2[i];
	}
	
	if (sumR + sumS != sumRS)
	{
		return sumRS / (sumR + sumS - sumRS);
	}
	else
	{
		return 0.0;
	}
}


double
SiMath::similarityByDice(unsigned int n, const double * d1, const double * d2, const double * weight)
{
	// Variables
	double sumRS(0.0);
	double sumR(0.0);
	double sumS(0.0);
	
	// Dice similarity
	for (unsigned int i = 0; i < n; ++i)
	{
		sumRS += d1[i] * d2[i] * weight[i];
		sumR  += d1[i] * d1[i] * weight[i];
		sumS  += d2[i] * d2[i] * weight[i];
	}
	
	if ( sqrt(sumR) != sqrt(sumS) )
	{
		return (2.0 * sumRS) / ( sqrt(sumR) + sqrt(sumS) );
	}
	else
	{
		return 0.0;
	}
}


double
SiMath::similarityByDice(unsigned int n, const double * d1, const double * d2)
{
	// Variables
	double sumRS(0.0);
	double sumR(0.0);
	double sumS(0.0);
	
	// Dice similarity
	for (unsigned int i = 0; i < n; ++i)
	{
		sumRS += d1[i] * d2[i];
		sumR  += d1[i] * d1[i];
		sumS  += d2[i] * d2[i];
	}
	
	if (sqrt(sumR) != sqrt(sumS))
	{
		return (2.0 * sumRS) / (sqrt(sumR) + sqrt(sumS));
	}
	else
	{
		return 0.0;
	}
}


double
SiMath::similarityByPearson(unsigned int n, const double * d1, const double * d2, const double * weight)
{
	// Variables
	double sumAB(0.0);
	double sumAA(0.0);
	double sumBB(0.0);
	double avgA(0.0);
	double avgB(0.0);
	
	// Pearson similarity
	// Averages
	for (unsigned int i = 0; i < n; ++i)
	{
		avgA += d1[i];
		avgB += d2[i];
	}
	avgA /= n;
	avgB /= n;
	
	for (unsigned int i = 0; i < n; ++i)
	{
		sumAB += (d1[i] - avgA) * (d2[i] - avgB) * weight[i];
		sumAA += (d1[i] - avgA) * (d1[i] - avgA) * weight[i];
		sumBB += (d2[i] - avgB) * (d2[i] - avgB) * weight[i];
	}
	
	if (sumAA * sumBB != 0.0)
	{
		return sumAB / sqrt(sumAA * sumBB);
	}
	else
	{
		return 0.0;
	}
}


double
SiMath::similarityByPearson(unsigned int n, const double * d1, const double * d2)
{
	// Variables
	// initial values
	double sumSqrtX = 0.0;
	double sumSqrtY = 0.0;
	double sumXY = 0.0;
	double mX = d1[0];
	double mY = d2[0];
		
	for (unsigned int i=1; i<n; ++i ){
		double sweep = (i - 1.0) / i;
		double dX = d1[i] - mX;
		double dY = d2[i] - mY;
				
		sumSqrtX += dX * dX * sweep;
		sumSqrtY += dY * dY * sweep;
		sumXY += dX * dY * sweep;
		
		mX += dY / i;
		mY += dY / i;
	}
	
	sumSqrtX = sqrt( sumSqrtX / n );
	sumSqrtY = sqrt( sumSqrtY / n );
	if ( sumSqrtX == 0 || sumSqrtY == 0 )
		return 0;
	else
		return sumXY/(n * sumSqrtX * sumSqrtY );	
}



