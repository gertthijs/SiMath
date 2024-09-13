/*
 *  KMeans.cpp
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

#include "KMeans.h"
#include "Similarity.h"
#include <map>

// constructors
SiMath::KMeans::KMeans(const KMeansParameters & p) : 
	_nbrClusters(p.nbrClusters),
	_nbrVar(p.dimension),
	_maxIterations(p.maxIterations), 
	_changeLevel(p.changeLevel),
	_isConverged(false), 
	_dataSeen(false), 
	_centers(p.nbrClusters)
{
	//std::cerr << _nbrClusters << std::endl;
	for ( unsigned int i=0; i<_nbrClusters; ++i )
	{
		_centers[i].coord.reset(_nbrVar);
		_centers[i].members = 0;
	}
}

// nothing to be done for the destructor
SiMath::KMeans::~KMeans()
{}

// main clustering procedure
std::vector<unsigned int>
SiMath::KMeans::cluster(const SiMath::Matrix & data)
{
	// get dimensions matrix=[mxn]
	if ( _nbrVar != data.nbrColumns() )
	{
		// create a new
		_nbrVar = data.nbrColumns();
	}
	
	// reset all cluster centra to 0.0
	for ( unsigned int i=0; i<_nbrClusters; ++i )
	{
		_centers[i].coord.reset(_nbrVar);
		_centers[i].members = 0;
	}
	
	unsigned int n = data.nbrRows();
	unsigned int c(0); // local variabel to indicate a cluster number
								 
	// initialisation by giving each data point a random cluster number;
	std::vector<unsigned int> labels(n);
	for ( unsigned int i=0; i<n ; ++i )
	{
		// random cluster index
		c = (unsigned int)(((double)rand()/RAND_MAX) * _nbrClusters);
		// add cluster number
		labels[i] = c;
	}

	// start the clustering procedure	
	unsigned int nIter = 0;
	bool bSame = false;
	
	// as long as the vector of cluster number changes and the number of iterations is smaller than the user defined number
	while ( bSame == false && nIter < _maxIterations )
	{	
		// augment number of iterations
		++nIter;		

		// STEP 1. Compute the cluster centers based on the current cluster assignment
		// reset the cluster centers
		unsigned int pruneCenters = 0;
		for ( c=0; c<_nbrClusters; ++c )
		{
			_centers[c].coord.reset(_nbrVar);
			_centers[c].members = 0;
		}
		
		for ( unsigned int i=0; i<n ; ++i )
		{
			// get the current cluster number
			c = labels[i];
			
			// add all columns
			for ( unsigned int j=0; j<_nbrVar; ++j)
			{
				_centers[c].coord[j] += data[i][j];
			}
			// augment the cluster counter
			_centers[c].members += 1;
		}
			
		// normalise center profiles
		for ( c=0; c<_nbrClusters; ++c )
		{
			if ( _centers[c].members != 0 )
			{
				for ( unsigned int j=0; j<_nbrVar; ++j)
				{
					_centers[c].coord[j] /= _centers[c].members;
				}
			}
			else
			{
				++pruneCenters;
			}
		}
		
		// STEP 2. Assign to each molecule the cluster number (=row index) of the closest cluster center 
		// keep track of possible changes in cluster numbers
		bSame = true;
		unsigned int nbrChanged(0);
		double mseDistance(0.0);
		double silhouette(0.0);
		for ( unsigned int i=0; i<n ; ++i )
		{
			c = 0;
			SiMath::Vector v = data.getRow(i);
			double dClosest(HUGE_VALF); // a large distance to start from
			double dSecond(HUGE_VALF); // next closest cluster 
 			for (unsigned int j=0; j<_nbrClusters; ++j)
			{
				if ( _centers[j].members == 0 ) continue; // skip empty cluster
				
				// compute distance between cluster center c and data point i
				double d = SiMath::similarityByEuclidean(_centers[j].coord, v);
				
				// update cluster number if distance is closer
				if ( d < dClosest )
				{
					c = j;
					dSecond = dClosest;
					dClosest = d;
				}
				else if ( d < dSecond )
				{
					dSecond = d;
				}
			}

			// check if the cluster number has changed
			if ( c != labels[i] )
			{
				labels[i] = c;
				bSame = false;
				++nbrChanged;
			}
			
			silhouette += (dSecond - dClosest)/dSecond;
			// mseDistance += (dClosest * dClosest);
			mseDistance += dClosest;
		}
		
		std::cerr << nIter << "\t=> labels changed " << nbrChanged;
		std::cerr << "\t=> Distance = " << mseDistance/n;
		std::cerr << "\t=> Silhouette = " << silhouette/n << std::endl;
		
		if ( nbrChanged < (_changeLevel * n) )
		{
			// fraction of compounds changed within this step is less than given level
			std::cerr << "Fraction of compounds that switched from cluster is less than " << _changeLevel << std::endl; 
			bSame = true;
		}
	} // end of while loop of main clustering 
	
	if ( nIter < _maxIterations )
	{	
		std::cerr << "Converged after " << nIter << " out of " << _maxIterations << " iterations." << std::endl << std::endl;
		_isConverged = true;
	}
	else
	{	
		std::cerr << "No convergence. Procedure stopped after maximal (" << _maxIterations << ") iterations was reached." << std::endl << std::endl;
		_isConverged = false;
	}
	
	// update cluster numbers based on the pruning 
	unsigned int nonEmpty(0);
	std::map<unsigned int, unsigned int> mapping;
	for (unsigned int i=0; i<_nbrClusters; ++i)
	{
		if ( _centers[i].members == 0 )
		{
			std::cerr << "pruning cluster: " << i << std::endl;
			// copy the last non-empty cluster to this row 
			for (unsigned int j=_nbrClusters-1; j>i; --j)
			{
				if ( _centers[j].members != 0 )
				{
					for ( unsigned int jj=0; jj<_nbrVar; ++jj)
					{
						_centers[i].coord[jj] = _centers[j].coord[jj];
						_centers[j].coord[jj] = 0.0;
					}
					_centers[i].members = _centers[j].members;
					_centers[j].members = 0;
					mapping[i] = j;
					++nonEmpty;
					break;
				}
			}
		}
		else
		{
			mapping[i] = i;
			++nonEmpty;
		}
	}
	
	if ( nonEmpty != _nbrClusters )
	{
		std::cerr << "Retained " << nonEmpty << " from " << _nbrClusters << " clusters." << std::endl;
		for ( unsigned int i=0; i<n; ++i )
			labels[i] = mapping[labels[i]];
		
		_nbrClusters = nonEmpty;
	}

	//  update labels of 
	return labels;
}


std::vector<unsigned int>
SiMath::KMeans::assign(const SiMath::Matrix & data)
{
	// get dimensions matrix=[mxn]
	unsigned int n = data.nbrRows();
	unsigned int m = data.nbrColumns();

	if ( _nbrVar != m )
		throw(KMeansError("Dimension of data and cluster centers do not match."));
	
	
	// initialisation by giving each data point a random cluster number;
	std::vector<unsigned int> labels(n);
	for ( unsigned int i=0; i<n ; ++i )
	{
		unsigned int c = 0;
		double dOld(HUGE_VALF); // a large distance to start from
		Vector v(data.getRow(i));
		for (unsigned int j=0; j<_nbrClusters; ++j)
		{
			// compute distance between cluster center c and data point i
			double d = SiMath::similarityByEuclidean(_centers[j].coord, v);
			
			// update cluster number if distance is closer
			if ( d < dOld )
			{
				c = j;
				dOld = d;
			}
		}
		// update label	
		labels[i] = c;
	}
	
	return labels;
}


unsigned int
SiMath::KMeans::assign(const SiMath::Vector & data)
{
	// get dimensions matrix=[mxn]
	unsigned int m = data.size();

	if ( _nbrVar != m )
		throw(KMeansError("Dimension of data and cluster centers do not match."));
	
	unsigned int label = 0;
	double dOld(HUGE_VALF); // a large distance to start from
	for (unsigned int j=0; j<_nbrClusters; ++j)
	{
		// compute distance between cluster center c and data point i
		double d = SiMath::similarityByEuclidean(_centers[j].coord, data);
		
		// update cluster number if distance is closer
		if ( d < dOld )
		{
			label = j;
			dOld = d;
		}
	}
	return label;
}


SiMath::Vector
SiMath::KMeans::getCenter(unsigned int i)
{
	if ( i < 0 || i >= _nbrClusters )
	{
		throw(KMeansError("Cluster index out of range."));
	}
	
	return _centers[i].coord;
}


void 
SiMath::KMeans::setNbrOfClusters(unsigned int n)
{
	if ( n != _nbrClusters )
		_centers.resize(n);
	_nbrClusters = n;
	for ( unsigned int i=0; i<_nbrClusters; ++i)
	{
		_centers[i].coord.reset(_nbrVar);
		_centers[i].members = 0;
	}
	_isConverged = false;
}

void
SiMath::KMeans::setChangeLevel(double d)
{
	if ( d < 0 || d > 1 )
		throw(KMeansError("Wrong change level defined. Should be between 0 and 1."));
	
	_changeLevel = d;
}


void 
SiMath::KMeans::setCenter(unsigned int i, const Vector & v)
{
	if ( i < 0 || i >= _nbrClusters )
		throw(KMeansError("Cluster index out of range."));
		
	if ( v.size() != _nbrVar )
		throw(KMeansError("Center vector is not compatible with the number of variables defined."));

	_centers[i].coord = v;	
	_centers[i].members = 0;
	return;
}

