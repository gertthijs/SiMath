/*
 *  LeaderClustering.cpp
 *
 *  Created by Gert Thijs on 09/08/06.
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

#include "LeaderClustering.h"
#include "Similarity.h"


SiMath::LeaderClustering::LeaderClustering(const LeaderParameters & p) :
	_leaders(0),
	_nbrClusters(0),
	_nbrVar(p.dimension),
	_threshold(p.threshold)
{
}


std::vector<unsigned int>
SiMath::LeaderClustering::cluster(const SiMath::Matrix & data)
{
	// get dimensions matrix=[mxn]
	if ( _nbrVar != data.nbrColumns() )
	{
		// create a new
		_nbrVar = data.nbrColumns();
	}
	
	// initialisation by giving each data point a random cluster number;
	unsigned int n = data.nbrRows();	
	std::vector<unsigned int> labels(n);
	std::vector<unsigned int> index(n);
	for ( unsigned int i=0; i<n; ++i )
	{	
		labels[i] = 0;
		index[i] = i;
	}

	// permutate the indices
	srand(rand());
	std::random_shuffle(index.begin(),index.end());
	
	// make the first element the first leader
	Leader newCenter;
	newCenter.coord = data.getRow(index[0]);
	newCenter.index = index[0];
	newCenter.members = 1;
	// first clear current set of leaders
	if ( !_leaders.empty() )
		_leaders.empty();

	// add first one
	_leaders.push_back(newCenter);
	_nbrClusters = 1; 
	labels[index[0]] = 0; // set label of first data point
	
	double mseDistance(0.0);	
	for ( unsigned int i=1; i<n; ++i )
	{
		// std::cerr << " " << index[i];
		SiMath::Vector r(data.getRow(index[i]));
		
		double dClosest(HUGE_VALF);
		unsigned int c(_leaders.size());
		for ( unsigned int j = 0; j<_leaders.size(); ++j )
		{
			double d = SiMath::similarityByEuclidean(_leaders[j].coord,r);
			if ( d < dClosest )
			{
				dClosest = d; 
				c = j; //store current cluster number
			}
		}
		
		if ( dClosest < _threshold ) // add to selected cluster
		{
			labels[index[i]] = c;
			_leaders[c].members += 1;
			mseDistance += (dClosest * dClosest);
		}
		else // create a new cluster
		{
			Leader newCenter;
			newCenter.coord = r;
			newCenter.index = index[i];
			newCenter.members = 1;
			_leaders.push_back(newCenter);
			_nbrClusters++;
			labels[index[i]] = _leaders.size() - 1;
			std::cerr << _nbrClusters << "\r";
		}
	}
	
	
	std::cerr << "MSE = " << mseDistance/n << std::endl;
	std::cerr << "Nbr. Clusters = " << _nbrClusters << std::endl;
	return labels;
}



std::vector<unsigned int>
SiMath::LeaderClustering::update(const SiMath::Matrix & data)
{
	// get dimensions matrix=[mxn]
	if ( _nbrVar != data.nbrColumns() )
	{
		throw(LeaderClusteringError("Unable to update clustering when dimension of data does not match."));
	}
	
	// initialisation by giving each data point cluster number 0
	unsigned int n = data.nbrRows();	
	std::vector<unsigned int> labels(n);
	std::vector<unsigned int> index(n);
	for ( unsigned int i=0; i<n; ++i )
	{	
		labels[i] = 0;
		index[i] = i;
	}
	
	// permutate the indices
	srand(rand());
	std::random_shuffle(index.begin(),index.end());
	
	double mseDistance(0.0);	
	for ( unsigned int i=0; i<n; ++i )
	{
		SiMath::Vector r(data.getRow(index[i]));
		
		double dClosest(HUGE_VALF);
		unsigned int c(_leaders.size());
		for ( unsigned int j = 0; j<_leaders.size(); ++j )
		{
			double d = SiMath::similarityByEuclidean(_leaders[j].coord,r);
			if ( d < dClosest )
			{
				dClosest = d; 
				c = j; //store current cluster number
			}
		}
		
		if ( dClosest < _threshold ) // add to selected cluster
		{
			labels[index[i]] = c;
			_leaders[c].members += 1;			
			mseDistance += (dClosest * dClosest);
		}
		else // create a new cluster
		{
			Leader newCenter;
			newCenter.coord = r;
			newCenter.members = 1;
			_leaders.push_back(newCenter);
			_nbrClusters++;
			labels[index[i]] = _leaders.size() - 1;			
		}
	}
	
	
	std::cerr << "MSE: " << mseDistance/n << std::endl;
	
	return labels;
}



std::vector<unsigned int>
SiMath::LeaderClustering::assign(const SiMath::Matrix & data)
{
	// get dimensions matrix=[mxn]
	unsigned int n = data.nbrRows();
	unsigned int m = data.nbrColumns();
	
	if ( _nbrVar != m )
		throw(LeaderClusteringError("Dimension of data and cluster centers do not match."));
	
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
			double d = SiMath::similarityByEuclidean(_leaders[j].coord, v);
			
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
SiMath::LeaderClustering::assign(const SiMath::Vector & data)
{
	// get dimensions matrix=[mxn]
	unsigned int m = data.size();
	
	if ( _nbrVar != m )
		throw(LeaderClusteringError("Dimension of data and cluster centers do not match."));
	
	// initialisation by giving each data point a random cluster number;
	unsigned int label = 0;
	double dOld(HUGE_VALF); // a large distance to start from
	for (unsigned int j=0; j<_nbrClusters; ++j)
	{

		double d = SiMath::similarityByEuclidean(_leaders[j].coord, data);
		
		// update cluster number if distance is closer
		if ( d < dOld )
		{
			label = j;
			dOld = d;
		}
	}
	// std::cerr << " = " << label << std::endl;
	return label;
}


SiMath::Vector
SiMath::LeaderClustering::getLeader(unsigned int i)
{
	if ( i < 0 || i >= _nbrClusters )
	{
		throw(LeaderClusteringError("Cluster index out of range."));
	}
	
	return _leaders[i].coord;
}

void
SiMath::LeaderClustering::reset()
{
	_leaders.clear();
	_nbrClusters = 0;
	return;
}


void 
SiMath::LeaderClustering::setNbrOfClusters(unsigned int n)
{
	if ( n != _nbrClusters )
		_leaders.resize(n);
		
	_nbrClusters = n;
	// set the leaders to an empty one
	Leader l;
	for ( unsigned int i=0; i<n; ++i )
	{
		_leaders[i] = l; 
	}	
	return;
}


void 
SiMath::LeaderClustering::setLeader(unsigned int i, const Vector & v)
{
	if ( i < 0 || i >= _nbrClusters )
		throw(LeaderClusteringError("Cluster index out of range."));
		
	if ( v.size() != _nbrVar )
		throw(LeaderClusteringError("Center vector is not compatible with the number of variables defined."));
	
	_leaders[i].coord = v;	
	_leaders[i].members = 0;
	
	return;
}


void 
SiMath::LeaderClustering::setThreshold(double d)
{
	if ( d < 0 ) throw(LeaderClusteringError("Threshold smaller than 0."));

	_threshold = d;
	return;
}



