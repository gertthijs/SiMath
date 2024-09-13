/*
 *  DivisiveKMeans.cpp
 *
 *  Created by Gert Thijs on 10/05/06.
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


#include "DivisiveKMeans.h"
#include "Similarity.h"
#include "MatrixOperations.h"


SiMath::DivisiveKMeans::DivisiveKMeans(const SiMath::DivisiveKMeansParameters & p) :
	_nbrClusters(0),
	_nbrVar(p.dimension),
	_maxMembers(p.maxMembers),
	_maxIterations(p.maxIterations), 
	_maxDepth(p.maxDepth),
	_latestNumber(0),
	_isConverged(false), 
	_dataSeen(false),
	_centers(),
	_processQueue(),
	_centerMap()
{
	// cluster number 0 is always NULL
	_centerMap[0] = NULL;
}


SiMath::DivisiveKMeans::~DivisiveKMeans()
{
	
}


std::vector<unsigned int> 
SiMath::DivisiveKMeans::cluster(const Matrix & data)
{
	// get dimensions matrix=[mxn]
	if ( _nbrVar != data.nbrColumns() )
		DivisiveKMeansError("Dimension mismatch between data matrix and model parameters.");
	
	unsigned int n = data.nbrRows();
	// create the top level center
	Center c0;
	c0.number = 1;
	c0.coord = SiMath::meanOfAllColumns(data);
	
	_latestNumber = 1;
	c0.members = n; // all datapoints belong to the first cluster
	_processQueue.push(c0);
		
	// new set of numbers
	unsigned int nbr1(0);
	unsigned int nbr2(0);
	
	// set all datapoints to have label 0 
	std::vector<unsigned int> labels(n,1);	
		
	while ( !_processQueue.empty() )
	{
		// get next center from the queue
		c0 = _processQueue.front();
		_processQueue.pop();

		// check if cluster can be splitted
		if ( c0.members > _maxMembers && c0.depth < _maxDepth )
		{
			// compute the next cluster numbers
			++_latestNumber;
			nbr1 = _latestNumber;
			++_latestNumber;
			nbr2 = _latestNumber;
			
			c0.child1 = nbr1;
			c0.child2 = nbr2;
		
			// split the selected cluster
			_splitCluster(data,labels,c0);
		}
		
		// add the center to list of centers
		_centers.push_back(c0);
		std::cerr << "center " << c0.number << " " << c0.members << " " << c0.depth << "  +" << c0.parent << " ->" << c0.child1 << " ->" << c0.child2 << " " << c0.coord;
	}

	// processing queue is empty, post process clustering results
	_nbrClusters = _centers.size();
	
	// create a map of cluster numbers and pointers to centers
	_createMap();
	
	_isConverged = true;
	_dataSeen = false;
	
	return labels;
}


std::vector<unsigned int> 
SiMath::DivisiveKMeans::assign(const Matrix & data)
{
	unsigned int n = data.nbrRows();
	std::vector<unsigned int> labels(n,0);

	if ( _centerMap.empty() )
		_createMap();
	
	for ( int i=0; i<n; ++i ) // loop over all datapoints
	{
		Vector v = data.getRow(i);
		Center * c = _centerMap[1];
		while ( c != NULL )
		{
			if ( c->child1 != 0 && c->child2 != 0 ) 
			{
				// assign point to either child1 or child2
				double d1 = SiMath::similarityByEuclidean(_centerMap[c->child1]->coord, v);
				double d2 = SiMath::similarityByEuclidean(_centerMap[c->child2]->coord, v);
				if ( d1 < d2 )
					c = _centerMap[c->child1];
				else
					c = _centerMap[c->child2];
			}
			else
			{
				// this center is not splitted any further
				// so assign its number to the data point
				labels[i] = c->number;
				c = NULL;
			}
		}
		c = NULL;
	}
	return labels;
}

unsigned int 
SiMath::DivisiveKMeans::assign(const Vector & data)
{
	unsigned int label(0);

	if ( _centerMap.empty() )
		_createMap();

	Center * c = _centerMap[1];
	while ( c != NULL )
	{
		if ( c->child1 != 0 && c->child2 != 0 ) 
		{
			// assign point to either child1 or child2
			double d1 = SiMath::similarityByEuclidean(_centerMap[c->child1]->coord, data);
			double d2 = SiMath::similarityByEuclidean(_centerMap[c->child2]->coord, data);
			if ( d1 < d2 )
				c = _centerMap[c->child1];
			else
				c = _centerMap[c->child2];
		}
		else
		{
			// this center is not splitted any further
			// so assign its number to the data point
			label = c->number;
			c = NULL;
		}
	}
		
	return label;
}



void
SiMath::DivisiveKMeans::setCenter(unsigned int i, const Vector & v, unsigned int p)
{
	// create a new center
	Center c;
	c.number = i;
	c.coord = v;

	if ( p == 0 ) 
		c.depth = 1;
	
	std::map<unsigned int, Center *>::iterator ai=_centerMap.find(p);
	
	if ( ai == _centerMap.end() )
	{
		// not found
		throw(DivisiveKMeansError("trying to add child node before parent node"));
	}
	else
	{
		c.parent = p;
		
		if ( (ai->second)->child1 == 0 ){
			(ai->second)->child1 = i;
		}else if ( (ai->second)->child2 == 0 ) {
			(ai->second)->child2 = i;
		}else{
			throw(DivisiveKMeansError("trying to add child node to parent node where both childs have already been defined."));
		}
		
		c.depth = (ai->second)->depth + 1;
	}
	
	// add centers 
	_centers.push_back(c);
	
	return;
}


//------------------------------------------------------------------------------------
// LOCAL HELPER FUNCTIONS
//------------------------------------------------------------------------------------
void
SiMath::DivisiveKMeans::_splitCluster(const Matrix & data, std::vector<unsigned int> & labels, Center & c)
{

	unsigned int n(labels.size());
	unsigned int m(data.nbrColumns());
	std::map<unsigned int, unsigned int> activeSet;
	std::map<unsigned int, unsigned int>::iterator ai;
	for ( unsigned int i=0; i<n ; ++i )
	{
		if ( labels[i] != c.number )
			continue;
		
		// define random label for active set based on child labels
		activeSet[i] = (rand()&01) ? c.child1 : c.child2;
	}
	
	// create two new cluster centers
	Center c1;
	c1.coord.reset(m);
	c1.number = c.child1;
	c1.depth = c.depth + 1;
	c1.members = 0;
	c1.parent = c.number;
	
	Center c2;
	c2.coord.reset(m);
	c2.number = c.child2;
	c2.depth = c.depth + 1;
	c2.members = 0;
	c2.parent = c.number;

	bool isConverged(false);
	unsigned int iter(0);
	while ( !isConverged && iter<_maxIterations )
	{
		// compute centers
		c1.coord = 0;
		c1.members = 0;
		c2.coord = 0;
		c2.members = 0;
		for ( ai=activeSet.begin(); ai != activeSet.end(); ++ai )
		{
			if ( ai->second == c1.number ) { 
				// add to center 1
				for ( int j=0; j<m; ++j )
					c1.coord[j] += data[ai->first][j];
				
				c1.members++;
			}else{
				// add to center 2
				for ( int j=0; j<m; ++j )
					c2.coord[j] += data[ai->first][j];
				
				c2.members++;
			}
		}
		
		// check if one of the centers does not contain any members;
		if ( c1.members == 0 || c2.members == 0 )
		{
			// nothing happens here, labels stay the parent label and no new centers are added to the queue
			c.child1 = 0;
			c.child2 = 0;
			// c.members = activeSet.size();
			return;
		}
			
		// normalise centers
		for ( int j=0; j<m; ++j )
		{
			c1.coord[j] /= c1.members;
			c2.coord[j] /= c2.members;
		}
		
		// assign each data point to closest center
		bool hasChanged(false);
		for ( ai=activeSet.begin(); ai != activeSet.end(); ++ai )
		{
			SiMath::Vector v = data.getRow(ai->first);
			double d1 = SiMath::similarityByEuclidean(c1.coord, v);
			double d2 = SiMath::similarityByEuclidean(c2.coord, v);
			
			// examine possible combinations of d1 and d2 together with labels
			if ( d1 < d2 && ai->second == c.child2 ) // point moves from c2 to c1
			{
				hasChanged = true;
				ai->second = c.child1;
			}
			else if ( d2 < d1 && ai->second == c.child1 ) // point moves from c1 to c2
			{
				hasChanged = true;
				ai->second = c.child2;		
			}
		}
		
		// check convergence (is converged if no point has changed cluster)
		isConverged = !hasChanged;
		++iter;
	}
	std::cerr << "converged after " << iter << std::endl;
	// assign labels to data points
	for ( ai=activeSet.begin(); ai != activeSet.end(); ++ai )
	{
		labels[ai->first] = ai->second;
	}
	
	// add create centers to either the queue or the set of final centers
	_processQueue.push(c1);
	_processQueue.push(c2);
	
	return;
}



void
SiMath::DivisiveKMeans::_createMap()
{
	_centerMap.clear();
	// cluster with index 0 
	_centerMap[0] = NULL;

	for ( unsigned int i=0; i<_centers.size(); ++i )
	{
		_centerMap[_centers[i].number] = &_centers[i];
	}
}
