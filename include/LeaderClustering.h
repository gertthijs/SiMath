/*
 *  LeaderClustering.h
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

#ifndef  __SIMATH_LEADERCLUSTERING__
#define  __SIMATH_LEADERCLUSTERING__

#include "Definitions.h"
#include "Matrix.h"
#include "Vector.h"

namespace SiMath
{
	
	/**
		\class  LeaderClusteringError
		\brief  Exception class
	 
		This exception class is derived from the std::runtime_error class and is used to 
		indicate that something went wrong during the leader clustering algorithm.
	 \ingroup _error_handling_
	 */	
	class SIMATH_EXCEPTION LeaderClusteringError : public SiMathError
	{
		public:
			LeaderClusteringError(const std::string & e) : SiMathError("LeaderClusteringError : " + e ) {};
	};


	/**
		\struct LeaderParameters
		\brief Struct to set the parameters of the leader clustering algorithm
		\ingroup _model_parameters_

		Usage:
		\code
		SiMath::LeaderParameters p;
		p.dimension = A.nbrColumns();
		p.threshold = 6.5;
		\endcode
	 */
	struct LeaderParameters {
		unsigned int dimension;   ///< Dimension of the data space
		double threshold;         ///< Threshold on the distance 
		
		/**
			\brief Default parameters:
			- dimension = 0,
			- threshold = 0.0
		 */
		LeaderParameters() : dimension(0), threshold(0.0) {};
	};

	/**
		\class LeaderClustering
		\brief Implenentation of leader clustering
	 \ingroup _clustering_
	 
		This class implements the leader clustering algorithm, which is a one pass 
		clustering algorithm. It starts by randomly selecting a data point and makes 
		it a cluster center (or leader). In the next steps a new data point is selected
		and the distance to the closest cluster is computed. If the distance to this 
		closest leader is smaller than a predefined threshold the point is added to this 
		cluster otherwise a new cluster is started of which this point will be the leader.
		This procedure is repeated until al points have been processed.
	 
		The main advantage of this method is its simplicity and its speed of computation. 
		The major drawbacks of this approach are the dependance on the order of the data points
		and the need for the cutoff value to be defined.
	 */
	class LeaderClustering
	{
		public:
	
		/**
			\brief Constructor
			\param p Reference to a predefined set of parameters
		 */
		LeaderClustering(const LeaderParameters & p);
	
		/**
			\brief Destructor
		 */
		~LeaderClustering() {};
	
		/**
			\brief Main clustering step
	 
			\param data Reference to a Matrix that holds the data points
			\return Vector with the cluster labels of all data points
		 */
		std::vector<unsigned int> cluster(const Matrix & data);
	
		/**
			\brief Clustering step that starts from the currently defined leaders
	 
			\param data Reference to a Matrix that holds the data points
			\return A Vector with the cluster labels of all data points
		 */
		std::vector<unsigned int> update(const Matrix & data);
	
		/**
			\brief Compute the cluster labels of new data
	 
			\param data Reference to a Matrix that holds the data points
			\return Vector with the cluster labels of all data points
		 */
		std::vector<unsigned int> assign(const Matrix & data);
	
		/**
			\brief Compute the cluster label of a new data point
	 
			\param data Reference to a Vector that holds the data point
			\return The cluster label of the data point
		 */
		unsigned int assign(const Vector & data);
	
		/**
			\brief Get the current number of clusters
		 */
		unsigned int getNbrOfClusters() { return _nbrClusters; };
	
		/**
			\brief Get the i-th leader as a SiMath::Vector 
		 */
		Vector getLeader(const unsigned int i);
	
		/** 
			\brief Get the index in the original matrix of the i-th leader
		 */
		unsigned int getLeaderIndex(const unsigned int i) { return _leaders[i].index;};
	
		/**
			\brief Reset the clustering to an empty set of leaders.
		 */
		void reset();

		/**
			\brief Set the i-th leader from a SiMath::Vector
		 
			This function is useful to load predefined leaders. It will only 
			work if the number of clusters has been preset.
		 
			\param i The index of the cluster
			\param v The leader 
			\throws LeaderClusteringError if i is out of bounds.
		 */
		void setLeader(unsigned int i, const Vector & v);

		/**
			\brief Set the number of clusters.
		 
		 This functions sets the expected number of clusters and 
		 also resets the current leaders to zero. This function is only
		 useful if the centers are going to be set directly.
		 
		 \param n Number of clusters
		 */
		void setNbrOfClusters(unsigned int n);

		/**
			\brief Set the threshold on the distance below which a data point is added to a cluster
		 
			\param d Threshold
			\throws LeaderClusteringError if the threshold is smaller than 0.
		 */
		void setThreshold(double d);
	
	private:
		
		/**
			\struct Leader 
			\brief Local representation of a cluster leader
			\param coord Vector with the coordinates of the leader
			\param index Cluster number of the leader
			\param members Number of data points assigned to this cluster
		 */
		struct Leader {
			Vector coord;
			unsigned int index;
			unsigned int members;
			
			Leader() : coord(0), index(0), members(0) {};
		};
	
	
		std::vector<Leader> _leaders;   ///< Storage of cluster centers or leaders
	
		unsigned int _nbrClusters;      ///< Number of clusters found
		unsigned int _nbrVar;           ///< dimension of the data matrix
		double _threshold;              ///< threshold on the distance to start a new cluster 
	
};


}; // end of namespace declaration


#endif   __SIMATH_LEADERCLUSTERING__
