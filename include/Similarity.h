/*
 *  Similarity.h
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

#ifndef __SIMATH_SIMILARITY__
#define __SIMATH_SIMILARITY__

#include "Definitions.h"
#include "Vector.h"

namespace SiMath
{
	
	/**
		\defgroup _similarity_ Similarity Measures 
		\brief Collection of similarity and distance measures
	*/
	///\name Distance and similarity measures
	//@{
	/**
		\brief Compute cosine similarity between two vectors
		\ingroup _similarity_
	  \param  d1  first data vector
	  \param  d2  second data vector
	  \return Value between -1 and 1.
	 */
	double similarityByCosine(const Vector & d1, const Vector & d2);
 
	/**
		\brief Compute weighted cosine similarity between two vectors 
	 \ingroup _similarity_
	  \param  d1  first data vector
	  \param  d2  second data vector
		\param weight   weighting vector 
	  \return Value between -1 and 1.
	 */
	double similarityByCosine(const Vector & d1, const Vector & d2, const Vector & weight);
	
	/**
		\brief Compute Tanimoto coefficient between two vectors.
	 \ingroup _similarity_
	  \param  d1  first data vector
	  \param  d2  second data vector
	  \return Value between 0 and 1.
	 */
	double similarityByTanimoto(const Vector & d1, const Vector & d2);
	
	/**
		\brief Compute weighted Tanimoto coefficient between two vectors. 
	 \ingroup _similarity_
	  \param  d1  first data vector
	  \param  d2  second data vector
	  \param weight   weighting vector 
	  \return double
	 */
	double similarityByTanimoto(const Vector & d1, const Vector & d2, const Vector & weight);

	/**
		\brief Compute Pearson correlation coefficient between two vectors. 
	 \ingroup _similarity_
		\param  d1  first data vector
		\param  d2  second data vector
		\return Value
	 */
	double similarityByPearson(const Vector & d1, const Vector & d2);
	
	/**
		\brief Compute weighted Pearson correlation coefficient between two vectors.
	 \ingroup _similarity_
		\param  d1  first data vector
 	  \param  d2  second data vector
	  \param weight   weighting vector 
	  \return double
	 */
	double similarityByPearson(const Vector & d1, const Vector & d2, const Vector & weight);
	
	/**
		\brief Compute Euclidean distance between two vectors.
	 \ingroup _similarity_
	  \param  d1  first data vector
	  \param  d2  second data vector
	  \return double
	 */
	double similarityByEuclidean(const Vector & d1, const Vector & d2);
	
	/**
		\brief  Compute weighted Euclidean distance between two vectors.
	 \ingroup _similarity_
		\param  d1  first data vector
		\param  d2  second data vector
		\param  weight   weighting vector 
		\return double
	 */
	double similarityByEuclidean(const Vector & d1, const Vector & d2, const Vector & weight);
	
	/**
		\brief Compute similarity 
	 \ingroup _similarity_
	  \param  d1  first data vector
	  \param  d2  second data vector
	  \return double
	 */
	double similarityByNormalisedEuclidean(const Vector & d1, const Vector & d2);
	
		/**
			\brief Compute similarity 
		 \ingroup _similarity_
		 \param  d1      first data vector
		 \param  d2      second data vector
		 \param weight   weighting vector 
		 \return double
		 */
	double similarityByNormalisedEuclidean(const Vector & d1, const Vector & d2, const Vector & weight);
	
	/**
		\brief Compute similarity with dice
	 \ingroup _similarity_
	 \param  d1  first data vector
	 \param  d2  second data vector
	 \return double
	 */
	double similarityByDice(const Vector & d1, const Vector & d2);
	
	/**
		\brief Compute similarity 
	 \ingroup _similarity_
	 \param  d1  first data vector
	 \param  d2  second data vector
	 \param weight   weighting vector 
	 \return double
	 */
	double similarityByDice(const Vector & d1, const Vector & d2, const Vector & weight);
	

	// --- array based implementation --
	
	/**
		\brief Compute Euclidean distance 
	 \ingroup _similarity_
 	  \param  n   number of variables 
	  \param  d1  first data vector
	  \param  d2  second data vector
	  \return double
	 */
	double similarityByEuclidean(const unsigned int n, const double * d1, const double * d2);
	
	/**
		\brief Compute weighted Euclidean distance
	 \ingroup _similarity_
		\param  n   number of variables 
		\param  d1  first data vector
		\param  d2  second data vector
		\param weight   weighting vector 
		\return double
	 */
	double similarityByEuclidean(unsigned int n, const double * d1, const double * d2, const double * weight);

	/**
		\brief Compute normalised Euclidean distance between two arrays of length n
	 \ingroup _similarity_
	  \param  n   number of variables 
	  \param  d1  first data vector
	  \param  d2  second data vector
	  \return double
	 */
	double similarityByNormalisedEuclidean(unsigned int n, const double * d1, const double * d2);

	/**
		\brief Compute normalised, weighted Euclidean distance between two arrays of length n
	 \ingroup _similarity_
	  \param  n       number of variables 
	  \param  d1      first data vector
	  \param  d2      second data vector
	  \param weight   weighting vector 
	  \return double
	 */
	double similarityByNormalisedEuclidean(unsigned int n, const double * d1, const double * d2, const double * weight);
		
	
	/**
		\brief Compute dice similarity between two arrays of length n 
	 \ingroup _similarity_
		\param  n   number of variables 
		\param  d1  first data vector
		\param  d2  second data vector
		\param weight   weighting vector 
		\return double
	 */
	double similarityByDice(unsigned int n, const double * d1, const double * d2, const double * weight);

	/**
		\brief Compute dice similarity between two arrays of length n
	 \ingroup _similarity_
		\param  n   number of variables 
		\param  d1  first data vector
		\param  d2  second data vector
		\return double
	 */
	double similarityByDice(unsigned int n, const double * d1, const double * d2);

	/**
		\brief Compute weighted Pearson correlation coefficient  between two arrays of length n
	 \ingroup _similarity_
		\param  n   number of variables 
		\param  d1  first data vector
		\param  d2  second data vector
		\param weight   weighting vector 
		\return double
	 */
	double similarityByPearson(unsigned int n, const double * d1, const double * d2, const double * weight);

	/**
		\brief Compute similarity with Pearson correlation  between two arrays of length n
	 \ingroup _similarity_
		\param  n   number of variables 
		\param  d1  first data vector
		\param  d2  second data vector
		\return double
	 */
	double similarityByPearson(unsigned int n, const double * d1, const double * d2);
	
		/**
		\brief Compute similarity  between two arrays of length n
		 \ingroup _similarity_
		 \param  n   number of variables 
		 \param  d1  first data vector
		 \param  d2  second data vector
		 \return double
		 */
	double similarityByTanimoto(unsigned int n, const double * d1, const double * d2);
	
	
		/**
		\brief Compute similarity  between two arrays of length n
		 \ingroup _similarity_
		 \param  n   number of variables 
		 \param  d1  first data vector
		 \param  d2  second data vector
		 \return double
		 */
	double similarityByCosine(unsigned int n, const double * d1, const double * d2);
	
	/**
		\brief Compute cosine similarity between two arrays of length n
	 \ingroup _similarity_
	 \param  n   number of variables 
	 \param  d1  first data vector
	 \param  d2  second data vector
	 \param weight   weighting vector 
	 \return Value between 0 and 1.0.
	 */
	double similarityByCosine(unsigned int n, const double * d1, const double * d2, const double * weight);
	
	
	/**
		\brief Compute similarity  between two arrays of length n
	 \ingroup _similarity_
	  \param  n   number of variables 
	  \param  d1  first data vector
	  \param  d2  second data vector
	 \param weight   weighting vector 
	 \return double
	 */
	double similarityByTanimoto(unsigned int n, const double * d1, const double * d2, const double * weight);
	
	
	//@}

}; // end of namespace declaration

#endif __SIMATH_SIMILARITY__
